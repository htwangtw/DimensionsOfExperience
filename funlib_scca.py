import numpy as np
import warnings

def quadratic1(a,b,c):
    '''
    Solve an quadratic equation
    '''
    if b**2-4*a*c < 0:
        x = np.nan
        # print('No solution') 
    elif b**2-4*a*c==0: 
        x = -b/(2*a)
        # print('Solution is',x) 
    else: 
        x = np.array(((-b+np.sqrt(b**2-4*a*c))/(2*a), (-b-np.sqrt(b**2-4*a*c))/(2*a)))
        # print('Solutions are', x)
    return x

def is_outliers(data, m=2.5):
    '''
    Check outliers.
    Checking if the giving number is 2.5 standard distributions away from the mean.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(data, axis=0)
        sd = np.sqrt(np.nanmean((data - mean)**2, axis=0))

    is_outliers = abs(data - mean) > m * sd
    return is_outliers

def imputedata(data, strategy='mean', missing=False):
    '''
    two impute strategys
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(data, axis=0)
        sd = np.sqrt(np.nanmean((data - mean)**2, axis=0))
    sign = np.sign(data - mean)
    is_out = is_outliers(data, m=2.5)
    data[is_out] = np.nan
    
    if strategy == '2sd':
        # impute as +-2sd m
        # reduce the change in distribution. 
        for i in range(data.shape[1]):
            if missing:
                sign[np.isnan(sign)] = 0 #missing data will be imputed as mean
            ind_nan = np.where(np.isnan(data[:,i]))
            data[ind_nan,i] = mean[i] + (sd[i] * 2 * sign[ind_nan,i])

    if strategy == 'mean':
        #impute as mean
        for i in range(data.shape[1]):
            ind_nan = np.where(np.isnan(data[:,i]))
            if missing: #missing data will be imputed as mean
                data[ind_nan,i] = mean[i]
            else: #missing data will be left as nan
                data[ind_nan,i] = mean[i] * abs(sign[ind_nan,i])
    return data

def demean(Y):
    '''
    Demean a 2-D data matrix.
    '''
    S = Y.sum(axis=0) / Y.shape[0]
    Y -= S[np.newaxis, :]
    var = (Y ** 2).sum(axis=0)
    var[var == 0] = 1
    Y /= var
    return Y


def mean_nonzero(data, axis):
    '''
    calculate the non zero elements in a 2-D matrix
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        temp_sum = np.sum(data!=0, axis=axis)
        temp_sum [temp_sum==0] = 1 
        output = np.sum(data, axis=axis)/temp_sum
    return output

from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.rinterface as ri
import pandas as pd

'''
Sparse Canonical Correlation Analysis (SCCA)
'''

def SCCA_r(X, Y, n_components, penX, penY):
    pandas2ri.activate()
    PMA = importr('PMA')
    df_X = pd.DataFrame(X)
    df_Y = pd.DataFrame(Y)
    rmat_X = pandas2ri.py2ri(df_X)
    rmat_Y = pandas2ri.py2ri(df_Y)
    ri.globalenv['X'] = rmat_X
    ri.globalenv['Y'] = rmat_Y

    out = PMA.CCA(x=X, z=Y, K=n_components, niter =100, standardize=False, penaltyx=penX, penaltyz=penY)
    u = pandas2ri.ri2py(out[0]) #the bridge function in python use 0-index now (R is 1-index)
    v = pandas2ri.ri2py(out[1])
    
    return u, v

def canonical_scores(X, Y, u, v):
    '''
    Calculate the canonical socres
    '''
    cano_scores_X = X.dot(u)
    cano_scores_Y = Y.dot(v)

    return cano_scores_X, cano_scores_Y

def csv_CanonicalScores(X, Y, u, v, subj_id, filename):
    '''
    Save the canonical socres as a csv file
    '''
    cano_scores_X, cano_scores_Y = canonical_scores(X, Y, u, v)
    n_components = cano_scores_X.shape[1]
    ID_CanoScores = np.column_stack((subj_id, cano_scores_X, cano_scores_Y))
    headers = (',').join(['IDNO']+['x_Factor_%i'%(i+1) for i in range(n_components)] + ['y_Factor_%i'%(i+1) for i in range(n_components)])
    with open(filename, 'wb') as f:
        f.write(headers+'\n')
        np.savetxt(f, ID_CanoScores, fmt='%10.8f', delimiter=',')
    return ID_CanoScores

from os.path import expanduser
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def SCCA_Output_Sheet(filename, X_roi_keys, Y_keys, u, v):
    '''
    Visualisation
    '''
    x_loadings = u
    y_loadings = v

    n_components = x_loadings.shape[1]
    n_connects = int(x_loadings.shape[0])
    # get the number of area from the number of connectivity network strength values we have.
    temp_areas = quadratic1(1, -1, -2*n_connects)
    n_areas = int(temp_areas[temp_areas>0])

    # Transform the flatten weights for network strength back to matrices for visualisation
    idx = np.triu_indices(n_areas, 1)
    corr_mat = np.zeros((n_areas, n_areas, n_components))
    for i in range(n_components):
        this_mat = np.zeros((n_areas, n_areas))
        this_mat[idx] = x_loadings[:, i]
        corr_mat[..., i] = this_mat + this_mat.T
    
    # load labels
    region_labels = np.load(expanduser(X_roi_keys))
    keys = np.load(expanduser(Y_keys))[1:]

    # plot the thing
    fig = plt.figure(figsize=(20, 10* n_components))
    # fig = plt.figure()
    fig.subplots_adjust(left=0.3, right=0.8, hspace = 0.2, wspace = 0.4)
    for i in range(n_components):

        ax = fig.add_subplot(n_components, 2, i*2 + 1)

        brain = ax.matshow(corr_mat[..., i], vmin=-0.9, vmax=0.9, cmap=plt.cm.RdBu_r)
        ax.set_xticks(np.arange(n_areas))
        ax.set_xticklabels(region_labels, rotation=90, fontsize='large')
        ax.set_yticks(np.arange(n_areas))
        ax.set_yticklabels(region_labels, fontsize='large')
        ax.plot([-0.5, n_areas-0.5], [-0.5, n_areas-0.5], ls='--', c='.3')
        # cb_brain = fig.colorbar(brain, fraction=0.046, pad=0.04)

        behav_ax = fig.add_subplot(n_components, 2, (i + 1)*2)
        behav_arr = np.zeros((len(keys),1))
        behav_arr.flat[:y_loadings.shape[0]] = y_loadings[:, i]
        behav = behav_ax.matshow(behav_arr, vmin=-0.9, vmax=0.9, cmap=plt.cm.RdBu_r)
        behav_ax.set_yticks(np.arange(len(keys)))
        behav_ax.set_yticklabels(keys, fontsize='large')
        behav_ax.set_xticklabels(' ')
        cb_behave = fig.colorbar(behav, fraction=0.046, pad=0.04)
        # fig.tight_layout()

    plt.savefig(filename + '_heatmaps.pdf')
    # plt.show()


from sklearn.linear_model import LinearRegression

def expVar(X, Y, penX, penY):
    '''
    Calculate the explained variable percentage by numbers of componet with a given penalty pair
    '''
    limit_exp_var = Y.shape[1] #save for later
    exp_var_X = []
    exp_var_Y = []
    for i in range(1, limit_exp_var+1):
        loadings = SCCA_r(X, Y, i, penX, penY)
        '''
        calculate the coefficent of determination (R square): 
        the proportion of the variance (fluctuation) of one variable that is predictable 
        from the other variable. In other words, the ratio of the explained variation to the total
        variation.
        '''
        P = loadings[0]
        lr = LinearRegression(fit_intercept=False)
        lr.fit(P, X.T)
        rec_X = lr.coef_.dot(P.T)  #a.dot(b) equals np.dot(a,b)
        exp_var_X.append(1 - (np.var(X - rec_X) / np.var(X)))
        Q = loadings[1]
        lr = LinearRegression(fit_intercept=False)
        lr.fit(Q, Y.T)
        rec_Y = lr.coef_.dot(Q.T)
        exp_var_Y.append(1 - np.var(Y - rec_Y) / np.var(Y))

    plt.close('all')
    plt.figure()
    plt.plot(np.arange(limit_exp_var) + 1, exp_var_X, label='Brain exp var')
    plt.plot(np.arange(limit_exp_var) + 1, exp_var_Y, label='Behavioral exp var')
    plt.ylim(-0.1, 1)
    plt.xlim(1, limit_exp_var)
    plt.legend(loc='lower right')

    np.set_printoptions(precision=3,suppress=True,linewidth=1000)
    x = np.transpose(np.array([range(1, limit_exp_var+1)] +[exp_var_X]+[exp_var_Y]))
    print ''
    print 'Explained data proportion'
    print '    n', '  exp_brain', '  exp_behaviour'
    print x
    plt.show()


import scikits.bootstrap as boots
def boots_SCCA_sk(X, Y, n_components, penX, penY):
    '''
    bootstrap v1
    '''
    data = (X, Y)
    def SCCA_boot(X,Y):
        '''
        The scikit bootstarp function can only accept data as varaibles, so we set the parameters here
        '''
        loadings = SCCA_r(X, Y, n_components, penX, penY)
        np.save('tmp/Random_bootSample', loadings)
        # calculate the canonical covarates: the diagonal of the matrix "u.T*X.T*Z*v"
        corvar = (np.corrcoef(np.dot(X,loadings[0]).T,np.dot(Y,loadings[1]).T)[n_components:, 0:n_components]).diagonal()
        return corvar

    
    #bootsrtap here
    ci = boots.ci(data, statfunction=SCCA_boot, method='pi') 
    print('Confident interval: 95%')
    print('High:', ci[1])
    print('Low:', ci[0])

    return ci


def boots_SCCA(X, Y, n_components, penX, penY):
    '''
    bootstrap v2
    '''
    n_samples =1000

    from numpy.random import randint
    # select a random index
    bootsInd = randint(X.shape[0],size=(n_samples, X.shape[0]))

    for i, I in enumerate(bootsInd):
        cur_X = X[I,:]
        cur_Y = Y[I,:]
        cur_loadings, cur_cors = SCCA_r(cur_X,cur_Y, n_components, penX, penY) # run SCCA
        cur_cors = (np.corrcoef(np.dot(cur_X,cur_loadings[0]).T,np.dot(cur_Y,cur_loadings[1]).T)[n_components:, 0:n_components]).diagonal()
        if i ==0:
            corr_master = np.zeros((cur_cors.shape[0], n_samples))
        corr_master[...,i] = cur_cors 

    #sort the canonical correlation score
    corr_master.sort(axis=1)

    # confidnece interval
    alpha = 0.05
    ind_low = int(n_samples*alpha/2)
    ind_high = int(n_samples - n_samples*alpha/2)
    ci_corr =  (corr_master[:,ind_low], corr_master[:,ind_high])
    np.set_printoptions(precision=3)
    print('Cannonical Correlation Coefficent: %s' %str(corr))
    print('Confident interval: 95%')
    print('Low: %s'%str(ci_corr[0]))
    print('High: %s'%str(ci_corr[1]))

    plt.hist(corr_master.sum(axis=0), bins='auto')  
    plt.plot([ci_corr[0].sum(), ci_corr[0].sum()], [0, 120], '--', lw=2)
    plt.plot([ci_corr[1].sum(), ci_corr[1].sum()], [0, 120], '--', lw=2)
    plt.plot([corr.sum(), corr.sum()], [0, 120], 'k-', lw=2)

    plt.title('Bootstrap distribution')
    plt.show()

    return ci_corr