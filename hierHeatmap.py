import pandas as pd

from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from mpl_toolkits.axes_grid1 import make_axes_locatable

def hier_clust(df):
    labels = list(df.index) #rows categories
    variables = list(df.columns) #column categories

    row_dist = pd.DataFrame(squareform(pdist(df, metric='euclidean')), columns=labels, index=labels)

    row_clusters = linkage(pdist(df, metric='euclidean'), method='complete')
    pd.DataFrame(row_clusters, 
                 columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
                 index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])

    # reorder rows with respect to the clustering
    row_dendr = dendrogram(row_clusters, labels=labels, orientation='left')
    df_rowclust = df.ix[row_dendr['leaves']]
    return df_rowclust
    
def heatmap(df):
    # plot
    # fig = plt.figure(figsize=(3,df.shape[1]*4))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.matshow(df, vmin=-1, vmax=1, interpolation='nearest', cmap=plt.cm.RdBu_r, aspect='equal')
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "15%", pad="20%")
    plt.colorbar(im, cax=cax)
    ax.set_yticks(range(df.shape[0]))
    ax.set_xticks(range(df.shape[1]))
    ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(list(df.columns), rotation='75', ha='left')
    ax.set_yticklabels(list(df.index))
    plt.show()
    return None