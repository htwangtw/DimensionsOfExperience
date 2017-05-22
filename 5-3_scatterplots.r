#load data

#load ggplot 2
require(ggplot2)
scatterplot <- function(data, X, Y, x_lable, y_lable, x_ticks_lims, y_ticks_lims, filename) {

#assign the plot to a variable
p2<- (
  # scatter plot main data
  # the data need to be in dataframe format. x and y can refer directly to the variables
  ggplot(data, aes(x=X, y=Y
 )) 
  # dot size, color and shape etc.
  + geom_point(shape=19, size = 2, colour="grey50") 
  # 95% confidence interval and regression line
  + geom_smooth(method=lm, colour="black", size = 1, fullrange = TRUE)
  # title and axis lable
  # + ggtitle("Blob")
  + xlab(x_lable)
  + ylab(y_lable)
  # limits of the axises
  + xlim(x_ticks_lims[1], x_ticks_lims[2]) 
  + ylim(y_ticks_lims[1], y_ticks_lims[2])
  
) 

# display the plot, this is the default look
p2

# add APA format theme
p2 + theme_classic() +theme(
  text = element_text(size=18),
  axis.line.x = element_line(colour = "black"),
  axis.line.y = element_line(colour = "black")
  )

#save the plot
ggsave(filename, width = 3.5, height = 3.3, dpi = 300)
}


data_MANOVA<-read.csv('./Results/MANOVA.csv', sep=',', header = TRUE)
data_COPE<-read.csv('./Results/DMN16_FSL_EV.csv', sep=',', header = TRUE)

scatterplot(data = data_COPE, X = data_COPE$SCCA_BOOTS_1_DMN16, Y = data_COPE$median,
            x_lable = "Positive Habitual", y_lable = "Functional Connectivity",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-0.75,0.75),
            filename = './Results/Graphs_paper/Fig4_DualReg_CCA1.png')

scatterplot(data = data_MANOVA, X = data_MANOVA$SCCA_BOOTS_1_DMN16, Y = data_MANOVA$FAC2_CogTask_SCCAsample,
            x_lable = "Positive Habitual", y_lable = "Executive Control",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-4,4),
            filename = './Results/Graphs_paper/Fig4_Task_EXE_CCA1.png')

scatterplot(data = data_MANOVA, X = data_MANOVA$SCCA_BOOTS_1_DMN16, Y = data_MANOVA$FAC1_QD_SCCAsample,
            x_lable = "Positive Habitual", y_lable = "Affective Disturbance",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-4,4),
            filename = './Results/Graphs_paper/Fig4_AD_CCA1.png')

scatterplot(data = data_COPE, X = data_COPE$SCCA_BOOTS_2_DMN16, Y = data_COPE$median,
            x_lable = "Spontaneous Off-Task", y_lable = "Functional Connectivity",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-0.75,0.75),
            filename = './Results/Graphs_paper/Fig4_DualReg_CCA2.png')

scatterplot(data = data_MANOVA, X = data_MANOVA$SCCA_BOOTS_2_DMN16, Y = data_MANOVA$WM_Minus_CRT_EFF,
            x_lable = "Spontaneous Off-Task", y_lable = "1-back - 0-back",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-2.5,2.5),
            filename = './Results/Graphs_paper/Fig4_NBack_CCA2.png')

scatterplot(data = data_MANOVA, X = data_MANOVA$SCCA_BOOTS_2_DMN16, Y = data_MANOVA$FAC3_CogTask_SCCAsample,
            x_lable = "Spontaneous Off-Task", y_lable = "Generation of Information",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-4,4),
            filename = './Results/Graphs_paper/Fig4_Task_GEN_CCA2.png')

scatterplot(data = data_MANOVA, X = data_MANOVA$SCCA_BOOTS_2_DMN16, Y = data_MANOVA$FAC1_QD_SCCAsample,
            x_lable = "Spontaneous Off-Task", y_lable = "Affective Disturbance",
            x_ticks_lims = c(-4,4), y_ticks_lims = c(-4,4),
            filename = './Results/Graphs_paper/Fig4_AD_CCA2.png')




