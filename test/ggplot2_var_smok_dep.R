#! /path/to/Rscript
args <- commandArgs(TRUE)
print(args)

variable=args[1]
folder=args[2]
var_detail=args[3]
library(stringr)	 
require("ggplot2")
require("igraph")

timesteps <- 2000
n.measures <- 3
no <- ""
runbin = 30

#dir(folder)


df_plot<- read.table(str_c(c(str_c(c(folder,"/",variable),collapse=""), "csv"),  collapse = "."),header=TRUE)

cols.scheme <- c("#3368A5","#55247A","#D33C3E")

myrunmean <- function(var,bin, ncond = 6) { 
	plen <- length(var)
	pcondlen <- plen / ncond
	olen <- plen - (ncond*(bin-1))
	ocondlen <- olen / ncond
	y <- c(1:olen)
	index <- array(rep(1:plen,times = ncond), c(pcondlen, ncond))
	for (x in 1:ncond) { y[1:ocondlen + ocondlen * (x-1)] <- running.mean( var[index[,x]] , bin )  }
	return (y)
}

#plotting_range=c(1,2000)#c(15:985)#c(runbin/2:dim(df_plot)[1]/3-runbin/2-1)#c(15:985)

myopts <- opts(
axis.title.x = theme_text(colour = "#626262", vjust = -0.2, face="bold", size=18),
axis.title.y = theme_text(colour = "#626262", vjust = 0.3, face="bold", size=18, angle=90),
axis.text.x=theme_text(colour = "#626262",size = 16),	
axis.text.y=theme_text(colour = "#626262",size = 16),	
legend.key.size = unit(2.0, "lines"),
legend.position = c(0.850,0.85),
legend.text=theme_text(colour = "#626262",face="bold", size=16),
legend.title= theme_blank(), #theme_text(hjust = -0,face="bold", size=16),
legend.background = theme_blank() # theme_rect(colour = "gray", fill="gray")
)

#f2.p3 <- ggplot(df_plot, aes( timestep[timestep %in% plotting_range ],mean ))
f2.p3 <- ggplot(df_plot, aes( timestep,mean ))

plot3 <- f2.p3 + geom_ribbon(fill = "#62626280",aes(
	ymin = lower, 
	ymax = upper )) + 
	geom_line(size=2.5, aes(color=condition)) +
	scale_color_manual("Model Dynamics", values =  c(cols.scheme),
	breaks = c("smoker_dyn", "smoker_infl", "smoker_full") , 
	labels = c("network","influence", "coupled")) + 
	labs(x ="Timesteps", y = str_c(c(var_detail, variable, folder),  collapse = " ")) +
	#coord_cartesian(ylim = c(-.1,1.2))+#,xlim = c(0,1000))+
	#coord_cartesian(xlim = c(0,1000))+
	myopts

pdf(str_c(c(str_c(c(folder,"/",variable,"_",folder),collapse=""), "pdf"),  collapse = "."), width = 9, height =6)
print(plot3)
dev.off()
