options(echo=TRUE)
library("grid")
library("ggplot2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
library("reshape2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
args <- commandArgs(trailingOnly = TRUE)
score <- function(x){
  x  = ((x-min(x))/(max(x)-min(x)))
  return(x)
}
frame <- read.table(args[1], header=FALSE)
ordered_frame <- frame#[order(-frame$V2),]
cleaned_ordered_frame <- ordered_frame#[,-2]
header <- cleaned_ordered_frame$V1
cleaned_ordered_frame_nohead <- cleaned_ordered_frame[ ,-1]
norm_cleaned_ordered_frame <- as.data.frame(apply(cleaned_ordered_frame_nohead, 1, score))
norm_cleaned_ordered_frame$x <- seq(-1000, 999, 100)
melted_norm_cleaned_ordered_frame <- melt(norm_cleaned_ordered_frame, id="x")
tiff(filename = paste0(args[1],".tiff"), width = 1800, height = 1600, res = 5)
#jpeg(filename = paste0(args[1],".jpg"))
ggplot(melted_norm_cleaned_ordered_frame, aes(x, variable)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high = "#000000") + xlab("Distance From Peak Center") + ylab("Increasing Peak Width") + theme(axis.text.x=element_text(size=1000, face="bold"), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_line(size=40), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", axis.ticks.length=unit(3, "cm"), panel.background=element_blank(), legend.box="horizontal", legend.key.size=unit(30, "cm"), legend.margin=unit(50, "cm"), plot.title=element_text(size=1500), plot.margin=unit(c(0,0,10,0), "cm"))
dev.off()
