options(echo=TRUE)
library("grid")
library("ggplot2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
library("reshape2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
args <- commandArgs(trailingOnly = TRUE)

frame <- read.table(args[1])
frame_nohead <- frame[ ,-1]
rev_frame_nohead <- as.data.frame(t(frame_nohead))
averaged <- rowMeans(rev_frame_nohead)
averaged_frame <- data.frame(x = seq(-2000, 1999, 10), Signal=averaged)
melted_averaged_frame <- melt(averaged_frame, id="x")
tiff(filename = paste0(args[1],".tiff"), width = 1800, height = 1600, res = 5)
ggplot(melted_averaged_frame, aes(x, value, colour=variable)) + geom_line(size=100) + scale_color_manual(values=c("#0072B2")) + xlab("Distance From Cluster Mode") + ylab("Average TSS Score") + theme(axis.line=element_blank(), axis.text.x=element_text(size=1000, face="bold"), axis.text.y=element_text(size=1000, face="bold"), axis.ticks.y=element_line(size=40), axis.ticks.x=element_line(size=40), axis.ticks.margin=unit(1, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.length=unit(3, "cm"), plot.margin=unit(c(20,10,10,10), "cm"), panel.border=element_blank(), panel.grid=element_blank(), panel.margin=element_blank(), panel.background=element_blank(), legend.position="none", legend.key.size=unit(50, "cm"), legend.text=element_text(size=1000), legend.key=element_blank(), legend.justification=c(1,1), strip.background=element_blank(), legend.box.just="right")
dev.off()


#coord_cartesian(xlim = c(-as.numeric(args[2]), as.numeric(args[3])), ylim = c(-as.numeric(args[4]), as.numeric(args[5])))
#, ylim = c(0,0.16)
