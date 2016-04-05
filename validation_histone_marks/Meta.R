options(echo=TRUE)
library("grid")
library("ggplot2")
library("reshape2")
library("digest")
args <- commandArgs(trailingOnly = TRUE)

frame <- read.table(args[1])
frame_nohead <- frame[ ,-1]
rev_frame_nohead <- as.data.frame(t(frame_nohead))
averaged <- rowMeans(rev_frame_nohead)
averaged_frame <- data.frame(x = seq(-1000, 999, 100), Signal=averaged)
melted_averaged_frame <- melt(averaged_frame, id="x")
tiff(filename = paste0(args[1],".tiff"), width = 1800, height = 1600, res = 5)
ggplot(melted_averaged_frame, aes(x, value, fill = variable)) + geom_bar(stat = "identity", position = "identity") + theme(axis.text.x=element_text(size=300), axis.text.y=element_text(size=300), axis.ticks.x=element_line(size=40), axis.title=element_blank(), legend.position = "none", axis.ticks.length=unit(3, "cm"))
dev.off()


#coord_cartesian(xlim = c(-as.numeric(args[2]), as.numeric(args[3])), ylim = c(-as.numeric(args[4]), as.numeric(args[5])))
