

options(echo=TRUE)
library("grid")
library("ggplot2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
library("reshape2", lib.loc="/data/ohler/Scott/Cluster_R_Packages/")
args <- commandArgs(trailingOnly = TRUE)
score <- function(x){
  x  = ((x-min(x))/(max(x)-min(x)))
  return(x)
}
frame <- read.table(args[1])
ordered_frame <- frame#[order(-frame$V2),]
ordered_frame_nohead <- ordered_frame[ ,-1]
rev_ordered_frame <- as.data.frame(t(ordered_frame_nohead))
averaged1 <- rowMeans(rev_ordered_frame)
#averaged_frame <- data.frame(x = seq(-1000, 999, 10), "Divergent TFIIB"=averaged1)
#melted_averaged_frame <- melt(averaged_frame, id="x")

#tail(melted_averaged_frame)
#head(melted_averaged_frame)

frame2 <- read.table(args[2])
ordered_frame2 <- frame2#[order(-frame$V2),]
ordered_frame2_nohead <- ordered_frame2[ ,-1]
rev_ordered_frame2 <- as.data.frame(t(ordered_frame2_nohead))
averaged2 <- rowMeans(rev_ordered_frame2)
#averaged_frame2 <- data.frame(x = seq(-1000, 999, 10), "Non-Divergent TFIIB"=averaged2)
#melted_averaged_frame2 <- melt(averaged_frame2, id="x")

frame3 <- read.table(args[3])
ordered_frame3 <- frame3#[order(-frame$V2),]
ordered_frame3_nohead <- ordered_frame3[ ,-1]
rev_ordered_frame3 <- as.data.frame(t(ordered_frame3_nohead))
averaged3 <- rowMeans(rev_ordered_frame3)
#averaged_frame3 <- data.frame(x = seq(-1000, 999, 5), "Divergent TFIIB"=averaged3)
#melted_averaged_frame3 <- melt(averaged_frame3, id="x")

frame4 <- read.table(args[4])
ordered_frame4 <- frame4#[order(-frame$V2),]
ordered_frame4_nohead <- ordered_frame4[ ,-1]
rev_ordered_frame4 <- as.data.frame(t(ordered_frame4_nohead))
averaged4 <- rowMeans(rev_ordered_frame4)
#averaged_frame4 <- data.frame(x = seq(-1000, 999, 5), "Non-Divergent TFIIB"=averaged4)
#melted_averaged_frame4 <- melt(averaged_frame4, id="x")


averaged_frame <- data.frame(x = seq(-2000, 1999, 10),  "all_other_DHS" = averaged3, "promoter"=averaged1,  "selected_DHS"=averaged2, "background"=averaged4, check.names=FALSE)
melted_averaged_frame <- melt(averaged_frame, id="x")
head(melted_averaged_frame)

tiff(filename = paste0(args[1],"_multi.tiff"), width = 1800, height = 1600, res = 5)
ggplot(melted_averaged_frame, aes(x, value, colour=variable)) + geom_line(size=100) + scale_color_manual(values=c("selected_DHS" = "#D55E00", "promoter" = "#009E73", "all_other_DHS" = "#56B4E9", "background" = "#F0E442")) + xlab("Distance From DHSs midpoints") + ylab("peak coverage score") + coord_cartesian(xlim = c(-1500,1500)) + theme(axis.line=element_blank(), axis.text.x=element_text(size=1000, face="bold"), axis.text.y=element_text(size=1000, face="bold"), axis.ticks.y=element_line(size=40), axis.ticks.x=element_line(size=40), axis.ticks.margin=unit(1, "cm"), axis.ticks.length=unit(3, "cm"), plot.margin=unit(c(20,10,10,10), "cm"), panel.border=element_blank(), panel.grid=element_blank(), panel.margin=element_blank(), panel.background=element_blank(), legend.position=c(1,1), legend.key.size=unit(50, "cm"), legend.text=element_text(size=750), legend.key=element_blank(), legend.justification=c(1,1), strip.background=element_blank(), legend.box.just="right")

#ggplot(melted_averaged_frame, aes(x, value, colour="Divergent TFIIB"), linetype="solid") + geom_line(size=100) + geom_line(aes(colour="Non-Divergent TFIIB"), linetype="solid", size=100, data=melted_averaged_frame2) + scale_color_manual(values=c("Divergent TFIIB" = "#0072B2", "Non-Divergent TFIIB" = "#E69F00")) + xlab("Distance From Cluster Mode") + ylab("Average TSS Score") + theme(axis.line=element_blank(), axis.text.x=element_text(size=1000, face="bold"), axis.text.y=element_text(size=1000, face="bold"), axis.ticks.y=element_line(size=40), axis.ticks.x=element_line(size=40), axis.ticks.margin=unit(1, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.length=unit(3, "cm"), plot.margin=unit(c(10,10,10,10), "cm"), panel.border=element_blank(), panel.grid=element_blank(), panel.margin=element_blank(), panel.background=element_blank(), legend.position=c(0,1), legend.key.size=unit(75, "cm"), legend.text=element_text(size=1000), legend.key=element_blank(), legend.justification=c(0,1), strip.background=element_blank(), legend.box.just="left")
dev.off()
# + coord_cartesian(ylim=c(0, 90)) axis.title.x=element_blank(), axis.title.y=element_blank()
#+ geom_line(aes(colour=variable), linetype="dotted", size=50, data=melted_averaged_frame4)
#+ scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted"))
#+ scale_color_manual(values=c("#D55E00", "#E69F00", "#0072B2", "#56B4E9"))
