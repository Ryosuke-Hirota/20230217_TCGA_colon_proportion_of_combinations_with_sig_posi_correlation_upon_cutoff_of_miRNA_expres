# This script is to draw line graph about percentage of combinations with significant positive correlation upon cutoff miRNA expression and number of samples without expression  
# 2023/02/17 made

# make new directory
setwd("C:/Rdata")
dir.create("20230217_TCGA_colon_proportion_of_combinations_with_sig_posi_correlation_upon_cutoff_of_miRNA_expression")

# import table about mean of miRNA and transcript expression
# this table is located at "https://github.com/Ryosuke-Hirota/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA"
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
mean.df <-read.table("TCGA_colon_confirmation_of_expression_level_of_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# delete unnecessary rows and columns, and do log2
mean.df <-mean.df[,c(1:5,9,13)]
mean.df <-mean.df[!is.na(mean.df[,3]),]
mean.df[,6] <-log2(mean.df[,6])
mean.df[,7] <-log2(mean.df[,7])

# set cutoff value about miRNA or transcript expression
mean.df <-mean.df[mean.df[,6]>=0,]
mean.df[mean.df[,3]>0&mean.df[,4]<0.05,8] <-"red"
mean.df[mean.df[,3]<=0|mean.df[,4]>=0.05,8] <-"black"

# import table about number of samples without miRNA or transcript expression
# this table is loated at "https://github.com/Ryosuke-Hirota/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_l"
setwd("C:/Rdata/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_level")
zero.sample.df <-read.table("TCGA_colon_table_about_number_of_sample_without_expression.txt",sep="\t",header = T,stringsAsFactors = F)
zero.sample.df <-zero.sample.df[,c(1,2,8)]

# merge table about mean and table about number of samples witout expression
mean.zero.df <-merge(mean.df,zero.sample.df,by=c("miRNA","transcript"))
mean.zero.df <-mean.zero.df[,c(1:7,9,8)]

# set cutoff value for number of sample without expression level
cutoff <-seq(0,280,20)

# make table for writing percentage
sm <-as.data.frame(matrix(nrow = length(cutoff),ncol = 5))

# calculate percentage upon cutoff miRNA expression and number of samples without expression
for (i in 1:length(cutoff)) {
  cutoff.df <-mean.zero.df[mean.zero.df[,8]<=cutoff[i],]
  sm[i,1] <-cutoff[i]
  sm[i,2] <-nrow(cutoff.df[cutoff.df[,3]>0&cutoff.df[,4]<0.05,])
  sm[i,3] <-nrow(cutoff.df)
  sm[i,4] <-sm[i,2]/sm[i,3]*100
  sm[i,5] <-min(cutoff.df[,5])
}

# draw line graph
setwd("C:/Rdata/20230217_TCGA_colon_proportion_of_combinations_with_sig_posi_correlation_upon_cutoff_of_miRNA_expression")
pdf("line_graph_about_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_miRNA_expression_and_number_of_sample_without_expression.pdf")
plot(sm[,1],sm[,4],xlab="cutoff for number of sample without expression",ylab="percentage of combinations with sig posi correlation",
     pch=19,type="b",ylim = c(60,75))
text(sm[,1],sm[,4], paste0(sm[,2],"/",sm[,3]), adj=c(0.5,-0.5),cex = 0.6)
text(sm[,1],sm[,4]+0.1, sm[,5], pos=3,cex = 0.6)
dev.off()

# draw volcano plot
for (i in 1:length(cutoff)) {
  cutoff.df <-mean.zero.df[mean.zero.df[,8]<=cutoff[i],]
  pdf(paste0("volcano_plot_about_TCGA_colon_correlation_between_expression_of_transcript_and_miRNA_cutoff_",cutoff[i],".pdf"))
  plot(cutoff.df[,3],-log10(cutoff.df[,4]),col=cutoff.df[,9],pch=19,xlab = "correlation coefficient",ylab = "-log10(p.value)",
       main = paste0("sig.posi.cor = ",nrow(cutoff.df[cutoff.df[,9]=="red",])," other = ",nrow(cutoff.df)-nrow(cutoff.df[cutoff.df[,9]=="red",])))
  abline(h=1.30103,v=0,lty=2)
  legend("topleft",legend =c("other","r>0, p<0.05"),col=unique(cutoff.df[,9]),pch=19)
  dev.off()
}

