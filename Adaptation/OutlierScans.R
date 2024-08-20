setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/")
dropbox_dir<-file.path("/Users/ruthrivkin/Library/CloudStorage/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

library(ggplot2)

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "black"),
            axis.line = element_line(linewidth=1), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=15, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,face="plain",size=17),
            axis.title.x=element_text(vjust=0.1,face="plain",size=17),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15))

library(LEA)
library(tidyverse)
pop.info <- read.csv(file.path(dropbox_dir, "Datasheets/24.05.24_master.csv"))
pop.info <- pop.info[order(pop.info$BEDorder),]


#PCAAdapt
library(pcadapt)

data.pca <- "PB_24.05.28_ldpruned.bed"
pca.ad <- read.pcadapt(data.pca, type = "bed")

#Find pop structure first, find K
findk <- pcadapt(input = pca.ad, K = 20) #k should be large to start
plot(findk, option = "screeplot", K =10) #looks like k=3 is the best fit

#Plot scores
plot(findk, option = "scores", pop = pop.info$Subpop) #lConsistent with plink pca. Going to stick with K = 2

plot(findk, option = "scores", i = 3, j = 4, pop = pop.info$Subpop) #NW splits out with PC3

#Detect outliers
outliers <- pcadapt(pca.ad, K = 3)
summary(outliers)

#definitely outliers
plot(outliers , option = "manhattan")
plot(outliers, option = "qqplot")
hist(outliers$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#Plot test stat D
plot(outliers, option = "stat.distribution")


#Choose cutoff point
library(qvalue)
#FDR 10%
qval <- qvalue(outliers$pvalues)$qvalues
alpha <- 0.01
outliers.qv <- which(qval < alpha)
length(outliers.qv) #724


#Benjamini-Hochberg Procedure 
padj <- p.adjust(outliers$pvalues, method="BH")
alpha <- 0.01
outliers.BH <- which(padj < alpha)
length(outliers.BH) #724, the same snps are id'd


#snmf outlier detection
project=load.snmfProject("PB_24.05.28_ldpruned.snmfProject")

p <- snmf.pvalues(project, entropy=TRUE, ploidy=2, K=4)
padj.snmf <- p.adjust(p$pvalues, method="BH")
alpha <- 0.01
snmf.BH <- which(padj.snmf < alpha)
length(snmf.BH) #1791

#create dataframe with contig information
# load map file
snp_info <- read.table("PB_24.05.28_ldpruned.bim")

# pull out only CHROM and POS
snp_info1 <- snp_info %>%
  dplyr::select(V1, V4) %>%
  rename(Contig = V1, 
         Pos = V4) %>%
  mutate(Order = c(1:464294))

head(snp_info1)

#Merge with candidate snps

#snmf
snmf.BH <- as.data.frame(snmf.BH)
names(snmf.BH) <- "Order"

snmf.outliers <- merge(snmf.BH, snp_info1)
write.csv(snmf.outliers, file.path(dropbox_dir, "snmfoutliers.csv"))
#pcadapt
outliers.BH <- as.data.frame(outliers.BH)
names(outliers.BH) <- "Order"

pcadapt.outliers <- merge(outliers.BH, snp_info1)
head(pcadapt.outliers)
write.csv(pcadapt.outliers, file.path(dropbox_dir, "pcadaptoutliers.csv"))


#Check for overlap
overlap <- merge(pcadapt.outliers, snmf.outliers)
length(overlap$Order) #67


