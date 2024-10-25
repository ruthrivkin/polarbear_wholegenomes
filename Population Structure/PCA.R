setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/")
dropbox_dir<-file.path("/Users/ruthrivkin/Library/CloudStorage/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

library(ggplot2)
library(tidyverse)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border=element_rect(fill = NA, color = "NA"),
            axis.line = element_line(linewidth = 0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=10, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=12),
            axis.title.x=element_text(vjust=0.1,size=12),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10))

#Generate pca files for pca



#Check for batch effects, min/maf filtering but nothing else
# read in data
pca <- read_table("PCA/PB_24.10.18_batch.eigenvec", col_names = TRUE)
eigenval <- scan("PCA/PB_24.10.18_batch.eigenval")
eigenval<- round(eigenval,2)

batch.data <- read.csv("PCA/pb_batch_effects.csv")

#set identifiers
pca$Subpop <- batch.data$PopulationCode
pca$SequenceYear <- as.factor(batch.data$SequenceYear)
names(pca)[1] <- "SampleName"
head(pca)

#Plot

cols = my.colors(length(unique(pca$Subpop)))

#subpops first
(pcaplot <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC1, y = PC2,color = Subpop), show.legend = TRUE ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",eigenval[1]," %)"),
         y = paste0("Principal component 3 (",eigenval[2]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)
ggsave(file.path(dropbox_dir, "Figures/PCA_plink.pdf"), width = 6.78, height = 5.3, dpi = 300)


#Check sequence year, no batch effects present
#PC1 & PC2
(pcaplot.batch <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC1, y = PC2,color = Subpop, shape = SequenceYear), show.legend = TRUE, size = 2 ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",eigenval[1]," %)"),
         y = paste0("Principal component 2 (",eigenval[2]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)

ggsave(file.path(dropbox_dir, "Figures/PCA_batch_downsampled.pdf"), width = 6.78, height = 5.3, dpi = 300)

#PC1 & PC3
(pcaplot.batch <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC1, y = PC3,color = Subpop, shape = SequenceYear), show.legend = TRUE, size = 2 ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",eigenval[1]," %)"),
         y = paste0("Principal component 3 (",eigenval[3]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)

#PC2 & PC3
(pcaplot.batch <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC2, y = PC3,color = Subpop, shape = SequenceYear), show.legend = TRUE, size = 2 ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 2 (",eigenval[2]," %)"),
         y = paste0("Principal component 3 (",eigenval[3]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)

#PC1 & PC4
(pcaplot.batch <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC1, y = PC4,color = Subpop, shape = SequenceYear), show.legend = TRUE, size = 2 ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",eigenval[1]," %)"),
         y = paste0("Principal component 4 (",eigenval[4]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)




#Check LD and HWE pruned samples next.. these look the same
# read in data
pca <- read_table("PB_24.05.28_ld_HWE.eigenvec", col_names = TRUE)
eigenval <- scan("PB_24.05.28_ld_HWE.eigenval")
eigenval<- round(eigenval,2)


# sort out the pca data

# set names
names(pca)[1] <- "Subpop"
names(pca)[2] <- "SampleName"


cols = my.colors(length(unique(pca$Subpop)))

(pcaplot <- ggplot(data = pca) +
    geom_point(mapping = aes(x = PC1, y = PC2,color = Subpop), show.legend = TRUE ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(x = paste0("Principal component 1 (",eigenval[1]," %)"),
         y = paste0("Principal component 2 (",eigenval[2]," %)")) + 
    scale_color_manual(values = cols) +
    ng
)
ggsave(file.path(dropbox_dir, "Figures/PCA_plink.pdf"), width = 6.78, height = 5.3, dpi = 300)


