setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/DemographicHistory/")
dropbox_dir<-file.path("~/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

getwd()

library(ggplot2)
c10 <- c(
  "dodgerblue2", "maroon", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6")# lt purple
ng <- theme(aspect.ratio=0.7,panel.background = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(linewidth=0.5), 
            axis.ticks=element_line(color="black"),
            axis.text=element_text(color="black",size=12, margin = 0.5), 
            axis.title=element_text(color="black",size=1), 
            axis.title.y=element_text(vjust=1,size=14),
            axis.title.x=element_text(vjust=0.1,size=14),
            axis.text.x=element_text(size=12),
            axis.text.y=element_text(size=12),
            # legend.position = "top", legend.direction="horizontal", 
            legend.text=element_text(size=10), 
            legend.title = element_text(size=12))



#Plot GONE results
#Investigate demographic history
#Plot GONE results
library(tidyverse)

# Load data
BB <- read.delim("GONE/36Contigs/Output_Ne_BB", header=T, skip = 1)
DS <- read.delim("GONE/36Contigs/Output_Ne_DS", header=T, skip = 1)
FB <- read.delim("GONE/36Contigs/Output_Ne_FB", header=T, skip = 1)
GB <- read.delim("GONE/36Contigs/Output_Ne_GB", header=T, skip = 1)
LS <- read.delim("GONE/36Contigs/Output_Ne_LS", header=T, skip = 1)
MC <- read.delim("GONE/36Contigs/Output_Ne_MC", header=T, skip = 1)
NB <- read.delim("GONE/36Contigs/Output_Ne_NB", header=T, skip = 1)
NW <- read.delim("GONE/36Contigs/Output_Ne_NW", header=T, skip = 1)
SB <- read.delim("GONE/36Contigs/Output_Ne_SB", header=T, skip = 1)
SH <- read.delim("GONE/36Contigs/Output_Ne_SH", header=T, skip = 1)


# Quick plot of the means (cut off at 200 generation time)
legend.colors <- c(
  "BB" = "dodgerblue2", 
  "DS" = "maroon", # red
  "FB" = "green4",
  "GB" = "#6A3D9A", # purple
  "LS" = "#FF7F00", # orange
  "MC" = "gold1",
  "NB" = "skyblue2",
  "NW" = "#FB9A99", # lt pink
  "SB" = "palegreen2",
  "SH" = "#CAB2D6")

ggplot()+
  geom_line(data=BB, aes(x=(Generation),y=Geometric_mean, color="BB"), lwd=1)+
  geom_line(data=DS, aes(x=(Generation),y=Geometric_mean, color="DS"), lwd=1)+
  geom_line(data=FB, aes(x=(Generation),y=Geometric_mean, color="FB"), lwd=1)+
  geom_line(data=GB, aes(x=(Generation),y=Geometric_mean, color="GB"), lwd=1)+
  geom_line(data=LS, aes(x=(Generation),y=Geometric_mean, color="LS"), lwd=1)+
  geom_line(data=MC, aes(x=(Generation),y=Geometric_mean, color="MC"), lwd=1)+
  geom_line(data=NB, aes(x=(Generation),y=Geometric_mean, color="NB"), lwd=1)+
  geom_line(data=NW, aes(x=(Generation),y=Geometric_mean, color="NW"), lwd=1)+
  geom_line(data=SB, aes(x=(Generation),y=Geometric_mean, color="SB"), lwd=1)+
  geom_line(data=SH, aes(x=(Generation),y=Geometric_mean, color="SH"), lwd=1)+
  ng +
  labs(y = "Effective Population Size", x = "Generation", color = "Subpopulation") +
  scale_color_manual(values = legend.colors) +
  xlim(0,200) +
  ylim(0,20000)
ggsave(file.path(dropbox_dir, "Figures/Gone_NE.pdf"), width = 6.78, height = 5.3, dpi = 300)


#Plot with generations
gen=11.5 
ggplot()+
  geom_line(data=BB, aes(x=(Generation*gen),y=Geometric_mean, color="BB"), lwd=1)+
  geom_line(data=DS, aes(x=(Generation*gen),y=Geometric_mean, color="DS"), lwd=1)+
  geom_line(data=FB, aes(x=(Generation*gen),y=Geometric_mean, color="FB"), lwd=1)+
  geom_line(data=GB, aes(x=(Generation*gen),y=Geometric_mean, color="GB"), lwd=1)+
  geom_line(data=LS, aes(x=(Generation*gen),y=Geometric_mean, color="LS"), lwd=1)+
  geom_line(data=MC, aes(x=(Generation*gen),y=Geometric_mean, color="MC"), lwd=1)+
  geom_line(data=NB, aes(x=(Generation*gen),y=Geometric_mean, color="NB"), lwd=1)+
  geom_line(data=NW, aes(x=(Generation*gen),y=Geometric_mean, color="NW"), lwd=1)+
  geom_line(data=SB, aes(x=(Generation*gen),y=Geometric_mean, color="SB"), lwd=1)+
  geom_line(data=SH, aes(x=(Generation*gen),y=Geometric_mean, color="SH"), lwd=1)+
  ng +
  labs(y = "Effective Population Size", x = "Years ago", color = "Subpopulation") +
  scale_color_manual(values = legend.colors) +
  xlim(0,1750) +
  ylim(0,8000)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_recentyears.pdf"), width = 6.78, height = 5.3, dpi = 300)



# Use the 100 iterations to make a 95% confidence interval (this is based off of Kardos et al. 2023's script: https://github.com/martykardos/KillerWhaleInbreeding/blob/main/FigureCode/rCode_Fig_ED_1.R)

library(scales)
library(matrixStats)

# BB first (high arctic)
# load all the iteration files and put it in a matrix
files <- paste("GONE/36Contigs/TEMPORARY_FILES_BB/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(BB <- ggplot()+
  # High Arctic
  geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
  geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
  # Cutting off at 150 generations
  xlim(0,150*gen)+
  ylim(0,15000)+
  theme_bw()+
  xlab("Years ago")+
  ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
  ggtitle("Baffin Bay") + 
  theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_BB_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)

#### DS
files <- paste("GONE/36Contigs/TEMPORARY_FILES_DS/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(DS <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,10000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Davis Strait") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_DS_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


#### Foxe Basin
files <- paste("GONE/36Contigs/TEMPORARY_FILES_FB/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(FB <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,7000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Foxe Basin") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_FB_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)



### GB
files <- paste("GONE/36Contigs/TEMPORARY_FILES_GB/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(GB <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,20000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Gulf of Boothia") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_GB_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


#### LS
files <- paste("GONE/36Contigs/TEMPORARY_FILES_LS/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(LS <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,18000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Lancaster Sound") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_LS_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


#### MC
files <- paste("GONE/36Contigs/TEMPORARY_FILES_MC/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(MC <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,40000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("M'Clintock Channel") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_MC_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


#### NB
files <- paste("GONE/36Contigs/TEMPORARY_FILES_NB/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(NB <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,30000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Northern Beaufort Sea") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_NB_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)

#### NW
files <- paste("GONE/36Contigs/TEMPORARY_FILES_NW/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(NW <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,7000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Norwegian Bay") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_NW_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


### SB

files <- paste("GONE/36Contigs/TEMPORARY_FILES_SB/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(SB <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,33000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Southern Beaufort Sea") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_SB_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)

#### SH
files <- paste("GONE/36Contigs/TEMPORARY_FILES_SH/outfileLD_TEMP/outfileLD_",1:100,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Cet CI for the recent 200 generations
NeCI <- matrix(NA,nrow=200,ncol=2)
for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

# Median
Ne_med <- as.data.frame(rowMedians(NeMat[1:200,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Define generation time b/c we want to plot in years
gen=11.5 #Taylor et al. 2007

(SH <- ggplot()+
    # High Arctic
    geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#92C5DE", alpha=0.5)+
    geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#4D9BD3", lwd=1.5)+
    # Cutting off at 150 generations
    xlim(0,150*gen)+
    ylim(0,10000)+
    theme_bw()+
    xlab("Years ago")+
    ylab(expression(paste("Effective population size (",italic(""*N*"")[e],")",sep="")))+
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(plot.margin = margin(0.2,0.6,0.2,0.3, "cm")) +  
    ggtitle("Southern Hudson Bay") + 
    theme(plot.title=element_text(margin=margin(t=10,b=-20, l = 10)))
)

ggsave(file.path(dropbox_dir, "Figures/Gone_NE_SH_withCI.pdf"), width = 6.78, height = 5.3, dpi = 300)


#Arrange all in a a tile
library(ggpubr)
#Update title and axis to make everything fit together
BB <- BB + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
DS <- DS + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
FB <- FB + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
GB <- GB + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
LS <- LS + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
MC <- MC + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
NB <- NB + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
NW <- NW + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
SB <- SB + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))
SH <- SH + theme(plot.title=element_text(margin=margin(t=10,b=-12, l = 5), size=10)) + ylab(expression(paste(italic(""*N*"")[e],sep="")))


ggarrange(BB, DS, FB, GB, LS, MC, NB, NW, SB, SH,
          ncol = 3, nrow = 4)
ggsave(file.path(dropbox_dir, "Figures/Gone_NE_allsubpops_withCI.pdf"))

