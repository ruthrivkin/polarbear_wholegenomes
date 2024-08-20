setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/")
dropbox_dir<-file.path("/Users/ruthrivkin/Library/CloudStorage/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

library(ggplot2)
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

#Load pixy files

library(dplyr)

pi<-read.table("pixy/withbed_pi.txt",sep="\t",header=T)
head(pi)
summary(pi$avg_pi)
hist(pi$avg_pi)

dxy<-read.table("pixy/withbed_dxy.txt",sep="\t",header=T)
head(dxy)

fst<-read.table("pixy/withbed_fst.txt",sep="\t",header=T)
head(fst)

##PI
# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
chroms <- unique(pi$chromosome)
chrOrder <- sort(chroms)
pi$chrOrder <- factor(pi$chromosome,levels=chrOrder)

length(chroms)

# Plot pi for each population found in the input file
# Saves a copy of each plot in the working directory
if("avg_pi" %in% colnames(pi)){
  pops <- unique(pi$pop)
  for (p in pops){
    thisPop <- subset(pi, pop == p)
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Pi for population", p))+
      labs(x="Position of window start", y="Pi")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("piplot_", p,".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

# Plot Dxy for each combination of populations found in the input file
# Saves a copy of each plot in the working directory
if("avg_dxy" %in% colnames(inp)){
  # Get each unique combination of populations
  pops <- unique(inp[c("pop1", "pop2")])
  for (p in 1:nrow(pops)){
    combo <- pops[p,]
    thisPop <- subset(inp, pop1 == combo$pop1[[1]] & pop2 == combo$pop2[[1]])
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_dxy, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Dxy for", combo$pop1[[1]], "&", combo$pop2[[1]]))+
      labs(x="Position of window start", y="Dxy")+
      theme(legend.position = "none")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("dxyplot_", combo$pop1[[1]], "_", combo$pop2[[1]],".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Dxy not found in this file")
}


#convert files to long format
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

pixy_folder <- "pixy"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
head(pixy_df)


# create a custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed

# plotting summary statistics across all chromosomes
pixy_df %>%
  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                       chromosome == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome)) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Chromsome")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))

