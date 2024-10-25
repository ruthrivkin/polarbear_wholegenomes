setwd("/Volumes/OneTouch5GB 1/PolarBearWholeGenomes/JGIRR00/")
dropbox_dir<-file.path("~/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

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

#snmf
library(LEA)
library(tidyverse)


# Create geno file
system("~/plink/plink --bfile PB_24.08.16_ld_HWE --allow-extra-chr -recode vcf --out PB_24.08.16_ld_HWE")

vcfin <- "PB_24.08.16_ld_HWE.vcf"
vcf2lfmm(vcfin)

#Sample information
pop.info <- read.csv(file.path(dropbox_dir, "Datasheets/24.08.16_master.csv"))
pop.info <- pop.info[order(pop.info$BedOrder),]
head(pop.info)
pop.info$Long <- as.numeric(pop.info$Long)
pop.info$Lat <- as.numeric(pop.info$Lat)
pop.info$LabAge <- as.numeric(pop.info$LabAge)


##identify best clusters
genoin <- "PB_24.08.16_ld_HWE.geno"
geno <- read.geno("PB_24.08.16_ld_HWE.geno")

obj.snmf = snmf(genoin, K = 1:10, ploidy = 2, entropy = T,
                alpha = 100, rep = 10, project = "new", seed = 42)

obj.snmf=load.snmfProject("PB_24.05.28_ldpruned.snmfProject")

summary(obj.snmf)

#obj.snmf = load.snmfProject("Subset.ld.hwe.snmfProject")

K <- summary(obj.snmf)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  dplyr::select(-temp)

# choose omptimal (K = 3)
ggplot(K, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 4, color = "black", fill = "blue") +
  scale_x_continuous(breaks = seq(0, 16, by = 1)) +
  labs(x = "Number of Clusters", y = "Cross-entropy criterion") +
  ng
ggsave(file.path(dropbox_dir, "Figures/Cross-entropy.pdf"), width = 6.78, height = 5.3, dpi=300)


#plot ancestry matrix
ce = cross.entropy(obj.snmf, K = 3)
ce
lowest.ce = which.min(ce)
lowest.ce

qmatrix = as.data.frame(Q(obj.snmf, K = 3, run = lowest.ce))
head(qmatrix)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

cols = my.colors(4)

pdf(file.path(dropbox_dir,"Figures/IndividualsAncestry.pdf"),  width = 6.78, height = 5.3)
barchart(obj.snmf, K = 3, run = lowest.ce,
         border = NA, space = 0,
         col = cols,
         xlab = "Individuals",
         ylab = "Ancestry proportions K = 5",
         main = "Ancestry matrix") -> bp
#create axis
axis(1, at = 1:length(bp$order), 
     labels = bp$subpop, las = 3, 
     cex.axis = .4)
dev.off()

# Label column names of qmatrix
ncol(qmatrix)
cluster_names = c()
for (i in 1:ncol(qmatrix)){
  cluster_names[i] = paste("Cluster", i)
}
cluster_names
colnames(qmatrix) = cluster_names
head(qmatrix)

# Add individual IDs
qmatrix$Ind = pop.info$SampleName

# Add site IDs
qmatrix$Site = pop.info$Subpop
head(qmatrix)

# Convert dataframe to long format

qlong = reshape2::melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

# Change order of sites by using the factor function
#site.order = c("BB","DS", "FB","GB","KB","LS","MC","NB","NW","SB","SH", "VM","WH")
qlong$Site= as.factor(qlong$Site)
levels(qlong$Site)

# Define colour palette
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

cols = my.colors(length(unique(qlong$variable)))


# Plot admixture barplot 
(admix.bar = ggplot(data=qlong, aes(x=Ind, y=value, fill=variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Site, scales = "free", ncol = 2)+
  scale_fill_manual(values = cols)+
  ylab("Admixture proportion")+
  # xlab("Individual")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(colour="black", size=12),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
)


# ----------------- #
#
# Prepare pie charts
#
# ----------------- #

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = reshape2::melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "FB", cols = cols)

# Apply function to all sites using for loop
site.ordered = sort(c("BB","DS", "FB","GB","LS","MC","NW", "NB","SB","SH"))

pies = list()
for (i in site.ordered){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}


# ----------------- #
#
# Prepare basemap
#
# ----------------- #

# Import csv file containing coordinates

(pop.coords <- pop.info  %>% 
  group_by(Subpop)  %>% 
  drop_na(Long, Lat) %>% 
  dplyr::summarize(N = n(),
         Mean.Lon = mean(Long),
         Mean.Lat = mean(Lat)
    )
)


# Order alphabetically by site]
coords = pop.coords[order(pop.coords$Subpop), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$Subpop)

# Set map boundary (xmin, xmax, ymin, ymax)
library(sf)
library(raster)
# Get map outlines from rworldmap package
library(rworldmap)
library(rworldxtra)
library(ggsn)
library(rgeos)
library(maps)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)

boundary = extent(-140,-55, 50, 90) 
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
(basemap = ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="lightgrey",
               colour="black", size=0.5)+
  coord_quickmap(expand=F)+
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 500, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", 
                 location = "bottomleft", anchor = c(x = -120, y = 50.4),
                 st.bottom = FALSE, st.size = 3, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  ng
)


# ----------------- #
#
# Add pie charts to basemap
#
# ----------------- #

# Extract coordinates for each site
coord.list = list()


for (i in site.ordered){
  coord.list[[i]] = c(subset(coords, Subpop == i)$Mean.Lon, subset(coords, Subpop == i)$Mean.Lat)
}

coord.list

# Define pie chart sizes
radius = 2.5

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(site.ordered)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
(pie.map = basemap + pies.ac)


ggsave(file.path(dropbox_dir, "Figures/snmf_pie_charts_map.pdf"), width = 6.78, height = 5.3, dpi = 300)



#Assign Individual cluster based on 75% ancestry proportion
names(qmatrix) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Ind", "Site")
cluster.assign <- qmatrix %>% 
  mutate(
    Assigned = case_when(
      Cluster1 > .75 ~ "Cluster 1",
      Cluster2 > .75 ~ "Cluster 2",
      Cluster3 > .75 ~ "Cluster 3",
      Cluster4 > .75 ~ "Cluster 4",
      .default = "Admixed"
    )
  )

write.csv(cluster.assign, file.path(dropbox_dir, "Datasheets/24.05.31_snmfcluster.csv"), row.names=FALSE)

clusternames <- dplyr::select(cluster.assign, Assigned, Ind)
colnames(clusternames) <- c("Cluster", "SampleName")

samples <- merge(pop.info, clusternames)
write.csv(samples, file.path(dropbox_dir, "Datasheets/24.05.31_master.csv"), row.names=FALSE)
