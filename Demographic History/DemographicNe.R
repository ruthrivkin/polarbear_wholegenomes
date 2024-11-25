setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/DemographicHistory/")
dropbox_dir<-file.path("~/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

getwd()


library(vcfR)
library(adegenet)
library(strataG)
library(dartR)
library(dartR.base)


# convert genpop
snpsR <- read.vcfR("CurrentNe/PB_24.10.30_ld_HWE.thin1.vcf", verbose = T)
# convert vcf to genpop
snps_genind <- vcfR2genind(snpsR)
pop.data <- read.csv("CurrentNe/IndMetadata.csv")

snps_genind@pop <- (as.factor(pop.data$pop))
snps_genind

class(snps_genind)

# convert genind to gtypes
snps_gtypes <- genind2gtypes(snps_genind)
#class(snps_gtypes)

Ne.1 <- ldNe(snps_gtypes, 
             maf.threshold=0, 
             by.strata=TRUE, 
             ci=0.95, 
             drop.missing=TRUE,
             num.cores = 2)


gl.1 <- gl.read.PLINK("CurrentNe/PB_24.10.30_ld_HWE.thin1", ind.metafile = "CurrentNe/IndMetadata.csv")
gl.1 <- gl.reassign.pop(gl.1, as.pop = "pop")
gl.1@loc.names <- gsub("\\.", "_", gl.1@loc.names)

#Subset by population
library(poppr)
gl.1@pop
gl.nw.1 <- popsub(gl.1, sublist = "NW" )

#Convert files to stratag req format
genind.NW.1 <- gl2gi(gl.nw.1)

# convert genind to gtypes
snps_gtypes.NW.1 <- genind2gtypes(genind.NW.1)
#class(snps_gtypes)

# Estimating Ne using ldNe (from https://github.com/jdalapicolla/Ne_StrataG.R/blob/master/Ne_Estimation.R)
Ne.1 <- ldNe(snps_gtypes.NW.1, 
           maf.threshold=0, 
           by.strata=TRUE, 
           ci=0.95, 
           drop.missing=TRUE,
           num.cores = 1)
Ne.1


#try NeEstimator
genepop.1 <- gl2genepop(gl.1)

ne.1 <- dartR.popgen::gl.LDNe(gl.1, outfile="cityLD.txt", 
                            outpath=getwd(),
                            neest.path ="~/Downloads/NeEstimator/", 
                            critical=c(0), singleton.rm=TRUE, mating='random',
                            Waples.correction = "nChromosomes", Waples.correction.value = 74,
                            naive = TRUE)


