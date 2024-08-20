setwd("/Volumes/OneTouch5GB/PolarBearWholeGenomes/JGIRR00/DemographicHistory/")
dropbox_dir<-file.path("~/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

getwd()

#devtools::install_github(repo="zakrobinson/RLDNe")

library(RLDNe)
#dartRverse::dartRverse_install()
library(dartRverse)
library(dartR.popgen)

#Load gl file that has been filtered for HWE and LD with pop level info
system("~/plink/plink --bfile PB_24.07.12_ldpruned.ne --allow-extra-chr -recode A --out PB_24.07.12_ldpruned.ne")

infile<- "PB_24.07.12_ldpruned.ne.raw"
myData <- read.PLINK(infile)
myData@loc.names
gl.save(myData, "PB_24.07.12_ldpruned.ne.gl")

myData <- gl.load("PB_24.07.12_ldpruned.ne.gl")

test <- myData[,1:100]
test@loc.names
gp.test <- gl2genepop(test)

nes <- dartR.popgen::gl.LDNe(myData, outfile="popsLD.txt", 
               outpath="DemographicHistory",
               neest.path ="~/Downloads/NeEstimator/", 
               critical=c(0,0.05), singleton.rm=TRUE, mating='random', 
               Waples.correction = "nChromosomes", Waples.correction.value = 37)

 #try snpR
#remotes::install_github("hemstrow/snpR")
library(snpR)
library(dartR)
gl <- gl.load("city.genlight.hw")
samples <- read.delim("SeedSamplesnpR.txt")

snpR <- import.snpR.data("ld.city.vcf", sample.meta = samples)
snpR@sample.meta
snpR@snp.meta

calc_ne(snpR, facets = "City", method = "LD", chr = "CHROM",
              NeEstimator_path = "/Users/ruthrivkin/Downloads/NeEstimator/Ne2-1M",
              outfile ="pop_ne.txt")
get.snpR.stats(x = snpR)
              