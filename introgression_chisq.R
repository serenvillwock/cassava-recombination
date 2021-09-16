#########################################
#########################################

# This script contains code for conducting chi square tests to test association between
# M. glaziovii introgression status and crossover number
# within and outside of the introgressed regions on Chr. 1 (25-end Mb) and 4 (5-25Mb)

# data used in this analysis is available on the Cassavabase FTP repository:
# https://www.cassavabase.org/ftp/manuscripts/Chan_et_al_2019/. 

# Originally written by Ariel W. Chan, 2018
# Revised by Seren S. Villwock, 2021

#########################################
#########################################

library(tidyverse, stringr)
setwd("~/Chan_et_al_2019/")


# read in the files
dosage <- readRDS("meanDoseGlaz_250Kwindows_AllPops_31519.rds")    # dosage for M. Glaz. alleles
key <- readRDS("parents_gbsID_50219.rds")    # list of parents that I asked Marnin to look up (n=196) plus their associated GBS ID

crossovers1.list <- readRDS("chr001_maxna0.30_t0.5_informative.rds")
crossovers1 <- as.data.frame(unlist(crossovers1.list))

crossovers4.list <- readRDS("chr004_maxna0.30_t0.5_informative.rds")
crossovers4 <- as.data.frame(unlist(crossovers4.list))


# Aggregate all of the meioses genome-wide:
meioses.gw.faircount = c()
filepath = "maxna/"
for(chrom in 1:18){
  chrompad <- str_pad(chrom, 3, "left", 0)
  meios_thischrom <- readRDS(paste(filepath,"chr", chrompad, "_maxna0.30_t0.5_informative.rds", sep=""))
  meioses.gw.faircount = c(meioses.gw.faircount, names(meios_thischrom))
}

# count crossovers within and outside of the introgression regions
# for each introgression class
COcountsum = matrix(nrow=3,ncol=4,
                    dimnames=list(c("dose0","dose1","dose2"),
                                  c("chr1_inside","chr1_outside","chr4_inside","chr4_outside")))

# get mean glaziovii dosage for each respective introgression region
for(chr in c(1,4)){ 
  chrom <- chr
  chrompad <- str_pad(chr, 3, "left", 0)
  
  target <- grep(paste("^", chrom, "$", sep=""), dosage$Chr)
  subset <- do.call(rbind, dosage$doseGlazInWindows[target])
  subset <- subset[subset$IID %in% key$FullSampleName, ]    # only include individuals listed in key

  if(chrom == 1){ subset <- subset[subset$Start >= 25000000,] } # subset for Chr 1 introgression region; 25 Mb to the end 
  if(chrom == 4){ subset <- subset[(subset$Start >= 5000000 | subset$Start == 5000000) & (subset$Stop < 25000000 | subset$Stop == 25000000),] }
  # subset for Chr 4 introgression region; 5-25 Mb to the end 
  
  probands <- unique(subset$IID)
  meandose <- sapply(probands, function(proband){ mean(subset$meanDoseGlaz[grep(paste("^", proband, "$", sep=""), subset$IID)]) })

  # round to either 0, 1, or 2
  meandose[meandose > 0 & meandose < 0.5] <- 0
  meandose[meandose > 0.5 & meandose < 1.5 | meandose == 0.5 & meandose < 1.5 ] <- 1
  meandose[meandose > 1.5 & meandose < 2 | meandose == 1.5 & meandose < 2] <- 2

  # rename meandose (currently using the FullSampleName; use germplasmName instead)
  newname <- c()
  for (i in names(meandose)){
    newname <- c(newname, key$germplasmName[grep(paste("^", i, "$", sep=""), key$FullSampleName)])
  }
  names(meandose) <- newname

  meandose = as.data.frame(meandose)
  meandose$parent = rownames(meandose)

  crossovers = get(paste0("crossovers",chrom))
  colnames(crossovers)[1] <- "interval"
  crossovers$meiosis <- rownames(crossovers)
  crossovers = separate(crossovers, col=meiosis, sep="\\.", into=c("offspring","parent"))
  crossovers = separate(crossovers, col=interval, sep="_", into=c("start","stop"))
  crossovers$start <- as.numeric(crossovers$start)
  crossovers$stop <- as.numeric(crossovers$stop)
  
  # pull glaziovii dosage of parents of each introgression class:
  
  if(chrom == 1){ k=0:2 }
  if(chrom == 4){ k=0:1 }
  
  for (j in k){ #iterate for each glaziovii status
    thisdose.parents <- meandose[meandose$meandose==j,2]
    
    totalgrep <- c() #which COs intervals (from unlisted dataframe) involve parents of a given dose?
    countgrep <- c() #which meioses (from original list) involve parents of a given dose? (needed for calculating expected p for chisq)
    
    for(m in 1:length(thisdose.parents)){
      totalgrep <- c(totalgrep, grep(thisdose.parents[m], crossovers$parent)) 
      pattern <- paste(".", thisdose.parents[m], sep="")
      countgrep <- c(countgrep, grep(pattern, meioses.gw.faircount)) #pull out the meioses involving parents of each introgression status
    }
    
    COs.thisdose <- crossovers[totalgrep,] #subset for the CO intervals from parents of this dose
    assign(paste0("chr",chrompad,"_j",j,"COs"), COs.thisdose) #name that subset to a new variable
    
    #meioses.thisdose <- get(paste0("crossovers",chrom,".list"))[countgrep] #subset for the meioses from parents of this dose
    meioses.thisdose <- meioses.gw.faircount[countgrep]
    
    assign(paste0("meio_count_chr",chrom,"_dose",j), meioses.thisdose) #name that subset to a new variable
    
    
    #write out files to save:
    #filepath = "/Chan_paper_data/"
    #saveRDS(COs.thisdose, paste(filepath,"chr", chrompad, "_j", j, "crossovers.rds", sep="")) 
    
    #reload RDS files if necessary:
    #chr001_j0COs <- readRDS("./Chan_paper_data/chr001_j0crossovers.rds")
    #chr001_j1COs <- readRDS("./Chan_paper_data/chr001_j1crossovers.rds")
    #chr001_j2COs <- readRDS("./Chan_paper_data/chr001_j2crossovers.rds")
    #chr004_j0COs <- readRDS("./Chan_paper_data/chr004_j0crossovers.rds")
    #chr004_j1COs <- readRDS("./Chan_paper_data/chr004_j1crossovers.rds")
    #COs.thisdose <- get(paste0("chr",chrompad,"_j",j,"COs"))
    
    if(chrom == 1){ #introgression is 25 Mb to the end of the chromosome 1
      inside.COs = subset(COs.thisdose, COs.thisdose$stop >= 25000000)
      outside.COs = subset(COs.thisdose, COs.thisdose$stop < 25000000)
      
      COcountsum[j+1, 1] = nrow(inside.COs) #count how many COs inside the region of interest
      COcountsum[j+1, 2] = nrow(outside.COs)
    }
    if(chrom == 4){ #introgression is 5 to 25 Mb on chromosome 4
      inside.COs = subset(COs.thisdose, COs.thisdose$stop >= 5000000 & COs.thisdose$start <= 25000000)
      outside.COs = subset(COs.thisdose, COs.thisdose$stop < 5000000 | COs.thisdose$start > 25000000)
      
      COcountsum[j+1, 3] = nrow(inside.COs)
      COcountsum[j+1, 4] = nrow(outside.COs)
    }
    
    #print(c("chrom=",chrompad," insideCOcount=",nrow(inside.COs), " outsideCOcount=", nrow(outside.COs)))
  }
}


COcountsum #a dataframe listing the number of meioses listed for parents of each dose


for (chrom in c(1,4)){

#calculate expected proportions for chi.square test:
  
  # This method ('length(meio_count_chrn_dosej)') counts the total number of meioses 
  # that are listed in the CO dataset for the introgressed chromosome
  # (a meiosis event may be listed up to 18 times if a crossover was detected on each chromosome)
  # that way, if there is no crossover observed for a certain chromosome within a certain meioses, 
  # it is not counted against the total, because we assume there was 
  # at least one crossover per chromosome per meiosis even if we couldn't detect it
  # so undected obligate crossovers are not counted in the diffence between observed & expected
  
  if(chrom == 1){ 
    totald0 = length(meio_count_chr1_dose0) #Ariel's text: 2645
    totald1 = length(meio_count_chr1_dose1) #Ariel's text: 4047
    totald2 = length(meio_count_chr1_dose2) #Ariel's text: 439
    total.all = sum(totald0 + totald1 + totald2)
    p = c((totald0/total.all), (totald1/total.all), (totald2/total.all))
    
    #check_expected = p * sum(COcountsum[1,1:2])
    #check_p = table(data$introgression)/sum(table(data$introgression))
  }
  
  if(chrom == 4){ 
    totald0 = length(meio_count_chr4_dose0)
    totald1 = length(meio_count_chr4_dose1)
    total.all = sum(totald0 + totald1)
    p = c((totald0/total.all), (totald1/total.all))
      }
  if(chrom == 1) {test = 1:2}
  if(chrom == 4) {test = 3:4}
for(t in test){
  Xsq = chisq.test(na.omit(COcountsum[,t]), p = p)
  assign(paste0("chrom",chrom,"chisq",t), Xsq) }
}

###Chi-sq results: Chr1 introgression, within the introgressed region
chrom1chisq1$observed
chrom1chisq1$expected
chrom1chisq1$p.value
#percent difference
(chrom1chisq1$observed-chrom1chisq1$expected) / (chrom1chisq1$expected) *100

###Chi-sq results: Chr1 introgression, outside the introgressed region
chrom1chisq2$observed
chrom1chisq2$expected
chrom1chisq2$p.value
#percent difference
(chrom1chisq2$observed-chrom1chisq2$expected) / (chrom1chisq2$expected) *100


###Chi-sq results: Chr4 introgression, within the introgressed region
chrom4chisq3$observed
chrom4chisq3$expected
chrom4chisq3$p.value
#percent difference
(chrom4chisq3$observed-chrom4chisq3$expected) / (chrom4chisq3$expected) *100


###Chi-sq results: Chr4 introgression, outside the introgressed region
chrom4chisq4$observed
chrom4chisq4$expected
chrom4chisq4$p.value
#percent difference
(chrom4chisq4$observed-chrom4chisq4$expected) / (chrom4chisq4$expected) *100
















####  Look Genome-wide: #######

#find genome-wide number of meioses and COs for each introgression status:
filepath="./Chan_et_al_2019/"


### Define meandose, mean glaziovii dosage, for chr1 and chr4 statuses: ###
#Chr 1 introgression meandose:
target1 <- grep(paste("^", 1, "$", sep=""), dosage$Chr)
subset1 <- do.call(rbind, dosage$doseGlazInWindows[target1])
subset1 <- subset1[subset1$IID %in% key$FullSampleName, ]         # only include individuals of interest (i.e. those listed in key)
subset1 <- subset1[subset1$Start >= 25000000,]
probands1 <- unique(subset1$IID)
meandose1 <- sapply(probands1, function(proband){ mean(subset1$meanDoseGlaz[grep(paste("^", proband, "$", sep=""), subset1$IID)]) })

meandose1[meandose1 > 0 & meandose1 < 0.5] <- 0
meandose1[meandose1 > 0.5 & meandose1 < 1.5 | meandose1 == 0.5 & meandose1 < 1.5 ] <- 1
meandose1[meandose1 > 1.5 & meandose1 < 2 | meandose1 == 1.5 & meandose1 < 2] <- 2

# rename (currently using the FullSampleName; use germplasmName instead)
newname <- c()
for (i in names(meandose1)){
  newname <- c(newname, key$germplasmName[grep(paste("^", i, "$", sep=""), key$FullSampleName)])
}
names(meandose1) <- newname

meandose1 = as.data.frame(meandose1)
meandose1$parent = rownames(meandose1)


#Chr 4 introgression meandose:
target4 <- grep(paste("^", 4, "$", sep=""), dosage$Chr)
subset4 <- do.call(rbind, dosage$doseGlazInWindows[target4])
subset4 <- subset4[subset4$IID %in% key$FullSampleName, ]         # only include individuals of interest (i.e. those listed in key)
subset4 <- subset4[(subset4$Start >= 5000000 | subset4$Start == 5000000) & (subset4$Stop < 25000000 | subset4$Stop == 25000000),] 
probands4 <- unique(subset4$IID)
meandose4 <- sapply(probands4, function(proband){ mean(subset4$meanDoseGlaz[grep(paste("^", proband, "$", sep=""), subset4$IID)]) })

meandose4[meandose4 > 0 & meandose4 < 0.5] <- 0
meandose4[meandose4 > 0.5 & meandose4 < 1.5 | meandose4 == 0.5 & meandose4 < 1.5 ] <- 1
meandose4[meandose4 > 1.5 & meandose4 < 2 | meandose4 == 1.5 & meandose4 < 2] <- 2

# rename (meandose4 (currently using the FullSampleName; use germplasmName instead)
newname <- c()
for (i in names(meandose4)){
  newname <- c(newname, key$germplasmName[grep(paste("^", i, "$", sep=""), key$FullSampleName)])
}
names(meandose4) <- newname

meandose4 = as.data.frame(meandose4)
meandose4$parent = rownames(meandose4)


meioses.genome.wide = data.frame()
COs.genome.wide = data.frame()

for (chrom in 1:18){
  
  chrompad <- str_pad(chrom, 3, "left", 0)

  COsthischrom <- readRDS(paste(filepath,"chr", chrompad, "_maxna0.30_t0.5_informative.rds", sep=""))
  
  for(intrg in c(1,4)){
  
  if(intrg == 1){ k=0:2 }
  if(intrg == 4){ k=0:1 }
  
  meandose = get(paste0("meandose",intrg))
  
  for (j in k){ #iterate for each glaziovii status
    thisdose.parents <- meandose[meandose$meandose==j,2] #pull out names of parents of this introgression status
    
    countgrep <- c() #for pulling out unique meioses involving parents of this introgression status
    for(m in 1:length(thisdose.parents)){
      pattern <- paste(".", thisdose.parents[m], sep="")
      countgrep <- c(countgrep, grep(pattern, names(COsthischrom))) #pull out any meioses involving parent m
    }
    
    # count total number of meioses in individuals of this introgression status for this chrom 
    meioses.thischrom.thisdose <- COsthischrom[countgrep] #subset for all of the meioses from a parent of this dose
    meioses.genome.wide <- rbind(meioses.genome.wide, c(chrom, intrg, j, length(meioses.thischrom.thisdose)))
    COs.genome.wide <- rbind(COs.genome.wide, c(chrom, intrg, j, length(unlist(meioses.thischrom.thisdose))))
  }
  }
}

colnames(meioses.genome.wide) = c("chrom","intrg","dose","meio_count")
colnames(COs.genome.wide) = c("chrom","intrg","dose","CO_count")

meio.genome.wide.sum = aggregate(meio_count~intrg+dose, data=meioses.genome.wide, sum)
COs.genome.wide.sum = aggregate(CO_count~intrg+dose, data=COs.genome.wide, sum)


# Genome-wide chi-square:

for(intrg in c(1,4) ){
  
    x = COs.genome.wide.sum$CO_count[COs.genome.wide.sum$intrg == intrg ]
    
    totalmeioses = sum(meio.genome.wide.sum$meio_count[meio.genome.wide.sum$intrg == intrg])
    
    p = meio.genome.wide.sum$meio_count[meio.genome.wide.sum$intrg == intrg] / totalmeioses
  
    Xsq <- chisq.test(x, p = p)
    
    assign(paste0("intrg",intrg,"Xsq.gw"), Xsq)
}


intrg4Xsq.gw$observed
intrg4Xsq.gw$expected
intrg4Xsq.gw$p.value
#percent difference (significant)
(intrg4Xsq.gw$observed - intrg4Xsq.gw$expected) / ( intrg4Xsq.gw$expected) *100


intrg1Xsq.gw$observed
intrg1Xsq.gw$expected
intrg1Xsq.gw$p.value
#percent difference (non-significant)
(intrg1Xsq.gw$observed - intrg1Xsq.gw$expected) / ( intrg1Xsq.gw$expected) *100


