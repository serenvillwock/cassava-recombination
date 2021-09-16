#########################################
#########################################

# This script contains code for Fig. 2 showing
# a comparison of the constructed genetic maps
# (broken down by introgression status)
# with the ICGMC genetic map

# data used in this analysis is available on the Cassavabase FTP repository:
# https://www.cassavabase.org/ftp/manuscripts/Chan_et_al_2019/. 

# Originally written by Ariel W. Chan, 2018
# Revised by Seren S. Villwock, 2021

#########################################
#########################################

library(tidyverse)
setwd("./Chan_et_al_2019")
includelegend <- "no"


# run first with chrom=1, then change to chrom=4
# to look at both chromosomes with introgression regions
chrom <- 4
backgroundchrom <- chrom
chrompad <- str_pad(chrom, 3, "left", 0)


dosage <- readRDS("meanDoseGlaz_250Kwindows_AllPops_31519.rds")    # dosage for M. Glaz. allele
parents <- readRDS("parents.rds")            # list of parents 
key <- readRDS("parents_gbsID_50219.rds")    # list of parents plus their associated GBS IDs

target <- grep(paste("^", chrom, "$", sep=""), dosage$Chr)
subset <- do.call(rbind, dosage$doseGlazInWindows[target])
subset <- subset[subset$IID %in% key$FullSampleName, ]         # only include individuals of interest (i.e. those listed in key)



# Wolfe et. al 2019 found long segments of M. glaziovii haplotypes in modern cassava 
# germplasm on chromosomes 1 and 4. The largest introgressions were detected on 
# chromosome 1, spanning from 25 Mb to the end of the chromosome and on chromosome 4 
# from 5Mb to 25Mb.

#subset is within the introgression region
if(chrom == 1){ subset <- subset[subset$Start > 25000000 | subset$Start == 25000000,] }
if(chrom == 4){ subset <- subset[(subset$Start > 5000000 | subset$Start == 5000000) & (subset$Stop < 25000000 | subset$Stop == 25000000),] }


probands <- unique(subset$IID)
#proband <- probands[1]
meandose <- sapply(probands, function(proband){ mean(subset$meanDoseGlaz[grep(paste("^", proband, "$", sep=""), subset$IID)]) })


##histogram of glaziovii allele dosages:
hist(meandose, breaks=25, xlim=c(0,2), 
     main=paste0("histogram of glaziovii allele \n dosage for chr ", chrom))
abline(v=0.7,lty=2)
abline(v=1.5,lty=2)

# round meandose to either 0, 1, or 2
meandose[meandose > 0 & meandose < 0.7] <- 0
meandose[meandose > 0.7 & meandose < 1.5 | meandose == 0.7 & meandose < 1.5 ] <- 1
length(meandose[meandose > 0.7 & meandose < 1.5 | meandose == 0.7 & meandose < 1.5 ])
meandose[meandose > 1.5 & meandose < 2 | meandose == 1.5 & meandose < 2] <- 2
length(meandose[meandose > 1.5 & meandose < 2 | meandose == 1.5 & meandose < 2])

# rename meandose (currently using the FullSampleName; use germplasmName instead)
newname <- c()
for (i in names(meandose)){
  newname <- c(newname, key$germplasmName[grep(paste("^", i, "$", sep=""), key$FullSampleName)])
}
names(meandose) <- newname


### Try to replicate the chi-square tests:
crossovers1 <- as.data.frame(unlist(
  readRDS("~/Desktop/Research/Chan_paper_revisions/Chan_paper_code_data/maxna/all_chroms/chr001_maxna0.30_t0.5_informative.rds")))
crossovers4 <- as.data.frame(unlist(
  readRDS("~/Desktop/Research/Chan_paper_revisions/Chan_paper_code_data/maxna/all_chroms/chr004_maxna0.30_t0.5_informative.rds")))

meandose = as.data.frame(meandose)
meandose$parent = rownames(meandose)

#get glaz dose information:
crossovers <- readRDS(paste("chr", chrompad, "_maxna0.30_t", t, "_informative.rds", sep=""))


# Ariel's interval method:

nperinterval = function(occupied, wall){       # left_right; left and right must be an entry of wall
  feet <- unlist(strsplit(occupied, "_"))
  mode(feet) <- "numeric"
  leftfoot <- feet[1]
  rightfoot <- feet[2]
  distance <- rightfoot-leftfoot                                             # how far did you travel?
  tempwall <- wall       # define rooms
  roomname <- c()
  while(length(tempwall)>1){
    roomname <- c(roomname, paste(tempwall[1], tempwall[2], sep="_"))
    tempwall <- tempwall[-1]
  }
  roomcount <- setNames(vector("numeric",length(roomname)), roomname)  # empty vector for tallying
  
  leftwall <- paste("^", leftfoot, "_", sep="")
  leftwall <- grep(leftwall, names(roomcount))
  
  rightwall <- paste("_", rightfoot, "$", sep="")
  rightwall <- grep(rightwall, names(roomcount))
  
  rooms <- names(roomcount[leftwall:rightwall])
  sizes <- sapply(strsplit(rooms, "_"), function(entry){
    a <- as.numeric(entry)
    return(a[2]-a[1])
  })
  roomcount[leftwall:rightwall] <- (sizes/sum(sizes))
  
  return(roomcount)
} 



t <- 0.5
maxna <- 0.3
require(ggplot2); require(stringr); require(parallel); require(plyr)


endpos <- c()
for (chr in 1:18){
  chrompad <- str_pad(chr, 3, "left", 0)
  fppos <- read.table(paste("/allchroms/chr", chrompad, "_familycorrected_", str_pad(maxna, 4, "right", 0), "_shapeit2.snp.mm", sep=""), header=T)
  fppos <- fppos[,c("id","position","missing")]
  endpos <- c(endpos, max(fppos$position))
}


ymax <- c()
for (i in 1:18){
  chrompad <- str_pad(i, 3, "left", 0)
  all <- readRDS(paste("chr", chrompad, "_introgressionco_JLscale-correctedm.rds", sep=""))
  ymax <- c(ymax, max(do.call(rbind, all)[,2]))
}



# look at each chromosome, separately with chr001 status (i.e. backgroundchrom <- 1)
for (chrom in setdiff(1:18, c(1,4))){

  #################
  
  chrompad <- str_pad(chrom, 3, "left", 0)
  fppos <- read.table(paste("/local/workdir/ac2278/shapeit2/maxna/chr", chrompad, "_familycorrected_", str_pad(maxna, 4, "right", 0), "_shapeit2.snp.mm", sep=""), header=T)
  fppos <- fppos[,c("id","position","missing")]
  pos <- fppos$position 
  crossovers <- readRDS(paste("chr", chrompad, "_maxna0.30_t", t, "_informative.rds", sep=""))
  
  # copies0 <- names(meandose[meandose==0])
  # copies0 <- sapply(copies0, function(proband){
  #   pattern <- paste(".", proband, "$", sep="")
  #   crossovers[grep(pattern, names(crossovers))]
  # })
  # copies0 <- unlist(copies0)
  # ncoperinterval0 <- mclapply(copies0, function(interval){
  #   nperinterval(occupied=interval, wall=pos)
  # }, mc.cores=50)
  # 
  # ncoperinterval0 <- do.call(rbind, ncoperinterval0)
  # ncoperinterval0 <- apply(ncoperinterval0, 2, sum)
  
  jessen <- read.delim("/workdir/ac2278/shapeit2/cassavaV6_0.chromosome.ICGMCv2.smoothspline.srt.norm.map", header=F)
  names(jessen) <- c("chromosome", "snpid", "genetic", "physical")
  subset <- jessen[grep(paste("Chromosome", str_pad(chrom, 2, "left", 0), sep=""), jessen$chromosome),]
  genetic <- subset$genetic
  physical <- subset$physical
  
  mono <- c(0, genetic) - c(genetic, 0)    # check to see that genetic positions are monotonically increasing
  mono <- mono[-1]
  mono <- mono[-length(mono)]
  start <- grep(T, mono>0)
  
  while(length(start)>0){
    from <- start[1]                         # physical position of the start of issue region 
    value <- genetic[from]                   # genetic value at this start position
    to <- grep(T, genetic>value)[1] - 1
    if(is.na(to)){ to <- nrow(subset) }
    genetic[from:to] <- rep(value, length(from:to))
    
    mono <- c(0, genetic) - c(genetic, 0)    # check to see that genetic positions are monotonically increasing
    mono <- mono[-1]
    mono <- mono[-length(mono)]
    start <- grep(T, mono>0)
  }
  
  jessenlength <- genetic[length(genetic)]
  
  #########################################
  
  
  # table(meandose) -- chr004
  # meandose
  # 0   1 
  # 153  43 
  
  # table(meandose) -- chr001
  # meandose
  # 0   1   2 
  # 67 114  15 
  
  
  if(length(table(meandose))==2){ k <- 0:1 } else k <- 0:2
  # m <- c()
  # for (j in k){
  #   copies <- names(meandose[meandose==j])
  #   copies <- sapply(copies, function(proband){
  #     pattern <- paste(".", proband, "$", sep="")
  #     meiosiscount <- length(grep(pattern, names(crossovers)))  # changed from number of parents to number of informative meioses
  #   })
  #   m <- c(m,unlist(sum(copies)))     
  # }
  # names(m) <- k
  # 
  # # for chr001: 1760 2514  252 (sum: 4526)
  if(backgroundchrom==1){ m <- setNames(c(2645, 4047, 439), k) }
  if(backgroundchrom==4){ m <- setNames(c(5318, 1813),k) }
  all <- list()
  
  for (j in k){
    copies <- names(meandose[meandose==j])
    copies <- sapply(copies, function(proband){
      pattern <- paste(".", proband, "$", sep="")
      crossovers[grep(pattern, names(crossovers))]
    })
    if(length(copies)==0){ print(paste("no ", j, " for chr", chrompad, sep="")) }
    
    m_i <- m[as.character(j)]          # changed from number of parents length(copies) to number of informative meioses
    copies <- unlist(copies)
    ncoperinterval <- mclapply(copies, function(interval){
      nperinterval(occupied=interval, wall=pos)
    }, mc.cores=10)
    if(length(copies)>0){ saveRDS(ncoperinterval, paste("chr", chrompad, "_maxna", maxna, "_t", t, "_j", j, "_ncoperinterval-corrected.rds", sep="")) }
    
    # ncoperinterval <- readRDS(paste("chr", chrompad, "_maxna", maxna, "_t", t, "_j", j, "_ncoperinterval-corrected.rds", sep=""))
    ncoperinterval <- do.call(rbind, ncoperinterval)
    if(length(copies)>0){ ncoperinterval <- apply(ncoperinterval, 2, sum) }
    
    # scale <- jessenlength/sum(ncoperinterval)
    scales <- readRDS("geneticmapscalingfactor.rds")
    scale <- scales[chrom]
    scaled <- ncoperinterval*scale
    print(scale)
    
    ariel <- cbind(pos, 0)
    mode(ariel) <- "numeric"
    for (i in 2:nrow(ariel)){
      previous <- pos[i-1]
      current <- pos[i]
      room <- paste(previous, "_", current, sep="")
      ariel[i,2] <- scaled[room]+ariel[(i-1),2]
    }
    ariel[,2] <- ariel[,2]*(sum(m)/m_i)    # JL's suggested scaling
    all[[as.character(j)]] <- ariel
    
  }


# saveRDS(all, paste("chr", chrompad, "_introgressionco_commonscale.rds", sep=""))
saveRDS(all, paste("chr", chrompad, "_introgressionco_JLscale-correctedm.rds", sep=""))
# all <- readRDS(paste("chr", chrompad, "_introgressionco.rds", sep=""))
# all <- readRDS(paste("chr", chrompad, "_introgressionco_commonscale.rds", sep=""))
# all <- readRDS(paste("chr", chrompad, "_introgressionco_JLscale.rds", sep=""))

AWC <- readRDS(paste("chr", chrompad, "_AWC.rds", sep=""))
all[['AWC']] <- AWC

if(chrom==4){ k <- 0:1 } else k <- 0:2

if(length(k)==2){
  x <- c(rep(all[[1]][,1], 3), subset$physical)/1000000                     # physical position 
  y <- c(all[[1]][,2], all[[2]][,2], all[[3]][,2],subset$genetic)
  category <- rep(c('0 copies', '1 copy', 'AWC', 'ICGMC'), c(rep(nrow(all[[1]]),3), nrow(subset)))
  data <- data.frame(x=x, y=y, category=category)
}
if(length(k)==3){
  x <- c(rep(all[[1]][,1], 4), subset$physical)/1000000                     # physical position 
  y <- c(all[[1]][,2], all[[2]][,2], all[[3]][,2], all[[4]][,2], subset$genetic)
  category <- rep(c('0 introgressions', '1 introgression', '2 introgressions', 'AWC', 'ICGMC'), c(rep(nrow(all[[1]]),4), nrow(subset)))
  data <- data.frame(x=x, y=y, category=category)
}
# scale <- jessenlength/sum(ncoperinterval0)
# scaled <- ncoperinterval0*scale
# 
# ariel <- cbind(pos, 0)
# mode(ariel) <- "numeric"
# for (i in 2:nrow(ariel)){
#   previous <- pos[i-1]
#   current <- pos[i]
#   room <- paste(previous, "_", current, sep="")
#   ariel[i,2] <- scaled[room]+ariel[(i-1),2]
# }

centromere <- read.delim('/workdir/ac2278/shapeit2/centromereV6.txt', header=T)
centromere_start <- centromere[,2]*1000000
centromere_end <- centromere[,3]*1000000
chrom_end <- centromere[,4]

# x <- c(ariel[,1], subset$physical)/1000000
# y <- c(ariel[,2], subset$genetic)
# category <- rep(c('AWC', 'ICGMC'), c(nrow(ariel), nrow(subset)))
# data <- data.frame(x=x, y=y, category=category)
# # geneticend <- c(geneticend, max(y))


# plot
if(length(k)==3){ colors <- c('orange', "blue", "darkgreen", "turquoise", "red2") }
if(length(k)==2){ colors <- c('orange', "blue", "turquoise", "red2") }
title <- paste("chromosome ", chrom, "; t=", t, sep="")
# if(includelegend=="no"){ pdf(paste("chr", chrompad, "_scaledgeneticmap-introgression-wAWC.pdf", sep="")) } else
#   pdf(paste("chr", chrompad, "_scaledgeneticmap-introgression-wAWC-wlegend.pdf", sep=""))

pdf(paste("chr", chrompad, "_wchr001introgressionstatus-correctedm_nolegend.pdf", sep=""))
# pdf(paste("chr", chrompad, "_scaledgeneticmap-introgression-wAWC-wlegend-wcommonscale.pdf", sep=""))
if(chrom==1){ introbegin <- 25 ; introend <- endpos[1]/1000000 }
if(chrom==4){ introbegin <- 5 ; introend <- 25 }
if(!chrom==1 & !chrom==4){ introbegin <- 0; introend <- 0 }
if(includelegend=="yes"){ print(ggplot(data, aes(x=x, y=y)) + geom_line(aes(colour=category)) + 
                                  scale_color_manual(values=colors) +
                                  ylab("genetic position (cM)") + 
                                  ggtitle(title) +
                                  xlab("physical position (Mb)") +
                                  theme(axis.text=element_text(size=20), 
                                        axis.title=element_text(size=18),
                                        plot.title = element_text(size=18)) + 
                                  xlim(0,max(endpos)/1000000) + 
                                  # ylim(0,max(jessen$genetic)+1) +
                                  ylim(0,max(data$y, na.rm=T)) +
                                  annotate("rect", xmin=centromere_start[chrom]/1000000, xmax=centromere_end[chrom]/1000000, ymin=0, ymax=Inf, alpha=0.1, fill="blue") +
                                  annotate("rect", xmin=introbegin, xmax=introend, ymin=0, ymax=Inf, alpha=0.1, fill="red"))
  
} else
  print(ggplot(data, aes(x=x, y=y)) + geom_line(aes(colour=category)) + 
          scale_color_manual(values=colors) +
          ylab("genetic position (cM)") + 
          ggtitle(title) +
          xlab("physical position (Mb)") +
          theme(legend.position="none", 
                axis.text=element_text(size=20), 
                axis.title=element_text(size=18),
                plot.title = element_text(size=18)) + 
          xlim(0,max(endpos)/1000000) + 
          # ylim(0,max(jessen$genetic)+1) +
          ylim(0, max(ymax)) +
          annotate("rect", xmin=centromere_start[chrom]/1000000, xmax=centromere_end[chrom]/1000000, ymin=0, ymax=Inf, alpha=0.1, fill="blue") +
          annotate("rect", xmin=introbegin, xmax=introend, ymin=0, ymax=Inf, alpha=0.1, fill="red"))



dev.off()
}

