
#########################################
#########################################

# This script contains code for conducting chi square tests to test
# for a difference in male and female crossover occurence
# and for generating Figure 3

# data used in this analysis is available on the Cassavabase FTP repository:
# https://www.cassavabase.org/ftp/manuscripts/Chan_et_al_2019/. 

# Originally written by Ariel W. Chan, 2018
# Revised by Seren S. Villwock, 2021

#########################################
#########################################


setwd("./Chan_et_al_2019/")

# do genome-wide test
t <- 0.5
require(ggplot2); require(stringr); require(parallel); require(plyr)
maxna <- 0.3
# calculate the total number of intervals (n)

nco_allmale <- c()
nco_allfemale <- c()
for (chrom in 1:18){
  
  chrompad <- str_pad(chrom, 3, "left", 0)
  
  male <- readRDS(paste("CObysex_data/chr", chrompad, "_maxna0.30_t", t, "_informative_male.rds", sep=""))
  female <- readRDS(paste("CObysex_data/chr", chrompad, "_maxna0.30_t", t, "_informative_female.rds", sep=""))
  
  selfed <- names(unlist(male, recursive=F))[names(unlist(male, recursive=F)) %in% names(unlist(female, recursive=F))]
  
  selfed <- gsub("TMEB419.", "", selfed)
  copy <- male[['TMEB419']]
  copy <- copy[names(copy) %in% selfed]
  
  male <- male[!names(male) %in% 'TMEB419']
  male <- c(male, 'TMEB419'=list(copy))
  
  female <- female[!names(female) %in% 'TMEB419']
  female <- c(female, 'TMEB419'=list(copy))
  
  
  male <- unlist(male)
  female <- unlist(female)
  
  nco_allmale <- c(nco_allmale, length(male))
  nco_allfemale <- c(nco_allfemale, length(female))
  
}
##################
nfemale <- sum(nco_allfemale) # total number of female COs
nmale <-sum(nco_allmale)

#chromosome-wide
for (i in 1:18){
  nmale_i <- nco_allmale[i]
  nfemale_i <- nco_allfemale[i]
  x <- c(nmale_i, nfemale_i)
  totaln <-  sum(x)
  
  malem <- 3486 
  femalem <- 3679 
  totalm <- malem + femalem
  pmale <-malem/totalm
  pfemale <- 1 - pmale
  Xsq <- chisq.test(x, p = c(pmale, pfemale))
  print(c(Xsq$observed, Xsq$expected, Xsq$p.value))
}

#genome-wide
x <- c(nmale, nfemale)
totaln <-  sum(x)

malem <- 3486 
femalem <- 3679 
totalm <- malem + femalem
pmale <-malem/totalm
pfemale <- 1 - pmale
Xsq <- chisq.test(x, p = c(pmale, pfemale))
Xsq$p.value




### BY INTERVAL
t <- 0.5
require(ggplot2); require(stringr); require(parallel); require(plyr)
maxna <- 0.3
# calculate the total number of intervals that we test (n)
pvalues_all <- list()

for (chrom in 1:18){
  chrompad <- str_pad(chrom, 3, "left", 0)
  windowsize <- 1000000
  #fppos <- read.table(paste("/local/workdir/ac2278/shapeit2/maxna/chr", chrompad, "_familycorrected_", str_pad(maxna, 4, "right", 0), "_shapeit2.snp.mm", sep=""), header=T)
  #fppos <- fppos[,c("id","position","missing")]
  #pos <- fppos$position
  
  # function to determine the number of co falling within artificial interval of equal size
  #ncointerval = function(occupied, propertylength, roomsize){
    # feet <- unlist(strsplit(occupied, "_"))
    # mode(feet) <- "numeric"
    # leftfoot <- feet[1]
    # rightfoot <- feet[2]
    # distance <- rightfoot-leftfoot                                             # how far did you travel?
    # 
    # tempwall <- wall <- seq(0, propertylength, roomsize)                    # define rooms
    # roomname <- c()
    # while(length(tempwall)>1){
    #   roomname <- c(roomname, paste(tempwall[1], tempwall[2], sep="_"))
    #   tempwall <- tempwall[-1]
    # }
    # roomcount <- setNames(vector("numeric",length(roomname)), roomname)  # empty vector for tallying
    # 
    # 
    # # does one or both foot/feet fall exactly ON the wall of a room? if yes, add/subtract a small perturbabtion to leftfoot/rightfoot
    # if(length(grep(T, leftfoot==wall))>0){
    #   leftfoot <- leftfoot + 0.0000001              # if start doesn't actually fall _within_ the boundaries of a room, add a small perturbation to start
    # }
    # if(length(grep(T, rightfoot==wall))>0){
    #   rightfoot <- rightfoot - 0.0000001                  # if end doesn't actually fall _within_ the boundaries of a room, subtract a small perturbation to end
    # }
    
    # do both feet fall within ONE room?
    # Which room does leftfoot fall within?
    # (1) what wall is to the left of leftfoot?
    # leftwall.lf <- wall[grep(T, leftfoot>wall)[length(grep(T, leftfoot>wall))]]     
    # # (2) what wall is to the right of leftfoot?
    # rightwall.lf <- wall[grep(T, leftfoot<wall)[1]]
    # # what room does the leftfoot fall within?
    # room.lf <- paste(leftwall.lf, rightwall.lf, sep="_")
    # 
    # # Which room does rightfoot fall within?
    # # what wall is to the left of rightfoot?
    # leftwall.rf <- wall[grep(T, rightfoot>wall)[length(grep(T, rightfoot>wall))]]
    # # what wall is to the right of rightfoot?
    # rightwall.rf <- wall[grep(T, rightfoot<wall)[1]]
    # # what room does the rightfoot fall within?
    # room.rf <- paste(leftwall.rf, rightwall.rf, sep="_")
    # 
    # # if both feet fall within ONE room, +1 to that room
    # if(room.lf==room.rf){ roomcount[room.lf] <- 1}
    # # if feet span multiple rooms, add [0-1] to rooms based on distance spanned
    # if(!room.lf==room.rf){
    #   roomcount[room.lf] <- roomcount[room.lf] + (rightwall.lf-leftfoot)/distance
    #   roomcount[room.rf] <- roomcount[room.rf] + (rightfoot-leftwall.rf)/distance
    #   sudoleft <- rightwall.lf
    #   sudoright <- leftwall.rf
    #   
    #   # what are the names of the rooms that exist between the two feet
    #   # (not including the rooms that the feet currently reside in)?
    #   sudoroomname <- c()
    #   sudotempwall <- seq(sudoleft, sudoright, roomsize)
    #   sudoroomname <- c()
    #   while(length(sudotempwall)>1){
    #     sudoroomname <- c(sudoroomname, paste(sudotempwall[1], sudotempwall[2], sep="_"))
    #     sudotempwall <- sudotempwall[-1]
    #   }
    #   roomcount[sudoroomname] <- rep(roomsize/distance, length(sudoroomname))
    # }
    # return(roomcount)
  #} 
  
  centromere <- read.delim('centromereV6.txt', header=T)
  centromere_start <- centromere[,2]*1000000
  centromere_end <- centromere[,3]*1000000
  chrom_end <- centromere[,4]
  
  male <- readRDS(paste("CObysex_data/chr", chrompad, "_maxna0.30_t", t, "_informative_male.rds", sep=""))
  female <- readRDS(paste("CObysex_data/chr", chrompad, "_maxna0.30_t", t, "_informative_female.rds", sep=""))
  

  both <- readRDS(paste("maxna/all_chroms/chr", chrompad, "_maxna0.30_t", t, "_informative.rds", sep=""))

  selfed <- names(unlist(male, recursive=F))[names(unlist(male, recursive=F)) %in% names(unlist(female, recursive=F))]
  
  selfed <- gsub("TMEB419.", "", selfed)
  copy <- male[['TMEB419']]
  copy <- copy[names(copy) %in% selfed]
  
  male <- male[!names(male) %in% 'TMEB419']
  male <- c(male, 'TMEB419'=list(copy))
  
  female <- female[!names(female) %in% 'TMEB419']
  female <- c(female, 'TMEB419'=list(copy))
  
  
  male <- unlist(male)
  female <- unlist(female)
  both <- unlist(both)
  
  # we use propertylength=round_any(chrom_end[chrom], windowsize) to avoid the scenario where some
  # interval entries produce different numbers of intervals (i.e., one more interval at the end of chromosome)
  # This happens because the chromosome is not perfectly divisable by windowsize, and some intervals 
  # may have a stop position at the very end of the chromosome (i.e., chromosome length is 34959721; 34959721/windowsize = 34.95972;
  # if we didn't use the round_any function, ncointerval would only account for 34 windows rather than 35 and an interval
  # could fall in window 35.
  maletallybyinterval <- mclapply(male, function(interval){
    ncointerval(occupied=interval, propertylength=round_any(chrom_end[chrom], windowsize, ceiling), roomsize=windowsize)
  }, mc.cores=10)
  maletallybyinterval <- do.call(rbind, maletallybyinterval)
  truemaletallybyinterval <- apply(maletallybyinterval, 2, sum)
  # maletallybyinterval <- truemaletallybyinterval*(7165/3486)      # scaled
  
  femaletallybyinterval <- mclapply(female, function(interval){
    ncointerval(occupied=interval, propertylength=round_any(chrom_end[chrom], windowsize, ceiling), roomsize=windowsize)
  }, mc.cores=10)
  femaletallybyinterval <- do.call(rbind, femaletallybyinterval)
  truefemaletallybyinterval <- apply(femaletallybyinterval, 2, sum)
  # femaletallybyinterval <- truefemaletallybyinterval*(7165/3679)   # scaled

  intervalend <- seq(windowsize, round_any(chrom_end[chrom], windowsize, ceiling), windowsize)
  
  ###########################
  #calculate expected number of crossovers
  male <- truemaletallybyinterval[-length(intervalend)]
  female <- truefemaletallybyinterval[-length(intervalend)]
  
  pvalues <- c()
  expectedmale <- c()
  expectedfemale <- c()
  for (i in 1:length(male)){
    x <- c(male[i], female[i])
    if (length(x[x==0])==2){ pvalues <- c(pvalues, NA); expectedmale <- c(expectedmale, NA); expectedfemale <- c(expectedfemale, NA) }   
    if (length(x[x==0])==0){ 
      totaln <-  male[i] + female[i]
      malem <- 3486 
      femalem <- 3679 
      totalm <- malem + femalem
      pmale <-malem/totalm
      pfemale <- 1 - pmale
      Xsq <- chisq.test(x, p = c(pmale, pfemale))
      pvalues <- c(pvalues, Xsq$p.value)
      expectedmale <- c(expectedmale, Xsq$expected[1])
      expectedfemale <- c(expectedfemale, Xsq$expected[2])
    }
  }
  pvalues_all[[paste("chr", chrompad, sep="")]] <- pvalues
  
  expectedmale <- c(expectedmale, NA)
  expectedfemale <- c(expectedfemale, NA)
  
  numberofintervals <- length(truemaletallybyinterval)
  x <- rep(intervalend/1000000, 4)     # include the last interval
  y <- c(truemaletallybyinterval, truefemaletallybyinterval, expectedmale, expectedfemale)
  category <- c(rep("observed male", numberofintervals), rep("observed female", numberofintervals),
                rep("expected male", numberofintervals), rep("expected female", numberofintervals))
  data <- data.frame(x=x, y=y, category=category)
  

  


  # title <- paste("chromosome ", chrom, "; t=", t, "; maxna=0.3", sep="")
  # 
  # end <- 3.5e+07/1000000
  # 
  # pdf('codensitybysex-expected.pdf', width=20, height=4)
  # ggplot(data, aes(x=x, y=y, group=category)) + 
  #   geom_line(aes(linetype=category, color=category)) + 
  #   geom_point(aes(colour=category)) +
  #   scale_color_manual(values=c('red', "blue", "red", "blue")) +
  #   scale_linetype_manual(values=c("dotted", "dotted", 'solid', "solid")) +
  #   ggtitle(title) +
  #   xlab("physical position (Mb)") +
  #   ylab("number of crossovers") + 
  #   theme(axis.text=element_text(size=18),
  #         axis.title=element_text(size=18),
  #         plot.title = element_text(size=18)) +
  #   annotate("rect", xmin=centromere_start[chrom]/1000000, xmax=centromere_end[chrom]/1000000, ymin=0, ymax=Inf, alpha=0.1, fill="blue")
  # 
 

  

print(paste("completed chr", chrompad, sep=""))
}
  
beta <- sum(sapply(pvalues_all, function(chrom){
  length(chrom)-length(grep(T, is.na(chrom)))
}))

significance <- sapply(pvalues_all, function(chrom){
  chrom < 0.05/beta
})

annotation <- sapply(significance, function(chrom){
  changed <- gsub(T, "*", chrom)
  changed <- gsub(F, "", changed)
  changed[is.na(changed)] <- "-"
  return(changed)
})

saveRDS(annotation, "chisquarewindowsex-new_annotation.rds")
length(grep(T, unlist(significance)))  # 45 of the 506 intervals tested were significance at the 0.05/beta threshold


########################
pdf('codensitybysex-expected-diffannoy.pdf', width=20, height=4)
for (chrom in 1:18){
  chrompad <- str_pad(chrom, 3, "left", 0)
  windowsize <- 1000000
  fppos <- read.table(paste("/local/workdir/ac2278/shapeit2/maxna/chr", chrompad, "_familycorrected_", str_pad(maxna, 4, "right", 0), "_shapeit2.snp.mm", sep=""), header=T)
  fppos <- fppos[,c("id","position","missing")]
  pos <- fppos$position
  
  # function to determine the number of co falling within artificial interval of equal size
  ncointerval = function(occupied, propertylength, roomsize){
    feet <- unlist(strsplit(occupied, "_"))
    mode(feet) <- "numeric"
    leftfoot <- feet[1]
    rightfoot <- feet[2]
    distance <- rightfoot-leftfoot                                             # how far did you travel?
    
    tempwall <- wall <- seq(0, propertylength, roomsize)                    # define rooms
    roomname <- c()
    while(length(tempwall)>1){
      roomname <- c(roomname, paste(tempwall[1], tempwall[2], sep="_"))
      tempwall <- tempwall[-1]
    }
    roomcount <- setNames(vector("numeric",length(roomname)), roomname)  # empty vector for tallying
    
    
    # does one or both foot/feet fall exactly ON the wall of a room? if yes, add/subtract a small perturbabtion to leftfoot/rightfoot
    if(length(grep(T, leftfoot==wall))>0){
      leftfoot <- leftfoot + 0.0000001              # if start doesn't actually fall _within_ the boundaries of a room, add a small perturbation to start
    }
    if(length(grep(T, rightfoot==wall))>0){
      rightfoot <- rightfoot - 0.0000001                  # if end doesn't actually fall _within_ the boundaries of a room, subtract a small perturbation to end
    }
    
    # do both feet fall within ONE room?
    # Which room does leftfoot fall within?
    # (1) what wall is to the left of leftfoot?
    leftwall.lf <- wall[grep(T, leftfoot>wall)[length(grep(T, leftfoot>wall))]]     
    # (2) what wall is to the right of leftfoot?
    rightwall.lf <- wall[grep(T, leftfoot<wall)[1]]
    # what room does the leftfoot fall within?
    room.lf <- paste(leftwall.lf, rightwall.lf, sep="_")
    
    # Which room does rightfoot fall within?
    # what wall is to the left of rightfoot?
    leftwall.rf <- wall[grep(T, rightfoot>wall)[length(grep(T, rightfoot>wall))]]
    # what wall is to the right of rightfoot?
    rightwall.rf <- wall[grep(T, rightfoot<wall)[1]]
    # what room does the rightfoot fall within?
    room.rf <- paste(leftwall.rf, rightwall.rf, sep="_")
    
    # if both feet fall within ONE room, +1 to that room
    if(room.lf==room.rf){ roomcount[room.lf] <- 1}
    # if feet span multiple rooms, add [0-1] to rooms based on distance spanned
    if(!room.lf==room.rf){
      roomcount[room.lf] <- roomcount[room.lf] + (rightwall.lf-leftfoot)/distance
      roomcount[room.rf] <- roomcount[room.rf] + (rightfoot-leftwall.rf)/distance
      sudoleft <- rightwall.lf
      sudoright <- leftwall.rf
      
      # what are the names of the rooms that exist between the two feet
      # (not including the rooms that the feet currently reside in)?
      sudoroomname <- c()
      sudotempwall <- seq(sudoleft, sudoright, roomsize)
      sudoroomname <- c()
      while(length(sudotempwall)>1){
        sudoroomname <- c(sudoroomname, paste(sudotempwall[1], sudotempwall[2], sep="_"))
        sudotempwall <- sudotempwall[-1]
      }
      roomcount[sudoroomname] <- rep(roomsize/distance, length(sudoroomname))
    }
    return(roomcount)
  } 
  
  centromere <- read.delim('centromereV6.txt', header=T)
  centromere_start <- centromere[,2]*1000000
  centromere_end <- centromere[,3]*1000000
  chrom_end <- centromere[,4]
  
  male <- readRDS(paste("chr", chrompad, "_maxna0.30_t", t, "_informative_male.rds", sep=""))
  female <- readRDS(paste("chr", chrompad, "_maxna0.30_t", t, "_informative_female.rds", sep=""))
  both <- readRDS(paste("chr", chrompad, "_maxna0.30_t", t, "_informative.rds", sep=""))
  
  selfed <- names(unlist(male, recursive=F))[names(unlist(male, recursive=F)) %in% names(unlist(female, recursive=F))]
  
  selfed <- gsub("TMEB419.", "", selfed)
  copy <- male[['TMEB419']]
  copy <- copy[names(copy) %in% selfed]
  
  male <- male[!names(male) %in% 'TMEB419']
  male <- c(male, 'TMEB419'=list(copy))
  
  female <- female[!names(female) %in% 'TMEB419']
  female <- c(female, 'TMEB419'=list(copy))
  
  
  male <- unlist(male)
  female <- unlist(female)
  both <- unlist(both)
  
  maletallybyinterval <- mclapply(male, function(interval){
    ncointerval(occupied=interval, propertylength=round_any(chrom_end[chrom], windowsize, ceiling), roomsize=windowsize)
  }, mc.cores=10)
  maletallybyinterval <- do.call(rbind, maletallybyinterval)
  truemaletallybyinterval <- apply(maletallybyinterval, 2, sum)
  # maletallybyinterval <- truemaletallybyinterval*(7165/3486)      # scaled
  
  femaletallybyinterval <- mclapply(female, function(interval){
    ncointerval(occupied=interval, propertylength=round_any(chrom_end[chrom], windowsize, ceiling), roomsize=windowsize)
  }, mc.cores=10)
  femaletallybyinterval <- do.call(rbind, femaletallybyinterval)
  truefemaletallybyinterval <- apply(femaletallybyinterval, 2, sum)
  # femaletallybyinterval <- truefemaletallybyinterval*(7165/3679)   # scaled
  
  # use tallybyinterval to set location for annotations
  tallybyinterval <- mclapply(both, function(interval){
    ncointerval(occupied=interval, propertylength=round_any(chrom_end[chrom], windowsize, ceiling), roomsize=windowsize)
  }, mc.cores=10)
  tallybyinterval <- do.call(rbind, tallybyinterval)
  tallybyinterval <- apply(tallybyinterval, 2, sum)
  
  intervalend <- seq(windowsize, round_any(chrom_end[chrom], windowsize, ceiling), windowsize)
  
  ###########################
  #calculate expected number of crossovers
  male <- truemaletallybyinterval[-length(intervalend)]
  female <- truefemaletallybyinterval[-length(intervalend)]
  
  pvalues <- c()
  expectedmale <- c()
  expectedfemale <- c()
  for (i in 1:length(male)){
    x <- c(male[i], female[i])
    if (length(x[x==0])==2){ pvalues <- c(pvalues, NA); expectedmale <- c(expectedmale, NA); expectedfemale <- c(expectedfemale, NA) }   
    if (length(x[x==0])==0){ 
      totaln <-  male[i] + female[i]
      malem <- 3486 
      femalem <- 3679 
      totalm <- malem + femalem
      pmale <-malem/totalm
      pfemale <- 1 - pmale
      Xsq <- chisq.test(x, p = c(pmale, pfemale))
      pvalues <- c(pvalues, Xsq$p.value)
      expectedmale <- c(expectedmale, Xsq$expected[1])
      expectedfemale <- c(expectedfemale, Xsq$expected[2])
    }
  }
  pvalues_all[[paste("chr", chrompad, sep="")]] <- pvalues
  
  expectedmale <- c(expectedmale, NA)
  expectedfemale <- c(expectedfemale, NA)
  
  numberofintervals <- length(truemaletallybyinterval)
  x <- rep(intervalend/1000000, 4)     # include the last interval
  y <- c(truemaletallybyinterval, truefemaletallybyinterval, expectedmale, expectedfemale)
  category <- c(rep("observed male", numberofintervals), rep("observed female", numberofintervals),
                rep("expected male", numberofintervals), rep("expected female", numberofintervals))
  data <- data.frame(x=x, y=y, category=category)
  
  

  title <- paste("chromosome ", chrom, "; t=", t, "; maxna=0.3", sep="")

  end <- 3.5e+07/1000000
  if(chrom==1){ introbegin <- 25 ; introend <- 35 }
  if(chrom==4){ introbegin <- 5 ; introend <- 25 }
  if(!chrom==1 & !chrom==4){ introbegin <- 0; introend <- 0 }
  
  print(ggplot(data, aes(x=x, y=y, group=category)) +
          geom_line(aes(linetype=category, color=category)) +
          geom_point(aes(colour=category)) +
          scale_color_manual(values=c('red', "blue", "red", "blue")) +
          scale_linetype_manual(values=c("dotted", "dotted", 'solid', "solid")) +
          ggtitle(title) +
          xlab("physical position (Mb)") +
          ylim(0,750) +
          ylab("number of crossovers") +
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=18),
                plot.title = element_text(size=18)) +
          annotate("rect", xmin=centromere_start[chrom]/1000000, xmax=centromere_end[chrom]/1000000, ymin=0, ymax=Inf, alpha=0.1, fill="blue") +
        annotate("rect", xmin=introbegin, xmax=introend, ymin=0, ymax=Inf, alpha=0.1, fill="red") + 
          annotate('text', size=15, label=c(annotation[[paste("chr", chrompad, sep="")]], "-"), x=intervalend/1000000, y=truefemaletallybyinterval+20) 
        
        )
  print(paste("completed chr", chrompad, sep=""))
}
dev.off()

