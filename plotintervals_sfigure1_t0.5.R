#########################################
#########################################

# This script contains code for creating supp. Fig. 1,
# mapping the crossover intervals detected for each chromosome

# data used in this analysis is available on the Cassavabase FTP repository:
# https://www.cassavabase.org/ftp/manuscripts/Chan_et_al_2019/. 

# Originally written by Ariel W. Chan, 2018
# Revised by Seren S. Villwock, 2021

#########################################
#########################################

#Read in files:

require(stringr); require(parallel); require(ggplot2); require(RColorBrewer)

directory <- "./Chan_et_al_2019"

t <- 0.5
centromere <- read.delim('centromereV6.txt', header=T)
centromere_start <- centromere[,2]*1000000
centromere_end <- centromere[,3]*1000000
thresholdna <- 0.3
pdf.fn <- paste("cointervals_familycorrected_nale",str_pad(thresholdna, 4, "right", 0), "_t", str_pad(t, 4, "right", 0), "_informative_merged.pdf", sep="")


pdf(pdf.fn)
for (chrom in 1:18){
  chrompad <- str_pad(chrom, 3, "left", 0)
  mapavg.fn <- paste("chr", chrompad, "_familycorrected_", str_pad(thresholdna, 4, "right", 0), ".mapavg", sep="")
  
  
  
  mapavg <- read.delim(paste(directory, mapavg.fn, sep="/"), header=T)
  passed <- mapavg[grep(T, mapavg$PROB_RECOMBINATION>t),]
  if(dim(passed)[1]==0){ return("no intervals") }
  pairname <- paste(passed$CHILD, passed$PARENT, sep=".")
  interval <- paste(passed$START, passed$END, sep="_")
  names(interval) <- pairname
  
  pairname <- unique(pairname)
  
  #   informative.fn <- paste("/workdir/ac2278/shapeit2/maxna/chr", chrompad, "_maxna0.30_t", t, "_informative.rds", sep="")
  #   interval <- readRDS(informative.fn)
  #   interval <- sapply(names(interval), function(pairname){
  #     names(interval[[pairname]]) <- rep(pairname, length(interval[[pairname]]))
  #     return(interval[[pairname]])
  #   })
  #   names(interval) <- NULL
  #   interval <- unlist(interval)
  
  # co intervals detected in the same parent-offspring pair share an pos value
  pos <- c(rep(NA, length(interval)))
  for (i in 1:length(pairname)){
    fm <- pairname[i]
    pos[grep(paste("^",fm,"$",sep=""), names(interval))] <- i
  }
  morethan1 <- table(pos)[table(pos)>1]  # parent-offspring pairs with more than one detected co
  names(morethan1) <- NULL
  # pos
  # 12  57  59 110 153 187 290 303 371 405 482 504 507 527 
  # 2   2   2   2   2   2   3   2   2   2   2   2   2   2 
  
  # color line based on parent of parent-offspring pair
  parents <- unique(do.call(rbind, strsplit(names(interval), "[.]"))[,2])
  col <- c(rep(NA,length(interval)))
  for(i in 1:length(parents)){
    parentname <- paste("[.]", parents[i], "$", sep="")
    col[grep(parentname, names(interval))] <- i
  }
  
  #   start <- c(as.numeric(do.call(rbind, strsplit(interval, "_"))[,1]), centromere_start[chrom]) # the last start point in the beginning of the centromere
  #   end <- c(as.numeric(do.call(rbind, strsplit(interval, "_"))[,2]), centromere_end[chrom])
  start <- c(as.numeric(do.call(rbind, strsplit(interval, "_"))[,1]))
  end <- c(as.numeric(do.call(rbind, strsplit(interval, "_"))[,2]))
  
  #   dat <- data.frame(
  #     pos = c(pos,pos[length(pos)]+1),
  #     start = start,
  #     end = end
  #   )
  
  dat <- data.frame(
    pos = pos,
    start = start,
    end = end
  )
  #     availcolors <- colors()[!1:length(colors()) %in% grep("white", colors())]
  # the centromere will always appear as red and will appear at the very bottom of the plot
  # generate n distinctive colors
  n <- length(unique(col))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  print(ggplot(dat) + 
          geom_segment(aes(x=start, y=pos, xend=end, yend=pos), color=c(col_vector[col]), size=0.5) +
          scale_y_reverse() + xlab("physical position (bp)") + ylab("parent-offspring pair") + 
          ggtitle(paste("chromosome ", chrom, sep="")) +
          #geom_rect(aes(xmin=centromere_start[chrom], xmax=centromere_end[chrom], ymin=0, ymax=Inf), alpha=0.005, fill="blue") 
          annotate("rect", xmin=centromere_start[chrom], xmax=centromere_end[chrom], ymin=0, ymax=Inf, alpha=0.1, fill="blue")
        #   return(table(morethan1))
  )
}
dev.off()