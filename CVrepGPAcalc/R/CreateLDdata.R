#' @title Function to generate low density SNP panels
#'
#' @description Selects SNPs to create low density SNP genotype files. SNPs can be selected folowing different methods and a number of SNP panel replicates can be generated for each density
#'
#' @param filename input file is a plink '.tped' data file
#' @param nreps the number of sampling replicates for each density
#' @param method the method of SNP selection (1 to select SNPs by random sampling across the entire genome, where the SNPs selected in the different replicates may overlap by chance; 2 to select SNPs with random sampling across the genome, but replicates are non-overlapping; 3 to select SNPs randomly within each chromosome but proportionally to the chromosome length; 4 to select SNPs based on their physical (or genetic) distance for pre-defined step sizes)
#' @param denstart the starting (lowest) SNP density to be sampled
#' @param denend the highest SNP density to be sampled
#' @param step the step by which the densities increase for method 3, or the SNP distance step size for method 4
#' @param gen_len the genome length if using method 3
#' @param chr_len_file for method 3, the input '.txt' file with chromosome, chromosome name, and chromosome length
#' @param stepstart for method 4, the lowest step size between selected SNPs
#' @param stepend for method 4, the highest step size between selected SNPs
#'
#' @return NULL
#'
#' @examples CreateLDdata(filename="genotypes.tped", nreps=3, method=1, denstart=200, denend=500, step=100)
#' @examples CreateLDdata(filename="genotypes.tped", nreps=3, method=3, denstart=200, denend=500, step=100, gen_len=2240.19, chr_len_file="Chr_len.txt")
#' @examples CreateLDdata(filename="genotypes.tped", method=4, stepstart=5000000, stepend=9000000, step=2000000)
#'
#' @export CreateLDdata

CreateLDdata <- function (filename, nreps, method=1, denstart, denend, step, gen_len, chr_len_file, stepstart, stepend) {

  options(scipen=999)
  df<- read.table(file=filename, header=F)

  #Method 1: creates LD snp genotype panels with random sampling across the genome, with overlap
  if (method==1) {  for (kk in seq(from=denstart, to=denend, by=step)) { print(noquote(paste("Density", kk)))
    for (jj in 1:nreps) { mydata <- df[sample(1:nrow(df), kk, replace=F),]

    filenameLD=paste0("LD", kk, "_rep", jj, ".tped")
    write.table(mydata, file=filenameLD, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "  ")
    }}
  }

  #Method 2: creates LD snp panels with random sampling across the genome, non-overlapping
  if (method==2) {  for (kk in seq(from=denstart, to=denend, by=step)) {  print(noquote(paste("Density", kk)))
    for (jj in 1:nreps) { mydata <- df[sample(1:nrow(df), kk, replace=F),]
    sel_snp <- mydata$V2

    filenameLD=paste0("LD", kk, "_rep", jj, ".tped")
    write.table(mydata, file=filenameLD, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "  ")

    df <- as.data.frame(df[-which(df$V2 %in% sel_snp),])
    }}
  }

  #Method 3: creates LD snp genotype panels with random sampling within each chromosome where the number of snps on each chromosome is proportional to chromosome length
  if (method==3) { chr1<-df$V1

  chr_len <- read.table(chr_len_file, header=T)
  Chr2<-as.numeric(chr_len$Chr)
  len<-as.numeric(chr_len$Size_Mb)

  for (kk in seq(from=denstart, to=denend, by=step)) { print(noquote(paste("Density", kk)))
    for (jj in 1:nreps) { chr_snp_list<-data.frame()
    snp_list<-data.frame()

    for (i in 1:length(unique(chr1))) { df_chr <- df[which(df$V1==i),]
    prop <- ceiling( kk * ( len[i] / gen_len ) )

    if (prop < nrow(df_chr)) { chr_snp_list <- df_chr[sample(1:nrow(df_chr), prop, replace=F),]
    snp_list<-rbind(snp_list, chr_snp_list)
    } else if (prop >= nrow(df_chr)) { chr_snp_list <- df_chr
    snp_list<-rbind(snp_list, chr_snp_list) }
    }
    filename=paste0("LD", kk, "_rep", jj, ".tped")
    write.table(snp_list, file=filename, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "    ")

    density <- kk
    rep <-jj
    num_snps <- nrow(snp_list)
    snp_count_df<-data.frame(density, rep, num_snps)
    colnames(snp_count_df) <- c("Density", "Replicate", "SNPs_selected")
    write.table(snp_count_df, file="snp_count.txt", append=TRUE, col.names=!file.exists("snp_count.txt"), row.names=FALSE, quote=FALSE, sep = "    ")
    }}
  }

  #Method 4: SNPs are selected to be on average equally spaced, based on their physical distance for varied step sizes
  if (method==4) { chr <- df$V1
  for (step in seq(from=stepstart, to=stepend, by=step)) { print(noquote(paste("Step", step)))
    snp_list <- list()
    for (i in 1:length(unique(chr))) { snps_on_chr <- vector()
    ichr_data<-df[which(chr==i), 1:4]

    snp <- ichr_data$V2
    Phys_pos <- ichr_data$V4
    startpos <- ichr_data$V4[1]
    endpos <- ichr_data$V4[nrow(ichr_data)]
    nsnps <- ceiling((endpos - startpos) / step )       #smallest integer not less than the value, to account for always selecting the last snp

    z <- which(Phys_pos == startpos)                    #always select the first snp
    snps_on_chr[1] <- as.character(snp[z])
    pos <- ( startpos + step )
    w <- 0

    for (j in 2:nsnps) { if ( (endpos - pos) < step) { snps_on_chr[j] <- as.character(snp[nrow(ichr_data)])     #always select the last snp, and if it's the last, stop
    break
    } else if (j > nrow(ichr_data)) { break                   #because of rounding with 'ceiling', if nsnps is > the number of snps in a chr, stop
    } else { z <- which( abs(Phys_pos-pos) == min(abs(Phys_pos-pos)) )    #find snp position which is closer to position of the starting snp + the step
    if (length(z)>1){ z <- z[2]                                   #for four instances of two snps equally spaced on either side of the step, chose the 2nd
    }
    if ( z == w ) { z <- (which( abs(Phys_pos-pos) == min(abs(Phys_pos-pos)))) +1  #if same snp selected twice, select next snp instead (always the next)
    }
    snps_on_chr[j]<- as.character(snp[z])
    pos <- ( pos + step )
    w <- z
    }
    }
    snp_list[[i]] <- snps_on_chr
    }

    Final_snp_list <- unlist (snp_list)
    Final_snp_list <- c("IID", "Chip", unique (Final_snp_list))
    df_final <- df[which(df$V2 %in% Final_snp_list),]

    total_snps<-length(Final_snp_list)                 #total number of snps selected
    snp_count_df<-data.frame(step, total_snps)
    colnames(snp_count_df) <- c("Step", "SNPs_selected")
    write.table(snp_count_df, file="snp_count.txt", append=TRUE, col.names=!file.exists("snp_count.txt"), row.names=FALSE, quote=FALSE, sep = "    ")

    ####Create outputfiles according to HD density origin, and step size
    filenameLD <- paste0("LDstep", step,".tped")
    write.table(df_final, file=filenameLD, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "    ")
  }

  ####Create outputfiles with the HD density panels only for the offspring - to be used later for CV
  #filenameHD<-paste0("HDoff.tped")
  #write.table(df, file=filenameHD, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep = "  ")
  }
}
