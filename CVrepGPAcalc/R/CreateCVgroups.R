#' @title Function to create cross-validation groups
#'
#' @description Randomises individuals and partitions into cross-validation groups. It repeats the process for a number of times, where each time different individuals are allocated into different cross-validation groups. It creates either one file with all groups, or separate files for each training set.
#'
#' @param datafile input phenotype file, the first column is the recoded id so that it starts from 1 and increases to the total number of individuals (rows) in the datafile; the second column is the phenotype trait value
#' @param ngroups number of cross-validation groups (folds), typically 5
#' @param nreps number of cross-validation replicates
#' @param multiPhenoFiles If F (default), a single phenotype file is produced for each replicate where the last column indicates the allocated cross-validation group number. If T multiple phenotype files are created, one for each cross-validation fold, containing the training sets.
#' @param PhenoHeader if T (default) output phenotypic files have headers
#'
#' @return NULL
#'
#' @examples createCVgroups("Phenotypes_for_CV.dat", nreps=10, multiPhenoFiles=F, PhenoHeader=T)
#'
#' @export createCVgroups

createCVgroups <- function(datafile, ngroups=5, nreps, multiPhenoFiles=F, PhenoHeader=T) {
  options(scipen=999)
  for (kk in 1:nreps){

    data<-read.table(datafile, header=TRUE)
    n <- nrow(data)
    rand <- runif(n)
    A<-cbind(data, rand)
    A <- A[order(rand),]

    fold <- floor(n/ngroups)
    extra <- n%% fold
    i<- 0
    Group <- rep(NA,n)

    for(igrp in 1:ngroups) {
      istart <- i+1
      iend   <- i+fold
      if(igrp <= extra) iend <- iend+1
      Group[istart:iend] <- igrp
      i <- iend
    }
    data2 <- cbind(A,Group)
    Rdatafile <- paste0("PartI_nrep", kk, ".RData")
    save(data2, file=Rdatafile)

    if (multiPhenoFiles==F) {
      if (PhenoHeader==T){ datafile2 <- paste0("Phenotypes_CV", kk, ".txt")
      write.table(data2, file=datafile2, row.names=FALSE, col.names=TRUE, quote=FALSE) } else if (PhenoHeader==F) {
        datafile2 <- paste0("Phenotypes_CV", kk, ".txt")
        write.table(data2, file=datafile2, row.names=FALSE, col.names= FALSE, quote=FALSE) }
    } else if (multiPhenoFiles==T) { if (PhenoHeader==T) { for (jj in 1:ngroups) { CV <- subset(data2, Group!=jj)
    datafileGroup <- paste0("Group", jj, "_CV", kk, ".txt")
    write.table(CV, file=datafileGroup, row.names=FALSE, quote=FALSE)
    }
    }
      if (PhenoHeader==F) { for (jj in 1:ngroups) { CV <- subset(data2, Group!=jj)
      datafileGroup <- paste0("Group", jj, "_CV", kk, ".txt")
      write.table(CV, file=datafileGroup, row.names=FALSE, col.names=FALSE, quote=FALSE)
      }
      }
    }
  }} 
