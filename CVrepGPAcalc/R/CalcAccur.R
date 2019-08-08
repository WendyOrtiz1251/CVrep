#' @title Function to calculate genomic prediction accuracy
#'
#' @description It calculates correlations between the trait value and the predicted breeding values, using the output files from ASReml analysis. It calculates the genomic prediction accuracies (see Options). Creates ASReml status report file.
#'
#' @param ngroups the number of cross validation groups
#' @param nreps the number of cross-validation replicates
#' @param skipnum the number of lines to skip when reading the ASReml “.sln” files
#' @param genomichh If F (default), it uses the heritability from each cross-validation fold to calculate the accuracy; If T, the user also needs to provide the genomic heritability (“h2”).
#' @param h2 the genomic heritability of the trait
#'
#' @return NULL
#'
#' @examples calcAccur(ngroups=5, nreps=10, skipnum=4, genomichh=F)
#'
#' @export calcAccur

calcAccur <- function(ngroups, nreps, skipnum, genomichh=F, h2 ) {

  for (kk in 1:nreps){

    Rdatafile <- paste0("PartI_nrep", kk, ".RData")
    load(Rdatafile)
    ordered <- data2[order(data2[,1]),]
    colnames(ordered)[1]<-"rec_id"
    colnames(ordered)[2]<-"trait_val"

    C <-vector()
    H <-vector()
    R <- vector()

    for (Grp in 1:ngroups) { filename = paste0("cvgroup", Grp, "_CV", kk, ".sln")
    G<-read.table(filename, header=FALSE, skip=skipnum)
    group_dat <- subset(G, select=c(2,3))
    names(group_dat) <- c("rec_id","Pred")
    GG<-merge(ordered, group_dat, by="rec_id")
    df<-subset(GG, Group==Grp)
    correl<-cor(df$trait_val, df$Pred, method = c("pearson"))

    if (genomichh==F) { filenamehh = paste0("cvgroup", Grp, "_CV", kk, ".pvc")           #Accuracies by dividing with sqrt of heritability wihin CrossValidation fold
    text <- readLines(filenamehh, encoding = "UTF-8")
    h2 <- as.numeric(strsplit(text[grep("hh", text)], "   ")[[1]][11])
    r<- correl/sqrt(h2) } else if (genomichh==T) { r<-correl/sqrt(h2) }  #Accuracies by dividing with the sqrt of one genomic heritability
    C<-c(C, correl)
    H<-c(H, h2)
    R<-c(R, r)

    filenameasr = paste0("cvgroup",Grp, "_CV", kk, ".asr")
    text2 <- readLines(filenameasr, encoding = "UTF-8")
    status <- text2[grep("Finished", text2)]
    df2 <- data.frame(kk, Grp, status)
    names(df2) <-c("Rep", "Group", "Status")
    reportFile <- paste0("asreml_status_report.txt")
    write.table(df2, file=reportFile, append=TRUE, col.names=!file.exists(reportFile), row.names=FALSE, quote=FALSE, sep = "    ")
    }

    ###########Mean correlation across folds
    meancor<-mean(C)

    ###Mean h2 across folds
    meanhh<-mean(H)

    ###Mean accuracy across folds
    meanr<-mean(R)

    ###
    table<-cbind(C, H, R)
    means<-c(meancor, meanhh, meanr)
    final<-rbind(table, means)
    rownames(final)<-c("group1", "group2", "group3", "group4", "group5", "means")
    colnames(final)<-c("Cor", "h2", "Accur")
    meansdf<-data.frame(kk, t(means))
    colnames(meansdf)<-c("Rep", "Cor", "h2", "Accur")

    write.table(final, file="CVResults.txt", append=TRUE, col.names=TRUE, row.names=TRUE, quote=FALSE, sep = "    ")
    write.table(meansdf, file="rr.txt", append=TRUE, col.names=!file.exists("rr.txt"), row.names=FALSE, quote=FALSE, sep = "    ")

  }
}
