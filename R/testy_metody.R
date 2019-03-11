suppressPackageStartupMessages(library('dendextend'))
suppressPackageStartupMessages(library('mclust'))

rysuj_zbior <- function(data, labels) {
    plot(data[, 1], data[, 2], col=labels, axes=FALSE, xlab="", ylab="", pch=4)
    box()
}


rysuj_wzorzec <- function(zbior) {
    data    <- as.matrix(read.table(paste0(zbior, ".data.gz")))
    labels  <- as.integer(read.table(paste0(zbior, ".labels0.gz"))[,1])
    
    # Rysowanie
    rysuj_zbior(data, labels+1*(min(labels)==0)) # jeśli indeksacja jest od zera, to przesuwamy
}


indeksy <- function(labs, labs_pred, indeks=NULL, print=TRUE) {
    if (print) 
        cat(paste(
            paste("FM:", FM_index(labs, labs_pred)),
            paste("AR:", adjustedRandIndex(labs, labs_pred)),
            sep="\n\n"))
    else {
        if (indeks == "FM")
            return(FM_index(labs, labs_pred))
        else
            return(adjustedRandIndex(labs, labs_pred))
    }
        # return(structure(
        #     list(FM_index(labs, labs_pred), adjustedRandIndex(labs, labs_pred)),
        #     names=c("FM", "AR")))
}


testuj <- function(zbior, ...) {
    
    data    <- as.matrix(read.table(paste0(zbior, ".data.gz")))
    labels  <- as.integer(read.table(paste0(zbior, ".labels0.gz"))[, 1])
    clust   <- max(labels) - min(labels) + 1 # zliczenie liczby skupień
    
    # Wyliczenie
    labels_predicted <- spectral_clustering(data, k=clust, ...)
    
    # Rysowanie
    rysuj_zbior(data, labels_predicted)
    
    # Wypisanie indeksów
    indeksy(labels, labels_predicted)
}