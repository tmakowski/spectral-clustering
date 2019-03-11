suppressPackageStartupMessages(library('mclust'))
suppressPackageStartupMessages(library('genie'))
suppressPackageStartupMessages(library('stringi'))
source("spectral.R")
source("testy_metody.R")


wyniki <- function(zbiory, etykiety, indeks, nazwa, m=8, std=FALSE) {
    stopifnot(indeks %in% c("FM", "AR"))
    stopifnot(length(zbiory) == length(etykiety))
    
    metody <- c(
        "Spectral Clustering",
        "Genie",
        "Hierarchical - ward.D",
        "Hierarchical - ward.D2",
        "Hierarchical - single",
        "Hierarchical - complete",
        "Hierarchical - average",
        "Hierarchical - mcquitty",
        "Hierarchical - median",
        "Hierarchical - centroid",
        "Mclust"
    )
    
    wynikowa <- structure(
        data.frame(matrix(ncol=length(metody), nrow=0)),
        names=metody)
    
    for (i in 1:length(zbiory)) {
        data   <- as.matrix(read.table(zbiory[i]))
        labels <- as.integer(read.table(etykiety[i])[, 1])
        clust  <- max(labels) - min(labels) + 1
        
        if (std) {
            data <- (data-mean(data))/sd(data)
        }
        
        set_name <- as.character(
            stri_match_first_regex(
                stri_match_first_regex(zbiory[i], "[^/]*$"), "^[^\\.]*"))
        print(set_name)
        
        curr <- lapply(list(
            spectral_clustering(data, k=clust, M=m),
            cutree(hclust2(objects=data), clust),
            cutree(hclust(dist(data), method="ward.D"), clust),
            cutree(hclust(dist(data), method="ward.D2"), clust),
            cutree(hclust(dist(data), method="single"), clust),
            cutree(hclust(dist(data), method="complete"), clust),
            cutree(hclust(dist(data), method="average"), clust),
            cutree(hclust(dist(data), method="mcquitty"), clust),
            cutree(hclust(dist(data), method="median"), clust),
            cutree(hclust(dist(data), method="centroid"), clust),
            Mclust(data, G=clust)$classification
        ), function(labs_pred) indeksy(labels, labs_pred, indeks, print=FALSE))
        
        wynikowa[set_name, ] <- unlist(curr)
        
        write.table(wynikowa[i, ], nazwa, sep=",", append=TRUE, col.names=(i==1))
    }
    
    return() #return(wynikowa)
}

lista_zbiorow <- list.files(".", pattern="data.gz", recursive=TRUE)
lista_etykiet <- list.files(".", pattern="labels0.gz", recursive=TRUE)


wyniki(lista_zbiorow, lista_etykiet, "FM", nazwa=file.path("indeksy", "FM_std.csv"))
wyniki(lista_zbiorow, lista_etykiet, "AR", nazwa=file.path("indeksy", "AR_std.csv"))
wyniki(lista_zbiorow, lista_etykiet, "FM", nazwa=file.path("indeksy", "FM_std.csv", std=TRUE))
wyniki(lista_zbiorow, lista_etykiet, "AR", nazwa=file.path("indeksy", "AR_std.csv", std=TRUE))
wyniki(lista_zbiorow, lista_etykiet, "FM", nazwa=file.path("indeksy", "FM_M15.csv", m=15))
wyniki(lista_zbiorow, lista_etykiet, "AR", nazwa=file.path("indeksy", "AR_M15.csv", m=15))
wyniki(lista_zbiorow, lista_etykiet, "FM", nazwa=file.path("indeksy", "FM_M25.csv", m=25))
wyniki(lista_zbiorow, lista_etykiet, "AR", nazwa=file.path("indeksy", "AR_M25.csv", m=25))

