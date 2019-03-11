suppressPackageStartupMessages(library('igraph'))
Rcpp::sourceCpp("spectral_aux.cpp")


uspojnij <- function(G) {
    visited <- rep(FALSE, nrow(G))  # lista odwiedzonych wierzchołków
    liczba_skladowych <- 0          # liczba składowych
    skladowe <- NULL
    
    iG <- graph_from_adjacency_matrix(G)
    
    while (!all(visited)) {
        liczba_skladowych <- liczba_skladowych + 1    # zaznaczamy, że weszliśmy do nowej składowej
        skladowe <- c(
            skladowe,
            list(subcomponent(iG, which(!visited)[1]))
        ) # zapisujemy bieżącą składową
        
        visited[
            unlist(
                skladowe[length(skladowe)])] <- TRUE # zapisujemy, że odwiedziliśmy te wierzchołki
    }
    
    # Nie ma co uspójniać
    if (liczba_skladowych == 1) return(G)

    mins  <- unlist(lapply(skladowe, min)[-1])
    maxes <- lapply(skladowe, max)
    maxes <- unlist(maxes[-length(maxes)])

    G[mins, maxes] <- 1
    G[maxes, mins] <- 1
    
    return(G)
}


Mnn_graph <- function(S) {
    G <- stworzG(S)
    # jeśli liczba składowych wynosi 1, to po prostu zwrócimy graf G,
    # a jeśli nie, to go uspójnimy
    return(uspojnij(G)) 
}


Laplacian_eigen <- function(G, k) {
    D <- diag(apply(G, 1, sum), nrow(G), ncol(G))
    
    ev <- eigen(D-G)
    
    # Permutacja sortująca wartości własne
    o <- insertionArgSort(ev$values)+1 # +1, bo indeksowanie od 1
    
    # Zwrócenie posortowanych kolumnowo wektorów własnych
    # zwracane jest k wektorów odpowiadające k najmniejszym wartościom własnym (z wykluczeniem najmniejszej)
    return(ev$vectors[, o][, 2:(k+1)])
}

spectral_clustering <- function(X, k, M) {
    stopifnot(k>=2)
    stopifnot(M<nrow(X))
    
    S <- Mnn(as.matrix(X), M)
    G <- Mnn_graph(S)
    E <- Laplacian_eigen(G, k)
    km <- kmeans(x=E, centers=k, iter.max=300, nstart=10)
    
    return(km$cluster)
}