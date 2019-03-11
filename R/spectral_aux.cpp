#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector insertionArgSort(NumericVector t) { 
    
    int curr, j, n = t.size();
    NumericVector o(n);
    
    // Wypełnienie tablicy porządku
    for (int i=0; i<n; i++) {
        o[i] = i;
    }
    
    for (int i=1; i<n; i++) { ; 
        curr= o[i];
        j = i-1; 
        
        while (j >= 0 && t[o[j]] > t[curr]) { 
            o[j+1] = o[j]; 
            j = j-1; 
        } 
        o[j+1] = curr; 
    } 
    
    return o;
} 

// [[Rcpp::export]]
NumericMatrix Mnn(NumericMatrix X, int M) {
    int n = X.nrow();
    int d = X.ncol();
    
    NumericMatrix dists(n, n);
    NumericMatrix S(n, M);
    NumericVector odl;
    NumericVector odleglosci(n);
    NumericVector posortowane;
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            
            // Obliczanie odległości
            odl = 0;
            for (int k=0; k<d; k++) {
                odl = odl + (X(i, k)-X(j, k))*(X(i, k)-X(j, k));
            }
            odl = sqrt(odl);
            
            // Zapisywanie odległości
            odleglosci[j] = odl[0];
        }
        
        posortowane = insertionArgSort(odleglosci);
        
        for (int l=0; l<M; l++) {
            S(i, l) = posortowane[l+1]+1; // +1, żeby indeksowanie z R się spięło
        }
    }
    
    return S;
}

// [[Rcpp::export]]
NumericMatrix stworzG(NumericMatrix S) {
    int n = S.nrow();
    int M = S.ncol();
    
    NumericMatrix G(n, n);
    
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int u=0; u<M; u++)
                if (S(i, u) == j || S(j, u) == i) {
                    G(i, j) = 1;
                    G(j, i) = 1;
                    break;
                }
        
    return G;
}
