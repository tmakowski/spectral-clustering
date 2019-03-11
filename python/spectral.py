from warnings import simplefilter
simplefilter("ignore", FutureWarning)

import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})

from sklearn.cluster import KMeans
from spectral_aux import Mnn, Mnn_graph


def Laplacian_eigen(G, k):
    """
    G -- (uspójniona) macierz sąsiedztwa M najbliższych sąsiadów,
    k -- liczba wektorów własnych
    Funkcja zwraca macierz rozmiaru (n, k), gdzie kolumnami są
    wektory własne odpowiadające 2., 3., ..., (k+1). co do wielkości wartościom własnym.
    """    
    
    D = np.zeros(G.shape, dtype=np.int)
    
    # Uzupełnianie przekątnej stopniami wierzchołków
    for i in range(D.shape[0]):
        D[i][i] = G[i].sum()
    
    eigenvalues, eigenvectors  = np.linalg.eigh(D-G)
    
    # transponujemy wartości własne, bo są kolumnowo
    # sortujemy wartości własne i odwracamy ich kolejność
    # wybieramy elementy 2, 3, ..., k+1 według porządku wyznaczonego przez te wartości własne
    # transponujemy wynik ponownie, ponieważ jest on wierszowy
    return eigenvectors.T[1:(k+1)].T


def spectral_clustering(X, k, M):
    """
    Funkcja spektralnej analizy skupień. 
    X -- zbiór danych,
    k -- liczba skupień,
    M -- parametr sterujący liczbą najliższych sąsiadów
    """
    assert k >= 2, M > 0
    
    S, _ = Mnn(X, M)
    G    = Mnn_graph(S)
    E    = Laplacian_eigen(G, k)
    km   = KMeans(n_clusters=k, random_state=58).fit_predict(E)
    
    return km


def spectral_clustering2(X, k, M):
    """ Funkcja wykorzystuje drugą wersję metody uspójniającej. """
    assert k >= 2, M > 0
    
    S, dists = Mnn(X, M)
    G        = Mnn_graph(S, dists)
    E        = Laplacian_eigen(G, k)
    km       = KMeans(n_clusters=k, random_state=58).fit_predict(E)
    
    return km