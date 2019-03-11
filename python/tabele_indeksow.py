# Algorytmy do testowania
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.cluster import SpectralClustering, MiniBatchKMeans, KMeans
from spectral import spectral_clustering
from genieclust.genie import Genie

# Analiza wyników
from sklearn.metrics import fowlkes_mallows_score, adjusted_mutual_info_score, adjusted_rand_score
from glob import glob
import re

import numpy as np
import pandas as pd


# Lista zbiorów
sciezki_etykiet = sorted(glob("pd4-zbiory-benchmarkowe/*/*labels0*") + glob("zbiory/*labels0*"))
sciezki_zbiorow = sorted(glob("pd4-zbiory-benchmarkowe/*/*data*") + glob("zbiory/*data*"))
metody = [
    'Spectral Clustering',
    'Linkage single',
    'Linkage complete',
    'Linkage average',
    'Linkage weighted',
    'Linkage centroid',
    'Linkage median',
    'Linkage ward',
    'Genie',
    'Spectral Clustering - sklearn',
    'MiniBatchKMeans',
    'KMeans'
]


def wyniki(sciezka_zbioru, sciezka_etykiet, index, std, M):
    """ Funkcja zwraca nazwę zbioru i wartości indeksów na danych zbiorach. """
    def __indeksy(d, labels, index):
        """ Funkcja zwraca wybrany indeks dla każdego elementu """
        if index == "FM":
            index_method = fowlkes_mallows_score
        elif index == "AM":
            index_method = adjusted_mutual_info_score
        elif index == "AR":
            index_method = adjusted_rand_score
        else:
            return
        
        return index, dict([(name, "" if labs is None else index_method(labels, labs)) for (name, labs) in d.items()])
        
        
    data   = np.loadtxt(sciezka_zbioru, ndmin=2)
    labels = np.loadtxt(sciezka_etykiet, dtype=np.int)
    
    # Standaryzacja zmiennych
    if std:
        tmp = data.T
        tmp = (tmp-tmp.mean())/tmp.std()
        data = tmp.T
    
    # Nazwa zbioru
    r1 = re.compile("[^/]*$")
    r2 = re.compile("[^\.]*")
    set_name = r2.search(r1.search(sciezka_etykiet).group()).group()
    print("\r%-15s" % set_name, end="")
    
    # Liczba skupień
    n_clust = int(labels.max()-labels.min()+1)
    
    return (set_name,) + __indeksy({
        "Spectral Clustering":           None if set_name == "s4" else spectral_clustering(data, n_clust, M), # na tym zbiorze algorytm się wywalał
        "Linkage single":                cut_tree(linkage(data, method='single'), n_clusters=n_clust).flatten(),
        "Linkage complete":              cut_tree(linkage(data, method='complete'), n_clusters=n_clust).flatten(),
        "Linkage average":               cut_tree(linkage(data, method='average'), n_clusters=n_clust).flatten(),
        "Linkage weighted":              cut_tree(linkage(data, method='weighted'), n_clusters=n_clust).flatten(),
        "Linkage centroid":              cut_tree(linkage(data, method='centroid'), n_clusters=n_clust).flatten(),
        "Linkage median":                cut_tree(linkage(data, method='median'), n_clusters=n_clust).flatten(),
        "Linkage ward":                  cut_tree(linkage(data, method='ward'), n_clusters=n_clust).flatten(),
        "Genie":                         Genie(n_clusters=n_clust).fit_predict(data),
        "Spectral Clustering - sklearn": SpectralClustering(n_clusters=n_clust, n_neighbors=M).fit_predict(data),
        "MiniBatchKMeans":               MiniBatchKMeans(n_clusters=n_clust).fit_predict(data),
        "KMeans":                        KMeans(n_clusters=n_clust).fit_predict(data)
    }, labels, index)


def wypisz_wyniki(nazwa, indeks, lista_zbiorow=sciezki_zbiorow, lista_etykiet=sciezki_etykiet, lista_metod=metody, std=False, M=8):
    """ Funkcja wypisuje wyniki działania wszystkich algorytmów na wszystkich zbiorach w kolejnych iteracjach. """
    with np.errstate(invalid='ignore'):
        df = pd.DataFrame(index=lista_metod)
        for (zbior, etykiety) in zip(lista_zbiorow, lista_etykiet):
            w = wyniki(zbior, etykiety, index=indeks, std=std, M=M)
            df[w[0]] = w[2].values()
    
        df.to_csv(nazwa)
    
    
# wypisz_wyniki("FM.csv", "FM")
# wypisz_wyniki("AM.csv", "AM")
# wypisz_wyniki("AR.csv", "AR")
# wypisz_wyniki("FM_std.csv", "FM", std=True)
# wypisz_wyniki("AM_std.csv", "AM", std=True)
# wypisz_wyniki("AR_std.csv", "AR", std=True)
