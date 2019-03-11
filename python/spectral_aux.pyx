import cython
import numpy as np
cimport numpy as np

@cython.boundscheck(False)


def __dist(np.ndarray[double] x, np.ndarray[double] y):
    """ Funkcja oblicza odległość między punktami x i y. """
    cdef int n = len(x)
    cdef double wynik=0.0
    
    for i in range(n):
        wynik = wynik + (x[i]-y[i])**2
        
    return wynik**0.5
    

def Mnn(np.ndarray[double, ndim=2] X, int M):
    """
    X -- macierz o rozmiarze (n, d), gdzie d to wymiar puntków, a n to ich liczba
    M -- liczba sąsiadów do wyznaczenia
    Funkcja wyznacza M najbliższych sąsiadów (względem metryki euklidesowej) każdego z punktów.
    """
    cdef int ix1=0, ix2=0, n=X.shape[0]
    assert M < n
    
    cdef np.ndarray[double, ndim=2] dists = np.zeros((n, n), dtype=np.double) # macierz odległości
    cdef np.ndarray[np.int64_t, ndim=2] S = np.zeros((n, M), dtype=np.int) # wynikowa macierz
    cdef float odl
    
    for ix1 in range(n): # mniejszy indeks
        for ix2 in range(ix1+1, n): # większy indeks
            odl = __dist(X[ix1], X[ix2])
            dists[ix1][ix2] = odl
            dists[ix2][ix1] = odl
        
        S[ix1] = np.argsort(dists[ix1])[1:M+1]
    
    return S, dists


def __dfs(np.ndarray[np.int64_t, ndim=2] G, np.ndarray[np.int64_t] vis, int v):
    """ Bazowa funkcja depth-first search, która przechodzi po kolejnych wierzchołkach. """
    cdef int n = G.shape[0]
    vis[v] = 1

    for w in range(n):
        if G[v][w] == 1:
            if vis[w] == 0: # nieodwiedzony sąsiad
                __dfs(G, vis, w) 

                
def dfs(np.ndarray[np.int64_t, ndim=2] G, int v):
    """ Funkcja, która wywołuje funkcję bazową. """
    cdef int n = G.shape[0]
    cdef np.ndarray[np.int64_t] visited = np.zeros(n, dtype=np.int)       
    
    __dfs(G, visited, v)
    
    return np.where(visited == 1)[0]


def __skladowe(np.ndarray[np.int64_t, ndim=2] G):
    """ Funkcja zwraza liczbę spójnych składowych i listy indeksów wierzchołków tych składowych. """
    to_visit = [i for i in range(G.shape[0])]
    skladowe = []
    cdef int liczba_skladowych = 0, ix
    
    while len(to_visit) > 0:
        liczba_skladowych += 1
        skladowe.append(dfs(G, to_visit[0]))
        
        for ix in skladowe[-1]:
            to_visit.remove(ix)
    
    return liczba_skladowych, skladowe


def __uspojnij(np.ndarray[np.int64_t, ndim=2] G):
    """ Funkcja uspójniająca graf G poprzez połączenie spójnych składowych krawędzią między
    wierzchołkiem o największym indeksie w jednej składowej a najmniejszym indeksem w drugiej składowej. """
    cdef int liczba_skladowych, i, j, s
    liczba_skladowych, lista_skladowych = __skladowe(G)
    
    # Uspójniamy tylko, gdy graf nie jest spójny.
    if liczba_skladowych > 1:
        
        # "Wiążemy" końce list reprezentujących spójne składowe
        for s in range(liczba_skladowych-1):
            i = lista_skladowych[s][-1]
            j = lista_skladowych[s+1][0]
            G[i][j] = 1
            G[j][i] = 1
        
    return G


def __uspojnij2(np.ndarray[np.int64_t, ndim=2] G, np.ndarray[double, ndim=2] dists):
    """ Funkcja uspójniająca graf G poprzez połączenie spójnych składowych krawędzią między najbliższymi wierzchołkami. """
    cdef int liczba_skladowych, najblizszy1, najblizszy2, v, w
    cdef double odleglosc_najblizszych
    liczba_skladowych, lista_skladowych = __skladowe(G)
    
    # Uspójniamy tylko, gdy graf nie jest spójny.
    if liczba_skladowych > 1:
        najblizszy1, najblizszy2 = -1, -1
        odleglosc_najblizszych = np.inf
        
        skladowa_glowna = []
        
        for skladowa in lista_skladowych:
            skladowa_glowna = np.concatenate((skladowa_glowna, skladowa))
            
            for v in skladowa_glowna:
                for w in skladowa:
                    if dists[v][w] < odleglosc_najblizszych and v != w:
                        najblizszy1 = v
                        najblizszy2 = w
                        odleglosc_najblizszych = dists[v][w]
            
            G[najblizszy1][najblizszy2] = 1
            G[najblizszy2][najblizszy1] = 1
            
    return G


def Mnn_graph(np.ndarray[np.int64_t, ndim=2] S, np.ndarray[double, ndim=2] dists=None):
    """
    S -- macierz o rozmiarze (n, M), gdzie n to liczba puntków, a M to liczba sąsiadów tego punktu
    Funkcja zwraca (uspójnioną, jeśli tego wymaga) macierz sąsiedztwa tych M najbliższych sąsiadów.
    """
    cdef int i, j, u, n=S.shape[0], M=S.shape[1]
    cdef np.ndarray[np.int64_t, ndim=2] G = np.zeros((n, n), dtype=np.int) 
    
    for i in range(n):
        for j in range(i+1, n):
            for u in range(M):
                if S[i][u] == j or S[j][u] == i:
                    G[i][j] = 1
                    G[j][i] = 1
                    break
    if dists is None:
        return __uspojnij(G)
    else:
        return __uspojnij2(G, dists)