# Algorytm spektralny analizy skupień 

W tym repozytorium zawarte są rozwiązania pracy domowych z kursów Programowania i analizy danych w R oraz Przetwarzania i analizy danych w języku Python prowadzonych przez prof. Marka Gągolewskiego na wydziale MiNI Politechniki Warszawskiej. Zadanie polegało na implementacji algorytmu spektralnego analizy skupień (ang. *spectral clustering*) oraz metod składających się na ten algorytm. Dodatkowym wymogiem było skorzystanie z, odpowiednio, Rcpp i Cythona przy realizacji pracy domowej.

## Zaimplementowane metody

 * `Mnn(X, M)` to funkcja, która dla danej macierzy X wielowymiarowych punktów wyznacza macierz zawierającą informację o *M* najbliższych sąsiadach każdego z punktów. **(Rcpp, Cython)**

 * `Mnn_graph(S)`, gdzie *S* jest macierzą wynikową powyższej metody. Funkcja ta generuje
symetryczną macierz *G*, gdzie *G*[i][j] = 1, jeśli *i*-ty punkt jest jednym z *M* najbliższych sąsiadów punktu *j*-tego lub vice versa. Dodatkowo graf reprezentowany przez macierz sąsiedztwa *G* jest uspójniany i taka zakutalizowana macierz *G* jest zwracana. **(Cython)** 

 * `Laplacian_eigen(G, k)`, gdzie *k* > 1, a *G* jest wynikiem powyższej metody. Jest to funkcja, która wyznacza laplasian *L* grafu *G* i zwraca macierz składającą się z wektorów własnych macierzy *L*, które odpowiadają od 2 do *k*+1 co do wielkości wartości własnej.

 * `spectral_clustering(X, k, M)` to funkcja wykorzystująca powyższe metody (oraz wbudowane implementacje algorytmu *k*-średnich), aby dokonać podziału danego zbioru *X* na *k* różnych skupień.

## Wyniki

Obie implementacje zawierają zwięzłe testy poprawności metod oraz raporty, które stanowią o skuteczności metody `spectral_clustering` w porównaniu do innych algorytmów analizy skupień.
