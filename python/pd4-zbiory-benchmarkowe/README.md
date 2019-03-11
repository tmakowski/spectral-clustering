The Grand Benchmark Suite for Clustering Algorithms [abridged version]
======================================================================

Maintained by [Marek Gagolewski](http://www.gagolewski.com)



The aim of this project is to **aggregate, polish, and standardize the existing
clustering benchmark suites** referred to across the machine learning
and data mining literature. Moreover, it introduces **new datasets**
of different dimensionalities, sizes, and cluster types.




Ground-truth/reference label vectors are provided alongside each dataset.
These have been created by experts.


Each label vector imposes the number of clusters, `k`, an algorithm
shall be detecting.

Cluster similarity measures (such as the adjusted Rand's or mutual information score)
should be used to compare the true partitions with the obtained ones,
see
[scikit-learn](https://scikit-learn.org/stable/modules/clustering.html#clustering-evaluation).





# Problem Suites


The following are provided **solely for research purposes**,
unless stated otherwise. Please cite the literature references mentioned
in the corresponding dataset description files in any publications
that make use of these.


### From External Sources

Note that there is some inherent overlap between the original databases.
We have tried to resolve any conflicts in the *best* possible manner.


1. fcps -
    the Fundamental Clustering Problem Suite proposed by A. Ultsch (2005)
    from the Marburg University, Germany

    Each dataset consists of 212-4096 observations in 2-3 dimensions.

    Source: https://www.uni-marburg.de/fb12/arbeitsgruppen/datenbionik/data


2. graves -
    *synthetic data sets* considered in the paper (Graves and Pedrycz, 2010)

    Each dataset consists of 200-1050 observations in 2 dimensions.


3. other includes:

    * `iris`, `iris5` - the (? - see Bezdek J.C. et al., 1999 for discussion)
        famous Iris dataset and its imbalanced version considered
        in (Gagolewski et al., 2016).

    as well as some datasets of unknown/unconfirmed origin (@HELP WANTED@)


4. sipu -
    datasets available at the SIPU (Speech and Image Processing Unit,
    School of Computing, University of Eastern Finland) website

    Many datasets were proposed by Fränti et al., see
    (Fränti, Sieranoja, 2018). However, some datasets gathered from other
    sources (see the referenced catalog for citations) but available
    for download via the SIPU website are also included.

    Source: https://cs.joensuu.fi/sipu/datasets/

    Note that the G2 sets are available as a separate package.
    Birch3 is not included as no ground-truth labels were provided.
    We excluded the `DIM-sets` as they turn out to be too easy
    for most algorithms.


### New Datasets

5. wut -
    authored by the fantastic students
    of Marek's [Python for Data Analysis course](http://www.gagolewski.com/teaching/padpy/) @
    [Warsaw University of Technology](https://ww4.mini.pw.edu.pl/):
    Przemysław Kosewski, Jędrzej Krauze, Eliza Kaczorek, and Anna Gierlak.






# How to Load and Use


### File Format Specification


> *Consistency matters. All the datasets are automatically tested if
they are in the format specified below.*


For each `dataset`, we have the following corresponding files:

* `dataset.txt` - provides the dataset description, gives copyright info, source, etc.

* `dataset.data.gz` - defines the `n`*`d` data matrix:

    * a gzipped text file storing data in tabular format (most software packages
    can decompress `.gz` inputs on the fly)
    * columns are whitespace-delimited
    * there are exactly `n` file lines (no column names, no headers, no comments)
    * possible values are in decimal or scientific notation (e.g, 1.0, 1.23e-8)

* `dataset.labels0.gz` - the ground truth label vector

    * a gzipped text file with exactly `n` integers, one per each line
    * the `i`-th label (line) corresponds to the `i`-th data point
    * class labels are consecutive integers: `1`, `2`, ..., `k`,
    where `k` is the total number of "meaningful" clusters




### Python

```python
import numpy as np
dataset = "..." # e.g., "h2mg/h2mg_1024_30"
data    = np.loadtxt(dataset+".data.gz", ndmin=2)
labels  = np.loadtxt(dataset+".labels0.gz", dtype=np.int) # or labels1, etc.
# recall that 0 denotes the noise class, 1 - 1st cluster, 2 - 2nd one, etc.
```


Note that some clustering similarity measures
are available in [scikit-learn](https://scikit-learn.org/stable/modules/clustering.html#clustering-evaluation),
e.g.,
`sklearn.metrics.adjusted_rand_score`,
`sklearn.metrics.adjusted_mutual_info_score`,
`sklearn.metrics.fowlkes_mallows_score`.


### R

```R
dataset <- "..." # e.g., "h2mg/h2mg_1024_30"
data    <- read.table(paste0(dataset, ".data.gz"))
labels  <- as.integer(read.table(paste0(dataset, ".labels0.gz"))[,1])  # or labels1, etc.
# recall that 0 denotes the noise class, 1 - 1st cluster, 2 - 2nd one, etc.
```

For some clustering similarity scores, see
`mclust::adjustedRandIndex` and `dendextend::FM_index`.



# Bibliography

Bezdek J.C. et al. (1999). Will the real iris data please
stand up?, *IEEE Transactions on Fuzzy Systems* **7**, pp. 368-369.
[doi:10.1109/91.771092](http://dx.doi.org/10.1109/91.771092)

Campello R., Moulavi D., Zimek A., Sander J. (2015).
Hierarchical Density Estimates for Data Clustering, Visualization, and Outlier Detection,
*ACM Transactions on Knowledge Discovery from Data* **10**, art. no. 5.
[doi:10.1145/2733381](http://dx.doi.org/10.1145/2733381)

Dasgupta S., Ng V. (2009). *Single Data, Multiple Clusterings*, In:
Proc. NIPS Workshop *Clustering: Science or Art? Towards Principled Approaches*.
Available at [clusteringtheory.org](http://clusteringtheory.org)

Dua D., Karra Taniskidou E. (2018). *UCI Machine Learning Repository*
[http://archive.ics.uci.edu/ml]. Irvine, CA: University of California,
School of Information and Computer Science.

Fränti P., Mariescu-Istodor R., Zhong C. (2016). XNN graph,
In: *Proc. IAPR Joint Int. Workshop on Structural, Syntactic,
and Statistical Pattern Recognition*, Merida, Mexico,
*Lecture Notes in Computer Science* **10029**, pp. 207-217.
[doi:10.1007/978-3-319-49055-7_19](http://dx.doi.org/10.1007/978-3-319-49055-7_19)

Fränti P., Sieranoja S. (2018).
K-means properties on six clustering benchmark datasets,
*Applied Intelligence* **48**, 2018, pp. 4743-4759.
[doi:10.1007/s10489-018-1238-7](http://dx.doi.org/10.1007/s10489-018-1238-7)

Gagolewski M., Bartoszuk M., Cena A. (2016).
Genie: A new, fast, and outlier-resistant hierarchical clustering algorithm,
*Information Sciences* **363**, pp. 8-23.
[doi:10.1016/j.ins.2016.05.003](http://dx.doi.org/10.1016/j.ins.2016.05.003)

Graves D., Pedrycz W. (2010).
Kernel-based fuzzy clustering and fuzzy clustering:
A comparative experimental study,
*Fuzzy Sets and Systems* **161**(4), pp. 522-543.
[doi:10.1016/j.fss.2009.10.021](http://dx.doi.org/10.1016/j.fss.2009.10.021)

Karypis G., Han E.H., Kumar V. (1999).
CHAMELEON: A hierarchical clustering algorithm using dynamic modeling,
*IEEE Transactions on Computers* **32**(8), pp. 68-75.
[doi:10.1109/2.781637](http://dx.doi.org/10.1109/2.781637)

Ultsch A. (2005). Clustering with SOM: U\*C,
In: *Proc. Workshop on Self-Organizing Maps*, Paris, France, pp. 75-82.
