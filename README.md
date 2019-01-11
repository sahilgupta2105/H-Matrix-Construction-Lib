# Hierarchical Matrix construction library

The project was a result of my Dual Degree Project at IIT Madras with Prof. S. Sundar. This library takes a sparse matrix as input and constructs a Hierarchical Matrix representation for the matrix. The library uses graph clustering algorithms to find sub-matrices inside the matrix which can be represented in compressed form using reduced-SVD. The library is based on algebraic construction process and uses no geometric information, whatsoever. This process has been described in this [paper](https://ir.uiowa.edu/cgi/viewcontent.cgi?article=1381&context=etd). Although the library gives out a correct H-matrix, it is slightly slow compared to other highly-optimized libaries avaliable on internet such as HLibPro. 

Note: the library uses a supporting library called Eigen, to handle all manipulations related to linear algebra.

For more info look [here](https://sahilgupta2105.github.io/project3/).
