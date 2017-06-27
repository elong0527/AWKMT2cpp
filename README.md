AWKMT2cpp
=========

Two-Sample Tests Based on Differences of Kaplan-Meier Curves

The package speed up the [survAWKMT2](https://CRAN.R-project.org/package=survAWKMT2) by rewriting the key function in C++. 

The paper to descrbied the method is:

[Uno H, Tian L, Claggett B, Wei LJ. A versatile test for equality of two survival functions based on weighted differences of Kaplan-Meier curves. Statistics in Medicine 2015, 34, 3680-3695.](http://onlinelibrary.wiley.com/doi/10.1002/sim.6591/abstract)

## Installation 

I assume that all auxillary packages are already installed

From an interactive R session:

```r
library(devtools)
install_github("elong0527/AWKMT2cpp")
library(AWKMT2cpp)
```

or directly download `AWKMT2cpp_1.0.1.tar.gz` and install the package from source

```r
install.packages("../AWKMT2cpp_1.0.1.tar.gz", repos = NULL, type = "source")
```
## Info

* The `AWKMT2cpp` package are 10X to 30X faster than `survAWKMT2`  
* The detail can be found in this [link](test_code/Validate_code.md)
