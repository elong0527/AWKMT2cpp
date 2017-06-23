Compare results and performances with AWKMT2cpp and survAWKMT2
================
2017-06-23

Evaluate Function `AWKMT2`
--------------------------

The source code to evaluate `AWKMT2` is saved in `Validate_AWKMT2cpp.R` with 100,000 replications.

``` r
load("Validate.Rdata")
res
```

    ##                             survAWKMT2 AWKMT2cpp    diff         type   test time_survAWKMT2 time_AWKMT2cpp
    ## crude_pvalue_T1_1_side         0.36679   0.36750 0.00071  permutation 1_side         1536.30         117.69
    ## crude_pvalue_T2_1_side         0.44692   0.44796 0.00104  permutation 1_side         1536.30         117.69
    ## bona_fide_pvalue_T1_1_side     0.40902   0.41025 0.00123  permutation 1_side         1536.30         117.69
    ## bona_fide_pvalue_T2_1_side     0.49318   0.49448 0.00130  permutation 1_side         1536.30         117.69
    ## crude_pvalue_T1_2_side         0.52866   0.52692 0.00174  permutation 2_side         1580.15         133.48
    ## crude_pvalue_T2_2_side         0.44588   0.44489 0.00099  permutation 2_side         1580.15         133.48
    ## bona_fide_pvalue_T1_2_side     0.55082   0.54916 0.00166  permutation 2_side         1580.15         133.48
    ## bona_fide_pvalue_T2_2_side     0.46766   0.46654 0.00112  permutation 2_side         1580.15         133.48
    ## crude_pvalue_T1_1_side1        0.37514   0.36074 0.01440 perturbation 1_side         1577.47          70.18
    ## crude_pvalue_T2_1_side1        0.46326   0.44570 0.01756 perturbation 1_side         1577.47          70.18
    ## bona_fide_pvalue_T1_1_side1    0.41796   0.40211 0.01585 perturbation 1_side         1577.47          70.18
    ## bona_fide_pvalue_T2_1_side1    0.51354   0.49350 0.02004 perturbation 1_side         1577.47          70.18
    ## crude_pvalue_T1_2_side1        0.55754   0.54376 0.01378 perturbation 2_side         1649.60          73.76
    ## crude_pvalue_T2_2_side1        0.47581   0.46180 0.01401 perturbation 2_side         1649.60          73.76
    ## bona_fide_pvalue_T1_2_side1    0.58443   0.56982 0.01461 perturbation 2_side         1649.60          73.76
    ## bona_fide_pvalue_T2_2_side1    0.50224   0.48656 0.01568 perturbation 2_side         1649.60          73.76

Evaluate Inner Function
-----------------------

### Load Data

``` r
library(AWKMT2cpp)
library(survAWKMT2)
## Test Data
D        = survival::pbc[1:312, c(2,3,4)] #The pbc data from 'survival' package
D$time   = as.numeric(D$time)
D$status = as.numeric(D$status==2)
D$trt    = as.numeric(D$trt==2)
names(D) = c("time", "status", "arm")
tau      = max(D[D[,2]==1,1])
nmethod  = 100 #Recommended to specify at least 10000 (default) or larger.
method   = "permutation"

indata=D
tau=tau
c_first=0
c_last=4
c_by=0.1
method="permutation"
nmethod=nmethod
seed=1
v1=TRUE
v2=TRUE
test="1_side"

#-- Get tau1 --
tau1 = 0

#-- Get tau2 --
tau2 = tau

#-- Get unique_time, n_times --
t_idx   = unique(sort(c(indata$time, 0, tau1, tau2, tau)))
n_times = length(t_idx)

#-- Get crange --
crange   = seq(c_first, c_last, by=c_by)

x0 <- subset(indata, arm == 0)$time
x1 <- subset(indata, arm == 1)$time

status0 <- subset(indata, arm == 0)$status
status1 <- subset(indata, arm == 1)$status
```

### Evaluate Inner function `funcKM2`

``` r
Rcpp::sourceCpp("funcKM2_eigen.cpp")

weight <- rep(1, length(x1) )
res <- funcKM2_eigen(t_idx, status1, x1, weight)
res1 <- do.call(cbind, res)
res2 <- survAWKMT2:::funcKM2(data = data.frame(time = x1, status = status1, arm = 1), t_idx, status1, x1, weight = NULL)

all( as.numeric(res1 - res2) < 1e-5)
```

    ## [1] TRUE

### Evaluate Inner function `funcAWKMT2`

``` r
Rcpp::sourceCpp("funcAWKMT2_eigen.cpp")

test_funcAWKMT2_eigen <- function(test, type){
  res1 <- funcAWKMT2_eigen(t_idx = t_idx, status0 = status0, status1 = status1, x0 = x0, x1 = x1,
                           tau1 = tau1, tau2 = tau2, crange = c(0.5,0.6,0.7), test = test, type = "observe", obs_survdiff = x1 )

  if(type != "observe"){
    res1 <- funcAWKMT2_eigen(t_idx = t_idx, status0 = status0, status1 = status1, x0 = x0, x1 = x1,
                       tau1 = tau1, tau2 = tau2, crange = c(0.5,0.6,0.7), test = test, type = type, obs_survdiff = res1$survdiff )
  }
  indata <- data.frame(time = c(x0,x1), status = c(status0, status1), arm = c(rep(0, length(status0)), rep(1, length(status1))) )
  res2 <- survAWKMT2:::funcAWKMT2_c1(indata = indata, t_idx = t_idx, tau1 = tau1, tau2 = tau2,
                             crange = c(0.5,0.6,0.7), test = test, type = type)

  if(type != "observe"){
    all(cbind(res1$V1, res1$V2) - res2 < 1e-5)
  }else{
    all(cbind(res1$V1, res1$V2) - res2$stats < 1e-5)
  }
}

test <- c("1_side", "2_side")
type <- c("observe","permutation")
par <- expand.grid(test = test, type = type, stringsAsFactors = FALSE)
par$res <- FALSE
for(i in 1:nrow(par)){
  par$res[i] <- test_funcAWKMT2_eigen(par$test[i], par$type[i])
}

par
```

    ##     test        type  res
    ## 1 1_side     observe TRUE
    ## 2 2_side     observe TRUE
    ## 3 1_side permutation TRUE
    ## 4 2_side permutation TRUE
