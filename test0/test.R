library(AWKMT2cpp)

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

## Test funcKM2

Rcpp::sourceCpp("funcKM2_eigen.cpp")

weight <- rep(1, length(x1) )
res <- funcKM2_eigen(t_idx, status1, x1, weight)
res1 <- do.call(cbind, res)
res2 <- survAWKMT2:::funcKM2(data = data.frame(time = x1, status = status1, arm = 1), t_idx, status1, x1, weight = NULL)

all( as.numeric(res1 - res2) < 1e-5)

## Test funcAWKMT2

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
  print(par[i,])
  par$res[i] <- test_funcAWKMT2_eigen(par$test[i], par$type[i])
}
## Test AWKMT2



library("survAWKMT2")
library("AWKMT2cpp")
# Rcpp::sourceCpp("foo.cpp")
# source("AWKMT2.R")


test <- "2_side"
type <- "permutation"
nmethod <- 2
nmethod <- 10

# res1 <- AWKMT2(indata=D, tau=tau, c_first=0, c_last=4, c_by=2, method=type,
#        nmethod=nmethod, seed=1, v1=TRUE, v2=TRUE, test=test)
#

res2 <- AWKMT2cpp(indata=D, tau=tau, c_first=0, c_last=4, c_by=2, method=type,
                   nmethod=nmethod, seed= NULL, test=test)

res2
