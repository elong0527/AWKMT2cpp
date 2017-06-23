
library("survAWKMT2")
library("AWKMT2cpp")

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

test <- "1_side"
type <- "perturbation"
nmethod <- 100000

time_res1 <- system.time(
  res1 <- AWKMT2(indata=D, tau=tau, c_first=0, c_last=4, c_by=2, method=type,
                 nmethod=nmethod, seed=1, v1=TRUE, v2=TRUE, test=test)
)

time_res2 <- system.time(
  res2 <- AWKMT2cpp(indata=D, tau=tau, c_first=0, c_last=4, c_by=2, method=type,
                    nmethod=nmethod, seed= NULL, test=test)
)

res <- data.frame( survAWKMT2 = unlist(res1[-1]), AWKMT2cpp = unlist(res2[-c(1,2)]) )
res$diff <- abs(res$survAWKMT2 - res$AWKMT2cpp)

save(res, time_res1, time_res2, file = paste(type, test,".Rdata", sep = "_"))


#---------------------------------------
#   Summary res information
#---------------------------------------

# load("permutation_1_side.Rdata")
# res$type <- c("permutation")
# res$test <- c("1_side")
# res$time_survAWKMT2 <- time_res1[3]
# res$time_AWKMT2cpp  <- time_res2[3]
# res1 <- res
# 
# load("permutation_2_side.Rdata")
# res$type <- c("permutation")
# res$test <- c("2_side")
# res$time_survAWKMT2 <- time_res1[3]
# res$time_AWKMT2cpp  <- time_res2[3]
# res2 <- res
# 
# load("perturbation_1_side.Rdata")
# res$type <- c("perturbation")
# res$test <- c("1_side")
# res$time_survAWKMT2 <- time_res1[3]
# res$time_AWKMT2cpp  <- time_res2[3]
# res3 <- res
# 
# load("perturbation_2_side.Rdata")
# res$type <- c("perturbation")
# res$test <- c("2_side")
# res$time_survAWKMT2 <- time_res1[3]
# res$time_AWKMT2cpp  <- time_res2[3]
# res4 <- res
# 
# res <- rbind(res1, res2, res3, res4)
# 
# save(res, file = "Validate.Rdata")
