rm(list=ls())

# go to the directory that stored the saved the dataset and codes.
setwd("...")

library(Rcpp)
library(survRM2)
library(survival)
sourceCpp("like_two_markers.cpp")
source("HelpFunc.R")

#----------------------- read in the dataset -----------------------------#
data <- read.csv("data.csv")

#----------------------- estimation --------------------------------------#
# use fractional polynomials to expand markers based on t
# And prepare index for estimation
data.ALL <- data_prep(data, borrow_hn=NULL) 

# fit the proposed kernel smoothed concordance-assisted objective function
df.fracpoly = 6 # make sure to match with the df used in the C function "fracpoly"
no_marker = 2
fit_ker_xtbt <- ker_xtbt(data.ALL$data_cb, data.ALL$index, df=df.fracpoly, no_marker=no_marker)

#-------------------  generate dynamic score ----------------------------#

######## At landmark time s = 0
s = 0 
score <- get.score(data=data, fit=fit_ker_xtbt, s=s, df=df.fracpoly)

# Obtain AUC and RMST
res_roc <- get.ROC.AUC_realdata(data=data,
                                fit=fit_ker_xtbt,
                                s=s,w=1.2,
                                logscale=1,offset=1
                                )

tau.rmst <- quantile(data[unique(data$id),"Y"], 0.95)
res_d <- get.D.s_realdata(data=data, fit=fit_ker_xtbt, s=s, logscale=1, 
                          offset=1, tau.rmst=tau.rmst)


######## At landmark time s = 0.6
s = 0.6
score <- get.score(data=data, fit=fit_ker_xtbt, s=s, df=df.fracpoly)

# Obtain AUC and RMST
res_roc <- get.ROC.AUC_realdata(data=data,
                                fit=fit_ker_xtbt,
                                s=s,w=1.2,
                                logscale=1,offset=1
)


res_d <- get.D.s_realdata(data=data, fit=fit_ker_xtbt, s=s, logscale=1, 
                          offset=1, tau.rmst=tau.rmst)
