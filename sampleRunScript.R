# The purpose of this R script is to demonstrate how to use the function
# getDD() to obtain the double description form of the single step stable 
# region for a given graph. 
#
# The algorithms implemented in this code were developed by Yue Yu, 
# Gianmarc Grazioli, and Carter T. Butts, and were introduced in a manuscript 
# titled: "Local Graph Stability in Exponential Family Random Graph Models"
# by Yue Yu, Gianmarc Grazioli, Nolan E. Phillips, and Carter T. Butts
# https://arxiv.org/abs/1908.09470 (this link will be updated)
#
# The example used is a network topology which has been observed in amyloid 
# fibril structures, called a 2-ribbon. For details on the significance
# of this fibril topology, please see the article titled:
# Network-Based Classification and Modeling of Amyloid Fibrils
# by Gianmarc Grazioli, Yue Yu, Megha H. Unhelkar, Rachel W. Martin, and Carter T. Butts
# https://doi.org/10.1021/acs.jpcb.9b03494
# 
# This code was written by Gianmarc Grazioli.
# 

source("getDD.R")
require(ergm.changestats)

# Load the network object:
load("2ribbonNetwork.RData")
plot(targetGraph)

# Define ERGM terms
myTerms<-"edges+kstar(2)+nsp(1:2)" #2-ribbon

# Calculate the double description that defines the 1-step stable region:
myStableRegionDD<-getDD(targetGraph, myTerms)
print(myStableRegionDD)
