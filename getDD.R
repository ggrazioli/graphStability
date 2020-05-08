# This file contains a collection of functions used to obtain the 
# double description form of the single step stable region for a 
# given graph. The primary function in this file is getDD(), with
# most other functions serving a supportive role to getDD(). More
# detailed descriptions of each function are included preceding 
# their declaration in this file. 
#
# The algorithms implemented in this code were developed by Yue Yu, 
# Gianmarc Grazioli, and Carter T. Butts, and were introduced in a manuscript 
# titled: "Local Graph Stability in Exponential Family Random Graph Models"
# by Yue Yu, Gianmarc Grazioli, Nolan E. Phillips, and Carter T. Butts
# https://arxiv.org/abs/1908.09470 
#
# This code was written by Gianmarc Grazioli.
#
library(ergm)
library(sna)
source("getVertices.R")


# This function returns a raw M matrix, i.e. a typically redundant set of halfspaces,
# the intersection of which defines the stable region. 
getRawOneStepMmatrix<-function(targetGraph, rhsAsChar, dropRepeats = TRUE){
  n<-targetGraph$gal$n
  edgeCombos<-t(combn(1:n, 2))
  refStats<-summary(as.formula(paste("targetGraph",rhsAsChar,sep="~")))
  cs<-list()
  for (i in 1:nrow(edgeCombos)) {
    cs[[i]]<-ergm.godfather(as.formula(paste("targetGraph",rhsAsChar,sep="~")),changes = matrix(edgeCombos[i,], ncol = 2)) - refStats
  }
  cs<-do.call(rbind, cs)
  testList<-vector("list", nrow(cs))
  if(dropRepeats){
    for(i in 1:nrow(cs)){
      testList[[i]]<-cs[i,]
    }
    out<-unique(testList)
    out<-as.matrix(do.call(rbind, out))
  }else out <- cs
  myRowNames<-vector()
  for (i in 1:nrow(out)) {
    myRowNames[i]<-i
  }
  rownames(out)<-myRowNames
  out
}

# This function initializes a double description object, created from
# rows of a raw M matrix. At minimum, pass the first dim - 1 rows 
# from the raw M matrix, where dim = number of sufficient statistics. 
# If it is tractable to calculate all n choose m potential vertices, where
# n is nrow(rawMmatrix) and m = dim - 1, the entire raw M matrix be passed 
# to this function, and it will return a complete double description. 
initializeDD<-function(m, radius = 100, cutoff = .00001){
  out<-vector("list", length = 4)
  out[[1]]<-m
  out[[2]]<-m
  out[[3]]<-getVertices(m, radius, cutoff)
  out[[3]]<-dropRedundantVertices(out[[3]])
  out[[4]]<-getHVmatFromVertexObject(out[[3]])
  names(out)<-c("rawM", "hRep", "vRep", "HsharesV")
  out<-updateHrepFromVrep(out)
  out
}

# This function leverages the full capability of the double description
# to calculate the one step stable region of a given network.  
getDD<-function(targetGraph, rhsAsChar, radius = 100, cutoff = .00001, vCountMax = 10000, feelingLucky = FALSE){
  #vCountMax is the most vertices allowed to calculate exhaustively
  mInit<-getRawOneStepMmatrix(targetGraph, rhsAsChar)
  vCount<-1
  rowMax<-1
  for (i in (ncol(mInit) - 1):nrow(mInit)) {
    vCount<-choose(i, (ncol(mInit) - 1))
    rowMax<-i
    if(vCount > vCountMax) break 
  }
  cat("rowMax =", rowMax,"and nrow(mInit) =",nrow(mInit),".\n", nrow(mInit)-rowMax, "more H-rep rows to include after initialization.\n")
  cat("Calculating",vCount,"initial vertices...\n")
  DDinit<-initializeDD(mInit[1:rowMax,], radius = radius, cutoff = cutoff)
  DDinit[[1]]<-mInit #put the full raw M back into the DD object
  names(DDinit)<-names(DDinit)
  cat("Double description initialized, checking for closed convex cone...\n")
  newDD<-DDinit
  if(feelingLucky==FALSE){
    coneIsClosed<-checkForClosed(DDinit)
  }else {
    cat("Feeling lucky, punk?\n", "Ok, fine, I'm skipping the check for a closed cone.\n", sep = "")
    coneIsClosed<-TRUE
  }
  if(!coneIsClosed){
    stop("Either initial H-rep rows were insufficient for establishing a closed
         convex cone, or a convex cone does not exist for this system. 
         Try increasing vCountMax.\n")
  }
  if(rowMax < nrow(mInit)){
    closedStart <- rowMax + 1
    for (j in c(closedStart:nrow(mInit))) {
      newDD <- addNewHalfspaceToClosedCone(newDD, j, radius = radius, cutoff = cutoff)
      vChangeData <- doesNewHchangeV(newDD, j)
      newDD <- updateVrepFromNewHrep(newDD, vChangeData)
    }
  }
  newDD<-updateHrepFromVrep(newDD)
  newDD
}

# This simple function can be used to calculate the double description
# by calculating all possible vertices and then eliminating the unstable
# ones. This is not the most efficient approach, but may be fast enough
# for models with few enough sufficient statistics.  
getDDexhaustive<-function(targetGraph, rhsAsChar, radius = 100, cutoff = .00001){
  mInit<-getRawOneStepMmatrix(targetGraph, rhsAsChar)
  ddInit<-initializeDD(mInit, radius = radius, cutoff = cutoff)
  ddInit
}

# This function takes a double description, and returns a new double description 
# where all unnecessary rows from the halfspace representation are removed. Only 
# the rows used to construct the vertex representation are kept. 
updateHrepFromVrep<-function(DDinput){
  allIDs<-as.vector(as.matrix(DDinput$vRep)[,((ncol(DDinput$vRep)+1)/2+1):ncol(DDinput$vRep)])
  allIDs<-sort(as.integer(unique(allIDs)))
  allIDs<-as.character(allIDs)
  newHrep<-DDinput$hRep[allIDs,]
  out<-list(DDinput$rawM, newHrep, DDinput$vRep, DDinput$HsharesV)
  names(out)<-names(DDinput)
  out
}

# This function tests whether or not adding a new row to the halfspace
# representation will change the vertex representation. If the new row
# has no effect on the vertices, the row is redundant. The inputs are
# an existing double description and the index indicating which row from
# the rawM element from the input will be incorporated into hRep. 
# It returns a list with two elements: a boolean indicating whether or not
# V is changed by the new H, and a vector of booleans that indicates which
# rows should be kept in the new version of V:
doesNewHchangeV<-function(DDinput, mIndex, cutoff = .00001){
  newH<-rbind(DDinput$hRep, DDinput$rawM[mIndex,])
  rownames(newH)<-c(rownames(DDinput$hRep), rownames(DDinput$rawM)[mIndex])
  vAsMatrix<-as.matrix(DDinput$vRep[,1:ncol(DDinput$hRep)])
  vertexTestResults<-vector(mode = "logical", length = nrow(DDinput$vRep))
  # Test which vertices are within the stable region:
  for (i in 1:nrow(DDinput$vRep)) {
    vertexTestResults[i] <- testWithMat(newH, vAsMatrix[i,], cutoff)
  }
  out<-list(!all(vertexTestResults), vertexTestResults)
  names(out)<-c("vChanged", "vRowsToKeep")
  out
}

# This function takes a double description, and returns....
# It takes the output from doesNewHchangeV() as second input...
# It assumes that the halfspace representation has already
# been changed, thus you should use doesNewHchangeV() to first
# confirm that the new row should be added to the halfspace
# representation, then add the row to the halfspace rep, then
# finally use this function to update the vertex representation
# using the the new halfspace representation. 
updateVrepFromNewHrep<-function(DDinput, changeData){
  out <- DDinput
  out$vRep <- DDinput$vRep[changeData$vRowsToKeep,]
  out$HsharesV <- getHVmatFromVertexObject(out$vRep) 
  out 
}

# This function is used to update the vertex representation
# using the h representation currently stored in the double
# description object used as input. 
updateVrepFromCurrentHrep<-function(DDinput, cutoff = .00001){
  vAsMatrix<-as.matrix(DDinput$vRep[,1:ncol(DDinput$hRep)])
  vertexTestResults<-vector(mode = "logical", length = nrow(DDinput$vRep))
  # Test which vertices are within the stable region:
  for (i in 1:nrow(DDinput$vRep)) {
    vertexTestResults[i] <- testWithMat(DDinput$hRep, vAsMatrix[i,], cutoff)
  }
  out<-DDinput
  out$vRep <- DDinput$vRep[vertexTestResults,]
  out
}

# This function tests whether the convex cone is closed by
# testing whether or not a cycle, of size >= the number of
# dimensions in parameter space (number of sufficient stats), 
# is present in the adjacency matrix of which which halfspaces 
# share a vertex. 
checkForClosed<-function(DDinput){
  cat("Testing DD for closed convex cone.\n")
  adjMat<-DDinput$HsharesV
  minCycleSize<-ncol(DDinput$rawM)
  if(nrow(DDinput$hRep) > ncol(DDinput$hRep)){
    maxCycleSize <- nrow(DDinput$rawM)
  }else maxCycleSize <- minCycleSize
  out<-FALSE
  for (i in minCycleSize:maxCycleSize) {
    cat("Testing for semicycle of size",i,"\n")
    out<-(summary(adjMat~cycle(i))>=1)
    if(out) break
  }
  out
}

# This function leverages the double description to add a new halfspace and the associated
# new vertices to the double description input object. Note that, if the new halfspace does
# not eliminate any of the vertices in the input double description, it is redundant. If the
# new halfspace is found to be redundant, addNewHalfspaceToClosedCone() will return the input
# double description unchanged.
addNewHalfspaceToClosedCone<-function(DDinput, rowIndex, radius = 100, cutoff = .00001){
  changeData<-doesNewHchangeV(DDinput, rowIndex)
  out<-DDinput 
  if(changeData$vChanged==TRUE){
    # Append the new row to the H-representation:
    newH<-rbind(DDinput$hRep, DDinput$rawM[rowIndex,])
    cat("Appending row",rowIndex, "i.e.", DDinput$rawM[rowIndex,],"to hRep\n")
    rownames(newH)<-c(rownames(DDinput$hRep), rownames(DDinput$rawM)[rowIndex])
    
    # Get all the vertices that will be eliminated, then extract the rows
    # from the halfspace representation that are affected by the deletions.
    # The new row is also appended onto the hRowsInvolved.
    vertsToDelete<-DDinput$vRep[!changeData$vRowsToKeep,]
    hRowsInvolvedIndices<-unique(as.vector(as.matrix(vertsToDelete[,(ncol(out$hRep)+1):ncol(out$vRep)])))
    hRowsInvolvedIndices<-as.character(sort(as.integer(hRowsInvolvedIndices)))
    hRowsInvolved<-rbind(newH[hRowsInvolvedIndices,], newH[nrow(newH),])
    rownames(hRowsInvolved)<-c(rownames(newH[hRowsInvolvedIndices,]), rownames(newH)[nrow(newH)])
    
    # Get all new vertices produced by affected halfspaces:
    newVertices<-getVertices(hRowsInvolved, radius, cutoff)
    # Update vertex representation by dropping eliminated vertices from vRep
    # then appending the new vertices onto vRep:
    newV<-rbind(DDinput$vRep[changeData$vRowsToKeep,], newVertices)
    out$hRep<-newH
    out$vRep<-newV
    out<-updateVrepFromCurrentHrep(out, cutoff = cutoff)
    out$HsharesV<-getHVmatFromVertexObject(out$vRep)
  }
  out
}

# This function adds a new halfspace to hRep and updates the vertices exhaustively 
# (creates and tests every possible vertex). Use this only until the hypercone is
# closed (possible combos grow factorially). Once the hypercone is closed, use the
# function addNewHalfspaceToClosedCone().
addNewHalfspaceExhaustive<-function(DDinput, rowIndex, radius = 100, cutoff = .00001){
  newH<-rbind(DDinput$hRep, DDinput$rawM[rowIndex,])
  rownames(newH)<-c(rownames(DDinput$hRep), rownames(DDinput$rawM)[rowIndex])
  newV<-getVertices(newH, radius, cutoff)
  newHsharesV<-getHVmatFromVertexObject(newV)
  out<-list(DDinput$rawM, newH, newV, newHsharesV)
  names(out)<-names(DDinput)
  out
}

# This function takes the vertex representation from a double description and outputs
# an adjacency matrix that indicates which halfspaces share at least one vertex. This 
# object is used to determine whether or not a closed hypercone is present.  
getHVmatFromVertexObject<-function(vertInput){
  if(ncol(vertInput) %% 2 == 0){
    cat("Error in getHVmatFromVertexObject(): vertex object must have odd number of columns")
    break
  }
  firstHindex = (ncol(vertInput)+1)/2 + 1
  indexMatrix <- vertInput[,firstHindex:ncol(vertInput)]
  uniqueIndices <- unique(as.vector(as.matrix(indexMatrix)))
  possibleIndexCombos <- combn(uniqueIndices, 2)
  uniqueRows <- unique(indexMatrix)
  out<-matrix(0, nrow = length(uniqueIndices), ncol = length(uniqueIndices))
  for (i in 1:ncol(possibleIndexCombos)) {
    for(j in 1:nrow(uniqueRows)){
      myVec <- vector()
      for (k in 1:ncol(uniqueRows)) {
        myVec[k] <- as.character(uniqueRows[j,k])
      }
      if(is.element(as.character(possibleIndexCombos[1,i]), myVec) && 
         is.element(as.character(possibleIndexCombos[2,i]), myVec)){
        out[which(uniqueIndices == as.character(possibleIndexCombos[1,i])), 
            which(uniqueIndices == as.character(possibleIndexCombos[2,i]))] = 1
        out[which(uniqueIndices == as.character(possibleIndexCombos[2,i])), 
            which(uniqueIndices == as.character(possibleIndexCombos[1,i]))] = 1
        break
      }
    }
  }
  out
}




