# getVertices() takes as input an m-Matrix small enough for an exhaustive search,
# (i.e. every possible vertex is calculated and then tested to determine whether
# or not it is within the stable region), a desired hypersphere radius, and a 
# cutoff that serves as a positive small number approximation for zero. Setting 
# the cutoff to zero can be a fine thing to do, but setting it to a small positive
# number is a more conservative approach for exploring parameter space, as the
# intersections between hyperplanes could produce very thin slices causing stable
# vertices to lie within the range of numerical error. The getVertices() function,
# and the other supporting functions in this file, are called by getDD(), and other
# functions within the getDD.R file. Descriptions of the supporting functions can
# be found below, preceding the functions they describe. 
#
# This code was written by Gianmarc Grazioli.
#
# The algorithms implemented in this code were developed by Yue Yu, 
# Gianmarc Grazioli, and Carter T. Butts, and were introduced in a manuscript 
# titled: "Local Graph Stability in Exponential Family Random Graph Models"
# by Yue Yu, Gianmarc Grazioli, Nolan E. Phillips, and Carter T. Butts
# https://arxiv.org/abs/1908.09470 
#

# This more straightforward approach is ideal for cases with only two
# sufficient statistics (e.g. the star graph representing the social
# network structure of a cult in section 3 of "Local Graph Stability 
# in Exponential Family Random Graph Models"). The getDD() function
# automatically checks for this. 
verticesForTwoStats<-function(mMatrix, radius, cutoff){
  solutions<-list()
  id<-diag(ncol(mMatrix))
  ctr<-1
  for (i in 1:nrow(mMatrix)) {
    for(j in 1:nrow(id)){
      hldr<-rbind(mMatrix[i,], id[j,]*radius)
      if(tryInvert(hldr)){
        solutions[[ctr]] <- solve(hldr, hldr[nrow(hldr),])
        ctr=ctr+1
        solutions[[ctr]] <- -solutions[[ctr-1]]
        ctr=ctr+1
      }
    }
  }
  out<-list()
  ctr=1
  for (i in 1:length(solutions)) {
    if(testWithMat(mMatrix, solutions[[i]], cutoff)){
      out[[ctr]] <- normIt(solutions[[i]])*radius
      ctr=ctr+1
    }
  }
  out<-do.call(rbind, out)
  out<-data.frame(out)
  out[!duplicated(out),]
}

# This function uses a more general approach to finding the vertices
# for network models with 3 or more sufficient statistics. 
getVertices<-function(mMatrix, radius, cutoff, dropRedundantVerts = TRUE){
  # Find all dim - 1 possible combinations of rows of M-matrix
  # that will form the rays constructed from halfspace intersections.
  planeCombos<-combn(1:nrow(mMatrix),ncol(mMatrix)-1)
  
  # Create identity matrix as source of one hot vectors that will
  # serve as hyperplanes that make up a hypercube for calculating
  # the intersection of rays with a half hypercube. These intersections 
  # between rays and the hypercube allow rays to be stored as points 
  # in parameter space, which can be rescaled to any length, so no 
  # information is lost. This code uses rescaling of the points so 
  # that the vertices outputted are intersections with a hypersphere. 
  cubeSideMat<-diag(ncol(mMatrix)) 
  # Construct matrices that store subsets of halfspaces whose intersections
  # are potential vertices in the vertex representation and store in matList:
  matList<-list()
  if(ncol(planeCombos) >= 100000) cat("Possible vertices =", ncol(planeCombos), "this could take a while...\n")
  for(i in 1:ncol(planeCombos)){
    if(i %% 100000 == 0) cat("Building matList",i, "of", ncol(planeCombos), " in getVertices()\n")
    matList[[i]]<-mMatrix[planeCombos[,i],]
  }
  # Here we solve for all vertices for which a solution exists. Since
  # the intersections are rays that go through the origin, a single 
  # point holds all information. Thus, if a solution exists (invertible),
  # a break statement is used to exit the loop calculating intersections
  # with every side of the half hypercube (we only need one point). 
  solnList<-list()
  planeComboListKeep<-list()
  ctr<-1
  for(i in 1:length(matList)){
    if(i %% 100000 == 0) cat("Solving vertex",i,"of",length(matList),"\n")
    # First test to see if the intersection of the ray with the
    # side of the hypercube is invertable. If it is, we use it
    # to find the intersection of the ray with the hypercube. If 
    # it is not, the intersection does not exist, and it will not
    # be included (e.g. two of the halfspaces are parallel).
    for (j in 1:nrow(cubeSideMat)) {
      matHldr <- rbind(matList[[i]], cubeSideMat[j,])
      if(tryInvert(matHldr)==TRUE){
        # Solution exists, so first get intersection with the
        # positive half cube
        solnList[[ctr]]<-solve(matHldr,matHldr[nrow(matHldr),]*radius)
        planeComboListKeep[[ctr]]<-planeCombos[,i]
        ctr<-ctr+1
        # Now get the solution with the negative half cube.
        # Need both because only one of these may represent
        # the side of the ray that is in the stable region.
        solnList[[ctr]] <- -solnList[[ctr-1]]
        planeComboListKeep[[ctr]] <- planeComboListKeep[[ctr-1]]
        ctr<-ctr+1
        break
      }
    }
  }
  # Normalize all solutions to a unit hypersphere:
  normList<-lapply(solnList,normIt)
  # Initialize some things before rescaling solutions
  # to the desired radius:
  outList<-list() # stores vertices that are in the stable region
  finalComboList<-list() # stores halfspace IDs that make up each stable vertex
  ctr<-1
  hldr<-vector()
  hldrIDs<-vector()
  # Now we test each vertex stored as a point in parameter
  # space to see if it lies within the stable region by 
  # multiplying it (dot product) with the M matrix. If all
  # elements in the resulting vector are less than or equal
  # to zero, then the vertex is stable, and we keep it.
  for(i in 1:length(normList)) {
    hldr<-normList[[i]]*radius
    hldrIDs<-planeComboListKeep[[i]]
    hldrIDs<-rownames(mMatrix[hldrIDs,])
      if(testWithMat(mMatrix, hldr, cutoff)==TRUE){
        outList[[ctr]]<-hldr
        finalComboList[[ctr]]<-hldrIDs
        ctr<-ctr+1
      }
  }
  myIDs<-do.call(rbind, finalComboList)
  IDnames<-vector()
  for (i in 1:ncol(myIDs)) {
    IDnames[i] <- paste("mRow.",i,sep="")
  }
  colnames(myIDs)<-IDnames
  out<-do.call(rbind, outList)
  out<-cbind.data.frame(out, myIDs)
  # Because the same vertex may appear multiple times due
  # to more than dim - 1 halfspaces intersecting at a single
  # ray, we use dropRedundantVertices() to remove the extra
  # copies of redundant vertices:
  if(dropRedundantVerts) out<-dropRedundantVertices(out)
  out
}

# Function to normalize a vector:
normIt<-function(vec){
  norm<-sqrt(sum(vec^2))
  vec/norm
}

# Function to see if a matrix is invertable:
tryInvert <- function(m) is.matrix(try(solve(m),silent=T))


# Function to test whether a vertex is in the stable region:
testWithMat<-function(mMatrix, vec, cutoff){
  if(is.vector(vec))! {
    myVec<-as.vector(vec)
  }else myVec<-vec
  all(mMatrix %*% myVec <= cutoff)
}

# We need a function to remove redundant vertices from the result. 
# This function does that for a single row of the V-representation.
# It is called by dropRedundantVertices():
dropRedundantWithOneVertex<-function(vertexDF, rowIndex, cutoff=.00001){
  dataColmax<-(ncol(vertexDF)+1)/2
  valueMat<-as.matrix(vertexDF[,1:dataColmax])
  refVec<-valueMat[rowIndex,]
  myBools<-vector(mode = "logical", length = nrow(vertexDF))
  for (i in 1:nrow(vertexDF)) {
    cutCurr <- sum((valueMat[i,] - refVec)^2)
    myBools[i] <- cutCurr >= cutoff
  }
  myBools[rowIndex] <- TRUE
  vertexDF[myBools,]
}

# This function removes extra vertices that appear more than
# once in the V-representation. This occurs when more than
# dim - 1 intersect at a particular ray:
dropRedundantVertices<-function(vertexDF, cutoff=.00001){
  ctr<-1
  mySize<-nrow(vertexDF)
  newDF<-vertexDF
  while (ctr <= mySize) {
    newDF<-dropRedundantWithOneVertex(newDF,ctr, cutoff = cutoff)
    mySize<-nrow(newDF)
    ctr<-ctr+1
  }
  newDF
}




