localMPR <-
function (baseData, allele.matrix = NULL, maxNStep = 5, 
            returnNumRecomEvents=FALSE, verbose = FALSE, strEND = "\n", ...) 
  {
    if (is.null(allele.matrix)) 
      allele.matrix <- base2Allele(baseData)
    if (nrow(baseData) == ncol(allele.matrix)) 
      allele.matrix <- t(allele.matrix)
    if (nrow(baseData) != nrow(allele.matrix)) 
      stop("nrow(baseData)!=nrow(allele.matrix), allele.matrix error!!!")
    ALLELE <- allele.matrix
    genoData <- base2Geno(baseData, allele.matrix)
    
    genoData <- matrix(t(genoData),nrow(genoData),ncol(genoData))
    mylog <- .C("core_localMPR",as.integer(genoData),as.integer(nrow(genoData)),as.integer(ncol(genoData)),
                                                as.integer(maxNStep),mylog=integer(nrow(genoData)))$mylog
    ALLELE[mylog==1, ] <- ALLELE[mylog==1, c(2, 1)]
    
    if (verbose) 
      cat("\tDone.", strEND)
    rownames(ALLELE) <- rownames(baseData)
    if (returnNumRecomEvents) {
      numR=NumRecomEvents(baseData,ALLELE)
      list(allele = ALLELE, numR)}
    else ALLELE
  }
