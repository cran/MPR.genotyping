NumRecomEvents <-
function (baseData, allele.matrix, genoData = NULL) 
  {
    if (is.null(genoData)) {
      genoData <- base2Geno(baseData, allele.matrix)
    }
    .C("core_NumRecomEvents",as.integer(genoData),
       as.integer(nrow(genoData)),as.integer(ncol(genoData)),R=integer(1))$R
  }
