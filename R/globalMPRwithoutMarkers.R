globalMPRwithoutMarkers <-
function (baseData)
    {
      alleleFinal <- matrix(0,nrow(baseData),2)
      rownames(alleleFinal) <- rownames(baseData)
      i=1
      initPart <- localMPR(baseData[i:(i+99),])
      alleleFinal[1:70,] <- initPart[1:70,]
      i=i+70
      while(i+99<=nrow(baseData))
      {
        nextPart <- localMPR(baseData[i:(i+99),])
        med=initPart[71:100,1]==nextPart[1:30,1]
        if(mean(med)>0.5)
        {
          alleleFinal[i:(i+69),] <- nextPart[1:70,]
        }
        else
        {
          alleleFinal[i:(i+69),] <- nextPart[1:70,c(2,1)]
        }
        initPart <- nextPart
        i=i+70
        cat(mean(med), sep = "\n")
      }
      nextPart <- localMPR(baseData[(nrow(baseData)-149):nrow(baseData),])
      med=alleleFinal[(nrow(alleleFinal)-149):(nrow(baseData)-120),1]==nextPart[1:30,1]
      if(mean(med)>0.5)
      {
        alleleFinal[(nrow(baseData)-149):nrow(baseData),] <- nextPart[1:150,]
      }
      else
      {
        alleleFinal[(nrow(baseData)-149):nrow(baseData),] <- nextPart[1:150,c(2,1)]
      }
      cat(mean(med), sep = "\n")
      alleleFinal
    }
