

library(doSNOW)
library(foreach)
library(SWTools)
library(lubridate)
library(tidyverse) 
library(animation)


#  FUNCTIONS ----------------------------------




fileSetup <- function(scenVars) {
  scenVars[["pFolder"]][2] <- dir.exists(scenVars[["pFolder"]][1])
  scenVars[["SourceVersion"]][2] <- dir.exists(scenVars[["SourceVersion"]][1])
  scenVars[["outputFolder"]][2] <- dir.exists(paste0(scenVars[["pFolder"]][1],
                                                     "/",
                                                     scenVars[["outputFolder"]][1]))
  scenVars[["sourceProjectName"]][2] <- file.exists(paste0(scenVars[["pFolder"]][1],
                                                           "/",
                                                           scenVars[["sourceProjectName"]][1]))
  scenVars[["inputSetFile"]][2] <- file.exists(paste0(scenVars[["pFolder"]][1],
                                                      "/",
                                                      scenVars[["inputSetFile"]][1]))
  # check the main directory exists and set wd if it is
  if(scenVars[["pFolder"]][2]){
    setwd(scenVars[["pFolder"]][1])
    print("Project Directory Exists.")
    # Check if the outputs dir exists.
    if(scenVars[["outputFolder"]][2]=="FALSE"){
      # if not, then create one
      dir.create(
        path = paste0(scenVars[["pFolder"]][1],
                      "/",
                      scenVars[["outputFolder"]][1])
        
      )
      scenVars[["outputFolder"]][2] <- dir.exists(paste0(scenVars[["pFolder"]][1],
                                                         "/",
                                                         scenVars[["outputFolder"]][1]))
      print(paste0("Outputs directory created at: ",
                   scenVars[["pFolder"]][1],
                   "/",
                   scenVars[["outputFolder"]][1]))
    } else {
      print("Output Directory Exists.")
      
    }
    
  }
}

writeInputSetFile3 <- function (Population,scenVars,sf) {
  f <- paste0(scenVars[["pFolder"]][1],"/run",1,"/",
              scenVars[["inputSetFile"]][1])
  outfile <- file(f,"w")
  
  for (z in 1:nrow(sf)){
    # newVal <- Population*sf$Value[z]
    if(sf$Units[z]=="none"){
      s1 <- paste0(sf$Expression[z],"=",Population)
    }else{
      s1 <- paste0(sf$Expression[z],"=",newVal," ",sf$Units[z])
    }
    writeLines(s1,outfile)
  }
  close(outfile)
}



getWaterYearAnnual <- function(df){
  
  for(i in 1:nrow(df)){
    if(df$Mon[i]<7){
      df$WY[i] <- paste0(df$Year[i]-1,"-",df$Year[i])
    } else {
      df$WY[i] <- paste0(df$Year[i],"-",df$Year[i]+1)
    }
  }
  return(df$WY)
  
}

getQuadraticRoot <- function(scenVars,x) {
  df <- read_csv(paste0(scenVars[["pFolder"]][[1]],"/",
                        scenVars[["summaryFilename"]][1])
  )
  
  # 2. Fit the quadratic
  
  d <- data.frame(Population=df$Population,Residual=df$Residual)
  
  d$Population2 <- d$Population^2
  # xPoll <- seq(10,200, by = 5)
  if(x > 10){
    q <- quantile(d$Residual,0.99)
    d <- d %>% 
      filter(Residual<q)
    # l <- x-4
    # rg <- nrow(d)
    # d <- d[l:rg,, drop=FALSE]
  }
  qm <- lm( Residual ~ Population + Population2, data = d)
  # summary(qm)
  # 3. Predict the root (phi=0)
  root <- -(qm[["coefficients"]][["Population"]])/(2*qm[["coefficients"]][["Population2"]])
  # 4. Run model with new estimate
  lastRoot <- root
  return(lastRoot)
}
getResRowVeneer <- function(scenVars, i, Population) {
  df <- fortify.zoo(VeneerGetTS(TSURL = "/runs/latest/location/Functions/element/Functions/variable/Functions@ACT_Net_Take@$f_Net_Take"))
  colnames(df) <- c("Date","NetTake")
  df$Date <- date(df$Date)
  # WIndow the data for analysis period
  df <- df %>% 
    filter(df$Date >= date(dmy(scenVars[["startWindow"]][1])) ) #%>% 
  df <- df %>% 
    filter(df$Date  <= date(dmy(scenVars[["endWindow"]][1]) ))
  df$Mon <- month(df$Date)
  df$Year <- year(df$Date) 
  df$WY <- 1
  df$WY <- getWaterYearAnnual(df)
  dfAnn <-  aggregate(df$NetTake, by = list(df$WY),sum)
  dfMean <- mean(dfAnn$x)  
  res <- abs(scenVars[["target"]][1]- dfMean /1000)^2
  
  nextRow <- data.frame("Run"=i,
                        "Population"=Population,
                        "NetTake" =dfMean/1000 ,
                        "Residual"=res)
  if(i==1){ 
    write_csv(nextRow,
              paste0(scenVars[["pFolder"]][[1]],"/",
                     scenVars[["summaryFilename"]][1]),
              append = FALSE
    )
  }else{
    write_csv(nextRow,
              paste0(scenVars[["pFolder"]][[1]],"/",
                     scenVars[["summaryFilename"]][1]),
              append = TRUE
    )
    
  }
  return(nextRow)
}

getLinearRoot <- function(scenVars,x){
  df <- read_csv(paste0(scenVars[["pFolder"]][[1]],"/",
                        scenVars[["summaryFilename"]][1])
  )
  d <- df  %>%
    filter(Run>=(x-2))
  qm <- lm(Population  ~ Residual, data = d)
  # root <- (scenVars$target[[1]]-qm[["coefficients"]][["(Intercept)"]])/qm[["coefficients"]][["Population"]]
  root=qm[["coefficients"]][["(Intercept)"]]
  # summary(qm)
  # 4. Run model with new estimate
  lastRoot <- round(root,0)
  return(lastRoot)
  
  
  
  
  
}
getLinearRootTarget <- function(scenVars,x){
  df <- read_csv(paste0(scenVars[["pFolder"]][[1]],"/",
                        scenVars[["summaryFilename"]][1])
  )
  d <- df  %>%
    filter(Run>=(x-2))
  fit <- lm(NetTake ~ Population, data = d)
  a <- coef(fit)[2]
  c <- coef(fit)[1]
  root <- (42.7 - c) / a
  # root <- (scenVars$target[[1]]-qm[["coefficients"]][["(Intercept)"]])/qm[["coefficients"]][["NetTake"]]
  # root=qm[["coefficients"]][["(Intercept)"]]
  # summary(qm)
  # 4. Run model with new estimate
  lastRoot <- round(root,0)
  return(lastRoot)

  
}





