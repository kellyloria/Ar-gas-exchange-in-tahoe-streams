# KAL made minor edits from Bob Hall et al. workflow 2024-08

#Assumes inHg pressure measurements in excel sheet.
#If mmHg:
###Go into the `gather_data` function
###Change the line: press <- mean(sub$Pressure) * 25.4
###To: press <- mean(sub$Pressure)

#Runs assumming an excel sheet with the following columns, but will run even if some are missing: 
###N2.Ar, He.Ar, Ar.Kr_83, Ar.Kr_84, X40, X32, X28, X4, X83, X84, X29.28
###Column names are not case sensitive
library(readxl)
library(readr)
library(dplyr)
setwd("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/")

######EDIT THESE, and be sure to comment out the pressure line that is wrong for you (whether your pressure is in inHg or mmHg). It currently runs for inHg of pressure.
MIMSdata <- "kelly_prelim_argon_31july24_KAL.xlsx"
rawFile <- "Kelly_rawdat_20240801.csv"


# rename for each save:
saveFile <- "2024_KAL_prelimdat_processed_09.csv"
pressure <- "inHg"
#pressure <- "mmHg"
source("/Users/kellyloria/Documents/UNR/Reaeration/AR_code_repo/mims_gas_functions_wHeKr.R")
##########################################################################################
#Formats and reads in the raw files in the best manner for the gather_data function
read_raw_file <- function(filepath){
  nCols <- max(count.fields(filepath, sep = ','))
  firstRead <- read.csv(filepath, header = FALSE, col.names = paste0("V",seq_len(nCols)), 
                        fill = TRUE, stringsAsFactors = FALSE)
  trimmedRead <- firstRead[cumsum(complete.cases(firstRead)) != 0,]
  header <- as.character(trimmedRead[1,])
  colnames(trimmedRead) <- header
  trimmedRead <- trimmedRead[-1, ]
  trimmedRead$Index <- as.numeric(as.character(trimmedRead$Index))
  trimmedRead$"N2/Ar" <- trimmedRead$"28" / trimmedRead$"40"
  trimmedRead$"O2/Ar" <- trimmedRead$"32" / trimmedRead$"40"
  trimmedRead$"29/28" <- as.numeric(trimmedRead$"29") / as.numeric(trimmedRead$"28")
  trimmedRead$"34/32" <- as.numeric(trimmedRead$"34") / as.numeric(trimmedRead$"32")
  trimmedRead$"30/28" <- as.numeric(trimmedRead$"30") / as.numeric(trimmedRead$"28")
  range <- c(min(trimmedRead$Index), max(trimmedRead$Index))
  return(list(trimmedRead, range))
}


# Averages data from the raw document
gather_data <- function(excelfile, rawFile) {
  name <- read_excel(excelfile)
  name <- name[!is.na(name$SampleID), ]
  name_range <- c(min(name$Index), max(name$Index))
  
  # Check for a raw file from the run date, then check contents to see if it includes all relevant indices
  pars <- read_raw_file(rawFile)
  rawname <- pars[[1]]
  raw_range <- pars[[2]]
  if (raw_range[1] > name_range[1] || raw_range[2] < name_range[2]) {
    print(paste0(rawFile, " does not contain all of the indices of your excel sheet. Check to make sure this is the raw file you'd like to use and either replace 'rawFile' with the correct file or run the other script to just use the excel sheet."))
  }
  datedf <- data.frame()
  for (j in unique(name$Samp)) {
    sub <- name[(name$Samp == j), ]
    if (nrow(sub) == 0) {
      next
    }
    indexes <- c(min(sub$Index), max(sub$Index))
    data <- rawname[(rawname$Index >= indexes[1]) & (rawname$Index <= indexes[2]), ]
    if (nrow(data) == 0) {
      print(paste("No data found for Samp:", j))
      next
    }
    data <- data[, -which(names(data) %in% c("Index"))]
    data$Time <- as.numeric(as.POSIXct(data$Time, format = "%m/%d/%Y %H:%M:%S"))
    
    data[] <- lapply(data, as.numeric)
    
    # Sneak in the inHg to mmHg correction
    newdat <- data.frame(lapply(data, function(x) type.convert(as.character(x), as.is = TRUE)), stringsAsFactors = FALSE)
    newdat$Time <- as.POSIXct(newdat$Time, origin = "1970-01-01")
    options(warn = -1)
    temp <- mean(sub$Temp, na.rm = TRUE)
    if (pressure == "inHg") {
      press <- mean(sub$Pressure, na.rm = TRUE) * 25.4
    } else if (pressure == "mmHg") {
      press <- mean(sub$Pressure, na.rm = TRUE)
    } else {
      print("Choose your pressure setting at the beginning of the code")
    }
    
    metadat <- data.frame(
      "Samp" = sub[["Samp"]][1], 
      "SampleID" = sub[["SampleID"]][1], 
      "Pressure" = press, 
      "Temp" = temp, 
      "Calibnum" = sub[["Calibnum"]][1],
      "Depth" = sub[["Depth"]][1],
      "WatDens" = watdens(temp), 
      "O2Sat" = osat1(temp, press), 
      "N2Sat" = nsat(temp, press), 
      "ArSat" = arsat(temp, press), 
      "O2.ArSat" = osat1(temp, press) / arsat(temp, press),
      "N2.ArSat" = nsat(temp, press) / arsat(temp, press),
      "HeSat" = hesat(temp, press),
      "KrSat_83" = krsat_83(temp, press),
      "KrSat_84" = krsat_84(temp, press),
      "He.ArSat" = hesat(temp, press) / arsat(temp, press),
      "Ar.KrSat_83" = arsat(temp, press) / krsat_83(temp, press),
      "Ar.KrSat_84" = arsat(temp, press) / krsat_84(temp, press)
    )
    
    options(warn = 0)
    
    metadat$Calibnum[is.na(metadat$Calibnum)] <- as.numeric(as.character(sub[["Sampnum"]][1][is.na(sub[["Calibnum"]][1])]))
    
    tempdf <- merge(metadat, newdat)
    
    datedf <- rbind(datedf, tempdf)
  }
  
  return(datedf)
}

avgdData <- gather_data(MIMSdata, rawFile)
str(avgdData)


targCols <- c("N2.Ar","He.Ar", "Ar.Kr_83", "Ar.Kr_84", "X40", "X32", "X28", "X4", "X83", "X84")
satCols <- c("N2.ArSat", "He.ArSat", "Ar.KrSat_83", "Ar.KrSat_84", "ArSat", "O2Sat", "N2Sat","HeSat","KrSat_83","KrSat_84") # where is the Ar/O2 
#Checks to make sure targCols and satCols are equal lengths (should correspond 1:1)
if(length(targCols) != length(satCols)){
  print("targCols and satCols need to be the same length.")
}

#Calculations concentrations for the targCols and satCols lists above
getRatioGeneral <- function(df, targCol, satCol){
  newcolname <- paste0(targCol, ".Conc")
  
  if(length(which(!is.na(match(tolower(colnames(df)), 
                               tolower(targCol))))) == 1){
    
    #Pull out differences of each calibnum
    for (i in 1:length(unique(df$Calibnum))){
      sub <- df[df$Calibnum == i | df$Calibnum == i+1,]
      
      calibs <- sub[grep("std", tolower(sub$SampleID)),]
      #Linear between low and high temps
      offsetfun <- lm(calibs[[satCol]] ~ 
                        calibs[, which(!is.na(match(tolower(colnames(calibs)), 
                                                    tolower(targCol))))])
      coeff1 <- coef(summary(offsetfun))[1]
      coeff2 <- coef(summary(offsetfun))[2]
      
      #Multiply currents
      sub[[newcolname]] <- 
        coeff1 + sub[, which(!is.na(match(tolower(colnames(sub)), 
                                          tolower(targCol))))]*coeff2
      
      df[[newcolname]][!is.na(base::match(df$Samp, sub$Samp))]<- sub[[newcolname]]
    }
    successTargs <- targCol
    failTargs <- NA
  }else{
    successTargs <- NA
    failTargs <- targCol
  }
  
  return(list(df, successTargs, failTargs))
}
#Calculates O2:Ar ratio
getRatioO2Ar <- function(df){
  newcolname <- "O2.Ar.Conc"
  
  if(length(which(!is.na(match(tolower(colnames(df)), 
                               tolower("O2.Ar"))))) == 1){
    
    for (i in 1:length(unique(df$Calibnum))){
      sub <- df[df$Calibnum == i | df$Calibnum == i+1,]
      
      calibs <- sub[grep("std", tolower(sub$SampleID)),]
      #Linear between 0 and all temps
      yvals <- c(0, calibs[, which(!is.na(match(tolower(colnames(calibs)), 
                                                tolower("O2.ArSat"))))])
      xvals <- c(0, calibs[, which(!is.na(match(tolower(colnames(calibs)), 
                                                tolower("O2.Ar"))))])
      offsetfun <- lm(yvals ~ xvals)
      coeff1 <- coef(summary(offsetfun))[1]
      coeff2 <- coef(summary(offsetfun))[2]
      
      #Multiply currents
      sub[[newcolname]] <- coeff1 + sub[, which(!is.na(match(tolower(colnames(sub)), 
                                                             tolower("O2.Ar"))))]*coeff2
      
      df[[newcolname]][!is.na(base::match(df$Samp, sub$Samp))]<- sub[[newcolname]]
    }
    successTargs <- "O2.Ar"
    failTargs <- NA
  }else{
    successTargs <- NA
    failTargs <- "O2.Ar"
  }
  return(list(df, successTargs, failTargs))
}
#Calculates 29:28 ratio
getRatio29.28 <- function(df){
  newcolname <- "del15N"
  
  if(length(which(!is.na(match(tolower(colnames(df)), 
                               tolower("X29.28"))))) == 1){
    
    for (i in 1:length(unique(df$Calibnum))){
      sub <- df[df$Calibnum == i | df$Calibnum == i+1,]
      
      calibs <- sub[grep("std", tolower(sub$SampleID)),]
      #Linear between 0 and all temps
      stdVals <- calibs[, which(!is.na(match(tolower(colnames(calibs)), 
                                             tolower("X29.28"))))]
      Rstd <- mean(stdVals)
      
      #Multiply currents
      sub[[newcolname]] <- (sub[, which(!is.na(match(tolower(colnames(sub)), 
                                                     tolower("X29.28"))))]/Rstd - 1)*1000
      
      df[[newcolname]][!is.na(base::match(df$Samp, sub$Samp))]<- sub[[newcolname]]
    }
    successTargs <- "X29.28"
    failTargs <- NA
  }else{
    successTargs <- NA
    failTargs <- "X29.28"
  }
  
  return(list(df, successTargs, failTargs))
}

successTargs <- c(NA)
failTargs <- c(NA)

#Runs the three functions above
if (length(targCols) == length(satCols)) {
  successTargs <- c()
  failTargs <- c()
  
  for (i in 1:length(targCols)) {
    allResults <- getRatioGeneral(avgdData, targCols[i], satCols[i])
    avgdData <- allResults[[1]]
    successTargs <- c(successTargs, allResults[[2]])
    failTargs <- c(failTargs, allResults[[3]])
  }
} else {
  stop("The number of target columns (targCols) and saturation columns (satCols) must be the same.")
}

allResults <- getRatio29.28(avgdData)
avgdData <- allResults[[1]]
successTargs <- c(successTargs, allResults[[2]])
failTargs <- c(failTargs, allResults[[3]])

#Reorders columns (Misc sample info, current values, saturation values, concentration values)
concentrations <- grep("Conc|del", colnames(avgdData))
saturations <- grep("Sat", colnames(avgdData))
others <- grep("Conc|Sat|del|Wat|Time", colnames(avgdData), invert = TRUE)
avgdData <- avgdData[, c(others, saturations, concentrations)]

## Saves data to a csv file
# write.csv(avgdData, saveFile, quote = FALSE)

failTargs <- failTargs[!is.na(failTargs)]
successTargs <- successTargs[!is.na(successTargs)]

print(paste0("Program ran successfully for ", 
             paste(unlist(successTargs), collapse = ", "), 
             "."), quote = FALSE)
print(paste0("Saved to ", saveFile, "."), quote = FALSE)
if(length(failTargs)>0){
  print(paste0("Program failed for ", failTargs, 
               ". Please check column names if you expected this to run."), quote = FALSE)
}

####
# end of temp edit 
####