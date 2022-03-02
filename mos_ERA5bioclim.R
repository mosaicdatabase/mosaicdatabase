#' ####################################################################### #
#' PROJECT: [VITAL-RATE TRADE-OFFS] Climate Data Retrieval
#' AUTHOR: [Erik Kusch]
#' PURPOSE: Loops over all records in COMADRE and COMPADRE and downloads queried climate variables from ERA5-Land for all unique combinations of locations and time-windows for which ERA5-Land provides data
#' ####################################################################### #
rm(list=ls()) # clear R environment
TRes <- "month"

# PREAMBLE ================================================================
## Functions ----------------------------------------
`%nin%` <- Negate(`%in%`) # negation of the %in% function

## CRAN packages ------------------------------------
install.load.package <- function(x) { # function which automatically installs/loads package
  if (!require(x, character.only = TRUE)){ # if the package is not installed
    install.packages(x, repos='http://cran.us.r-project.org') # install the package from CRAN
  }
  require(x, character.only = TRUE) # load package
}
package_vec <- c( # vector of package names to be installed/loaded from CRAN
  "Rcompadre", # required for access to COMADRE and COMPADRE
  "Rage", # needed for generation time computation
  "devtools" # required for GitHub installation of KrigR
)
sapply(package_vec, install.load.package) # apply package function to all package names
rm(package_vec) # remove package names from environment

## KrigR --------------------------------------------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # If KrigR is not yet installed
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") # disable benign warnings as stop-flags in installation process
  devtools::install_github("https://github.com/ErikKusch/KrigR/tree/Develop", force = TRUE) # install KrigR, required for climate data retrieval, from GitHub
}
library(KrigR) # load KrigR package

## API credentials ----------------------------------
try(source("PersonalSettings.R")) # This file contains my API credentials and numberofCore specification. It is not shared on GitHub, due to privacy reasons
# CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- as.numeric(readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER."))
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if(!exists("Cores")){ # Core check: if number of cores for parallel processing has not been set yet
  Cores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

## Directories --------------------------------------
Dir.Base <- getwd() # read the current directory, this is our project directory
Dir.Data <- file.path(Dir.Base, TRes) # register Data directory as subdirectory to project directory
if(!dir.exists(Dir.Data)){dir.create(Dir.Data)} # create Data subdirectory if it doesn't exist yet

## Land Mask ----------------------------------------
## LandMask is needed for checking whether COM(P)ADRE locations fall onto land or not
if(!file.exists(file.path(Dir.Base, "LandMask.nc"))) { # If landmask has not been downloaded yet
  ## a generic function call to KrigR to obtain ERA5-Land data at a global scale
  LandMask <- download_ERA(
    Variable = "2m_temperature",
    DataSet = "era5-land",
    DateStart = "1995-01-01",
    DateStop = "1995-01-31",
    TResolution = "month",
    TStep = 1,
    FileName = "LandMask",
    Dir = Dir.Base,
    API_User = API_User,
    API_Key = API_Key,
    Cores = Cores
  )
}else{ # If landmask has already been downloaded previously
  LandMask <- raster(file.path(Dir.Base, "LandMask.nc")) # load LandMask
}

# DATA DOWNLOAD FUNCTION ==================================================
FUN.Clim <- function(Variable = NULL, # which climate variable to retrieve data for
                     force = FALSE, # whether to check for which MPMs we already have the relevant data
                     MPMFrames = TRUE, # whether to export time-series of variable data for each MPM
                     Cores = 1, # over how many cores to parallelise download calls
                     TRes = "month", # temporal resolution of data from which to calculate aggregate metrics ("month" or "day" recommended)
                     database = NULL # for which database to obtain data
){
  print(paste("Starting data retrieval for", Variable))
  
  ### Year-Limitation ----
  if(TRes == "month" | TRes == "year"){
    LimYear <- 1981
  }else{
    LimYear <- 1950
  }
  
  ### MPM Download  ----
  data_db <-  database
  database <- tolower(database@version[["Database"]])
  message(paste("Figuring out which MPMs in", database, "require climate data."))
  
  data_LatLonYr <- data.frame(matrix(NA, nrow=length(data_db$SpeciesAccepted), ncol=6)) # create data frame in which to store the important information
  # Loop extracting from R
  for(db_iter in 1:length(data_db$SpeciesAccepted)){ # entry loop: extract relevant information for each MPM
    data_LatLonYr[db_iter,1] <- data_db$MatrixID[[db_iter]] # MatrixID
    data_LatLonYr[db_iter,2] <- data_db$Lat[[db_iter]] # Latitude
    data_LatLonYr[db_iter,3] <- data_db$Lon[[db_iter]] # Longitude
    data_LatLonYr[db_iter,4] <- suppressWarnings(as.numeric(data_db$StudyStart[[db_iter]])) # Study Start, suppressing warnings due to faulty entries
    data_LatLonYr[db_iter,5] <- suppressWarnings(as.numeric(data_db$StudyEnd[[db_iter]])) # Study End, suppressing warnings due to faulty entries
    data_LatLonYr[db_iter,6] <- try(suppressWarnings(Rage::gen_time(data_db$mat[[db_iter]]@matU, data_db$mat[[db_iter]]@matF)), silent = TRUE) # Study End, suppressing warnings due to faulty entries
  } # end of entry loop
  colnames(data_LatLonYr) <- c("MatrixID", "Lat", "Lon", "StudyStart", "StudyEnd", "GenTime") # set column names
  data_LatLonYr$db <- database # add data base identifier to the data frame
  MPMs_df <- data_LatLonYr # save data as MPMs_df
  rm("data_LatLonYr", "data_db", "db_iter") # remove superfluous objects from environment
  
  ### Generation Time Fixes ----
  MPMs_df$GenTime <- as.numeric(MPMs_df$GenTime)
  
  ### Data merging ----
  if(file.exists(paste0(database, "_MPMs_climate.csv"))){ # if we already have some version of the data
    data_file <- read.csv(paste0(database, "_MPMs_climate.csv"))[-1] # load existing data
    MPMs_df <- rbind(data_file, MPMs_df[which(MPMs_df$MatrixID %nin% data_file$MatrixID),]) # append new records
    rm(data_file) # remove data_file from environment
  }
  write.csv(MPMs_df, paste0(database, "_MPMs_climate.csv")) # write .csv to project directory
  
  ### Data Checks ----
  ## Here, we check for which MPMs ERA5-Land data can be obtained
  MPMs_df$Download <- 2 # create download query column, 2 indicates insufficient metadata data
  MPMs_df$Download[rowSums(!is.na(MPMs_df[,c(2:5)])) == 4] <- 1 # 1 indicates sufficient metadata availability (lat/lon, years)
  MPMs_df$Download[which(as.numeric(as.character(MPMs_df$StudyEnd)) < LimYear & MPMs_df$Download == 1)] <- 3 # 3 indicates study end date pre-dating first available year of ERA5-land data availability
  DownloadQueries <- which(MPMs_df$Download == 1) # store rownumbers of MPMs with sufficient metadata
  CoordinateCheck_df <- MPMs_df[DownloadQueries,2:3] # extract Lat/Lot for MPMs with sufficient metadata
  OverCheck_df <- extract(x = stack(LandMask), y = data.frame(CoordinateCheck_df$Lon, CoordinateCheck_df$Lat), method = "bilinear") # extract values of LandMask for all MPMs with sufficient metadata
  OutofBounds <- which(is.na(OverCheck_df)) # these are the points which have download query 1, but fall outside of ERA5-land shape
  MPMs_df$Download[DownloadQueries][OutofBounds] <- 4 # 4 indicates locations outside of land mask
  rm("DownloadQueries", "CoordinateCheck_df", "OverCheck_df", "OutofBounds") # remove superfluous objects from environment
  message(paste(database, "contains", nrow(MPMs_df), "MPM entries. Out of these",  table(MPMs_df$Download)[4], "fall outside of terrestrial regions. Data for", table(MPMs_df$Download)[3], "was collected before", LimYear, ", while", table(MPMs_df$Download)[2], "do not provide sufficient metadata for the climate data retrieval. For", table(MPMs_df$Download)[1], "MPM entries, climate data can be retrieved."))
  
  ### Column/Variable check ----
  if(Variable == "2m_temperature"){ # if targeted Variable is 2m_temperature
    VariableName <- "air_temperature" # assign variable name air_temperature (names starting with numbers lead to .csv artefacts)
  }else{
    VariableName <- Variable
  }
  if(VariableName %nin% colnames(MPMs_df)){ # If currently targeted variable is not stored in local .csv file yet
    MPMs_df$XYZ <- NA # create column in last position
    colnames(MPMs_df)[ncol(MPMs_df)] <- VariableName # assign variable name to column name
    MPMs_df$XYZ <- NA # create another column in last position
    colnames(MPMs_df)[ncol(MPMs_df)] <- paste0(VariableName, "_SD") # assign variable name for standard deviation of data to column name
    MPMs_df$XYZ <- NA # create column in last position
    colnames(MPMs_df)[ncol(MPMs_df)] <- paste0("T_", VariableName) # assign variable name to column name
    MPMs_df$XYZ <- NA # create another column in last position
    colnames(MPMs_df)[ncol(MPMs_df)] <- paste0("T_", VariableName, "_SD") # assign variable name for standard deviation of data to column name
  }
  
  ### force check ----
  if(isTRUE(force)){ # if force argument is TRUE
    DownsNeeded <- MPMs_df$MatrixID[MPMs_df$Download == 1] # query download for all MPMs for which we can obtain ERA5-Land data
  }else{  # if force argument is not TRUE
    DownsNeeded <- MPMs_df$MatrixID[is.na(MPMs_df[,which(colnames(MPMs_df) == VariableName)]) & MPMs_df$Download == 1] # identify which of the entries already have the necessary data and query downloads for all other MPMs for which we can obtain ERA5-Land data
  }
  
  ### duplicate check ---
  ## isolate all metadata necessary for ERA5-Land data retrieval for all MPMs for which downloads will be queried
  Duplicate_df <- with(MPMs_df[MPMs_df$MatrixID %in% DownsNeeded, ], 
                       data.frame(MatrixID = MatrixID,
                                  Lat = Lat,
                                  Lon = Lon,
                                  StudyStart = StudyStart,
                                  StudyEnd = StudyEnd)
  )
  Duplication_df <- unique(Duplicate_df[,-1]) # find the unique combinations of metadata
  
  ### Download loop ----
  message(paste("There are",  nrow(Duplication_df), "unique combinations of locations and years for which climate data will be retrieved."))
  
  for(Down_Iter in 1:nrow(Duplication_df)){ # download loop: stage ERA5-downloads
    if(nrow(Duplication_df) == 0){break}
    message(paste("############### DOWNLOAD ", Down_Iter, "/", nrow(Duplication_df)))
    ## Duplictaes
    ID_vec <- c() # Indicator vector for position of duplicated MatrixIDs
    for(ID_Iter in 1:nrow(Duplicate_df)){ # figuring out which MatrixIDs belong to this unique download specification
      ID_vec <- c(ID_vec, 
                  ifelse(sum(Duplicate_df[ID_Iter,-1] == Duplication_df[Down_Iter, ]) == 4, TRUE, FALSE) # is TRUE if Down_Iter-row in Duplication_df is the exat same as ID_Iter row in Duplicate_df
      )
    }
    ID_Down <- Duplicate_df$MatrixID[ID_vec] # extract IDs of duplicates
    
    ## Metadata
    Iter_df <- Duplication_df[Down_Iter,] # extract relevant metadata
    Iter_df$MatrixID <- ID_Down[1] # set arbitrary MatrixID for download_ERA() function call
    Y_1 <- as.numeric(as.character(Iter_df$StudyStart)) # extract start year for time-window
    if(Y_1 < LimYear){Y_1 <- LimYear} # set to LimYear (currently earliest year for which ERA5-Land data is available) if it predates LimYear
    Y_2 <- as.numeric(as.character(Iter_df$StudyEnd)) # extract study end year for time window
    
    ## ERA5-Download
    if(startsWith(x = Variable_Iter, prefix = "total_") | Variable_Iter == "runoff" | Variable_Iter == "snowfall"){PrecipFix <- TRUE}else{PrecipFix <- FALSE} # apply PrecipFix of download_ERA() if target variable is total_precipitation
    
    ## DATE SEQUENCE for checking whether singular download can be executed
    if(Variable=="total_precipitation" & TRes == "day" | TRes == "hour"){ # date sequence for total precipitation at sub-monthly temporal resolution would stage one more year of data than asked for originally in Singular download
      DateMod <- 366
    }else{
      DateMod <- 0
    }
    # date sequence between start of data availability and end of MPM time-window
    if(TRes == "day" | TRes == "hour"){ # repeat each day by 24 (to represent hours)
      Dates_seq <- seq.Date(from = as.Date(paste0(LimYear, "-01-01")), to = as.Date(paste0(Y_2, "-12-31"))+DateMod, by = "day")
      Dates_seq <- rep(Dates_seq, each = 24)
    }else{
      Dates_seq <- seq.Date(from = as.Date(paste0(LimYear, "-01-01")), to = as.Date(paste0(Y_2, "-12-31"))+DateMod, by = TRes)
    }
    if(length(Dates_seq) > 1e5){ # check if download is too big or not for singular download
      print("Parallel Download")
      SingularDL <- FALSE
    }else{
      print("Singular Download")
      SingularDL <- TRUE
    }
    Clim_ras <- download_ERA(
      Variable = Variable, 
      PrecipFix = PrecipFix,
      DateStart =  paste0(LimYear, "-01-01"),
      DateStop = paste0(Y_2, "-12-31"),
      TResolution = TRes,
      Extent = Iter_df,
      Buffer = 0.1,
      ID = "MatrixID",
      Dir = Dir.Data,
      FileName = as.character(Down_Iter),
      API_User = API_User,
      API_Key = API_Key,
      TryDown = Inf,
      Cores = 1,
      SingularDL = SingularDL,
      verbose = FALSE
    )
    
    ## compilation of raw data if queried
    if(isTRUE(MPMFrames)){ # if raw data storage is queried
      setwd(Dir.Data) # set working directory to temporary directory
      Temps_vec <- colMeans(extract(x = stack(Clim_ras), y = data.frame(Iter_df$Lon, Iter_df$Lat), buffer = 2e4)[[1]], na.rm = TRUE) # extract climate data for each layer
      ## figuring out intervals
      Down_start <- lubridate::date(paste0(LimYear, "-01-01"))
      Down_end <- lubridate::date(paste0(Y_2, "-12-31"))
      T_seq <- seq(Down_start, Down_end, by = TRes)
      if(TRes == "hour"){
        T_seq <- paste(rep(T_seq, each = 24), str_pad(0:23, 2, 0, side = "left"), sep = "_")
      }
      ## data frame and saving
      Raw_df <- data.frame(
        Date = T_seq,
        Mean = as.numeric(Temps_vec)
      )
      colnames(Raw_df)[2] <- VariableName
      FNames <- file.path(Dir.Data, paste0(ID_Down, "_Raw.csv")) # file names for which to save the above information
      if(file.exists(FNames[1])){ # if any of these already exists
        FName1 <- read.csv(FNames[1])[-1] # read the already existing data frame
        Raw_df <- cbind(FName1, Raw_df)[,-ncol(FName1)-1] # append current data as columns
      }
      for(write_Iter in 1:length(ID_Down)){ # write data frame to filenames
        write.csv(x = Raw_df, file = FNames[write_Iter]) 
      }
      setwd(Dir.Base)
    }
    
    ## figuring out intervals
    Down_start <- lubridate::date(paste0(LimYear, "-01-01"))
    Down_end <- lubridate::date(paste0(Y_2, "-12-31"))
    TFull_seq <- seq(Down_start, Down_end, by = TRes)
    Down_start <- lubridate::date(paste0(Y_1, "-01-01"))
    TStudy_seq <- seq(Down_start, Down_end, by = TRes) #!!! this is where generation time can be worked in
    
    ## Calculation of Aggregate Metrics
    Mean_ras <- mean(Clim_ras[[which(TFull_seq %in% TStudy_seq)]]) # calculate mean value
    SD_ras <- stackApply(Clim_ras[[which(TFull_seq %in% TStudy_seq)]], 
                         rep(1, nlayers(Clim_ras[[which(TFull_seq %in% TStudy_seq)]])), sd) # calculate standard deviation
    if(Y_1 == LimYear){
      T_1Mean <- NA
      T_1SD <- NA
    }else{
      T_1_Mean_ras <- mean(Clim_ras[[which(TFull_seq %nin% TStudy_seq)]]) # calculate mean value
      T_1_SD_ras <- stackApply(Clim_ras[[which(TFull_seq %nin% TStudy_seq)]], 
                               rep(1, nlayers(Clim_ras[[which(TFull_seq %nin% TStudy_seq)]])), sd) # calculate standard deviation
      T_1Mean <- mean(extract(x = T_1_Mean_ras, y = data.frame(Iter_df$Lon, Iter_df$Lat), buffer = 2e4)[[1]], na.rm = TRUE)
      T_1SD <- mean(extract(x = T_1_SD_ras, y = data.frame(Iter_df$Lon, Iter_df$Lat), buffer = 2e4)[[1]], na.rm = TRUE)
    }
    
    ## data extraction and adding to metadata dataframe
    Clim_vec <- c(mean(extract(x = Mean_ras, y = data.frame(Iter_df$Lon, Iter_df$Lat), buffer = 2e4)[[1]], na.rm = TRUE),
                  mean(extract(x = SD_ras, y = data.frame(Iter_df$Lon, Iter_df$Lat), buffer = 2e4)[[1]], na.rm = TRUE), 
                  T_1Mean, T_1SD
    )
    
    MPMs_df[MPMs_df$MatrixID %in% ID_Down, 
            which(colnames(MPMs_df) == VariableName):(which(colnames(MPMs_df) == VariableName)+3)
    ] <- rep(Clim_vec, each = length(ID_Down))
    
    
    ## writing metadata frame
    write.csv(MPMs_df, file = file.path(Dir.Base,  paste0(database, "_MPMs_climate.csv")))
    ## cleaning of local files
    unlink(file.path(Dir.Data, list.files(Dir.Data, pattern = ".nc")), recursive = TRUE)
  } # end of download loop
  # sink(paste0(Variable, "FINISHED.txt"))
  # print("This is just a hack placeholder to enable a while loop further down.")
  # sink()
  message(paste(Variable, "data is up-to-date for all MPMS in for which climate data can be obtained."))
}

# FUNCTION CALLS ==========================================================
Variable_vec <- c("2m_temperature", "total_precipitation", "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2", "snowfall", "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4", "runoff", "total_evaporation")

## the following contains a while loop to restart the climate data retrieval function upon error with the server connection
db <- cdb_fetch("compadre")
for(Variable_Iter in Variable_vec){ # loop over all variables for which we want data
  # while(!file.exists(file.path(Dir.Base, paste0(Variable_Iter, "FINISHED.txt")))){ # as long as indicator for finished retrieval is not present
    closeAllConnections() # close parallel core connections
    try( # run data retrieval function
      FUN.Clim(Variable = Variable_Iter,
               force = FALSE,
               MPMFrames = TRUE,
               Cores = Cores,
               TRes = TRes,
               database = db)
    )
  # }
}

## the following contains a while loop to restart the climate data retrieval function upon error with the server connection
db <- cdb_fetch("comadre")
for(Variable_Iter in Variable_vec){ # loop over all variables for which we want data
  # while(!file.exists(file.path(Dir.Base, paste0(Variable_Iter, "", "FINISHED.txt")))){ # as long as indicator for finished retrieval is not present
    closeAllConnections() # close parallel core connections
    try( # run data retrieval function
      FUN.Clim(Variable = Variable_Iter, 
               force = FALSE, 
               MPMFrames = TRUE, 
               Cores = Cores,
               TRes = TRes,
               database = db)
    )
  # } 
}
