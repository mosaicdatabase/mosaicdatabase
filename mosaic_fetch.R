#' #                                                                                                        -----
#'*   PROJECT NAME    [Fetch_ATP]
#'*   AUTHOR          [Connor Bernard. connor.bernard@zoo.ox.ac.uk.]
#'*   LAST UPDATE     [16 Feb 2022]
#'*   PURPOSE         [Downloading and Coercing the data in r]
#'*   DEPENDENCIES    [NA]
#' #                                                                                                        -----
#-------------------------------
#---------------------

#  Libraries
install.load.package <- function(x) { # Automate installs & load packages from CRAN where needed
  if (!require(x, character.only = TRUE)){ # If package is not yet installed, then install package from CRAN
    install.packages(x, repos='http://cran.us.r-project.org')
  }
  require(x, character.only = TRUE) # Require loads packages that are downloaded
}

package_vec <- c( # vector of package/library names - note: CRAN-dependent (no GitHub, local, &c.)
  "RCurl"
)


url <- "https://github.com/mosaicdatabase/mosaicdatabase/blob/main/MatrixIDs.rds?raw=true"
Indices <- readRDS(url(url, method="libcurl"))
Indices

sapply(package_vec, install.load.package) # Load the vector of libraries based on the install function - T & Fs returned
rm(package_vec,install.load.package) # Once loaded, kill the load function and name vector


metaIndexer <- function(i){
  return(1:8+(8*i)-8)
}


mos_fetch <- function(id_key){
  id_key <- id_key
  url_link <- paste("https://raw.githubusercontent.com/mosaicdatabase/mosaicdatabase/main/MOSAIC_", id_key, ".csv", sep="")
  mosaic <- read.csv(url(url_link))
  setClass("mosaic_meta",
           representation(value = "character",
                          author = "character",
                          year = "character",
                          journal = "character",
                          doi = "character", 
                          database = "character",
                          mosaic = "character",
                          notes = "logical")
  )
  
  theMetMet <- list()
  
  for(i in 1:14){
    theMetMet[[i]] <- new("mosaic_meta",
                          value = unlist(mosaic[,metaIndexer(i)[1]], use.names = FALSE),
                          author = unlist(mosaic[,metaIndexer(i)[2]], use.names = FALSE),
                          year = unlist(mosaic[,metaIndexer(i)[3]], use.names = FALSE),
                          journal = unlist(mosaic[,metaIndexer(i)[4]], use.names = FALSE),
                          doi = unlist(mosaic[,metaIndexer(i)[5]], use.names = FALSE),
                          database = unlist(mosaic[,metaIndexer(i)[6]], use.names = FALSE),
                          mosaic = unlist(mosaic[,metaIndexer(i)[7]], use.names = FALSE),
                          notes = unlist(mosaic[,metaIndexer(i)[8]], use.names = FALSE)
    )
  }
  
  biomass <- theMetMet[[1]]
  height <- theMetMet[[2]]
  growthdet <- theMetMet[[3]]
  regen <- theMetMet[[4]]
  dimorph <- theMetMet[[5]]
  matsyst <- theMetMet[[6]]
  hermaph <- theMetMet[[7]]
  seqherm <- theMetMet[[8]]
  dispcap <- theMetMet[[9]]
  disptype <- theMetMet[[10]]
  modedisp <- theMetMet[[11]]
  dispclass <- theMetMet[[12]]
  volancy <- theMetMet[[13]]
  aquadep <- theMetMet[[14]]
  
  class(biomass)
  
  setClass("mosaic_base",
           representation(species = "character",
                          index = "list",
                          biomass = "mosaic_meta",
                          height = "mosaic_meta",
                          growthdet = "mosaic_meta",
                          regen = "mosaic_meta",
                          dimorph = "mosaic_meta", 
                          matsyst = "mosaic_meta",
                          hermaph = "mosaic_meta",
                          seqherm = "mosaic_meta",
                          dispcap = "mosaic_meta",
                          disptype = "mosaic_meta",
                          modedisp = "mosaic_meta",
                          dispclass = "mosaic_meta",
                          volancy = "mosaic_meta",
                          aquadep = "mosaic_meta")
  )
  
  mosiac_main <- new("mosaic_base",
                     species = unlist(mosaic[,113], use.names = FALSE),
                     index = Indices,
                     biomass = biomass,
                     height = height,
                     growthdet = growthdet,
                     regen = regen,
                     dimorph = dimorph, 
                     matsyst = matsyst,
                     hermaph = hermaph,
                     seqherm = seqherm,
                     dispcap = dispcap,
                     disptype = disptype,
                     modedisp = modedisp,
                     dispclass = dispclass,
                     volancy = volancy,
                     aquadep = aquadep
  )
  
  
  return(mosiac_main)
}