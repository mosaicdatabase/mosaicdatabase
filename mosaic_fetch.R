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


url <- "https://github.com/mosaicdatabase/mosaicdatabase/blob/main/new_IDs.rds?raw=true"
Indices <- readRDS(url(url, method="libcurl"))
Indices

sapply(package_vec, install.load.package) # Load the vector of libraries based on the install function - T & Fs returned
rm(package_vec,install.load.package) # Once loaded, kill the load function and name vector

comp <- cdb_fetch("compadre")
coma <- cdb_fetch("comadre")
cpa <- comp
cma <- coma




metaIndexer <- function(i){
  return(1:8+(8*i)-8)
}


mos_fetch <- function(id_key){
  id_key <- id_key
  url_link <- paste("https://raw.githubusercontent.com/mosaicdatabase/mosaicdatabase/main/MOSAIC_", id_key, ".csv", sep="")
  mosaic <- read.csv(url(url_link))
  
  spp <- mosaic$specieslist
  
  compadre <- gsub("_", " ", cpa@data$SpeciesAuthor)
  compadre <- sub('^(\\w+\\s+\\w+).*', '\\1', compadre)
  
  comadre <- gsub("_", " ", cma@data$SpeciesAuthor)
  comadre <- sub('^(\\w+\\s+\\w+).*', '\\1', comadre)
  
  taxMeta <- list()
  
  Kingdom_a <- list()
  Phylum_a <- list()
  Class_a <- list()
  Order_a <- list()
  Family_a <- list()
  Genus_a <- list()
  Species_a <- list()
  
  for(i in 1:length(spp)){
    Kingdom_a[[i]] <- try(cma@data$Kingdom[match(spp[[i]], comadre)])
    Phylum_a[[i]] <- try(cma@data$Phylum[match(spp[[i]], comadre)])
    Class_a[[i]] <- try(cma@data$Class[match(spp[[i]], comadre)])
    Order_a[[i]] <- try(cma@data$Order[match(spp[[i]], comadre)])
    Family_a[[i]] <- try(cma@data$Family[match(spp[[i]], comadre)])
    Genus_a[[i]] <- try(cma@data$Genus[match(spp[[i]], comadre)])
    Species_a[[i]] <- try(cma@data$Species[match(spp[[i]], comadre)])
  }
  
  Kingdom_a <- unlist(Kingdom_a)
  Phylum_a <- unlist(Phylum_a)
  Class_a <- unlist(Class_a)
  Order_a <- unlist(Order_a)
  Family_a <- unlist(Family_a)
  Genus_a <- unlist(Genus_a)
  Species_a <- unlist(Species_a)
  
  Kingdom_p <- list()
  Phylum_p <- list()
  Class_p <- list()
  Order_p <- list()
  Family_p <- list()
  Genus_p <- list()
  Species_p <- list()
  
  for(i in 1:length(spp)){
    Kingdom_p[[i]] <- try(cpa@data$Kingdom[match(spp[[i]], compadre)])
    Phylum_p[[i]] <- try(cpa@data$Phylum[match(spp[[i]], compadre)])
    Class_p[[i]] <- try(cpa@data$Class[match(spp[[i]], compadre)])
    Order_p[[i]] <- try(cpa@data$Order[match(spp[[i]], compadre)])
    Family_p[[i]] <- try(cpa@data$Family[match(spp[[i]], compadre)])
    Genus_p[[i]] <- try(cpa@data$Genus[match(spp[[i]], compadre)])
    Species_p[[i]] <- try(cpa@data$Species[match(spp[[i]], compadre)])
  }
  cpa@data$Kingdom[match(spp[[1003]], compadre)]
  
  Kingdom_p <- unlist(Kingdom_p)
  Phylum_p <- unlist(Phylum_p)
  Class_p <- unlist(Class_p)
  Order_p <- unlist(Order_p)
  Family_p <- unlist(Family_p)
  Genus_p <- unlist(Genus_p)
  Species_p <- unlist(Species_p)
  
  
  TaxaMeta <- data.frame(matrix(NA, ncol=7, nrow=length(mosaic$specieslist)))
  rownames(TaxaMeta) <- mosaic$specieslist
  colnames(TaxaMeta) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for(i in 1:length(mosaic$specieslist)){
    if(is.na(Kingdom_p[i]==FALSE)){
      TaxaMeta[i,1]<- Kingdom_a[i]
      TaxaMeta[i,2]<- Phylum_a[i]
      TaxaMeta[i,3]<- Class_a[i]
      TaxaMeta[i,4]<- Order_a[i]
      TaxaMeta[i,5]<- Family_a[i]
      TaxaMeta[i,6]<- Genus_a[i]
      TaxaMeta[i,7]<- Species_a[i]
    }else if(is.na(Kingdom_p[2]==TRUE)){
      TaxaMeta[i,1]<- Kingdom_p[i]
      TaxaMeta[i,2]<- Phylum_p[i]
      TaxaMeta[i,3]<- Class_p[i]
      TaxaMeta[i,4]<- Order_p[i]
      TaxaMeta[i,5]<- Family_p[i]
      TaxaMeta[i,6]<- Genus_p[i]
      TaxaMeta[i,7]<- Species_p[i]
    }
  }
  
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
                          taxaMetadat = "data.frame",
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
                     taxaMetadat = TaxaMeta,
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