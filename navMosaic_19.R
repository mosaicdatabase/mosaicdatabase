#' #                                                                                                        -----
#'*   PROJECT NAME    [Fetch_ATP]
#'*   AUTHOR          [Connor Bernard. connor.bernard@zoo.ox.ac.uk.]
#'*   LAST UPDATE     [14 Feb 2022]
#'*   PURPOSE         [Downloading and Coercing the data in r]
#'*   DEPENDENCIES    [NA]
#' #                                                                                                        -----
#-------------------------------
#---------------------
rm(list = ls()) # Clear workspace
if(!is.null(dev.list())) dev.off() # Clear plots
cat("\014") # Clear console


library(devtools)
source_url("https://raw.githubusercontent.com/mosaicdatabase/mosaicdatabase/main/navMosaic_19.R")



library(devtools)
source_url("https://raw.githubusercontent.com/mosaicdatabase/mosaicdatabase/main/mosaic_fetch.R")
mosaic <- mos_fetch("v1.0.0")
rm(metaIndexer, mos_fetch,url,Indices)


cpa <- cdb_fetch("compadre")
cma <- cdb_fetch("comadre")



# Linking MOSAIC to COMADRE, COMPADRE, + PADRINO records
equivalency_A <- function(x, y){
  any(x==y)
}

IndexToRecord <- function(indexValue){#Index-to-Index 
  return(which(unlist(lapply(mosaic@index, equivalency_A, indexValue))))
}

IndexToRecord(242472)

RecordToIndex <- function(mosaicIndex){#Index-to-Index 
  return(mosaic@index[[mosaicIndex]])
}

RecordToIndex(885)



# Promping options for the different routes
match(RecordToIndex(885), cpa$MatrixID)
cpa$mat[[match(RecordToIndex(885), cpa$MatrixID)]]

RecordToIndex(5)
RecordToIndex(2)



mosaic@biomass


# Simple trait search

traitAllRecords <- function(trait){
  return(slot(slot(mosaic, trait),"value"))
}

traitAllRecords("biomass")
traitAllRecords("volancy")

plot(as.numeric(traitAllRecords("biomass")))


# Note the extreme outlier
max(as.numeric(traitAllRecords("biomass")), na.rm = TRUE)
maxOL_index <- which(as.numeric(traitAllRecords("biomass"))==max(as.numeric(traitAllRecords("biomass")), na.rm = TRUE))


plot(as.numeric(traitAllRecords("biomass"))[-maxOL_index])
plot(log(as.numeric(traitAllRecords("biomass"))[-maxOL_index]))
hist(log(as.numeric(traitAllRecords("biomass"))[-maxOL_index]), breaks=15)




# Species Trait Summary
# NOTE THE SENSITIVITY OF THIS MEASURE TO CHANGES IN THE METADATA -- the 1:2 element that is removed
speciesSummary  <- function(speciesNameLatin){
  speciesIndex <- which(mosaic@species==speciesNameLatin)
  metadataNames <- slotNames(mosaic)[-c(1:2)]
  meta <- data.frame(matrix(NA, nrow=1, ncol=length(metadataNames)))
  colnames(meta) <- metadataNames
  for(i in 1:length(metadataNames)){
    meta[,i] <- slot(slot(mosaic, metadataNames[[i]]), "value")[[speciesIndex]]
  }
  return(meta)
}

speciesSummary("Aepyceros melampus")
speciesSummary("Aepyceros melampus")
sppListTemp <- mosaic@species[2:20]


# Queerying multiple species:

multiSpp <- function(sppList){
  summary <- lapply(sppList, speciesSummary)
  framesummary <- data.frame(matrix(NA, nrow=length(summary), ncol=length(summary[[1]])))
  colnames(framesummary) <- colnames(summary[[1]])
  for(i in 1:length(summary)){
    framesummary[i,] <- summary[[i]]
  }
  return(framesummary)
}

multiSpecOut <- multiSpp(sppListTemp)
multiSpecOut[,1]

hist(as.numeric(multiSpecOut[,1]), breaks=10)




lapply(sppListTemp, speciesSummary)

?lapply()




# CONNECT TO COMADRE/COMPADRE

RecordToIndex(4)

# From MOSAIC ID to a fetching of matrix population models

# We will need to modify the language around the compadre elements
# We will need to look @ a command to consider comadre v compadre

allMat <- function(RecordToIndex){
  allmatrices <- list()
  for(i in 1:length(match(RecordToIndex, cma$MatrixID))){
    allmatrices[[i]] <- cma$mat[[match(RecordToIndex, cma$MatrixID)[[i]]]]
    }
  return(allmatrices)
}

allMat(RecordToIndex(4))



MatA <- function(RecordToIndex){
  allmatrices <- list()
  for(i in 1:length(match(RecordToIndex, cma$MatrixID))){
    allmatrices[[i]] <- cma$mat[[match(RecordToIndex, cma$MatrixID)[[i]]]]@matA
  }
  return(allmatrices)
}

MatA(RecordToIndex(4))




indexToModel <- function(indexList){
  elements <- list()
  for(i in 1:length(indexList)){
    elements[[i]] <- paste("model", i)
  }
  elements[[length(elements)+1]] <- "all models"
  return(unlist(elements))
}

indexToModel(RecordToIndex(2))

listPull <- function(n, list){
  return(list[[n]])
}


# Package this function into its own function
# Need to fix the all models feature

listPull(menu(indexToModel(RecordToIndex(4))), MatA(RecordToIndex(4)))

link_comadre <- function(mosaicIndex){
  print(paste("Matrix Population Models available for:", mosaic@species[mosaicIndex]))
  print("More than one model in database - please make selection")
  return(listPull(menu(indexToModel(RecordToIndex(mosaicIndex))),
                  c(MatA(RecordToIndex(mosaicIndex)))))
}

link_comadre(2)
2





##########

mylist1 <- list(cpa$mat[[2]]@matA, cpa$mat[[4]]@matA)

match.arg(c(my1,my2), mylist1)
?match.arg



menu(c("model 1", "model 2"))

1


listPull <- function(n, list){
  return(mylist1[[n]])
}

listPull(menu(c("model 1", "model 2")))
indexToModel(RecordToIndex(10))


#########














#-----------------

switch(menu(c("model 1", "model 2")) + 1,
       cat("Nothing done\n"))

2


list1 <- list(cpa$mat[[2]]@matA, cpa$mat[[4]]@matA)


?switch()



RecordToIndex(2)
switch(menu(c(indexToModel(RecordToIndex(2)))) + 1,
       cat("Nothing done\n"),
       cpa$mat[[2]]@matA, cpa$mat[[4]]@matA)



RecordToIndex(4)


modelIndexing <- function(indexSeries){
  listofmodels <- list()
  for(i in 1:length(indexSeries)){
    listofmodels[[i]] <- paste(cpa$mat[[indexSeries[[i]]]]@matA)
    }
    return(listofmodels)
}


paste(cpa$mat[[match(RecordToIndex(4), cma$MatrixID)]]@matA)


cma$mat[[RecordToIndex(4)[[1]]]]
modelIndexing(RecordToIndex(4))
RecordToIndex(4)[[1]]
cpa$mat[[2]]@matA


cpa$MatrixID
match(RecordToIndex(4), cma$MatrixID)

cma$mat[[match(RecordToIndex(4), cma$MatrixID)[[1]]]]


# Set this up in a for loop where the indexing gets broken down with subpaste and subsequent re-matching

cue <- list()
for(k in 1:5){
  cue[[k]] <- try(paste("cma$mat[[match(RecordToIndex(4), cma$MatrixID)[[", k, "]]]]", sep=""))
}

unlist(cue)



RecordToIndex(4)



rownames(mosaic@biomass) <- rnorm()


# Checking whether there is a mosaic record for a given species (T/F)

spp_check <- function(binomialName){ # is any given species in the database
  return(length(which(mosaic@species==binomialName))>0)
}

spp_check("Fritillaria biflora")
spp_check("Ursus meritimus")










































