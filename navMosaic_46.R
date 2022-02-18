#' #                                                                                                        -----
#'*   PROJECT NAME    [Fetch_ATP]
#'*   AUTHOR          [Connor Bernard. connor.bernard@zoo.ox.ac.uk.]
#'*   LAST UPDATE     [14 Feb 2022]
#'*   PURPOSE         [Downloading and Coercing the data in r]
#'*   DEPENDENCIES    [NA]
#' #                                                                                                        -----
#-------------------------------
#---------------------

IndexToRecord <- function(indexValue){#Index-to-Index 
  return(which(unlist(lapply(mosaic@index, equivalency_A, indexValue))))
}

RecordToIndex <- function(mosaicIndex){#Index-to-Index 
  return(mosaic@index[[mosaicIndex]])
}

# Simple trait search - trait All Records
traitAllRecords <- function(trait){
  return(slot(slot(mosaic, trait),"value"))
}


# For pulling all metadata for a specific trait-individual record
metaRecords <- function(trait, index){
  metas <- data.frame(matrix(NA, nrow=1, ncol=5))
  colnames(metas) <- c("author", "year", "journal", "database", "mosaic")
  metas[1,1] <- slot(slot(mosaic, trait), "author")[[index]]
  metas[1,2] <- slot(slot(mosaic, trait), "year")[[index]]
  metas[1,3] <- slot(slot(mosaic, trait), "journal")[[index]]
  metas[1,4] <- slot(slot(mosaic, trait), "database")[[index]]
  metas[1,5] <- slot(slot(mosaic, trait), "mosaic")[[index]]
  return(metas)
}


# Species Trait Summary
speciesSummary  <- function(speciesNameLatin){
  speciesIndex <- which(mosaic@species==speciesNameLatin)
  metadataNames <- slotNames(mosaic)[-c(1:3)]
  meta <- data.frame(matrix(NA, nrow=1, ncol=length(metadataNames)))
  colnames(meta) <- metadataNames
  for(i in 1:length(metadataNames)){
    meta[,i] <- slot(slot(mosaic, metadataNames[[i]]), "value")[[speciesIndex]]
  }
  return(meta)
}


# Is a species included in MOSAIC - T/F

spp_check <- function(binomialName){ # is any given species in the database
  return(length(which(mosaic@species==binomialName))>0)
}


# Summary of multiple species

multiSummary <- function(sppList){
  summary <- lapply(sppList, speciesSummary)
  framesummary <- data.frame(matrix(NA, nrow=length(summary), ncol=length(summary[[1]])+1))
  colnames(framesummary) <- c("sppnames", colnames(summary[[1]]))
  for(i in 1:length(summary)){
    framesummary[i,] <- c(sppList[[i]],summary[[i]])
  }
  return(framesummary)
}

# From MOSAIC ID to a fetching of matrix population models

allMat <- function(RecordToIndex){
  allmatrices <- list()
  for(i in 1:length(match(RecordToIndex, cma$MatrixID))){
    allmatrices[[i]] <- cma$mat[[match(RecordToIndex, cma$MatrixID)[[i]]]]
    }
  return(allmatrices)
}

allMat(RecordToIndex(4))

# Spp count by trait

recordCount <- function(trait, numeric=FALSE){
  allTrait <- traitAllRecords(trait)
  if(numeric==FALSE){
    allTrait <- allTrait[!is.na(allTrait)]
    allTrait <- allTrait[!allTrait=="NDY"]
    b <- length(allTrait)
  }else if(numeric==TRUE){
    suppressWarnings(allTrait <- as.numeric(allTrait))
    allTrait <- allTrait[!is.na(allTrait)]
    b <- length(allTrait)
  }
  return(b)
}


# Proportions

TraitFrequency <- function(trait){
  focalTrait <- slot(slot(mosaic, trait), "value")
  focalTrait <- as.factor(focalTrait)
  hold <- matrix(NA,
                 nrow=length(levels(focalTrait)),
                 ncol=2)
  
  rownames(hold) <- levels(focalTrait)
  colnames(hold) <- c("counts", "freq")
  
  for(i in 1:length(levels(focalTrait))){
    hold[i,1] <- sum(focalTrait==levels(focalTrait)[[i]])
  }
  
  hold[which(rownames(hold)=="NDY"),2] <- NA
  hold[-which(rownames(hold)=="NDY"),2] <- round(hold[-which(rownames(hold)=="NDY"),1]/sum(hold[-which(rownames(hold)=="NDY"),1]), 3)
  return(hold)
}

cue <- list()
for(k in 1:5){
  cue[[k]] <- try(paste("cma$mat[[match(RecordToIndex(4), cma$MatrixID)[[", k, "]]]]", sep=""))
}

# Metadata
metaRecords <- function(trait, index){
  metas <- data.frame(matrix(NA, nrow=1, ncol=5))
  colnames(metas) <- c("author", "year", "journal", "database", "mosaic")
  metas[1,1] <- slot(slot(mosaic, trait), "author")[[index]]
  metas[1,2] <- slot(slot(mosaic, trait), "year")[[index]]
  metas[1,3] <- slot(slot(mosaic, trait), "journal")[[index]]
  metas[1,4] <- slot(slot(mosaic, trait), "database")[[index]]
  metas[1,5] <- slot(slot(mosaic, trait), "mosaic")[[index]]
  return(metas)
}

multiMetaRecords <- function(trait, sequenceOfIds){
  formerMeta <- metaRecords("volancy", sequenceOfIds[[1]])
  for(i in 2:length(sequenceOfIds)){
    formerMeta <- rbind(metaRecords("volancy", sequenceOfIds[[i]]), formerMeta)
  }
  return(formerMeta)
}




# Package this function into its own function

listPull <- function(n, list){
  return(list[[n]])
}

equivalency_A <- function(x, y){
  any(x==y)
}

is.integer0 <- function(x){is.integer(x) && length(x) == 0L}

simpleExtract <- function(n){
  if(is.integer0(intersect(cma@data$MatrixID, mosaic@index[[n]][[1]]))==FALSE){
    durat <- cma@data$StudyDuration[match(mosaic@index[[n]], cma@data$MatrixID)]
    aut <- cma@data$Authors[match(mosaic@index[[n]], cma@data$MatrixID)]
    year <- cma@data$YearPublication[match(mosaic@index[[n]], cma@data$MatrixID)]
    dim <- cma@data$MatrixDimension[match(mosaic@index[[n]], cma@data$MatrixID)]
  }else if(is.integer0(intersect(cpa@data$MatrixID, mosaic@index[[n]][[1]]))==FALSE){
    durat <- cpa@data$StudyDuration[match(mosaic@index[[n]], cpa@data$MatrixID)]
    aut <- cpa@data$Authors[match(mosaic@index[[n]], cpa@data$MatrixID)]
    year <- cpa@data$YearPublication[match(mosaic@index[[n]], cpa@data$MatrixID)]
    dim <- cpa@data$MatrixDimension[match(mosaic@index[[n]], cpa@data$MatrixID)]
  }
  modelMeta <- list()
  for(i in 1:length(durat)){
    d1 <- NA
    d2 <- NA
    d3 <- NA
    d4 <- NA
    d1 <- try(durat[[i]])
    d2 <- try(gsub("\\;.*","", aut[[i]]))
    d3 <- try(year[[i]])
    d4 <- try(dim[[i]])
    modelMeta[[i]] <- paste(d1, "yr", "_", d2,"_", d3, "_", d4,"d", sep="")
  }
  return(modelMeta)
}


# if then commant re: integration

indexToModel <- function(indexList, n){
  elements <- list()
  for(i in 1:length(indexList)){
    elements[[i]] <- paste("model", i)
  }
  for(i in 1:length(indexList)){
    elements[[i]] <- simpleExtract(n)[[i]]
  }
  elements[[length(elements)+1]] <- "all models"
  return(unlist(elements))
}


# Mat As only
MatA_ca <- function(RecordToIndex){
  allmatrices <- list()
  for(i in 1:length(match(RecordToIndex, cma$MatrixID))){
    allmatrices[[i]] <- cma$mat[[match(RecordToIndex, cma$MatrixID)[[i]]]]@matA
  }
  return(allmatrices)
}

MatA_cp <- function(RecordToIndex){
  allmatrices <- list()
  for(i in 1:length(match(RecordToIndex, cpa$MatrixID))){
    allmatrices[[i]] <- cpa$mat[[match(RecordToIndex, cpa$MatrixID)[[i]]]]@matA
  }
  return(allmatrices)
}

MatA <- function(Index){
  if(is.integer0(intersect(cpa@data$MatrixID, RecordToIndex(Index)))==FALSE){
    b <- MatA_cp(RecordToIndex(Index))
  }else if(is.integer0(intersect(cma@data$MatrixID, RecordToIndex(Index)))==FALSE){
    b <- MatA_ca(RecordToIndex(Index))
  }
  return(b)
}

link_comadre <- function(mosaicIndex){
  print(paste("Matrix Population Models available for:", mosaic@species[mosaicIndex]))
  print("More than one model in database - please make selection")
  return(listPull(menu(indexToModel(RecordToIndex(mosaicIndex), mosaicIndex)),
                  c(MatA(mosaicIndex))))
}
