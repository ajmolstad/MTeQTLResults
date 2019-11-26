getMTeQTLResults <- function(filepath, gene = NULL, chr = NULL, pos = NULL, extract = NULL){
  
  # -----------------------------------------
  # Access manifest
  # -----------------------------------------
  manifest <- readRDS(paste(filepath,"/MTeQTLResults/Data/Manifest.RDS", sep=""))
  
  if(is.null(gene) & is.null(chr) | is.null(gene) & is.null(pos)){
    stop("The user must supply either a gene list or a list of chromosomes/positions!")
  }
  outList <- NULL

  if(!is.null(gene)){
    
    genes <- unlist(gene)
    getGenes <- match(genes, manifest[,1])
    manifest_trunc <- manifest[getGenes,,drop=FALSE]
    if(dim(manifest_trunc)[1] > 1){
      getChr <- unique(manifest_trunc[,2])
    } else {
      getChr <- manifest_trunc[1,2]
    }
    getChr <- as.numeric(as.character(getChr))
    for(j in 1:length(getChr)){
      temp <- readRDS(paste(filepath,"/MTeQTLResults/Data/Chr", getChr[j],"_MTeQTLResults.RDS", sep=""))
      genesChr <- manifest_trunc[which(manifest_trunc[,2] == as.character(getChr[j])),1]
      keepInds <- match(genesChr, unlist(lapply(temp, "[[", 1)))
      outListtemp <- NULL
      for(k in 1:length(keepInds)){
        outListtemp[[k]] <- temp[[keepInds[k]]]
        outListtemp[[k]]$chr <- getChr[j]
        outListtemp[[k]]$TSS <- as.numeric(as.character(manifest[match( outListtemp[[k]]$gene, manifest[,1]),3]))
        outListtemp[[k]]$TES <- as.numeric(as.character(manifest[match( outListtemp[[k]]$gene, manifest[,1]),4]))
      }
      outList <- append(outList, outListtemp)
    }
    
    if(!any(extract == "TSS")){
      for(j in 1:length(outList)){
        outList[[j]]$TSS <- NULL
      }
    }
    if(!any(extract == "TES")){
      for(j in 1:length(outList)){
        outList[[j]]$TES <- NULL
      }
    }
    if(!any(extract == "chr")){
      for(j in 1:length(outList)){
        outList[[j]]$chr <- NULL
      }
    }
    if(!any(extract == "Omega")){
      for(j in 1:length(outList)){
        outList[[j]]$Omega <- NULL
      }
    }
    if(!any(extract == "R2")){
      for(j in 1:length(outList)){
        outList[[j]]$R2 <- NULL
      }
    }
    if(!any(extract == "beta") & !is.null(extract)){
      for(j in 1:length(outList)){
        outList[[j]]$beta <- NULL
      }
    }
    
    
    
  } else {
    if(!is.null(chr) & !is.null(pos)){
      
      chrUnlist <- unlist(chr)
      getChr <- as.numeric(as.character(unique(chrUnlist)))
      outListtemp <- NULL
      
      for(j in 1:length(getChr)){
        
        temp <- readRDS(paste(filepath,"/MTeQTLResults/Data/Chr", getChr[j],"_MTeQTLResults.RDS", sep=""))
        getPos <- unlist(pos)[which(chr == getChr[j])]
        outListtemp <- NULL
        
        for(k in 1:length(getPos)){
          getInds <- which((manifest[,2] == as.character(getChr[j])) & (as.numeric(as.character(manifest[,3])) - 1e-6 <= getPos[k]) & (as.numeric(as.character(manifest[,4])) + 1e-6 >= getPos[k]))  
          if(length(getInds) != 0){
            geneName <- manifest[getInds,1]
            keepInds <- match(geneName, lapply(temp,"[[", 1))
            outListtemp[[k]] <- temp[[keepInds[k]]]
            outListtemp[[k]]$chr <- getChr[j]
            outListtemp[[k]]$TSS <- as.numeric(as.character(manifest[match(outListtemp[[k]]$gene, manifest[,1]),3]))
            outListtemp[[k]]$TES <- as.numeric(as.character(manifest[match(outListtemp[[k]]$gene, manifest[,1]),4]))
          }
        }
        
        outList <- append(outList, outListtemp)
            
      }
      
      if(!any(extract == "TSS")){
        for(j in 1:length(outList)){
          outList[[j]]$TSS <- NULL
        }
      }
      if(!any(extract == "TES")){
        for(j in 1:length(outList)){
          outList[[j]]$TES <- NULL
        }
      }
      if(!any(extract == "chr")){
        for(j in 1:length(outList)){
          outList[[j]]$chr <- NULL
        }
      }
      if(!any(extract == "Omega")){
        for(j in 1:length(outList)){
          outList[[j]]$Omega <- NULL
        }
      }
      if(!any(extract == "R2")){
        for(j in 1:length(outList)){
          outList[[j]]$R2 <- NULL
        }
      }
      if(!any(extract == "beta") & !is.null(extract)){
        for(j in 1:length(outList)){
          outList[[j]]$beta <- NULL
        }
      }
    }
  }
    
    return(outList)
}