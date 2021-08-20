#---------------------------------------------------------------------------------------------#
# LOADDATA FUNCTION                                                                           #
#---------------------------------------------------------------------------------------------#

# Borrows Dr. Sun's loadData function to read in the BEDFiles and Count Data files.
# Source: isoform.0.99.1 R Library
#         Dr. Wei Sun (UNC)

loadData_djEdit <- function (countFile, bedFile, readLen, lmax = 500) {
  
  dat = countFile
  
  groupIDs = strsplit(dat$exons, split = ";", fixed = TRUE)
  
  splitFun <- function(vx) {
    s1 = strsplit(vx, split = "|", fixed = TRUE)
    if (any(sapply(s1, length) != 3)) {
      stop("expected three items\n")
    }
    matrix(unlist(s1), ncol = 3, byrow = TRUE)
  }
  
  groupIDs = lapply(groupIDs, splitFun)
  
  clustIDs = lapply(groupIDs, function(xm) {
    unique(xm[, 1])
  })
  clustIDs[1:2]
  
  geneIDs = lapply(groupIDs, function(xm) {
    unique(xm[, 2])
  })
  geneIDs[1:2]
  
  ncIDs = sapply(clustIDs, length)
  table(ncIDs)
  
  cID2rm = which(ncIDs > 1)
  
  ngIDs = sapply(geneIDs, length)
  table(ngIDs)
  w2check = which(ngIDs > 1)
  chkGIDs <- function(g1) {
    g2s = strsplit(g1, split = ":", fixed = TRUE)
    gus = unique(unlist(g2s))
    foundOne = FALSE
    for (gu1 in gus) {
      if (all(sapply(g2s, function(x, y) {
        y %in% x
      }, y = gu1))) {
        foundOne = TRUE
        break
      }
    }
    foundOne
  }
  gIDchk = sapply(geneIDs[w2check], chkGIDs)
  
  gID2rm = w2check[which(!gIDchk)]
  if (!all(cID2rm %in% gID2rm)) {
    stop("hm... cID2rm should be a subset of gID2rm :[\n")
  }
  str1 = "because they belong to different transcript clusters"
  message(length(cID2rm), " exon combinations are skipped ", 
          str1)
  str1 = "combinations are skipped because they belong to different genes"
  message("additional ", length(gID2rm) - length(cID2rm), " exon ", 
          str1)
  dim(dat)
  
  if (length(gID2rm) > 0) {
    dat = dat[-gID2rm, ]
    cID1 = unlist(clustIDs[-gID2rm])
    groupIDs = groupIDs[-gID2rm]
  }
  else {
    cID1 = unlist(clustIDs)
  }
  
  dim(dat)  
  
  # info = read.table(bedFile, sep = "\t", as.is = TRUE)
  # colNames = c("chr", "start", "end", "exon", "score", "strand")
  # cN = paste(colNames, collapse = ", ")
  # 
  # if (ncol(info) != 6) {
  #   stop(bedFile, " should have 6 columns: ", cN, "\n")
  # }
  # names(info) = colNames
  # dim(info)
  # info[1:2, ]
  
  info = bedFile
  
  infoIDs = strsplit(info$exon, split = "|", fixed = TRUE)
  if (any(sapply(infoIDs, length) != 3)) {
    stop("expected three items in infoIDs\n")
  }
  infoIDs = matrix(unlist(infoIDs), ncol = 3, byrow = TRUE)
  dim(infoIDs)
  infoIDs[1:2, ]
  
  infoN = info[, 1:3]
  infoN$gene = infoIDs[, 2]
  infoN$exon = as.numeric(infoIDs[, 3])
  dim(infoN)
  infoN[1:2, ]
  
  infoChr = split(infoN$chr, infoIDs[, 1])
  infoStart = split(infoN$start, infoIDs[, 1])
  infoEnd = split(infoN$end, infoIDs[, 1])
  infoGene = split(infoN$gene, infoIDs[, 1])
  infoExon = split(infoN$exon, infoIDs[, 1])
  
  exonIDs = lapply(groupIDs, function(Mx) {
    sort(as.numeric(Mx[, 3]))
  })
  length(exonIDs)
  exonIDs[1:5]
  
  dim(dat)
  dat[1:2, ]
  dat[["exons"]] = sapply(exonIDs, paste, collapse = ";")
  dim(dat)
  dat[1:2, ]
  cat("\n-----------------------------------------------------\n")
  cat("Extract information per gene")
  cat("\n-----------------------------------------------------\n")
  datCount = split(dat$count, cID1)
  datExons = split(dat$exons, cID1)
  geneMod = list()
  uIds = unique(cID1)
  for (k in 1:length(uIds)) {
    if (k%%2000 == 0) {
      cat(k, date(), "\n")
    }
    uId = uIds[k]
    gene1 = data.frame(count = datCount[[uId]], exons = datExons[[uId]], 
                       stringsAsFactors = FALSE)
    exonInfo = data.frame(chr = infoChr[[uId]], start = infoStart[[uId]], 
                          end = infoEnd[[uId]], gene = infoGene[[uId]], exon = infoExon[[uId]], 
                          stringsAsFactors = FALSE)
    geneMod[[uId]] = list(info = exonInfo, count = gene1)
  }
  
  geneMod
  
}