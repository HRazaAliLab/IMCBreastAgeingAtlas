## =============================================================================
## generalUtilities.R
##
## Purpose:
##   Centralized utility functions for the IMC Breast Ageing Atlas project.
##   Includes data accessors, preprocessing helpers, clustering wrappers,
##   image/mask utilities, plotting helpers, and statistical models.
##
## Usage:
##   source(here::here("code/generalUtilities.R"))
##   listGeneralUtilities()
##
## Notes:
##   - Functions are ordered from lightweight helpers → data access → analysis
##     → image processing → plotting.
##   - Heavy dependencies (EBImage, terra, ComplexHeatmap) are used only in
##     later sections to minimize sourcing failures.
## =============================================================================

suppressPackageStartupMessages({
  # Core deps that show up across functions
  library(data.table)
  library(here)
})

## =============================================================================
## Script introspection
## =============================================================================

#' List utility functions currently loaded
#'
#' By default this inspects the *calling* environment (parent.frame()).
listGeneralUtilities <- function(
    env = parent.frame(),
    pattern = NULL,
    include_internal = FALSE,
    print = TRUE
) {
  nms <- ls(envir = env, all.names = TRUE)
  is_fun <- vapply(nms, function(nm) {is.function(get(nm, envir = env, inherits = FALSE))}, logical(1))
  nms <- nms[is_fun]
  if (!include_internal) nms <- nms[!grepl("^\\.", nms)]
  if (!is.null(pattern)) nms <- nms[grepl(pattern, nms)]
  nms <- sort(nms)
  if (print) {
    cat(sprintf("Utility functions available (%d):\n", length(nms)))
    if (length(nms)) cat(paste0("  - ", nms, collapse = "\n"), "\n") else cat("  - \n")
    return(invisible(nms))
  }
  nms
}

## =============================================================================
## Basic helpers (pure utilities)
## =============================================================================

#' shorthand for as.data.table
adt <- function(x, ...) data.table::as.data.table(x, ...)

#' shorthand for summary
su <- function(x, ...) summary(x, ...)

#' shorthand for sort unique elements
sortunique <- function(x, ...) sort(unique(x,...))

#' Capitalise first letter; similar to Stata's -proper-
#' @param string
#' @return string in proper format
capStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  out <- paste(toupper(substring(c, 1,1)), substring(c, 2),
               sep="", collapse=" ")
  return(out)
}

#' Get matched substrings using regular expressions
getMatched <- function(pattern, vector) regmatches(vector, regexpr(pattern, vector))

#' Add zeros to make a character a set length (usually to standardise IDs)
padZeros <- function(number, nCharacters = 0L) {
  Zeros <- nCharacters - nchar(number)
  Zeros <- paste0(rep(0, Zeros), collapse = '')
  return(paste0(Zeros, number))
}

#' Convert scientific notation into neat LaTeX strings
mkEnumPower <- Vectorize(function(smallNum){
  smallNum <- paste0(smallNum)
  if(!grepl('e', smallNum)) return(paste0('$', smallNum, '$'))
  
  split <- strsplit(smallNum, split = 'e')
  mkTeXNum <- function(x) paste0("$", x[1], " \\times 10^{", x[2], "}$")
  out <- mkTeXNum(split[[1]])
  return(out)
})

## =============================================================================
## File system, IDs, and paths
## =============================================================================

#' Strips filepaths and fileendings to give ImageIDs of files within a folder
mkImgID <- function(filePath) gsub('_a0.*', '_a0', basename(filePath))

#' returns current cluster partition for passing to slurm requests
getPartition <- function(){
  partition <- gsub('\\..*', '', Sys.info()['nodename'])
  partition <- strsplit(partition, split = '-')[[1]][2]
  if(partition == 'node') return('general')
  return(partition)
}

#' Create a directory if it doesn't exist
#' @param path to new directory
mkDir <- function(path){
  if(dir.exists(path)) message(path, ' exists')
  else {
    dir.create(path, recursive = TRUE)
    message('Dir created')
  }
  return(fs::path(path))
}

#' Write a compressed parquet file (gzip, max compression).
writeOut <- function(x, file) {
  ifhasIndexRemove <- function(x){
    if(is.null(data.table::indices(x))) return(x)
    else{
      data.table::setindex(x, NULL)
      message('The data.table had a secondary index. This was removed to enable safe writing as a parquet file. It has not been restored.')
      return(x)
    }
  }
  arrow::write_parquet(ifhasIndexRemove(x), file, compression = 'gzip', compression_level = 9)
}

#' Write a temporary fst file from a data.table for later aggregation with aggregateTmpfst
writeTmpfst <- function(d, id, outdir){
  outFile <- tempfile(pattern = paste0('TmpFST', id), tmpdir = outdir, fileext = '.fst')
  fst::write_fst(d, outFile)
  stopifnot('file writing failed' = file.exists(outFile))
  return(outFile)
}

#' Aggregate the contents of (temporary) fst files and return aggregated data.table
aggregateTmpfst <- function(indir, delete = TRUE, ncores = 2L){
  stopifnot('indir does not exist' = dir.exists(indir))
  allFiles <- list.files(indir, pattern = '^TmpFST.*\\.fst$', full = TRUE)
  stopifnot('no relevant files fitting the \'^TmpFST.*\\.fst$\' pattern' = length(allFiles) > 0L)
  readTmpfst <- function(path) fst::read_fst(path, as.data.table = TRUE)
  allDat <- if(ncores > 1L) parallel::mclapply(allFiles, readTmpfst, mc.cores = ncores) else
    lapply(allFiles, readTmpfst)
  out <- data.table::rbindlist(allDat, fill = TRUE)
  if(delete){
    fs::file_delete(allFiles)
    if(length(list.files(indir)) == 0L) unlink(indir)
    else warning(paste(indir, 'contains additional files; directory not removed'))
  }
  return(out)
}

#' Write indexed dataset for fast slicing by ImageNumber
mkIndexedDataset <- function(dt, outDir = here('scratch')){
  require('fst')
  data.table::setkey(dt, 'ImageNumber')
  rowIndices <- data.table::copy(dt[, .(ImageNumber)])
  rowIndices[, index := seq_len(.N)]
  rowIndices[, min := min(index), by = ImageNumber]
  rowIndices[, max := max(index), by = ImageNumber]
  rowIndices <- rowIndices[, .SD[1], by = ImageNumber]
  TmpFile <- function(x) tempfile(pattern = x, tmpdir = outDir, fileext = '.fst')
  indexFile <- TmpFile('indexFile')
  cellFile <- TmpFile('cellsIndexed')
  fst::write_fst(rowIndices, indexFile)
  fst::write_fst(dt, cellFile)
  return(list(cellFile = cellFile, indexFile = indexFile))
}

#' Get data for a single ImageNumber using indices created by mkIndexedDataset()
getCellsImgNo <- function(imgNo, indexFiles){
  require('fst')
  rowIndices <- fst::read_fst(indexFiles$indexFile, as.data.table = TRUE)
  rowIndices <- rowIndices[ImageNumber == imgNo]
  dat <- fst::read_fst(indexFiles$cellFile, as.data.table = TRUE,
                       from = rowIndices[, min], to = rowIndices[, max])
  stopifnot(dat[, data.table::uniqueN(ImageNumber)] == 1L)
  stopifnot(dat[1L, ImageNumber] == imgNo)
  return(dat)
}

#' save a ggplot as pdf named ImageNumber#####.pdf
savePlotByImgNo <- function(plot, imgNo, h = 7, w = 7, outDir){
  imgNo <- padZeros(imgNo, nCharacters = 5)
  outfile <- fs::path(outDir, paste0('ImageNumber', imgNo, '.pdf'))
  ggplot2::ggsave(outfile, plot, height = h, width = w)
  return(outfile)
}

#' List files in a folder in a data.table with an ImageID column
fetchFileListWithImgID <- function(pathToFolder, pattern = '.tiff', recursive = FALSE){
  d <- data.table(files = list.files(pathToFolder, full = TRUE, pattern = pattern, recursive = recursive))
  d[, ImageID := mkImgID(files)]
  d <- d[grep('delete|ROI|orientation|check', invert = TRUE, x = ImageID)]
  return(d)
}

## =============================================================================
## Dataset accessors
## =============================================================================

#' Load a parquet data set as a data.table
read_parquet_adt <- function(file, ...) adt(arrow::read_parquet(file, ...))

#' Gets single-cell data
getCells <- function(){
  cells <- read_parquet_adt(here('data/derived/cells.parquet'))
  return(cells)
}

#' Gets clinical data
getClinical <- function(){
  clinical <- read_parquet_adt(here('data/derived/clinical.parquet'))
  clinical <- clinical[, c("ImageID", "PatientID", "Age", "TissueArea", "TotalArea", "ImageWidth", "ImageHeight")]
  return(clinical)
}

#' Gets cell context / histology context data
#' @param path Path to cellContext parquet
#' @return data.table
getCellContext <- function(path = here("data/derived/cellContext.parquet")){
  ctx <- read_parquet_adt(path)
  ctx <- data.table::as.data.table(ctx)
  data.table::setkey(ctx, ImageID, CellID)
  return(ctx)
}

#' Load cell phenotype annotation
getCellClusters <- function(annotationFile = here('data/derived/cellClusterAnnotation.csv')){
  stopifnot('The annotation file path is not valid' = file.exists(annotationFile))
  annotation <- data.table::fread(annotationFile)
  data.table::setkey(annotation, 'ClusterID')
  minimalColumns <- c('ClusterID', 'ClusterLabel', 'PrintOrder', 'Colour', 'isEpithelial')
  
  columnError <- paste0('Columns are missing from the annotation file.\n
    Required columns are: ', paste0(minimalColumns, collapse = '; '), ' > But the following are missing: ',
                        minimalColumns[!(minimalColumns %in% names(annotation))])
  if(!all(minimalColumns %in% names(annotation))) stop(columnError)
  
  annotation[, TeXClusterLabel := mkTeX(ClusterLabel, ignoreChar = '/')]
  annotation[, BackupTeXClusterLabel := mkTeX(BackupLabel, ignoreChar = '/')]
  
  epithelial <- annotation[isEpithelial == 1L]
  tme <- annotation[isEpithelial == 0L]
  
  Layers <- grep('isLayer[0-9]', names(annotation), value = TRUE)
  stopifnot('isLayer columns must have sequential numeric suffixes, starting at 1.'
            = seq_along(Layers) == as.integer(gsub('[A-Za-z]', '', Layers)))
  
  mkLayerVector <- function(n, d){
    layerVar <- paste0('isLayer', n)
    labD <- d[get(layerVar) == 1L][order(Group, PrintOrder)]
    outVec <- setNames(labD[, Colour], labD[, TeXClusterLabel])
    return(outVec)
  }
  
  mkNestedLayerVectors <- function(d){
    if(d[, .N] == 0L) return(NULL)
    out <- lapply(seq_along(Layers), mkLayerVector, d)
    names(out) <- Layers
    return(out)
  }
  
  out <- lapply(list(epithelial, tme), mkNestedLayerVectors)
  out <- append(out, list(annotation), after = 0)
  names(out) <- c('table', 'epi', 'tme')
  out <- stats::na.omit(out$table)
  return(out)
}

getPanel <- function(pathToPanel = here('data/raw/panel/panel.csv')){
  panel <- data.table::fread(pathToPanel)
  panel <- panel[Full == 1L]
  
  imctoolsReadOrder <- data.table::fread(here('data/raw/panel/metalReadOrder.csv'), header = FALSE)
  imctoolsReadOrder <- imctoolsReadOrder[[1]]  # take first (and only) column
  print(panel)
  print(imctoolsReadOrder)
  panelOut <- panel[, FormattedName]
  names(panelOut) <- panel[, Metal]
  
  stopifnot(
    'Your panel csv does not account for all of the metals in fullstack images.' =
      length(panelOut) == length(imctoolsReadOrder)
  )
  
  panelOut <- panelOut[imctoolsReadOrder]
  message('Panel in imctools read order')
  return(list(panel = panelOut))
}

getData <- function(type, curated = TRUE, justColNames = FALSE, selectCols = NULL){
  types <- list.files(here('data/derived'), pattern = '*.parquet$')
  types <- gsub('\\.parquet', '', types)
  types <- grep('Curated', types, invert = TRUE, value = TRUE)
  datName <- match.arg(type, choices = types)
  fileToRead <- getFilePath(datName, curated = curated)
  if(!is.null(justColNames) | !is.null(selectCols)){
    f <- arrow::ParquetFileReader$create(fileToRead)
    colNames <- f$GetSchema()$names
    if(justColNames) return(colNames)
    stopifnot(all(selectCols %in% colNames))
  }
  d <- if(!is.null(selectCols)) read_parquet_adt(fileToRead, col_select = selectCols) else read_parquet_adt(fileToRead)
  return(d)
}

printDatasetColumns <- function() {
  files <- list(
    "cells.parquet"                   = names(getCells()),
    "clinical.parquet"                = names(getClinical()),
    "cellContext.parquet"             = names(getCellContext()),
    "cellAdjacency.parquet"           = names(read_parquet_adt(here("data/derived/cellAdjacency.parquet"))),
    "compartmentSummary.parquet"      = names(read_parquet_adt(here("data/derived/compartmentSummary.parquet"))),
    "cellCompartmentDistances.parquet"= names(read_parquet_adt(here("data/derived/cellCompartmentDistances.parquet"))),
    "cellClusterAnnotation.csv"       = names(fread(here("data/derived/cellClusterAnnotation.csv"))),
    "spatialClusterAnnotation.csv"    = names(fread(here("data/derived/spatialClusterAnnotation.csv")))
  )
  for (nm in names(files)) {
    cat(sprintf("\n%s (%d columns)\n", nm, length(files[[nm]])))
    cat(paste(files[[nm]], collapse = ", "), "\n")
  }
}

## =============================================================================
## Numeric transforms and preprocessing
## =============================================================================

#' Rescale vector to 0 1
rescale01 <- function(x){
  cat('rescale01 is deprecated and will eventually be removed, use scales::rescale instead.\n')
  if(all(x == 0L)) return(x)
  x <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  return(x)
}

#' Arcsine transformation asin(sqrt(x))
arcsine <- function(x) asin(sqrt(x))

#' Scale (z-score) and clip
scale_clip <- function(x, limit = 2){
  x <- scale(x)
  x[x > limit] <- limit
  x[x < -limit] <- -limit
  return(x)
}

#' Clip vector at 99th centile if it's >0
clip99 <- function(vec, centile = 0.99){
  if(any(is.na(vec))) {
    print(paste("NA values found in the input vector of length:", length(vec)))
    na_indices <- which(is.na(vec))
    print(paste("Indices of NA values:", toString(na_indices)))
  }
  c99 <- quantile(vec, probs = centile)
  if(c99 > 0L) vec[vec > c99] <- c99
  return(vec)
}

## =============================================================================
## Cell counts, fractions, and densities
## =============================================================================

#' Make a data.table of cell counts and proportions by specified columns
getCellCounts <- function(d, clusterColumn = 'ClusterID',
                          idCols = 'PatientID', sepBy = 'isEpithelial') {
  
  mkCellCounts <- function(d, clusterColumn, idCols){
    allCols <- c(clusterColumn, idCols)
    d <- d[, .SD, .SDcols = allCols]
    d[, totalCells := .N, by = idCols]
    d[, nCellsPerType := .N, by = allCols]
    counts <- d[, .SD[1], by = allCols]
    
    uniqueElements <- lapply(allCols, function(x) counts[, unique(get(x))])
    allCombo <- adt(expand.grid(uniqueElements))
    setnames(allCombo, allCols)
    
    counts <- merge(x = allCombo, y = counts, by = allCols, all.x = TRUE)
    counts <- counts[order(totalCells)][, totalCells := totalCells[1], by = idCols]
    counts <- counts[!is.na(totalCells)]
    counts[is.na(nCellsPerType), nCellsPerType := 0L]
    counts[, proportion := nCellsPerType / totalCells]
    stopifnot(all.equal(counts[, sum(proportion), by = idCols][,V1],
                        rep(1L, counts[, uniqueN(.SD), .SDcols = idCols])))
    setkeyv(counts, allCols)
    return(counts)
  }
  
  return(d[, mkCellCounts(.SD, clusterColumn = clusterColumn, idCols = idCols), by = sepBy])
}

#' Make a data.table of fractions of cell phenotype positive for a specified marker
getPositiveFractions <- function(d,
                                 clusterColumn = 'ClusterID',
                                 idCols = 'PatientID',
                                 statusColumn = 'isKi67Pos',
                                 sepBy = 'isEpithelial'){
  
  posFrac <- getCellCounts(d, clusterColumn = c(clusterColumn, statusColumn),
                           idCols = idCols, sepBy = sepBy)
  setkeyv(posFrac, c(idCols, sepBy, clusterColumn, statusColumn))
  posFrac[, totalCellsPerType := sum(nCellsPerType), by = c(idCols, clusterColumn)]
  posFrac[, positiveFraction := nCellsPerType / totalCellsPerType]
  posFrac <- posFrac[get(statusColumn) == TRUE]
  return(posFrac)
}

#' Get a table of cell densities by a specified phenotype/cluster column
#'
#' @param d data.table/data.frame of cells
#' @param clusterColumn string column name giving phenotype/cluster labels to count
#' @param imgID string ImageID column name
#' @param ptID string PatientID column name
#' @param byPt if TRUE, aggregate counts and area per patient
#' @param sepBy optional string or character vector to stratify counts
#' @param areaCol area column in clinical to use (default: "TissueArea")
#'
#' @return data.table with Density, nCellsPerType, totalArea
getDensities <- function(d,
                         clusterColumn,
                         imgID = "ImageID",
                         ptID  = "PatientID",
                         byPt  = TRUE,
                         sepBy = "isEpithelial",
                         areaCol = "TissueArea") {
  
  if (missing(clusterColumn) || !is.character(clusterColumn) || length(clusterColumn) != 1) {
    stop("clusterColumn must be a single string column name.")
  }
  dt <- as.data.table(d)
  need_cols <- c(imgID, clusterColumn)
  if (!is.null(sepBy)) need_cols <- c(need_cols, sepBy)
  
  missing_cols <- setdiff(need_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `d`: ", paste(missing_cols, collapse = ", "))
  }
  
  clinical <- getClinical()
  
  if (!(imgID %in% names(clinical))) {stop("clinical is missing column: ", imgID)}
  if (!(areaCol %in% names(clinical))) {stop("clinical is missing area column: ", areaCol)}
  
  if (byPt && !(ptID %in% names(dt))) {
    if (!(ptID %in% names(clinical))) {stop("PatientID required in `d` or clinical when byPt = TRUE.")}
    dt <- merge(
      dt,
      clinical[, .(ImageID = get(imgID), PatientID = get(ptID))],
      by.x = imgID,
      by.y = "ImageID",
      all.x = TRUE
    )
  }

  area_dt <- clinical[, c(imgID, ptID, areaCol), with = FALSE]
  area_dt <- unique(area_dt, by = imgID)
  area_dt[, totalArea := get(areaCol)]
  
  idCol <- if (byPt) ptID else imgID
  
  if (byPt) {area_dt <- area_dt[, .(totalArea = sum(totalArea, na.rm = TRUE)),
                       by = idCol]
  } else {
    area_dt <- area_dt[, .(totalArea = totalArea[1]),
                       by = imgID]
  }

  count_by <- c(idCol, clusterColumn, sepBy)
  count_by <- count_by[!is.null(count_by)]
  
  counts <- dt[
    !is.na(get(clusterColumn)),
    .(nCellsPerType = .N),
    by = count_by
  ]
  
  out <- merge(counts, area_dt, by = idCol, all.x = TRUE)
  out[, Density := nCellsPerType / totalArea]
  
  keep <- c(idCol, clusterColumn, sepBy, "Density", "nCellsPerType", "totalArea")
  keep <- keep[keep %in% names(out)]
  setcolorder(out, keep)
  
  out[]
}


## =============================================================================
## Statistical models and summaries
## =============================================================================

format_custom_pval <- function(p_values) {
  sapply(p_values, function(p) {
    if (p < 1e-3) {
      sprintf("%.0e", p)
    } else if (p > 0.1) {
      format(p, scientific = FALSE, digits = 1, nsmall = 0)
    } else {
      format(p, scientific = FALSE, digits = 1, nsmall = 1)
    }
  })
}

#' Retrieve point estimate, 95% CI, pval and sample size from a model
grabEstimates <- function(modelOut, lineN = 2L, CIByProfiling = FALSE, CIBylme4 = FALSE){
  model_tidy <- adt(broom::tidy(modelOut))
  estimate <- model_tidy[lineN, .(estimate, p.value, std.error)]
  if(CIByProfiling) statsFun <- stats::confint
  else statsFun <- stats::confint.default
  if(CIBylme4) {
    statsFun <- lme4::confint.merMod
    CI <- adt(broom::confint_tidy(modelOut, func = statsFun, method = 'Wald')[lineN,])
  }
  else CI <- adt(broom::confint_tidy(modelOut, func = statsFun)[lineN,])
  modelout <- cbind(estimate, CI)
  modelout[, samplesize := nobs(modelOut)]
  return(modelout)
}

#' Get conventional hierarchical clustering order
getHC <- function(X, distance = "euclidean", clustering = "ward.D") {
  d <- dist(X, method = distance)
  hc <- hclust(d, method = clustering)
  return(hc)
}

#' Compare the proportions of clusters by a variable
doClusterComparison <- function(countDat, otherLinearResponse = NULL,
                                imgID = 'ImageID', ptID = 'PatientID', clusterCol = 'ClusterID', compareCol,
                                countColumn = 'nCellsPerType', totalColumn = 'totalCells',
                                addPredictors = NULL, lme4Control = NULL){
  
  stopifnot(data.table::is.data.table(countDat))
  d <- data.table::copy(countDat)
  
  stopifnot('compareCol is not a valid column in countDat' = compareCol %in% names(d))
  if(!is.null(otherLinearResponse)) stopifnot('otherLinearResponse is not a valid column in countDat' =
                                                otherLinearResponse %in% names(d))
  if(!is.null(imgID)){
    stopifnot('imgID and/or ptID are not valid column(s) in countDat' =
                c(imgID, ptID) %in% names(d))
  }else{
    stopifnot('ptID is not a valid column in countDat' = ptID %in% names(d))
  }
  
  stopifnot('countDat does not contain nCellsPerType and totalCells columns.
		Both are needed. countDat should be output of getCellCounts' = c(countColumn, totalColumn) %in% names(d))
  
  startN <- d[, .N]
  d <- d[!is.na(get(clusterCol))][!is.na(get(compareCol))]
  endN <- d[, .N]
  diff <- startN - endN
  if(startN != endN) message('Removed ', diff, ' NAs')
  predictors <- if(is.null(addPredictors)) compareCol else paste(compareCol, addPredictors, sep = ' + ')
  predictors <- paste(predictors, '+')
  
  doComparisonBylme <- function(cluster, d, imgID, ptID) {
    d[, indicator := get(clusterCol) == cluster]
    hasImgID <- if(any(names(d) == imgID)) paste('+ (1|', imgID, ')') else NULL
    
    if(is.null(otherLinearResponse)){
      lme4Control <- if(is.null(lme4Control)) lme4::glmerControl() else lme4Control
      fmla <- paste0(countColumn, '/', totalColumn, ' ~ ', predictors, '(1|', ptID, ')', hasImgID)
      lmeOut <- lme4::glmer(as.formula(fmla), family = binomial, weights = get(totalColumn),
                            data = d[(indicator)], control = lme4Control )
    }else{
      lme4Control <- if(is.null(lme4Control)) lme4::lmerControl() else lme4Control
      fmla <- paste0(otherLinearResponse, ' ~ ', predictors, '(1|', ptID, ')', hasImgID)
      lmeOut <- lme4::lmer(as.formula(fmla), data = d[(indicator)], control = lme4Control)
    }
    
    mklinfctMatrix <- function(fitted){
      m <- rep(0, nrow(coef(summary(fitted))))
      m[2] <- 1
      return(matrix(m, 1))
    }
    
    ht <- multcomp::glht(lmeOut, linfct = mklinfctMatrix(lmeOut))
    result <- summary(ht)
    tolog2OR <- function(x) log2(exp(x))
    effect <- if(is.null(otherLinearResponse)) tolog2OR(coef(result)) else coef(result)
    pV <- result$test$pvalues[1]
    ci  <- confint(ht)
    ci <- c(ci$confint[,"lwr"], ci$confint[,"upr"])
    ci <- if(is.null(otherLinearResponse)) sapply(ci, tolog2OR) else ci
    
    out <- data.table(ClusterID = cluster,
                      nPatients = d[(indicator), uniqueN(get(ptID))],
                      effect = effect, lci = ci[1], uci = ci[2], p = pV)
    if(is.null(otherLinearResponse)) setnames(out, 'effect', 'log2OR')
    d[, indicator := NULL]
    return(out)
  }
  
  clusters <- d[, unique(get(clusterCol))]
  result <- rbindlist(lapply(clusters, doComparisonBylme, d = d, imgID = imgID, ptID = ptID))
  result[, adjP := p.adjust(p, method = 'BH')]
  return(result)
}

## =============================================================================
## Clustering and embeddings
## =============================================================================

#' FlowSOM wrapper for SOM and (optionally) metaclusters
doFlowSOM <- function(
    dat,
    xdim = 10,
    ydim = 10,
    metacluster = FALSE,
    k = 30,
    seed = 357951
){
  data_FlowSOM <- flowCore::flowFrame(as.matrix(dat))
  out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out, xdim = xdim, ydim = ydim)
  som_nodes <- out$map$mapping[, 1]
  
  if(metacluster){
    out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)
    fsom_cluster <- out[som_nodes]
    FlowSOM_out <- list(
      som_nodes = som_nodes,
      fsom_clusters = fsom_cluster
    )
  } else{
    FlowSOM_out <- list(
      som_nodes = som_nodes
    )
  }
  return(FlowSOM_out)
}

#' Run EmbedSOM::SOM using all available threads for speed
doSOM <- function(dat, dimn, nThreads = parallel::detectCores(), epochs = 50, ...){
  som <- EmbedSOM::SOM(as.matrix(dat),
                       parallel = TRUE, threads = nThreads,
                       batch = TRUE, xdim = dimn, ydim = dimn, rlen = epochs, ...)
  return(som)
}

#' Runs phenograph using Rphenoannoy::Rphenoannoy
doPhenograph <- function(dat, k = 30, id){
  
  chk_no_ids_dat <- grep('^id$|^ID$|^som_nodes$|^apc_cluster$', names(dat))
  if(length(chk_no_ids_dat) != 0L) {
    stop('You have included an id variable with data to pass to phenograph')
  }
  
  pgclust <- Rphenoannoy::Rphenoannoy(as.matrix(dat), k = k)
  idsout <- igraph::membership(pgclust[[2]])
  idsout <- as.numeric(names(idsout))
  membership_out <- as.numeric(igraph::membership(pgclust[[2]]))
  idsout <- id[idsout]
  
  return(
    data.table(
      id = idsout,
      pg_membership = membership_out
    )
  )
}

#' Utilises Rphenoannoy to generate cell clusters
mkPhenoannoyClusters <- function(dat, channels, k = 30, seed){
  set.seed(seed)
  medians <- dat[, lapply(.SD, median), .SDcols = channels, by = somNodes]
  pgclust <- Rphenoannoy::Rphenoannoy(as.matrix(medians), k = k)$community$membership
  return(setNames(pgclust, medians$somNodes))
}

## =============================================================================
## Geometry and mask utilities
## =============================================================================

#' Compute area of a polygon based on its convex hull
convexHullArea <- function(xCoords, yCoords){
  coordDat <- as.matrix(cbind(xCoords, yCoords))
  ChullOrder <- grDevices::chull(x = xCoords, y = yCoords)
  polgyonOrder <- coordDat[ChullOrder,]
  polgyonOrder <- rbind(polgyonOrder, polgyonOrder[1,])
  polygonArea <- sp::Polygon(polgyonOrder, hole = FALSE)@area
  return(polygonArea)
}

#' Convert a set of polygon coordinates to a mask matrix
coords2Mask <- function(coords, idVar, xmax = 0L, ymax = 0L){
  stopifnot('coords should contain three columns in this order: x, y and idVar' = ncol(coords) == 3L)
  stopifnot('coords must be a data.table' = data.table::is.data.table(coords))
  polys <- split(coords, by = idVar, keep.by = FALSE)
  polys <- lapply(polys, as.matrix)
  v <- terra::vect(polys, type = 'polygons')
  r <- terra::rast(v, res = 1)
  ymx <- ifelse(nrow(r) > ymax, nrow(r), ymax)
  xmx <- ifelse(ncol(r) > xmax, ncol(r), xmax)
  r <- terra::extend(r, c(0, xmx, 0, ymx))
  rasterMat <- terra::rasterize(v, r, field = coords[, unique(get(idVar))])
  rasterMat <- terra::as.matrix(rasterMat, wide = TRUE)
  rasterMat[is.na(rasterMat)] <- 0
  return(rasterMat)
}

#' Get a data.table of oriented contours for plotting with geom_polygon
getContour <- function(pathToMask){
  stopifnot(file.exists(pathToMask))
  mkContour <- function(pathToMask){
    EBImage::ocontour(EBImage::rotate(EBImage::flip(tiff::readTIFF(pathToMask, as.is = TRUE)), angle = 90))
  }
  contours <- mkContour(pathToMask)
  mkContourDT <- function(no) {
    d <- adt(contours[[no]])
    data.table::setnames(d, names(d), c('x', 'y'))
    d[, CellID := as.integer(no)]
    d <- d[, x := x + 1L][, y := y + 1L]
    return(d)
  }
  imgContours <- rbindlist(lapply(names(contours), mkContourDT))
  return(imgContours)
}

#' Overlap two image masks
maskOverlap <- function(parentMask, targetMask, parallel = TRUE){
  chkMask <- function(x) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    return(isTRUE(is.matrix(x) & all(is.wholenumber(x))))
  }
  
  stopifnot('parentMask is not a matrix or contains non-integers' = chkMask(parentMask))
  stopifnot('targetMask is not a matrix or contains non-integers' = chkMask(targetMask))
  stopifnot('The masks are not of the same dimensions' = dim(parentMask) == dim(targetMask))
  
  allObjects <- data.table(ObjectNumber = as.integer(parentMask[parentMask != 0]))
  allObjects <- allObjects[, area := .N, by = ObjectNumber][, .SD[1], by = ObjectNumber]
  setkey(allObjects, 'ObjectNumber')
  
  mkRelatedToTable <- function(object){
    d <- data.table(relatedTo = targetMask[parentMask == object])[relatedTo != 0L]
    if(d[, .N] == 0L) return(NULL)
    d[, overlap := .N, by = relatedTo]
    d <- d[, .SD[1], by = relatedTo]
    d[, ObjectNumber := object]
    return(d)
  }
  
  out <- parallel::mclapply(allObjects[, ObjectNumber], mkRelatedToTable,
                            mc.cores = parallel::detectCores())
  if(!parallel) out <- lapply(allObjects[, ObjectNumber], mkRelatedToTable)
  out <- rbindlist(out)
  
  if(out[, .N == 0L]) {
    message('No objects are overlapping')
    return(NULL)
  }
  
  allObjects <- merge(x = allObjects, y = out, by = 'ObjectNumber')
  allObjects[, proportionOverlap := overlap / area]
  return(allObjects)
}

## =============================================================================
## Image processing and spatial analysis
## =============================================================================

getAllImagePaths <- function(dir = here("data/raw/images")) {
  files <- list.files(dir, pattern = "\\.tiff$", full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) {
    return(data.table(ImageID = character(), Type = character(), path = character()))
  }
  bn <- basename(files)
  
  img_id <- sub("(FullStack|CellMask|NuclearMask|NucleusMask|MyoepiMask|VesselMask|LymphaticMask)\\.tiff$",
                "", bn, ignore.case = TRUE)
  known_suffix <- grepl("(FullStack|CellMask|NuclearMask|NucleusMask|MyoepiMask|VesselMask|LymphaticMask)\\.tiff$",
                        bn, ignore.case = TRUE)
  img_id[!known_suffix] <- NA_character_
  
  type <- rep(NA_character_, length(bn))
  type[grepl("FullStack\\.tiff$",    bn, ignore.case = TRUE)] <- "FullStack"
  type[grepl("CellMask\\.tiff$",     bn, ignore.case = TRUE)] <- "Cell"
  type[grepl("(NuclearMask|NucleusMask)\\.tiff$", bn, ignore.case = TRUE)] <- "Nucleus"
  type[grepl("MyoepiMask\\.tiff$",   bn, ignore.case = TRUE)] <- "Myoepi"
  type[grepl("VesselMask\\.tiff$",   bn, ignore.case = TRUE)] <- "Vessel"
  type[grepl("LymphaticMask\\.tiff$",bn, ignore.case = TRUE)] <- "Lymphatic"
  images <- data.table(
    ImageID = img_id,
    Type   = type,
    path   = files
  )
  images <- images[!is.na(ImageID) & !is.na(Type)]
  dup <- images[, .N, by = .(ImageID, Type)][N > 1]
  if (nrow(dup) > 0) {
    stop(
      "Found multiple files for the same ImageID+Type:\n",
      paste(sprintf("  %s / %s (n=%d)", dup$ImageID, dup$Type, dup$N), collapse = "\n")
    )
  }
  
  setkeyv(images, c("ImageID", "Type"))
  images
}


getImagePaths <- function(imgID, dir = here("data/raw/images")) {
  files <- list.files(dir, pattern = paste0("^", imgID, ".*\\.tiff$"), full.names = TRUE)
  pick  <- function(p) files[grepl(p, basename(files), ignore.case = TRUE)]
  out <- list(FullStack = pick("FullStack"), Cell = pick("CellMask"), Nucleus = pick("NucleusMask"))
  bad <- names(out)[lengths(out) != 1]
  if (length(bad)) stop("Expected exactly 1 file for: ", paste(bad, collapse = ", "), " (", imgID, ")")
  out
}

#' Get IMC data for plotting
getRGBLong <- function(TiffPath, RGB, medianFilter = FALSE, size = 1L,
                       noRed = FALSE, noBlue = FALSE, noMarkers = FALSE) {
  stopifnot(file.exists(TiffPath))
  panel <- getPanel()
  imcData <- tiff::readTIFF(TiffPath, as.is = TRUE, all = TRUE)
  names(imcData) <- panel$panel
  
  ChannelsFound <- panel$panel[panel$panel %in% RGB]
  if(length(ChannelsFound) != 3L){
    stop(cat(ChannelsFound, 'are present in the formatted panel. There should be 3. \n
      Channels found were', ChannelsFound))
  }
  imcData <- imcData[RGB]
  
  if(medianFilter){
    doMedianFilter <- function(x) EBImage::medianFilter(EBImage::normalize(x), size = size)
    imcData <- parallel::mclapply(imcData, doMedianFilter, mc.cores = 3L)
  }
  
  mkLong <- function(dat){
    dat <- data.table(dat)
    longDat <- suppressWarnings(data.table::melt(dat))
    longDat[, x := as.integer(gsub('V', '', variable))][, variable := NULL]
    data.table::setnames(longDat, 'value', 'counts')
    data.table::setkey(longDat, 'x')
    longDat[, y := seq_len(.N), by = x]
    return(longDat)
  }
  
  imcToPlot <- lapply(imcData, mkLong)
  addMarkerCol <- function(marker){ imcToPlot[[marker]][, Marker := marker] }
  imcToPlot <- rbindlist(lapply(RGB, addMarkerCol))
  
  rescaleAndClip <- function(x){
    x <- scales::rescale(clip99(x))
    return(x)
  }
  
  imcToPlot[, scaledCounts := rescaleAndClip(counts), by = Marker]
  imcToPlot <- imcToPlot[, .(x, y, Marker, scaledCounts)]
  imcToPlot <- data.table::dcast(imcToPlot, x + y ~ Marker, value.var = 'scaledCounts')
  
  if (noRed) {
    remove <- RGB[1]
    imcToPlot[, (remove) := 0]
  }
  if (noBlue) {
    remove <- RGB[3]
    imcToPlot[, (remove) := 0]
  }
  if (noMarkers) {
    imcToPlot[, (RGB[1]) := 0]
    imcToPlot[, (RGB[2]) := 0]
    imcToPlot[, (RGB[3]) := 0]
  }
  
  data.table::setnames(imcToPlot, RGB, c('r', 'g', 'b'))
  imcToPlot[, rgb := grDevices::rgb(r, g, b)]
  return(list(data = imcToPlot[, .(x, y, rgb)], rgb = RGB))
}

#' Get a ggplot object of rgb raster image
mkRGBplot <- function(imgID, RGB = c('CD45', 'panCK', 'DNA1'),
                      noRed = FALSE, noBlue = FALSE, ...) {
  ImgPaths <- getImagePaths(imgID)
  fullStack <- ImgPaths$FullStack
  rgb <- getRGBLong(fullStack, RGB = RGB, noRed = noRed, noBlue = noBlue, ...)
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = rgb$data, ggplot2::aes(x = x, y = y, fill = rgb)) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
  return(list(plot = p, ImagePaths = ImgPaths))
}

#' Make a mask of 'holes' in tissue from an imc tiff
findHoles <- function(fullStackPath, kDilate = 5, kClose = 3,
                      minObjectSize = 100, maxObjectSize = NULL, removeCHull = FALSE,
                      mF = 0, markers = NULL){
  
  img <- tiff::readTIFF(fullStackPath, all = TRUE)
  N <- function(x) EBImage::normalize(x)
  
  clipNorm <- function(x) {
    x <- N(clip99(x))
    x[x > 0] <- mean(x)
    return(x)
  }
  
  if(!is.null(markers)) {
    names(img) <- suppressMessages(getPanel()$panel)
    img <- img[markers]
  }
  
  img <- lapply(img, clipNorm)
  sumImg <- Reduce('+', img)
  
  if(mF >=1){sumImg <- EBImage::medianFilter(sumImg, mF)}
  
  otsuMask <- function(x) N(x) > EBImage::otsu(N(x))
  
  sumImg <- otsuMask(sumImg)
  kBrush <- function(k) EBImage::makeBrush(k, shape='disc')
  sumImg <- EBImage::dilate(sumImg, kBrush(kDilate))
  sumImg <- EBImage::closing(sumImg, kBrush(kClose))
  objects <- EBImage::fillHull(EBImage::bwlabel(!sumImg))
  
  raster_chull_mask <- function(xy, r) {
    idchull <- grDevices::chull(xy)
    polygon <- sf::st_polygon(list(as.matrix(xy[c(idchull, idchull[1]), ])))
    polygon <- terra::vect(polygon)
    dummy <- terra::rasterize(polygon, r)
    dummy[is.na(dummy)] <- 0
    terra::crs(dummy) <- terra::crs(r)
    dummy
  }
  
  objectsPass <- objects
  
  if(removeCHull){
    mat <- which(sumImg == TRUE, arr.ind = TRUE)
    ch <- raster_chull_mask(xy = mat, terra::rast(sumImg))
    isTissue <- as.matrix(ch, wide = TRUE)
    objectsPass[isTissue==0] <- 0L
    objectsPass <- EBImage::reenumerate(objectsPass)
  }
  
  area <- EBImage::computeFeatures.shape(objectsPass)
  maxObjectSize <- if(!is.null(maxObjectSize)) maxObjectSize else max(area[, 's.area'])
  objectsToRemove <- which(area[, 's.area'] < minObjectSize | area[, 's.area'] > maxObjectSize)
  objectsPass <- EBImage::rmObjects(objectsPass, objectsToRemove, reenumerate = TRUE)
  
  return(objectsPass)
}

#' Run two class kmeans to determine which pixels are 'positive'
mkKmeansCounts <- function(imgID, marker, msize = 1, osize = 1,
                           addMarker = NULL, ignoreEdges = TRUE, returnImage = FALSE){
  set.seed(13)
  panel <- suppressMessages(getPanel()$panel)
  if(!is.null(addMarker)) stopifnot(all(c(marker, addMarker) %in% panel))
  else stopifnot(all(marker %in% panel))
  
  imgPaths <- getImagePaths(imgID)
  d <- tiff::readTIFF(imgPaths$fullStack, all = TRUE, as.is = TRUE)
  names(d) <- panel
  nmask <- tiff::readTIFF(imgPaths$Nucleus, as.is = TRUE)
  mask <- tiff::readTIFF(imgPaths$Cell, as.is = TRUE)
  
  mkLong <- function(channel){
    clipNorm  <- function(x) EBImage::normalize(clip99(x))
    img <- d[channel][[1]]
    img[nmask == 0L] <- 0L
    if(msize > 0L) img <- EBImage::medianFilter(img, msize)
    img <- clipNorm(img)
    long <- suppressWarnings(data.table::melt(adt(img)))
    long[, x := as.integer(gsub('^V', '', variable))]
    long[, y := seq_len(.N), x]
    setnames(long, 'value', channel)
    return(long[, .SD, .SDcols = c('x', 'y', channel)])
  }
  
  channels <- if(!is.null(addMarker)) c(marker, addMarker) else marker
  long <- lapply(channels, mkLong)
  long <- if(length(long) > 1L) Reduce(function(x, y) merge(x, y), long) else long[[1]]
  
  if(long[get(marker) == 0L, .N] == long[, .N]) return(NULL)
  km <- tryCatch(long[, kmeans(.SD, 2), .SDcols = channels],
                 error = function(e){NULL})
  if(is.null(km)){
    message('Kmeans failed. All marker values are probably 0.')
    return(km)
  }
  
  if(length(channels) == 1L) colnames(km$centers) <- marker
  posVal <- ifelse(km$centers[1, marker] > km$centers[2, marker], 1, 2)
  
  mkPositivePxImg <- function(posVal){
    long[, isPositive := km$cluster == posVal]
    long <- long[(isPositive)][, isPositive := as.integer(isPositive)][, .(x, y, isPositive)]
    counts <- Matrix::sparseMatrix(j = long$x, i = long$y, x = 1L, dims = dim(d[[1]]))
    img <- as.matrix(counts)
    return(img)
  }
  
  img <- mkPositivePxImg(posVal)
  
  findEdges <- function(pathToMask){
    edges <- getContour(pathToMask)
    mask <- tiff::readTIFF(pathToMask, as.is = TRUE)
    edges <- Matrix::sparseMatrix(j = edges$x, i = edges$y, x = 1L, dims = dim(mask))
    return(as.matrix(edges))
  }
  
  if(ignoreEdges){
    edges <- findEdges(imgPaths$Cell)
    edges[edges == 1L] <- -1
    img[edges == -1] <- 0L
  }
  
  newVar <- paste0(marker, 'KmeansAvg')
  
  mkKMImg <- function(img, osize){
    objects <- EBImage::bwlabel(img)
    if(all(objects == 0L)) {
      message('No segmented objects detected. Marker expression is likely very sparse.')
      return(NULL)
    }
    area <- EBImage::computeFeatures.shape(objects)
    area <- adt(area)[, ObjectNumber := .I][, .(s.area, ObjectNumber)]
    toRemove <- area[s.area <= osize][, ObjectNumber]
    objects <- EBImage::rmObjects(objects, toRemove)
    img[objects == 0L] <- 0L
    return(img)
  }
  
  measureObjects <- function(img, mask){
    counts <- EBImage::computeFeatures.basic(mask, img, basic.quantiles = c(0))
    counts <- adt(counts)
    counts[, ObjectNumber := .I]
    counts[, ImageID := imgID]
    setnames(counts, 'b.mean', newVar)
    return(counts)
  }
  
  kmimg <- mkKMImg(img = img, osize = osize)
  if(is.null(kmimg)) return(NULL)
  counts <- measureObjects(img = kmimg, mask = mask)
  counts <- counts[, .SD, .SDcols = c('ImageID', 'ObjectNumber', newVar)]
  if(returnImage) list(counts = counts, kmimg = kmimg)
  else return(list(counts = counts))
}

#' Measure object neighbours (cellprofiler-like)
measureObjectNeighbours <- function(pathToMask, dilateBy, includeDiagonal = TRUE, parallel = FALSE){
  require('data.table')
  stopifnot('Mask file does not exist; pass a full, valid path' = file.exists(pathToMask))
  Rcpp::cppFunction('IntegerMatrix getMN(int y, int x)
  {
    int P1y = y;     int P1x = x - 1;
    int P2y = y + 1; int P2x = x - 1;
    int P3y = y + 1; int P3x = x;
    int P4y = y + 1; int P4x = x + 1;
    int P5y = y;     int P5x = x + 1;
    int P6y = y - 1; int P6x = x + 1;
    int P7y = y - 1; int P7x = x;
    int P8y = y - 1; int P8x = x - 1;

    IntegerVector mooreneighbours  =
      {
       P1y,P2y,P3y,P4y,P5y,P6y,P7y,P8y,
       P1x,P2x,P3x,P4x,P5x,P6x,P7x,P8x
      };

    IntegerMatrix mn(8, 2, mooreneighbours.begin());
    colnames(mn) = CharacterVector::create("y", "x");
    return(mn);
  }')
  
  mask <- tiff::readTIFF(pathToMask, as.is = TRUE)
  if(dilateBy == 0L) dilated <- mask
  else{
    kern <- EBImage::makeBrush(dilateBy, shape = 'diamond')
    dilated  <- EBImage::dilate(mask, kern)
  }
  bounds <- dim(dilated)
  contours <- EBImage::ocontour(dilated)
  
  getNeighbours <- function(object, includeDiagonal = FALSE) {
    objectContour <- contours[[as.character(object)]] + 1L
    if(includeDiagonal) getMoore <- function(m) getMN(m[1], m[2])
    else getMoore <- function(m) getMN(m[1], m[2])[c(1, 3, 5, 7),]
    MNs <- apply(objectContour, 1, getMoore, simplify = FALSE)
    MNs <- do.call(rbind, MNs)
    
    yCol <- (MNs[, 1] <= bounds[1]) & (MNs[, 1] > 0L)
    xCol <- (MNs[, 2] <= bounds[2]) & (MNs[, 2] > 0L)
    MNs <- MNs[rowSums(cbind(yCol, xCol)) == 2L, ]
    
    Ns <- setdiff(unique(dilated[MNs]), c(0, object))
    if(length(Ns) == 0L) return(NULL)
    Ns <- data.table(from = object, to = Ns)
    return(Ns)
  }
  
  objects <- sort(unique(dilated[dilated != 0]))
  if(parallel) Neighbours <- rbindlist(parallel::mclapply(objects,
                                                          getNeighbours, includeDiagonal = includeDiagonal,
                                                          mc.cores = parallel::detectCores()))
  else Neighbours <- rbindlist(lapply(objects, getNeighbours, includeDiagonal = includeDiagonal))
  return(Neighbours)
}

#' Find nearest cell neighbour of a particular type
DistanceToNearestNeighbour <- function(imgDT, class){
  suppressPackageStartupMessages(require('spatstat'))
  
  getDistances <- function(imgD, Xcoord, Ycoord, class){
    xrange <- range(imgD[, get(Xcoord)])
    yrange <- range(imgD[, get(Ycoord)])
    boundary <- spatstat.geom::owin(xrange = xrange, yrange = yrange)
    imgPP <- spatstat.geom::ppp(imgD[, get(Xcoord)], imgD[, get(Ycoord)],
                                marks = imgD[, as.factor(get(class))], window = boundary)
    out <- adt(spatstat.geom::nnwhich(imgPP, by = spatstat.geom::marks(imgPP), k = 1))
    out[, eval(names(out)) := lapply(.SD, function(x) imgD$ObjectNumber[x])]
    outDist <- adt(spatstat.geom::nndist(imgPP, by = spatstat.geom::marks(imgPP), k = 1))
    setnames(outDist, paste0('DistTo', names(outDist)))
    out <- cbind(out, outDist, ObjectNumber = imgD$ObjectNumber)
  }
  
  distances <- getDistances(imgD = imgDT, Xcoord = 'Location_Center_X' ,
                            Ycoord = 'Location_Center_Y', class = class)
  return(distances)
}

## =============================================================================
## Plotting, themes, and visualization helpers
## =============================================================================

#' Utility function to add '$' to LaTeX expressions for latex2exp
mkTeX <- function(x, ignoreChar = NULL){
  x <- paste0('$', x, '$')
  notTeX <- function(s, x) gsub(s, paste0('$', s, '$'), x)
  dontParse <- c(' ', '&', '-', ignoreChar)
  for(i in paste0('\\', dontParse)) x <- notTeX(i, x)
  x <- gsub("(\\\\phi)\\$$", "\\1$ ", x, perl = TRUE)
  return(x)
}

#' Utility to parse latex expressions as axis labels
#' To use e.g. scale_y_discrete(labels = label_parse)
#' @import ggplot2 latex2exp
label_parse <- function(breaks) parse(text = latex2exp::TeX(breaks))

#' Convert ggplot into two objects - plot and legend
#' Thanks to https://github.com/taylor-lab/composite-mutations/blob/master/r/prerequisites.R
#' @param p a ggplot2 object
#' @return named list of 2 - plot and legend
#' @import cowplot ggplot2 grid gtable
extract_gglegend <- function(p){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) > 0) leg <- tmp$grobs[[leg]]
  else leg <- NULL
  
  legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
  plot <- p + ggplot2::theme(legend.position='none')
  list(plot=plot,legend=legend)
}

#' Returns key colours for analyses
getNormalBreastProjectColours <- function(){
  colsOut <- list(
    Compartment  = c('Epithelial' = '#6C7ACE', 'TME' = '#E8DDA6'),
    BroadPhenotype = setNames(c('#6C7ACE', '#7CA876', '#EFAA9C'), c('Epithelial', 'Immune', 'Stromal')),
    DuctLobule = c("duct" = "#4169E1", "lobule" = "#DC143C"),
    AgeGroup = c("lessthan50" = "gray", "greaterthan50" = "dimgray")
  )
  return(colsOut)
}