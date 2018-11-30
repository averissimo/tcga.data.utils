#' Create MultiAssayExperiment object from data
#'
#' @param clinical use custom clinical (that can be pre-processed)
#'
#' @return a MultiAssayExperiment object
#' @export
build.assay <- function(clinical.custom, 
                        gdc.custom, 
                        mutation.custom) {
  
  #
  # Expression data (FPKM)
  
  # map expression data with clinical
  es.map.fpkm <- data.frame(primary = strtrim(colnames(gdc.custom$rnaseq.fpkm), 12), 
                            colname  = colnames(gdc.custom$rnaseq.fpkm), 
                            stringsAsFactors = FALSE)
  
  # filter only valid date.. i.e expression that have clinical data
  valid.ix.fpkm  <- es.map.fpkm$primary %in%  clinical.custom$bcr_patient_barcode
  valid.dat.fpkm <- gdc.custom$rnaseq.fpkm[, valid.ix.fpkm]@assays@.xData$data[['HTSeq - FPKM']]
  colnames(valid.dat.fpkm) <- colnames(gdc.custom$rnaseq.fpkm[, valid.ix.fpkm])
  rownames(valid.dat.fpkm) <- rownames(gdc.custom$rnaseq.fpkm[, valid.ix.fpkm])
  
  #sample.barcode.fpkm <- strtrim(colnames(valid.dat.fpkm), 16)
  #valid.codes.fpkm  <- sample.barcode.fpkm[sample.barcode.fpkm %in% gdc.custom$bio.sample$bcr_sample_barcode]
  
  #fpkm.df       <- Biobase::AnnotatedDataFrame(gdc.custom$bio.sample[valid.codes.fpkm,])
  fpkm.df <- Biobase::AnnotatedDataFrame(gdc.custom$rnaseq.fpkm[, valid.ix.fpkm]@colData %>% as.data.frame)
  #rownames(fpkm.df) <- colnames(valid.dat.fpkm)
  
  # build expression set
  es.fpkm <- Biobase::ExpressionSet(assayData = valid.dat.fpkm, phenoData = fpkm.df)
  
  #
  # Expression data (Counts)
  
  # map expression data with clinical
  es.map.counts <- data.frame(primary = strtrim(colnames(gdc.custom$rnaseq.counts), 12), 
                              colname = colnames(gdc.custom$rnaseq.counts), 
                              stringsAsFactors = FALSE)
  
  # filter only valid date.. i.e expression that have clinical data
  valid.ix.counts  <- es.map.counts$primary %in%  clinical.custom$bcr_patient_barcode
  valid.dat.counts <- gdc.custom$rnaseq.counts[, valid.ix.counts]@assays@.xData$data[['HTSeq - Counts']]
  colnames(valid.dat.counts) <- colnames(gdc.custom$rnaseq.counts[, valid.ix.counts])
  rownames(valid.dat.counts) <- rownames(gdc.custom$rnaseq.counts[, valid.ix.counts])
  
  
  #sample.barcode.counts <- strtrim(colnames(valid.dat.counts), 16)
  #valid.codes.counts  <- sample.barcode.counts[sample.barcode.counts %in% gdc.custom$bio.sample$bcr_sample_barcode]
  
  counts.df <- Biobase::AnnotatedDataFrame(gdc.custom$rnaseq.counts[, valid.ix.counts]@colData %>% as.data.frame)
  
  # build expression set
  es.counts  <- Biobase::ExpressionSet(assayData = valid.dat.counts, phenoData = counts.df)
  
  #
  # Mutation data
  mutation.colnames <- colnames(mutation.custom)  
  valid.ix          <- colnames(mutation.custom) %in% clinical$bcr_patient_barcode
  
  mut.map <- data.frame(primary = mutation.colnames[valid.ix], colname = mutation.colnames[valid.ix])
  mut     <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mutation.custom))
  
  #
  # Setup to create MultiAssayExperiment object
  
  futile.logger::flog.info('Building Assay...')
  listmap        <- list(es.map.fpkm, es.map.counts, mut.map)
  names(listmap) <- c("RNASeqFPKM", "RNASeqCounts", "Mutation")
  
  dfmap    <- MultiAssayExperiment::listToMap(listmap)
  objlist  <- list("RNASeqFPKM" = es.fpkm, "RNASeqCounts" = es.fpkm, "Mutation" = mut)
  my.assay <- MultiAssayExperiment::MultiAssayExperiment(objlist, clinical.custom, dfmap)
  
  return(my.assay)
}