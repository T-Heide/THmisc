#' Creates a data.frame of named tiles from a BSgenome object
#'
#' @param geno A BSgenome object
#' @param chr_filter Optional vector of chromosomes to subset to.
#' @param width The width of the tiles (Default: 1e6)
#'
#' @return A data.frame containing the positions of the tiles.
#' @export
window_list_from_genome = function(geno, chr_filter=NULL, width=1e6) {
  
  checkmate::assertClass(geno, "BSgenome")
  checkmate::assertCharacter(chr_filter, null.ok = TRUE)
  
  res = GenomicRanges::GRanges(seqnames(geno), IRanges::IRanges(0, seqlengths(geno))) %>%
    tile(width=width) %>% 
    lapply(data.frame) %>%
    do.call(what=rbind) %>% 
    data.frame() %>% 
    dplyr::mutate(name=paste0("win", seq_along(seqnames))) %>% 
    dplyr::select(seqnames, start, end, name) %>% 
    magrittr::set_colnames(c("chr","start","end","name"))
  
  if (!is.null(chr_filter)) {
    res = res %>% dplyr::filter(chr %in% chr_filter)
  }
  
  return(res)
}
#' Function to look-up (average) gene positions.
#'
#' @param x 
#' @param txdb 
#' @param egdb 
#' @param return_intervals 
#'
#' @return
#' @export
#'
#' @examples
get_gene_pos = function(x, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, egdb=org.Hs.eg.db::org.Hs.eg.db, return_intervals=FALSE) {
  
  # check inputs
  ec = checkmate::makeAssertCollection()
  checkmate::assertCharacter(x, any.missing = FALSE, null.ok=TRUE, add = ec)
  checkmate::assertClass(txdb, "TxDb", null.ok = TRUE, add = ec)
  checkmate::assertClass(egdb, "OrgDb", null.ok = TRUE, add = ec)
  checkmate::assertFlag(return_intervals, add = ec)
  checkmate::reportAssertions(ec)
  
  ##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  gene_to_id = 
    suppressMessages(
      AnnotationDbi::select(
        egdb,
        x,
        columns = c("SYMBOL", "ENTREZID"),
        keytype = "SYMBOL"
      ))
  
  pos_data = 
    suppressMessages(
      AnnotationDbi::select(
        txdb,
        gene_to_id$ENTREZID,
        columns = c("CDSCHROM", "CDSSTART", "CDSEND", "GENEID"),
        "GENEID"
      ) %>% split(.$GENEID)
    )
  
  pos_data = pos_data[gene_to_id$ENTREZID]
  
  if (return_intervals) {
    
    res = pos_data
    names(res) = gene_to_id$SYMBOL
    
  } else {
    
    pos_data_avg = 
      data.frame(
        ENTREZID=names(pos_data),
        chr=gsub("chr", "", sapply(pos_data, function(x) unique(x$CDSCHROM))),
        pos=sapply(pos_data, function(x) mean((x$CDSSTART + x$CDSEND)/2)), 
        row.names = NULL
      )
    
    res = merge(gene_to_id, pos_data_avg, by="ENTREZID", all=TRUE)
    
  }
  
  return(res)
}
