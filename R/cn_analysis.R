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
