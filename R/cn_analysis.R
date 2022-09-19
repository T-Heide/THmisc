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

  res = GenomicRanges::GRanges(GenomeInfoDb::seqnames(geno), IRanges::IRanges(0, GenomeInfoDb::seqlengths(geno))) %>%
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

#' Calculate average copy-number values in bins.
#'
#' @param cn_data A list of data frames or a single data frame that contains interval definitions (start/start.pos, end/end.pos, chr/chromosome columns) and copy-number values (CN/CNt).
#' @param windows A data.frame that can be converted to a GRanges object defining the window intervals.
#'
#' @return A data.frame containing all columns of the 'windows' argument and average copy-number values of all entries in 'cn_data'.
#' @export
#'
calc_avg_cn_per_bin = function(cn_data, windows) {

  .try_fix_input = function(x) {
    alt_enc = c(start="start.pos", end="end.pos", chr="chromosome", CN="CNt")
    for (i in seq_along(alt_enc)) {
      cn = colnames(x)
      if (!names(alt_enc)[i] %in% cn & alt_enc[i] %in% cn) {
        x[,names(alt_enc)[i]] = x[,alt_enc[i]]
      }
    }

    checkmate::assertDataFrame(x)
    checkmate::checkSubset(c("chr","start","end","CN"), colnames(x))

    return(x)
  }

  .calc_avg_cn_in_bin = function(x, y) {
    dg = as(x[, c("chr","start","end","CN")], "GRanges")
    cb = dplyr::mutate(y, value=NA)
    pairs = findOverlaps(as(cb, "GRanges"), dg , minoverlap=1)
    average_value = tapply(dg$CN[to(pairs)], from(pairs), mean, na.rm=TRUE)
    cb$value[as.numeric(names(average_value))] = average_value
    return(cb$value)
  }

  checkmate::assertTRUE(is.data.frame(cn_data) | is.list(cn_data))
  checkmate::assertDataFrame(windows, null.ok = FALSE)

  if (is.data.frame(cn_data)) {
    checkmate::assertSubset("sample_barcode", colnames(cn_data))
    cn_data = split(cn_data, cn_data$sample_barcode)
  }

  cn_data = lapply(cn_data, .try_fix_input)
  checkmate::assertNamed(cn_data)
  checkmate::assertList(cn_data)

  cn_list = lapply(cn_data, .calc_avg_cn_in_bin, y=windows)
  cn_mat = cbind(windows, do.call(what=cbind, cn_list))

  return(cn_mat)
}

findAmplificationBins = function(dpw, m_cn, amp_factor = 2, max_amp_size = 10*10^6) {

  # Expects ordered by sample, then chromosome, then position
  stopifnot(!any(duplicated(rle(dpw$sample_barcode)$values)))
  stopifnot(sum(diff(as.numeric(gsub("win", "", dpw$name))) < 0) == length(unique(dpw$sample_barcode))-1)

  # drop na and determine amplified bins
  dpw = dpw[!is.na(dpw$CN),]
  crit_val = m_cn[dpw$sample] # sample ids are/have to be unique
  dpw$amped = dpw$CN > amp_factor * crit_val

  ids = paste0(dpw$sample, ".", dpw$chr)
  d_per_chr_and_bc = split(dpw, ids)

  amplified_bins =
    lapply(d_per_chr_and_bc, function(d) {
      # change to rle
      df_rle = rle(d$amped)
      df_rle$idx_s = cumsum(c(1, df_rle$lengths[-length(df_rle$lengths)]))
      df_rle$idx_e = df_rle$idx_s + df_rle$lengths - 1

      # check length of segments
      length = d$end[df_rle$idx_e] - d$start[df_rle$idx_s]
      wh_not = which(length > max_amp_size)
      idx = lapply(wh_not, function(i) seq(df_rle$idx_e[i], df_rle$idx_s[i]))
      d$amped[unlist(idx)] = FALSE
      return(d[, c("amped","patient","name")])
    }) %>% do.call(what=rbind)

  amp_frac_mat = with(amplified_bins, tapply(amped, list(name, patient), mean))

  return(amp_frac_mat)
}

#' Calculate copy-number alteration summary stats
#'
#' @param d A data.frame containing CN data. Must contain a 'sample' or 'sample_barcode', a 'patient' , 'CN' or 'value' (the copy-number value), and 'name' (the CN window name).
#' @param window_infos A optional interval annotation data frame of windows. Must contain "chr","start","end","name","window" columns.
#' @param ... arguements passed to `findAmplificationBins`
#'
#' @return A list of  data.frame containing summary statics about how often a genomic window ('name') is lost, gained or amplified in each 'patient'.
#' @export
#'
calc_cn_summary_stats = function(d, window_infos=NULL, ...) {

  .getmode = function(v) {
    uniqv = unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  .get_median_cn = function(x) with(x, tapply(CN, sample, median, na.rm=TRUE))


  bcs = colnames(d)
  if ("sample" %in% bcs) {
  } else if ("sample_barcode" %in% colnames(d)) {
    d$sample = d$sample_barcode
  } else {
    stop("Missing 'sample' column in data.\n")
  }

  if (!"patient" %in% bcs) {
    stop("Missing 'patient' column in data.\n")
  }

  if ("CN" %in% bcs) {
  } else if ("value" %in% bcs) {
    d$CN = d$value
  } else {
    stop("Missing 'CN' column in data.\n")
  }

  d_per_pat = split(d, d$patient)
  median_cn = lapply(d_per_pat, .get_median_cn)
  mc_ploidy = lapply(median_cn, .getmode)

  # use most common ploidy as critical cn value
  crit_vals = unlist(mc_ploidy[rep(names(median_cn), sapply(median_cn, length))])
  names(crit_vals) = unlist(sapply(median_cn, names))

  # determine gains, losses and amplifications for each sample
  d$gained = d$CN > crit_vals[d$sample_barcode]
  frac_gains = with(d, t(tapply(gained, list(patient, name), mean)))

  d$loss = d$CN < crit_vals[d$sample_barcode]
  frac_losses = with(d, t(tapply(loss, list(patient, name), mean)))

  amped = data.frame(findAmplificationBins(d, crit_vals, ...))
  amped = amped[rownames(frac_gains), colnames(frac_gains)]
  dimnames(amped) = list(rownames(frac_gains), colnames(frac_gains))


  # bin results into a list
  result_list =
    list(
      gains = frac_gains,
      losses = frac_losses,
      focal_amplifications = amped
    )

  # mask any interval with missing data
  wh_mask = apply(is.na(do.call(what=cbind, result_list)), 1, any)
  for (i in seq_along(result_list)) result_list[[i]][wh_mask,] = NA

  # check dimnames are all equal
  dn = lapply(result_list, dimnames)
  for (i in seq_along(dn)) stopifnot(all.equal(dn[[1]], dn[[i]]))

  # extract window data from input and add these to results
  if (!is.null(window_infos)) {
    wh = colnames(window_infos) %in% c("chr","start","end","name","window")
    windows = unique(window_infos[,wh])
    stopifnot(!duplicated(windows[, colnames(windows) %in% c("chr","start","end")]))
  } else {
    names = unique(d$name)
    windows = data.frame(name=names, 'chr'='unknown', start=seq_along(names), end=seq_along(names)+1)
  }

  result_list = c(list(windows=windows), result_list)
  return(result_list)
}

#' Plot per-case copy-number alteration summary stats
#'
#' @param data A list of summary stats returned by the 'calc_cn_summary_stats' function.
#'
#' @return A plot showing the fraction of samples showing copy-number alterations in each case.
#' @export
#'
#'
plot_cn_altered_per_patient = function(data) {

  data_per_window_patient_summary_gains =
    cbind(data$windows, data$gains[data$windows$name,]) %>%
    reshape2::melt(id.vars=colnames(data$windows))

  data_per_window_patient_summary_losses =
    cbind(data$windows, data$losses[data$windows$name,]) %>%
    reshape2::melt(id.vars=colnames(data$windows)) %>%
    dplyr::mutate(value=-value)

  plot =
    data_per_window_patient_summary_gains %>%
    ggplot(aes(x=(start+end)/2, y=value)) +
    geom_line(color="red") +
    geom_line(data=data_per_window_patient_summary_losses, color="blue") +
    facet_grid(variable~chr, scales="free_x", space="free_x") +
    background_grid(major = "y", minor = "y") +
    cowplot::theme_cowplot() +
    xlab("Position on Chromosome") +
    ylab("Fraction") +
    ylim(c(-1,1)) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    scale_x_continuous(breaks=vector()) +
    theme(strip.text.x=element_text(angle=90)) +
    theme(strip.text.y=element_text(angle=0))

  return(plot)
}

#' Plot copy-number alteration summary stats across cases
#'
#' @param data A list of summary stats returned by the 'calc_cn_summary_stats' function.
#' @param cutoff_clonal The fraction of samples that have to show a copy-number alteration in a case to consider it to be 'clonal' (i.e., present in all samples).
#' @param cutoff_subclonal The fraction of samples that have to show a copy-number alteration in a case to consider it to be 'subclonal' (i.e., present in some samples).
#' @param bands Optional information on chromosome banding. The format should be like those in the 'D3GB' package (e.g., D3GB::GRCh38.bands'), that is a data.frame containing 'chr', 'start', 'end' 'name' and 'score' columns. Rows with a 'score' of 'acen' will be used to calculate actual centromer positions.
#' @param gene_labels Optional gene names to label. Positions will be retrieved with the 'get_gene_pos' position.
#' @param ... Optional arguments that will be passed to the 'get_gene_pos' function (see '?'get_gene_pos' for details).
#'
#' @return A plot showing the frequency of clonal and subclonal copy-number alterations across cases.
#' @export
#'
plot_cn_freq = function(data, cutoff_clonal=1, cutoff_subclonal=0, bands=NULL, gene_labels=NULL, ...) {

  alpha_vals = c("clonal"=1, "sub-clonal"=0.45)
  col_vals = c("Lost"="#3c73a8","Gained"="#ec2d01", "Amplified"="black")

  # get data for chr banding if available
  if (!is.null(bands)) {

    chr_centers =
      dplyr::filter(bands) %>%
      split(., .$chr) %>%
      sapply(function(x) max(x$end)/2)

    centromer_position =
      dplyr::filter(bands, score=="acen") %>%
      split(., .$chr) %>%
      sapply(function(x) mean(range(c(x$end,x$start))))

    chr_order = as.character(unique(D3GB::GRCh38.bands$chr))

  } else {

    chr_centers =
      with(data$windows,
           tapply(
             c(start, end),
             gsub("chr", "", c(chr, chr)),
             function(x) mean(range(x), na.rm = TRUE)
           )
      )

    centromer_position = NULL
    chr_order = as.character(unique(data$windows$chr))
  }


  # get data for plotting
  gains_w = data$gains[data$windows$name,]
  loss_w = data$losses[data$windows$name,]
  foc_gains_w = data$focal_amplifications[data$windows$name,]

  data_per_window_summary =
    data$windows %>%
    magrittr::set_rownames(data$windows$name) %>%
    cbind(clonal_gains=apply(gains_w >= cutoff_clonal, 1, sum, na.rm=1)) %>%
    cbind(subclonal_gains=apply(gains_w > cutoff_subclonal, 1, sum, na.rm=1)) %>%
    cbind(clonal_losses=apply(loss_w >= cutoff_clonal, 1, sum, na.rm=1)) %>%
    cbind(subclonal_losses=apply(loss_w > cutoff_subclonal, 1, sum, na.rm=1)) %>%
    cbind(amps=apply(foc_gains_w > cutoff_subclonal, 1, sum, na.rm=1)) %>%
    dplyr::mutate(chr=factor(gsub("chr", "", chr))) %>%
    dplyr::filter(!chr %in% c("X","Y")) %>%
    dplyr::mutate(center = (start + end) / 2) %>%
    dplyr::mutate(pos = center - chr_centers[as.character(chr)]) %>%
    dplyr::mutate(chr = factor(chr, chr_order, ordered=1))


  # variables for plotting
  y_range = c(-ncol(gains_w), ncol(gains_w))
  y_breaks = pretty(y_range)


  # create plot
  plot_cna_recurence =
    data_per_window_summary  %>%
    ggplot(aes(x=pos, y=subclonal_gains))  +
    cowplot::theme_cowplot() +
    geom_area(aes(alpha="sub-clonal", fill="Gained")) +
    geom_area(aes(y=clonal_gains, alpha="clonal", fill="Gained")) +
    geom_area(aes(y=-subclonal_losses, alpha="sub-clonal", fill="Lost")) +
    geom_area(aes(y=-clonal_losses, alpha="clonal", fill="Lost")) +
    geom_area(aes(y=amps, alpha="clonal", fill="Amplified")) +
    geom_hline(yintercept=0) +
    facet_grid(.~chr, scales="free_x", space="free_x") +
    theme(strip.text.x=element_blank()) +
    theme(strip.background.x=element_blank()) +
    scale_alpha_manual(values=alpha_vals) +
    scale_fill_manual(values=col_vals) +
    scale_x_continuous(breaks=0, labels="") +
    scale_y_continuous(breaks=y_breaks, labels=abs(y_breaks), limits=y_range) +
    xlab("") + ylab("Number of cases") + labs(fill="") +
    guides(alpha="none") +
    theme(legend.position = "bottom")


  # add manual labels of chr below plot
  chr_pos_data =
    data.frame(
      chr=factor(names(chr_centers), chr_order),
      pos=chr_centers,
      chr_lab_y_pos = y_range[1] - diff(y_range)*0.02
    ) %>% dplyr::filter(!chr %in% c("X","Y"))

  plot_cna_recurence =
    plot_cna_recurence +
    geom_text(data=chr_pos_data, aes(x=0, y=-Inf, label=chr), vjust=2, size=3.8) +
    coord_cartesian(ylim = y_range, clip = 'off')


  # add mark of centromer position
  if (!is.null(bands)) {
    plot_cna_recurence =
      plot_cna_recurence +
      geom_vline(
        data = data.frame(
          chr = factor(names(centromer_position), chr_order, ordered=T),
          pos = centromer_position - chr_centers
        ) %>% dplyr::filter(!chr %in% c("X", "Y")),
        aes(xintercept = pos),
        linetype = 3,
        size = 0.5,
        alpha = 0.2
      )
  }

  # add gene labels
  if (!is.null(gene_labels)) {

    if (is.character(gene_labels)) {
      gene_labels = get_gene_pos(gene_labels, ...)
    }

    gene_labels$win_name = NA
    for (i in seq_along(gene_labels$win_name)) {
      wh_c = as.character(gene_labels$chr[i]) == data_per_window_summary$chr
      wh_d = which.min(abs(data_per_window_summary$center[wh_c] - gene_labels$pos[i]))
      gene_labels$win_name[i] = data_per_window_summary$name[wh_c][wh_d]
    }

    gene_labels$y_pos = NA
    for (i in seq_along(gene_labels$y_pos)) {
      win = gene_labels$win_name[i]
      gain = data_per_window_summary[win, "subclonal_gains"]
      loss = data_per_window_summary[win, "subclonal_losses"]
      opts = c(-loss, gain)
      gene_labels$y_pos[i] = opts[nnet::which.is.max(abs(opts))]
    }

    gene_labels$chr = factor(gene_labels$chr, chr_order)
    gene_labels$pos = gene_labels$pos - chr_centers[as.character(gene_labels$chr)]

    for (c_sign in c(-1, 1)) {

      plot_cna_recurence =
        plot_cna_recurence +
        ggrepel::geom_text_repel(
          data = gene_labels[sign(gene_labels$y_pos) == c_sign, ],
          aes(y=y_pos, label=SYMBOL),
          min.segment.length = 0,
          nudge_y = c_sign * diff(y_range) * 0.1,
          force = 100,
          size=2.5
        )

    }

  }

  # switch of clipping for the strip texts:
  plot_cna_recurence = ggplotGrob(plot_cna_recurence)
  for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
    plot_cna_recurence$grobs[[i]]$layout$clip = "off"
  }

  return(plot_cna_recurence)
}

#' Function to look-up (average) gene positions.
#'
#' @param x A vector of gene symbols.
#' @param txdb Optional TxDb annotation object to use (default: TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene).
#' @param egdb Optional Bioconductor annotation package to use to lookup gene ids (default: egdb=org.Hs.eg.db::org.Hs.eg.db).
#' @param return_intervals Argument defining if all gene intervals should be returned (default: FALSE, i.e., a data.frame of gene center).
#'
#' @return
#' @export
#'
#' @examples get_gene_pos("KRAS") # data.frame that contains the position of the KRAS gene (around position 25220000 on chr12).
#' @examples get_gene_pos("KRAS", TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene) # data.frame that contains the position of the KRAS gene in the hg19 reference genome (around position 25380000 on chr12).
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
        ENTREZID = names(pos_data),
        chr = gsub("chr", "", sapply(pos_data, function(x) unique(x$CDSCHROM))),
        pos = sapply(pos_data, function(x) mean((x$CDSSTART + x$CDSEND)/2)),
        start = sapply(pos_data, function(x) min(x$CDSSTART)),
        end = sapply(pos_data, function(x) max(x$CDSEND)),

        row.names = NULL
      )

    res = merge(gene_to_id, pos_data_avg, by="ENTREZID", all=TRUE)

  }

  return(res)
}


#' Calculate PGA stats
#'
#' @param cna_data Copy-number summary stats, must contain the following columns: sample_barcode, CN (the copy-numner), start (start of the CN interval), end (end of the CN interval), analyte_name (the type of the sequencing the data came from),  tissue_type (the type of the analysed tissue, e.g., cancer and adenoma), msi_status (MSI status or any additional tumour property).
#'
#' @return A list that contains PGA summary stats. Can be passed to plot_pga_stats to create a plot of these.
#' @export
#'
calc_pga_stats = function(cna_data) {

  names(cna_data$analyte_name) = "NULL"

  cna_data =
    cna_data %>%
    dplyr::mutate(win_size = (end - start) / 10^6) %>%
    # If there is no CN data, also make the window equal to NA to avoid an incorrect sum
    dplyr::mutate(win_size = ifelse(is.na(value), NA, win_size))

  # Collect averages
  epicc_pl = plyr::ddply(
    cna_data,
    plyr::.(sample_barcode),
    summarize,
    Assay = unique(analyte_name),
    Tissue = unique(tissue_type),
    MSI = unique(msi_status),
    Ploidy = median(round(CN), na.rm = T),
    mPloidy = mean(round(CN), na.rm = T)
  )

  # ploidy lookup table
  ploidy = epicc_pl$Ploidy
  names(ploidy) = epicc_pl$sample_barcode
  ploidy = ploidy[as.character(cna_data$sample_barcode)]

  # Round the ploidy for it is 1n, 2n etc
  cna_data$bl_ploidy = as.numeric(round(cna_data$value) != ploidy)

  # Calculate the percentage not equal to the baseline ploidy whilst normalising for the window size
  epicc_pga = plyr::ddply(
    cna_data,
    plyr::.(sample_barcode),
    summarize,
    PGA = (sum(bl_ploidy * win_size, na.rm = T) / sum(win_size, na.rm = T))*100
  )

  # Gather what we need
  epicc_plot = merge(epicc_pl, epicc_pga)

  # Get specific tumour type annotation
  epicc_plot$Tumour_MSI = paste(epicc_plot$Tissue, epicc_plot$MSI)
  epicc_plot$Tumour_MSI = gsub("adenoma", "Adenoma", epicc_plot$Tumour_MSI)
  epicc_plot$Tumour_MSI = gsub("cancer", "Cancer", epicc_plot$Tumour_MSI)

  return(epicc_plot)
}


#' Plot PGA stats
#'
#' @param pga_stats data.frame containing data returned from `calc_pga_stats` (see `?calc_pga_stats` for details).
#' @param layers list of additional ggplot layers added to each subplot.
#'
#' @return A plot showing PGA stats.
#' @export
#'
plot_pga_stats = function(pga_stats, layers=list()) {

  # This plot looks at spread of _mean_ ploidy (not rounded)
  p1 =
    ggplot(
      pga_stats,
      aes(x = Tumour_MSI, y = mPloidy, fill = Tumour_MSI)
    ) +
    ylab("Ploidy") +
    xlab("") +
    geom_boxplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1)) +
    layers

  # Here we look at the spread of PGAs
  p2 =
    ggplot(
      pga_stats,
      aes(x = Tumour_MSI, y = PGA, fill = Tumour_MSI)
    ) +
    ylab("mean PGA") +
    xlab("") +
    geom_boxplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1)) +
    layers

  # This compares the two
  p3 =
    ggplot(
      pga_stats[order(pga_stats$Tumour_MSI, decreasing = T),],
      aes(x = mPloidy, y = PGA, col = Tumour_MSI, shape = Assay)
    ) +
    xlab("Ploidy") +
    ylab("PGA") +
    geom_point() +
    layers

  # ggarrange makes it look neat
  ggpubr::ggarrange(
    p1, p2, p3,
    labels = LETTERS[1:3],
    widths = c(1,1,2),
    ncol = 3,
    nrow = 1
  )


}
