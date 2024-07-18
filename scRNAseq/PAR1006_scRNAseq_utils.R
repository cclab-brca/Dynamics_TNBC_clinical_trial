require(R.utils)
require(data.table)
require(ggpubr)
require(rstatix)
require(tgutil)
require(tgstat)
require(tgconfig)
require(glue)
require(dplyr)
require(reticulate)
require(officer)

# wrap for opening a plot (png, ps or pdf). Expecting fn extension to be .png, and it will be chanhed according to the device parameter
.plot_start = function(fn, w, h, device='pdf', res=72, pointsize=12)
{
  if (device == "png") {
    png(filename=sub("ps$", "png", fn), width=w, height=h, res=res, pointsize = pointsize)
  }
  else if (device == "ps") {
    postscript(file=sub("png$", "ps", fn), width=w/res, height=h/res)
  }
  else if (device == "pdf") {
    pdf(file=sub("png$", "pdf", fn), width=w/res, height=h/res)
  }
  else {
    stop(sprintf("unknown output device type: %s", device))
  }
}

####
# reload packages and re-init parameters and directories
rl = function(scdb_dir="scdb", scfigs_dir="figs", config_fn=NULL, force_init=T)
{
  if (!is.null(config_fn) && file.exists(config_fn)) {
    tgconfig::override_params(config_fn, package="metacell")
  }
  
  dir.create(scdb_dir, recursive = T, showWarnings = F)
  dir.create(scfigs_dir, recursive = T, showWarnings = F)
  if (!exists(".scdb_base") || force_init) {
    scdb_init(scdb_dir, force_reinit=force_init)
  }
  scfigs_init(scfigs_dir)
}

####
# Test if object exists
scdb_obj_exists = function(obj_type, obj_id) 
{
  if(!exists(".scdb")) {
    message("scdb not initialized")
    invisible(F)
  } else {
    invisible(file.exists(sprintf("%s/%s.%s.Rda", .scdb_base, obj_type, obj_id)))
  }
}


###
# return a named vector of colors to names (from the mc color_key table)
get_mc_col2group = function(mc, white_is_undet=T) {
  col2group = as.character(mc@color_key$group)
  names(col2group) = as.character(mc@color_key$color)
  col2group = col2group[ names(col2group) %in% unique(mc@colors)] 
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      col2group = c(col2group, c('white'='Undet'))
    }
  }
  col2group[unique(names(col2group))]
}

####
# return a named vector of colors to names (from the mc color_key table)
get_mc_group2col = function(mc, white_is_undet=T) {
  group2col = as.character(mc@color_key$color)
  names(group2col) = as.character(mc@color_key$group)
  group2col = group2col[ group2col %in% unique(mc@colors)] 
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      group2col = c(group2col, c('Undet'='white'))
    }
  }
  
  group2col[unique(names(group2col))]
}


###
# return a objects for mc annotation using pheatmap
get_mc_pheatmap_ann = function(mc, name) {
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  
  df = data.frame(row.names=colnames(mc@mc_fp))
  df[, name] = col2grp[mc@colors]
  
  list(ann=df, ann_cols=grp2col[unique(names(grp2col))])
}

###
#
mcell_mc_plot_legend = function(mc_id, ofn=NULL, ncols=1, p_cex=1, p_pch=19, l_bty='n', names_ord=NULL, sort_names=T, filter_names=NULL, name=NULL)
{
  mc = scdb_mc(mc_id)
  if (is.null(ofn)) {
    ofn = scfigs_fn(mc_id, sprintf("%sannot_legend", ifelse(is.null(name), "", paste0(name, "_"))))
  }
  
  grp2col = get_mc_group2col(mc)
  nms = names(grp2col)
  if (!is.null(filter_names)) {
    nms = filter_names
  }
  if (!is.null(names_ord)) {
    stopifnot(all(names_ord %in% nms))
    nms = names_ord
  }
  else if (sort_names) {
    nms = sort(nms)
  }
  plot.new
  sh = max(strheight(nms, unit='in', cex=p_cex))
  sw = max(strwidth(nms, unit='in', cex=p_cex))
  dev.off()
  
  fig_w = ncols * (sw + 0.2) + 0.5
  fig_h = ceiling(length(nms) / ncols) * sh * 2 + 0.5
  .plot_start(ofn, fig_w * 72, fig_h * 72)
  par(mar=c(0, 0, 0, 0))
  plot.new()
  plot.window(0:1, 0:1)
  legend("topleft", legend=nms, col=grp2col[nms], pch=p_pch, cex=p_cex, ncol=ncols, bty=l_bty)
  dev.off()
  
  invisible(c(fig_w, fig_h))
}

####
# Scatter plot lfp values of 2 genes (or any other value per mc if supplier in x and/or y)
plt = function(nm1, nm2, lfp, cols, ofn=NULL, x=NULL, y=NULL, show_mc_ids=T, cex=3, cex.lab=1, add_grid=F, main="", xlim=NULL, ylim=NULL, plot_legend=F, leg_ncol=1, leg_exp=1.5, leg_cex=0.8, leg_ord=NULL, mc_obj=NULL, h_lines=NULL, v_lines=NULL, hlines_lty=2, vlines_lty=2) {
  stopifnot(!plot_legend | !is.null(mc_obj))
  
  if (is.null(x)) {
    x = lfp[nm1, ]
  }
  if (is.null(y)) {
    y = lfp[nm2, ]
  }
  if (!is.null(ofn)) {
    .plot_start(ofn, 450 * ifelse(plot_legend, leg_exp, 1), 450)
  }
  if (is.null(xlim)) {
    xlim = c(min(x), max(x) +  diff(range(x)) * ifelse(plot_legend, leg_exp-1, 0))
  }
  if (is.null(ylim)) {
    ylim = range(y)
  }
  plot(x, y, pch=21, cex=cex, bg=cols, xlab=nm1, ylab=nm2, cex.lab=cex.lab, main=main, xlim=xlim, ylim=ylim)
  if (show_mc_ids) {
    text(x, y, colnames(lfp), cex=cex/4)
  }
  if (add_grid) {
    grid(col='black', lwd=0.5)
  }
  if (!is.null(h_lines)) {
    abline(h=h_lines, lty=hlines_lty)
  }
  if (!is.null(v_lines)) {
    abline(v=v_lines, lty=vlines_lty)
  }
  
  if (plot_legend) {
    grp2col = get_mc_group2col(mc_obj)
    col2grp = get_mc_col2group(mc_obj)
    if (is.null(leg_ord)) {
      grp2col = grp2col[sort(unique(col2grp[mc_obj@colors[as.numeric(colnames(lfp))]]))]
    } else {
      stopifnot(all(leg_ord %in% names(grp2col)))
      grp2col = grp2col[leg_ord]
    }
    legend("topright", legend=names(grp2col), col = grp2col, pch=19, bty='n', cex=leg_cex, ncol = leg_ncol)
  }
  if (!is.null(ofn)) {
    dev.off()
  }
}

# 
#
####
mc_plot_e_gc_barplots = function(mc_id, name, genes=NULL, gene_groups=NULL, ncolumns=2, panel_height=60, panel_width=400, mc_ord=NULL, ord_first_by_color=T, grp_ord=NULL, intra_grp_ord_by=NULL, n_ideal_umi=1000, shared_ylim=F) 
{
  mc = scdb_mc(mc_id)
  lfp = log2(mc@mc_fp)
  col2group = get_mc_col2group(mc)
  
  if (is.null(gene_groups)) {
    if (is.null(genes)) {
      marks_gset = scdb_gset(mc_id)
      genes = names(marks_gset@gene_set)
    }
    e_gc = mc@e_gc[genes, ] * n_ideal_umi
    nms = genes
  } else {
    e_gc = t(sapply(names(gene_groups), function(nm) colMeans(mc@e_gc[gene_groups[[nm]], ]))) * n_ideal_umi
    nms = names(gene_groups)
  }
  
  if (!is.null(mc_ord)) {
    e_gc = e_gc[, mc_ord]
  }
  else if (ord_first_by_color) {
    if (is.null(grp_ord)) {
      grp_ord = unique(col2group)
    }
    intra_grp_ord = 1:ncol(lfp)
    if (!is.null(intra_grp_ord_by)) {
      stopifnot(all(names(intra_grp_ord_by) %in% col2group) && all(unlist(intra_grp_ord_by) %in% rownames(e_gc)))
      gene_ord_vals = sapply(names(intra_grp_ord_by), function(nm) colMeans(e_gc[intra_grp_ord_by[[nm]], , drop=F]) )
      intra_grp_ord = gene_ord_vals[cbind(1:ncol(e_gc), col2group[mc@colors])]
    }
    
    e_gc = e_gc[, order(as.numeric(ordered(col2group[mc@colors], levels=grp_ord)) + intra_grp_ord * 1e-6)]
  }
  
  .plot_start(scfigs_fn(mc_id, sprintf("mc_geom_mean_%s", name)), ncolumns * panel_width, panel_height * ceiling(length(nms) / ncolumns))
  layout(matrix(1:(length(nms) + length(nms) %% 2), ncol=ncolumns))
  par(mar=c(0.5, 12, 0.5, 1))
  
  for (g in nms) {
    barplot(e_gc[g, ], border=NA, col=mc@colors[as.numeric(colnames(e_gc))], xaxt='n', yaxt='n', space=0, ylim=if (shared_ylim) { range(e_gc) } else { range(e_gc[g,]) })
    yaxp = par("yaxp")
    axis(2, yaxp=c(yaxp[1], yaxp[2], 1), las=2, cex=1)
    mtext(g, 2, line=1.5, cex=1.2, las=2)
  }
  dev.off()
  
}


####
# Calculate %UMIs from genes in the given gene set
get_gset_f_umis = function(mat_id, gset_id, update_mat=F)
{
  
  mat = scdb_mat(mat_id)
  gset = scdb_gset(gset_id)
  
  stopifnot(!is.null(mat) && !is.null(gset))
  
  genes = names(gset@gene_set)
  missed_genes = setdiff(genes, mat@genes)
  if (length(missed_genes) > 0) {
    message(sprintf("UMI matrix miss %d genes from %s (%s ...)", length(missed_genes), gset_id, paste0(missed_genes[1:min(length(missed_genes), 10)], collapse=", ")))
    genes = intersect(genes, mat@genes)
  }
  f = colSums(mat@mat[genes, ]) / colSums(mat@mat)
  
  if (update_mat) {
    mat@cell_metadata[colnames(mat@mat), paste0('f_', gset_id)] = f
    scdb_add(mat_id, mat)
  }
  
  invisible(f)
}

###
# differential expression between 2 groups of metacell (if mcs1/2 supplied) or cells (if nms1/2 supplier)
# loo_grps - perform the comparison with leave-one-out and only keep genes that appear in all runs (and either enriched or depleted in all). loo_grps is a list of the loo groups mapping to cell names
diff_expr = function(mc, mat_ds, mcs1=NULL, mcs2=NULL, reg=5, min_max_umi=0, nms1=NULL, nms2=NULL, filter_outlier_genes=F, compare_top_mc_to_n_highest=3, max_top_to_n_highest_ratio=3, verbose=T, geo_mean=F, geo_mean_per_cell=F, enr_min_mean_umi=0.1, enr_min_f_pos=0.5, calculate_p_value=F, p_adjust_method='fdr', loo_grps=NULL)
{
  if (is.null(nms1)) {
    nms1 = names(mc@mc)[mc@mc %in% mcs1]
  }
  if (is.null(nms2)) {
    nms2 = names(mc@mc)[mc@mc %in% mcs2]
  }
  nms1 = intersect(colnames(mat_ds), nms1)
  nms2 = intersect(colnames(mat_ds), nms2)
  if (verbose) {
    message(sprintf("comparing %d vs %d cells", length(nms1), length(nms2)))
  }
  
  loo_de_run = function(c_loo_grp) {
    c_nms1 = nms1
    c_nms2 = nms2
    if (!is.null(c_loo_grp)) {
      c_nms1 = setdiff(nms1, loo_grps[[c_loo_grp]])
      c_nms2 = setdiff(nms2, loo_grps[[c_loo_grp]])
    }
    
    if (geo_mean) {
      df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), mu1=apply(mat_ds[, c_nms1], 1, function(y) {exp(mean(log(1+y)))-1}), mu2=apply(mat_ds[, c_nms2], 1, function(y) {exp(mean(log(1+y)))-1}), stringsAsFactors = F)
      df$tot1 = df$mu1 * length(c_nms1)
      df$tot2 = df$mu2 * length(c_nms2)
    } else {
      df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), tot1=Matrix::rowSums(mat_ds[, c_nms1]), tot2=Matrix::rowSums(mat_ds[, c_nms2]), stringsAsFactors = F)
    }
    df$mean1 = Matrix::rowMeans(mat_ds[, c_nms1])
    df$mean2 = Matrix::rowMeans(mat_ds[, c_nms2])
    df$f_pos1 = Matrix::rowMeans(mat_ds[, c_nms1] > 0)
    df$f_pos2 = Matrix::rowMeans(mat_ds[, c_nms2] > 0)
    
    norm_by = min(sum(df$tot1), sum(df$tot2))
    df$tot1 = df$tot1 / sum(df$tot1) * norm_by
    df$tot2 = df$tot2 / sum(df$tot2) * norm_by
    
    if (geo_mean && geo_mean_per_cell) {
      df$enr = log2( (df$mu1 + reg) / (df$mu2 + reg))
    } else {
      df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
    }
    
    
    # Filtering
    df = df[pmax(df$tot1, df$tot2) >= min_max_umi &
              ifelse(df$enr > 0, df$mean1, df$mean2) >= enr_min_mean_umi &
              ifelse(df$enr > 0, df$f_pos1, df$f_pos2) >= enr_min_f_pos, ]
    
    df = df[order(df$enr, decreasing=T), ]
    
    if (filter_outlier_genes) {
      fp = mc@mc_fp[intersect(rownames(mc@mc_fp), df$gene), ]
      if (!is.null(mcs1)) {
        fp = fp[, mcs1, drop=F]
      }
      gmax = apply(fp, 1, max)
      gnext = apply(fp, 1, function(v) head(tail(sort(v), n=compare_top_mc_to_n_highest), n=1) )
      df[rownames(fp), 'out_r'] =  gmax/gnext
      to_filt = !is.na(df$out_r) & df$out_r > max_top_to_n_highest_ratio & df$enr > 0
      if (sum(to_filt) > 0) {
        if (verbose) {
          message(sprintf("filtering %d outlier genes:", sum(to_filt)))
        }
        print(df[to_filt, ])
        df = df[!to_filt, ]
      }
    }
    return(df)
  }
  
  n_cores = tgconfig::get_param("mc_cores", "metacell")
  if (n_cores > 1) {
    doMC::registerDoMC(n_cores)
  }
  
  df = loo_de_run(NULL)
  if (!is.null(loo_grps)) {
    loo_retained_genes = do.call("rbind", alply(names(loo_grps), 1, loo_de_run, .parallel = n_cores > 1)) %>%
      group_by(gene) %>%
      summarise(n=n(), all_enr = all(enr > 0), all_dep = all(enr < 0)) %>%
      filter(n == length(loo_grps) & (all_enr | all_dep)) %>% 
      pull(gene)
    df = df[loo_retained_genes, ]
  }
  
  # Add p-value if requested (wilcox test per gene, correct p-value for multiple tests)
  if (calculate_p_value) {
    
    #all_pvals = setNames(unlist(plyr::alply(rownames(mat_ds), 1, function(g) { wilcox.test(mat_ds[g, c_nms1], mat_ds[g, c_nms2])$p.value }, .parallel = n_cores > 1)),	                       rownames(mat_ds))
    all_pvals = setNames(unlist(plyr::alply(df$gene, 1, function(g) { wilcox.test(mat_ds[g, nms1], mat_ds[g, nms2])$p.value }, .parallel = n_cores > 1)), df$gene)
    adj_pvals = p.adjust(all_pvals, method=p_adjust_method)
    
    df$pval = all_pvals[df$gene]
    df$pval_adjust = adj_pvals[df$gene]
  }
  
  invisible(df)	        
}

###
# General function to select feature genes, auto-compute parameres by default
select_feature_genes = function(gstat_id, gset_id=gstat_id, lateral_gset_id = NULL, T_vm = NULL, T_tot = NULL, T_top3 = NULL, verbose=T) 
{
  gs = scdb_gstat(gstat_id)
  if (is.null(T_vm)) {
    T_vm = quantile(gs$ds_vm_norm, 0.75)
  }
  if (is.null(T_tot)) {
    T_tot = round(median(gs$tot))
  }
  if (is.null(T_top3)) {
    T_top3 = round(median(gs$ds_top3))
  }
  if (verbose) {
    message(sprintf("marker selection filters for %s: \ntotal umis:\t%d (%.2f) > %.2f\nds top3:\t%d (%.2f) > %.2f\nds vm norm:\t%d (%.2f) > %.2f", gstat_id, sum(gs$tot > T_tot), mean(gs$tot > T_tot), T_tot, sum(gs$ds_top3 > T_top3), mean(gs$ds_top3 > T_top3), T_top3, sum(gs$ds_vm_norm > T_vm), mean(gs$ds_vm_norm > T_vm), T_vm))
  }	
  
  mcell_gset_filter_varmean(gstat_id, gset_id, T_vm=T_vm, force_new=T)
  mcell_gset_filter_cov(gstat_id, gset_id, T_tot=T_tot, T_top3=T_top3)
  mcell_plot_gstats(gstat_id, gset_id, max_vm=NULL)
  
  if (!is.null(lateral_gset_id)) {
    marker_gset = scdb_gset(gset_id)	
    lateral_gset = scdb_gset(lateral_gset_id)
    if (verbose){
      message(sprintf("removing %d lateral genes from markers, left with %d", length(intersect(names(marker_gset@gene_set), names(lateral_gset@gene_set))), length(setdiff(names(marker_gset@gene_set), names(lateral_gset@gene_set)))))
    }
    marker_gset = gset_new_restrict_gset(marker_gset, lateral_gset, inverse=T, "cgraph markers w/o lat genes")
    scdb_add_gset(gset_id, marker_gset)
  }
}

####
# Build blaclisted gene sets by gene names (mitochondrial, IG, ncRNA genes, RP-# genes and snoXXX genes)
meta_build_blacklist_gsets_by_gene_nms = function(all_id, ds_nm, gene_pref=NULL, mito_mode='minimal') 
{
  full_m = scdb_mat(all_id)
  
  # mitochondrial gene set
  mito_gset_id = sprintf("%s_mito", ds_nm)
  
  mt_cands = grep("^MT-", full_m@genes, v=T, perl=T, ignore.case = T)
  mt_genes = setNames(rep('neto MT', length(mt_cands)), mt_cands)
  if (mito_mode != 'minimal') {
    mt_cands = grep(paste0(paste0("^", gene_pref, c("MT-", "MTRN", "MTAT", "MTND", "MRP")), collapse="|"), full_m@genes, v=T, perl=T, ignore.case = T)
    mito = data.table::fread("/Users/eyal-l01/proj/mc_common/data/MitoCarta2_human.txt", header=T, sep="\t", stringsAsFactors=F)
    mito_genes = paste0(gene_pref, mito$Symbol)
    mt_both = intersect(mt_cands, mito_genes)
    mt_cands = setdiff(mt_cands, mt_both)
    mitocarta = setdiff(mito_genes, mt_both)
    mt_genes = c(rep('regexp MT', length(mt_cands)), rep('MitoCarta2', length(mitocarta)), rep('MitoCarta2 and regexp MT', length(mt_both)))
    names(mt_genes) = c(mt_cands, mitocarta, mt_both)
  }
  scdb_add_gset(mito_gset_id, gset_new_gset(mt_genes, 'mitochondrial genes'))
  
  # IG gene set
  ig_gset_id = sprintf("%s_ig", ds_nm)
  ig_nms = grep(paste0(paste0("^", gene_pref, c("IGK", "IGL", "IGJ", "IGH", "IGBP", "IGSF")), collapse="|"), full_m@genes, v=T, perl=T, ignore.case = T)
  ig_genes = rep("IG", length(ig_nms))
  names(ig_genes) = ig_nms
  scdb_add_gset(ig_gset_id, gset_new_gset(ig_genes, 'IG genes'))
  
  # ncRNA gene set
  ncrna_gset_id = sprintf("%s_ncrna", ds_nm)
  ncrna_nms = paste0(gene_pref, c('MALAT1', 'XIST', 'NEAT1', 'hsa-mir-6723'))
  ncrna_genes = rep("ncRNA", length(ncrna_nms))
  names(ncrna_genes) = ncrna_nms
  scdb_add_gset(ncrna_gset_id, gset_new_gset(ncrna_genes, "ncRNA genes"))
  
  # RP pseudo genes set
  ncrp_gset_id = sprintf("%s_ncrp", ds_nm)
  ncrp_nms = grep(paste0("^", gene_pref, "sRP[0-9]+-"), full_m@genes, v=T, perl=T, ignore.case = T)
  ncrp_genes = rep("ncRP", length(ncrp_nms))
  names(ncrp_genes) = ncrp_nms
  scdb_add_gset(ncrp_gset_id, gset_new_gset(ncrp_genes, "RP##- genes"))
  
  # sno Genes
  sno_gset_id = sprintf("%s_sno", ds_nm)
  sno_nms = grep(paste0("^", gene_pref, "SNOR[AD][0-9]+"), full_m@genes, v=T, perl=T, ignore.case = T)
  sno_genes = rep("sno", length(sno_nms))
  names(sno_genes) = sno_nms
  scdb_add_gset(sno_gset_id, gset_new_gset(sno_genes, "SNORA/D genes"))
  
}

#
# build clean mat: remove blaclisted genes (mitochondrial, ncRNA, RP[0-9], IG) and small cells
meta_build_blist_filtered_master_mat = function(all_id, filt_id, ds_nm, min_umis_post_gene_ignore = 500, min_umis_pre_gene_ignore=500, max_mito_f = 0.6, filt_mat_by_column=NULL, sample_field=NULL, gsets_to_filter=c("mito", "ig", "ncrna", "ncrp", "sno"))
{
  full_m = scdb_mat(all_id)
  stopifnot(is.null(filt_mat_by_column) || filt_mat_by_column %in% colnames(full_m@cell_metadata) || all(is.logical(full_m@cell_metadata[, filt_mat_by_column])))
  
  blist_gsets = sapply(gsets_to_filter, function(gnm) scdb_gset(sprintf("%s_%s", ds_nm, gnm)))
  blist_genes = unlist(lapply(blist_gsets, function(gs) names(gs@gene_set)))
  
  uc = Matrix::colSums(full_m@mat)
  
  mt_genes = intersect(names(blist_gsets[["mito"]]@gene_set), full_m@genes)
  mito_f = Matrix::colSums(full_m@mat[mt_genes, ]) / uc
  
  full_m@cell_metadata$f_mito = mito_f
  scdb_add_mat(all_id, full_m)
  
  # plot %mito vs log umis
  .plot_start(scfigs_fn(all_id, "fMito_vs_logUMIs"), 300, 300)
  valid_cols = ifelse(mito_f >= max_mito_f, 'red', 'black')
  if (!is.null(filt_mat_by_column)) {
    valid_cols = ifelse(full_m@cell_metadata[, filt_mat_by_column], 'black', 'red')
  }
  plot(log2(uc), mito_f, pch=19, cex=0.5, col=valid_cols , xlab="UMIs (log2)", ylab="%mito")
  dev.off()
  
  # plot %mito by sample (if field was supplied)
  if (!is.null(sample_field)) {
    .plot_start(scfigs_fn(all_id, "fMito_by_sample"), 500, 500)
    par(mar=c(16, 4, 1, 1))
    boxplot(mito_f ~ full_m@cell_metadata[, sample_field], las=2, col='navyblue', notch=T, pch=19, cex=0.5, xlab='', ylab="% mito")
    abline(h=max_mito_f, col='red', lty=2)
    dev.off()
    
  }
  # filter cells with low counts, large MT fraction, ignore MT genes, RP##- and mega-strong RNA genes
  mcell_mat_ignore_genes(filt_id, all_id, blist_genes)
  
  filt_mat = scdb_mat(filt_id)
  to_ignore = filt_mat@cells[mito_f >= max_mito_f | Matrix::colSums(filt_mat@mat) <= min_umis_post_gene_ignore | uc <= min_umis_pre_gene_ignore]
  if (!is.null(filt_mat_by_column)) {
    to_ignore = filt_mat@cells[!filt_mat@cell_metadata[, filt_mat_by_column]]
  }
  
  mcell_mat_ignore_cells(filt_id, filt_id, union(filt_mat@ignore_cells, to_ignore))
  
  # working on the filtered mat
  mcell_add_gene_stat(filt_id, filt_id)
}

####
# Wrapper to cluster genes based on the mat given the input anchor genes
select_gene_modules_by_anchor_genes = function(mat_id, gene_anchors, gset_nm, cor_thresh = 0.1, gene_anti=c(), sz_cor_thresh=0.7, nclusts=20, downsample_n=NA) 
{
  tab_fn = sprintf("%s/%s.txt", scfigs_dir(mat_id, "gmods_by_anchors"), gset_nm)
  message("mcell_mat_rpt_cor_anchors")
  
  mcell_mat_rpt_cor_anchors(mat_id=mat_id,
                            gene_anchors = gene_anchors,
                            cor_thresh = cor_thresh,
                            gene_anti = gene_anti,
                            tab_fn = tab_fn,
                            sz_cor_thresh=sz_cor_thresh)
  
  
  foc_gcor = read.table(tab_fn, sep="\t", h=T, stringsAsFactors=F, check.names=F)
  foc_neto = foc_gcor[,setdiff(colnames(foc_gcor),c("sz_cor","max","neg_max"))]
  foc_genes = if (is.null(ncol(foc_neto))) { setNames(rep(1, length(foc_neto)), rownames(foc_gcor)) } else { apply(foc_neto, 1, which.max) }
  gset = tgGeneSets(foc_genes, gset_nm)
  
  scdb_add_gset(gset_nm, gset)
  
  sub_mat_id = paste(mat_id, gset_nm, sep="_")
  
  mcell_mat_ignore_genes(sub_mat_id, mat_id, names(foc_genes), reverse=T)
  
  message("mcell_gset_split_by_dsmat")
  mcell_gset_split_by_dsmat(gset_nm, sub_mat_id, nclusts)
  mcell_plot_gset_cor_mats(gset_nm, sub_mat_id)
  
}

####
# Generate lateral gene sets (cluster by anchor genes and selects clusters containing the anchror genes) - cell cycle, IFN response and stress
meta_generate_lateral_gene_sets = function(filt_id, ds_nm, specie="human", anchors=NULL, nclusts=20, add_by_name=NULL, cor_thresh=0.1, gene_pref=NULL, selected_sets=NULL, rebuild=T, n_cores=1)
{
  if (is.null(anchors)) {
    anchors = switch(specie,
                     human = list(cell_cycle=c('MKI67', 'HIST1H1D', 'PCNA', 'SMC4', 'MCM3', 'TYMS', 'TOP2A', 'TUBB'), ifn=c('ISG15', 'OAS1', 'WARS', 'IFIT1'), stress=c("TXN", "HSP90AB1", "HSPA1A", "FOS", "JUN", "HIF1A"), hypoxia=c("VEGFA", "NDRG1", 'BNIP3')),
                     mouse = list(cell_cycle=c('Mki67', 'Hist1h1d', 'Pcna', 'Smc4', 'Mcm3', 'Tyms', 'Top2a'), ifn=c('Isg15', 'Wars', 'Ifit1'), stress=c("Hsp90ab1", "Hspa1a", "Fos", "Jun", "Hif1a")))
  }
  if (is.null(add_by_name)) {
    add_by_name = switch(specie, 
                         human = list(cell_cycle='^HIST|^CENP|^SMC[0-9]', ifn="^IFI"),
                         mouse = list(cell_cycle='^Hist|^Cenp|^Smc[0-9]', ifn="^Ifi"))
  }
  
  mat_f = scdb_mat(filt_id)
  
  if (is.null(selected_sets)) {
    selected_sets = names(anchors)
  }
  
  if (n_cores > 1 & length(selected_sets) > 1) {
    doMC::registerDoMC(cores=n_cores)
  }	
  
  inner_gset_by_anchor = function(gnm) {
    gs_id = paste(ds_nm, gnm, sep="_")
    gs_filt_id = paste(gs_id, "filt", sep="_")
    if (rebuild | !scdb_obj_exists("gset", gs_filt_id)) {
      message(sprintf("Processing %s...", gnm))
      
      curr_anchors = paste0(gene_pref, anchors[[gnm]])
      select_gene_modules_by_anchor_genes(filt_id, curr_anchors, gset_nm=gs_id, cor_thresh=cor_thresh, sz_cor_thresh=cor_thresh, nclusts=nclusts)
      
      gset = scdb_gset(gs_id)
      selected_clusts = unique(gset@gene_set[curr_anchors])
      for (cl in selected_clusts) {
        file.rename(sprintf("%s/%s_%s.gset_cors/%d.png", .scfigs_base, ds_nm, gnm, cl), sprintf("%s/%s_%s.gset_cors/filt_%d.png", .scfigs_base, ds_nm, gnm, cl))
      }
      mcell_gset_remove_clusts(gs_id, filt_clusts=selected_clusts, new_id = gs_filt_id, reverse=T)
      
      gset = scdb_gset(gs_filt_id)
      add_re = add_by_name[[gnm]]
      if (!is.null(add_re)) {
        gset = gset_add_genes(gset, setdiff(grep(add_re, mat_f@genes, v=T, perl=T), names(gset@gene_set)), max(gset@gene_set)+1)
      }
      scdb_add_gset(gs_filt_id, gset)		
    }
    gset = scdb_gset(gs_filt_id)
    names(gset@gene_set)
  }
  
  lat_genes = plyr::alply(selected_sets, 1, inner_gset_by_anchor, .parallel=n_cores > 1 & length(selected_sets) > 1)
  names(lat_genes) = selected_sets
  lat_genes_memb = unlist(lapply(seq_along(lat_genes), function(i) { v = lat_genes[[i]]; r = rep(i, length(v)); names(r) = v; r }))
  scdb_add_gset(paste0(ds_nm, "_lateral"), gset_new_gset(lat_genes_memb, sprintf('lateral: %s', paste0(names(anchors), collapse=", "))))
  
}

####
# Pipline generating mc from mat
meta_mat2mc = function(mat_id, cells=NULL, name="", lateral_gset_id=NULL, 
                       T_vm = NULL, T_tot = NULL, T_top3 = NULL, T_lfc = 3000, 
                       cgraph_knn = NULL, cgraph_downsamp=T, 
                       bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
                       mc_K=30, min_mc_size=30, mc_alpha=2, feat_gset_id=NULL,
                       split_mc_with_dbscan_and_filt_outliers=T, rebuild=T) 
{
  set.seed(42)
  
  # create submat if required
  if (!is.null(cells)) {
    stopifnot(nchar(name) > 0)
    new_id = paste(mat_id, name, sep="_")
    mcell_mat_ignore_cells(new_id, mat_id, cells, reverse=T)
    mcell_add_gene_stat(new_id, new_id)
  } else {
    new_id = mat_id
  }
  
  # Add gstats if it doesnt exist
  mcell_add_gene_stat(new_id, new_id, force=rebuild)
  
  # select genes to affect graph creation
  if (is.null(feat_gset_id)) {
    select_feature_genes(gstat_id = new_id, gset_id=new_id, lateral_gset_id = lateral_gset_id, T_vm = T_vm, T_tot = T_tot, T_top3 = T_top3) 
    feat_gset_id = new_id
  }
  
  # create cgraph
  if (is.null(cgraph_knn)) {
    mat = scdb_mat(new_id)	
    cgraph_knn = min(max(ceiling(mat@ncells / 200), 20), 250)
  }
  if (rebuild || !scdb_obj_exists('cgraph', new_id)) {
    message(sprintf("creating balanced Knn graph, K = %d", cgraph_knn))
    mcell_add_cgraph_from_mat_bknn(mat_id=new_id, gset_id=feat_gset_id, graph_id=new_id, K=cgraph_knn, dsamp=cgraph_downsamp)
  }
  message(new_id, " cgraph done")
  
  # bootstrap
  if (rebuild || !scdb_obj_exists('coclust', new_id)) {
    mcell_coclust_from_graph_resamp(new_id, new_id,
                                    min_mc_size=round(cgraph_knn/5),
                                    p_resamp=bootstrap_p_resamp,
                                    n_resamp=bootstrap_n_resamp)
  }
  message(new_id, " bootstrap done")
  
  # mc from coclust matrix
  if (rebuild || !scdb_obj_exists('mc', new_id)) {
    mcell_mc_from_coclust_balanced(new_id, new_id, new_id, K=mc_K, min_mc_size=min_mc_size, alpha=mc_alpha)
    
    # find and remove outliers
    # breaking down heterogenous metacells and removing outliers by extreme gene expression
    if (split_mc_with_dbscan_and_filt_outliers) {
      mcell_mc_split_filt(new_id, new_id, new_id, T_lfc=T_lfc, plot_mats=F)
      mcell_plot_outlier_heatmap(mc_id=new_id, mat_id = new_id, T_lfc  = T_lfc)
    }	
  }
}

####
#
#' Wrapper function for assigning groups (and colors) to metacells based on their confusion matrix clustering. 
#' Run first to create the clustered confusion matrix, use it and the mc_sup object to classify sup_ids and create the supmc_file (and optionally the marks_file) and then rerun it to assign colors to the metacell object
#' 
#'
#' @param mc_id id of the metacell object (both used as the input and as the output object)
#' @param graph_id id of the cell knn graph object
#' @param supmc_file file name assigning name and colors to clusters of metacells (optional, if null then no coloring is done. tab delimited with these columns: supid, color, name)
#' @param marks_file file name assigning name and colors to metacells by thresholding expression of genes (optional, will override colors assigned by the supmc_file, tab-delimited with these columns: name, gene, color, T_fold. T_fold is the threshold on the gene log2 fp value).
#' @param col_by_cutree_k Color supmcs by cutree col_by_cutree_k nodes on the hc (using default colors and C# cluster names. Ignored if marks_file != NULL)
#' @param res return value of this function (can be supplied after the first run to shorten run-time)
#' @param show_mc_ids in heatmap plot
#'
#' @return list containining mc_hc (hierrarchical clustering of the metacells) and mc_sup (info on derived clusters of metacells) 
#' 
colorize_by_confusion_mat = function(mc_id, graph_id=mc_id, supmc_file=NULL, marks_file=NULL, col_by_cutree_k=0, annot_ref="hg_Blueprint", res=NULL, show_mc_ids=F, T_gap=0.04, min_nmc=2, width=1600, height=3000) 
{
  message(sprintf("colorize by confusion mat: mc_id=%s\tgraph_id=%s", mc_id, graph_id))
  mc = scdb_mc(mc_id)
  
  # Cluster metacells by confusion matrix
  if (is.null(res)) {	
    mc_hc = mcell_mc_hclust_confu(mc_id=mc_id,
                                  graph_id=graph_id)
  }
  else {
    mc_hc = res$mc_hc
  }
  
  # Annotate metacell clusters (globally and locally enriched genes)
  if (is.null(res)) {
    mc_sup = mcell_mc_hierarchy(mc_id=mc_id,
                                mc_hc=mc_hc, T_gap=T_gap, n_min_outcells=min(300, mean(table(mc@mc))))
  }
  else {
    mc_sup = res$mc_sup
  }
  
  # Colorize metacells based on the manually created supmc_file (and optionally by the marks_file)
  if (is.null(supmc_file)) {
    if (col_by_cutree_k > 0) {
      message(sprintf("Color mcs by cutree=%d", col_by_cutree_k))
      cls = cutree(mc_hc, k=col_by_cutree_k)
      col_pal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075')
      if (col_by_cutree_k > length(col_pal)) { stop(sprintf("In colorize_by_confusion_mat, coloring by cutree k=%d > maximum supported (%d)", col_by_cutree_k, length(col_pal)))	  }
      
      mc = scdb_mc(mc_id)
      mc@colors = col_pal[cls]
      mc@color_key = data.frame(gene="", group=sprintf("C%d", 1:col_by_cutree_k), color=col_pal[1:col_by_cutree_k])
      
      if (!is.null(annot_ref)) {
        annot_df = switch(annot_ref,
                          hg_HPCA = HumanPrimaryCellAtlasData(),
                          hg_Blueprint = celldex::BlueprintEncodeData(),
                          hg_ImmCED = DatabaseImmuneCellExpressionData(),
                          hg_NovHemD = NovershternHematopoieticData(),
                          hg_MonacoImm = MonacoImmuneData(),
                          mm_ImmGen = ImmGenData(),
                          mm_RNAseq = MouseRNAseqData())
        
        sr_annot = SingleR(test = log2(mc@mc_fp), ref = annot_df, labels = annot_df$label.main)
        sr_supmc_nms = tapply(sr_annot$labels, cls, function(v) { nm = names(which.max(table(v))); sprintf("%s (%d/%d)", nm, sum(v == nm), length(v))})
        mc@color_key = data.frame(gene="", group=sprintf("C%d: %s", 1:col_by_cutree_k, sr_supmc_nms), color=col_pal[1:col_by_cutree_k])
      }
      
      scdb_add_mc(mc_id, mc)
    }
  } else {
    mc_colorize_sup_hierarchy(mc_id=mc_id,
                              supmc = mc_sup,
                              supmc_key = supmc_file,
                              gene_key= marks_file)
  }
  
  # generate metacell clusters heatmap
  mcell_mc_plot_hierarchy(mc_id=mc_id,
                          graph_id=graph_id,
                          mc_order=mc_hc$order,
                          sup_mc = mc_sup,
                          width=width, height=height, min_nmc=min_nmc, show_mc_ids=show_mc_ids)
  
  list(mc_hc=mc_hc, mc_sup=mc_sup)
}

####
# plots wrapper
meta_mc_plots = function(mc_id, mc2d_id=mc_id, graph_id=mc_id, mat_id=mc_id, lateral_gset_id=NULL, metadata_fields_to_export=NULL, mc_ord=NULL, plot_2d=T, color_by_conf=F, col_by_cutree_k=0, T_gap=0.04, max_gene_len=8, add_metadata=NULL, fp_T_fold=2, min_confu_nmc=2, annot_ref='hg_Blueprint')
{
  set.seed(42)
  
  mc = scdb_mc(mc_id)
  if (!color_by_conf && length(mc@colors) == 1 && mc@colors == "white") {
    mc_colorize_default(mc_id)
  }
  
  conf = colorize_by_confusion_mat(mc_id=mc_id, graph_id=graph_id, show_mc_ids=T, col_by_cutree_k=col_by_cutree_k, T_gap=T_gap, min_nmc=min_confu_nmc, annot_ref=annot_ref)
  if (color_by_conf) {
    supmc_fn = sprintf("config/%s_supmc.txt", mc_id)
    stopifnot(file.exists(supmc_fn))
    
    marks_fn = sprintf("config/%s_marks.txt", mc_id)
    
    conf = colorize_by_confusion_mat(mc_id=mc_id, graph_id=graph_id, res=conf, show_mc_ids=T, supmc_file=supmc_fn, marks_file=if (file.exists(marks_fn)) { marks_fn } else { NULL}, annot_ref=annot_ref)	
  }
  
  mc_marks_gset_id = paste0(mc_id, "_mc_marks")
  
  if (is.null(lateral_gset_id)) {
    mcell_gset_from_mc_markers(gset_id=mc_marks_gset_id, mc_id=mc_id)
    
    for (plot_cells in c(T,F)) {
      mcell_mc_plot_marks(mc_id=mc_id, gset_id=mc_marks_gset_id, mat_id=mat_id, mc_ord=mc_ord, plot_cells=plot_cells, max_gene_len=max_gene_len, add_metadata=if (plot_cells) { add_metadata } else { NULL })
    }
  } else {
    mcell_gset_from_mc_markers(gset_id=mc_marks_gset_id, mc_id=mc_id, blacklist_gset_id=lateral_gset_id)
    mcell_gset_from_mc_markers(gset_id=paste0(mc_marks_gset_id, "_lateral"), mc_id=mc_id, filt_gset_id=lateral_gset_id)
    
    for (plot_cells in c(T,F)) {
      mcell_mc_plot_marks(mc_id=mc_id, gset_id=mc_marks_gset_id, mat_id=mat_id, lateral_gset_id=paste0(mc_marks_gset_id, "_lateral"), mc_ord=mc_ord, plot_cells=plot_cells, max_gene_len=max_gene_len, add_metadata=if (plot_cells) { add_metadata } else { NULL })
    }
  }
  
  if (plot_2d) {
    mcell_mc2d_force_knn(mc2d_id, mc_id, graph_id)
    mcell_mc2d_plot(mc2d_id=mc2d_id)	
    
    if (!is.null(metadata_fields_to_export)) {
      for (mf in metadata_fields_to_export) {
        message(sprintf("plot 2d by %s...", mf))
        for (single_plot in c(T, F)) {
          mcell_mc2d_plot_by_factor(mc2d_id, mat_id=mat_id, meta_field=mf, single_plot=single_plot)
        }
      }	
    }
  }
  if (!color_by_conf && col_by_cutree_k == 0) {
    # return original mc (erase default coloring if done)
    scdb_add_mc(mc_id, mc)
  }
  
  mcell_mc_export_tab(mc_id=mc_id, gstat_id=mat_id, mat_id=mat_id, T_fold=fp_T_fold, metadata_fields=metadata_fields_to_export)
  
  invisible(conf)
}

###
# common mat to mc and plots pipline
common_pipe = function(ds_nm, build_mat_func, max_f_mit=0.6, min_umis_pre_gene_ignore=500, min_umis_post_gene_ignore=500, gsets_to_filter=c("mito", "ig", "ncrna", "ncrp", "sno"), filt_mat_by_column=NULL,
                       gset_cor_thresh=0.1, mc_cells=NULL, mc_name=NULL, 
                       T_vm = NULL, T_tot = NULL, T_top3 = NULL, T_lfc = 3000, 
                       cgraph_knn = NULL, cgraph_downsamp=T, 
                       bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
                       mc_K=30, min_mc_size=30, mc_alpha=2, 
                       metadata_fields = NULL, color_by_conf = F, col_by_cutree_k=0,
                       gene_pref=NULL, selected_sets=NULL,
                       sample_field=NULL, specie="human", annot_dfs=c(human="hg_Blueprint", mouse="mm_RNAseq"),
                       rebuild=T, out_base=ds_nm, filt_lateral_genes=T, fp_T_fold=2, min_confu_nmc=2, seed=42,
                       split_mc_with_dbscan_and_filt_outliers=T) 
{
  set.seed(seed)
  
  all_id <<- sprintf("%s_raw", ds_nm)
  filt_id <<- sprintf("%s_filt", ds_nm)
  lateral_gset_id = if (filt_lateral_genes) sprintf("%s_lateral", ds_nm) else NULL
  dir.create(ds_nm, showWarnings = F)
  
  if (!is.null(out_base)) {
    rl(scdb_dir=paste0(out_base, "/scrna_db"), scfigs_dir=paste0(out_base, "/figs"), config_fn=sprintf("config/%s.yaml", ds_nm))
  }
  
  # build raw mat
  if (!is.null(build_mat_func) & (rebuild | !(paste0("mat.", all_id) %in% scdb_ls("mat")))) {
    message(sprintf("%s: building raw mat...", ds_nm))
    do.call(build_mat_func, args=list(mat_id=all_id))
    scdb_init(.scdb_base, T)
  }
  
  if (!scdb_obj_exists("gstat", all_id)) {
    mcell_add_gene_stat(all_id, all_id)
  }
  min_umis_cutoff = mcell_plot_umis_per_cell(all_id)
  
  # filter cells and genes (similar to melanoma - mito, IG, ncRNA)
  if (rebuild | !(paste0("mat.", filt_id) %in% scdb_ls("mat"))) {
    message(sprintf("%s: filtering mat...", ds_nm))
    meta_build_blacklist_gsets_by_gene_nms(all_id, ds_nm, gene_pref=gene_pref)	
    meta_build_blist_filtered_master_mat(all_id, filt_id, ds_nm, max_mito_f=max_f_mit, min_umis_pre_gene_ignore = min_umis_pre_gene_ignore, min_umis_post_gene_ignore = min_umis_post_gene_ignore, filt_mat_by_column=filt_mat_by_column, sample_field=sample_field, gsets_to_filter = gsets_to_filter)
  }
  
  # build lateral gene sets
  if (filt_lateral_genes && (rebuild || !(paste0("gset.", lateral_gset_id) %in% scdb_ls("gset")))) {
    message(sprintf("%s: building lateral gene sets...", ds_nm))
    meta_generate_lateral_gene_sets(filt_id, ds_nm, specie=specie, cor_thresh = gset_cor_thresh, gene_pref=gene_pref, selected_sets=selected_sets)
  }
  
  # partition to metacells
  mc_id = sprintf("%s%s", filt_id, ifelse(is.null(mc_name), "", paste0("_", mc_name)))
  new_mat_id = mc_id
  if (rebuild | !(paste0("mc.", mc_id) %in% scdb_ls("mc"))) {
    message(sprintf("%s: partitioning to metacells...", ds_nm))
    meta_mat2mc(filt_id, lateral_gset_id = lateral_gset_id, cells=mc_cells, name=mc_name, 
                T_vm = T_vm, T_tot = T_tot, T_top3 = T_top3, T_lfc = T_lfc, 
                cgraph_knn = cgraph_knn, cgraph_downsamp = cgraph_downsamp, 
                bootstrap_n_resamp = bootstrap_n_resamp, bootstrap_p_resamp = bootstrap_p_resamp,
                mc_K = mc_K, min_mc_size = min_mc_size, mc_alpha = mc_alpha, split_mc_with_dbscan_and_filt_outliers= split_mc_with_dbscan_and_filt_outliers,
                rebuild=rebuild)
  }
  
}



####
# Mark suspected doublet cells (currently only by Scrublet)
mcell_mat_filter_doublets = function(mat_id, new_mat_id=mat_id, method="Scrublet", break_by="amp_batch_id", break_by_cutoffs=NULL, cells=NULL, genes=NULL, name="", plot_hist=T, plot_umap=F, update_mat_metadata=T)
{
  mat = scdb_mat(mat_id)
  stopifnot(!is.null(mat))
  
  m = mat@mat
  if (!is.null(cells)) {
    stopifnot(all(cells %in% mat@cells))
    m = m[, cells]
  }
  
  if (!is.null(genes)) {
    stopifnot(all(genes %in% mat@genes))
    m = m[genes, ]
  }
  
  doub = c()
  doub_score = c()
  
  if (method == "Scrublet") {
    reticulate::use_condaenv(condaenv="r-scrublet", required=T)
    scrub = reticulate::import("scrublet")
    pyplot = reticulate::import("matplotlib.pyplot")
    
    if (is.null(break_by)) {
      ms = list(all=colnames(m))
    } else {
      stopifnot(break_by %in% colnames(mat@cell_metadata))
      ms = split(colnames(m), mat@cell_metadata[colnames(m), break_by])
    }
    
    for (nm in names(ms)) {
      message(sprintf("Finding doublets for %s...", nm))
      nms = ms[[nm]]
      
      res = scrub$Scrublet(t(m[, nms]))
      res2 = res$scrub_doublets()
      
      doub_score[nms] = res2[[1]]
      new_co = NULL
      if (is.null(break_by_cutoffs) || !(nm %in% names(break_by_cutoffs))) {
        doub[nms] = res2[[2]]
      } else {
        doub[nms] = res2[[1]] >= break_by_cutoffs[nm]
      }
      
      if (plot_hist) {
        res$plot_histogram()
        pyplot$savefig(scfigs_fn(mat_id, sprintf("scrublet_doublet_scores%s_%s", name, nm)))
      }
      if (plot_umap) {
        uc = scrub$get_umap(res$manifold_obs_)
        rownames(uc) = nms
        .plot_start(scfigs_fn(mat_id, sprintf("scrublet_doublet_umap%s_%s", name, nm)), 600, 600)
        plot(uc[,1], uc[,2], pch=19, cex=0.2, col='lightgray', main=sprintf("%s %s (%d/%d)", name, nm, sum(doub[nms]), length(nms)), xlab='umap1', ylab='umap2')
        cdoubs = nms[doub[nms]]
        points(uc[cdoubs,1], uc[cdoubs,2], pch=19, cex=0.2)
        dev.off()
      }
    }
  }
  
  if (update_mat_metadata) {
    mat@cell_metadata[names(doub), paste0("doublet_by_", method)] = doub
    mat@cell_metadata[names(doub), paste0(method, "_doublet_score")] = doub_score
    
    scdb_add_mat(mat_id, mat)
  }
  
  if (!is.null(new_mat_id)) {
    mcell_mat_ignore_cells(new_mat_id, mat_id, c(mat@ignore_cells, names(doub[doub])))
  }
  
  invisible(list(doub_score=doub_score, doub=doub))
}


###
#
mcell_mat_f_umi_by_gene_biotype = function(mat_id, genes_ifn, gene_pref="GRCh38_", biotypes=c("protein_coding", "lncRNA"))
{
  mat = scdb_mat(mat_id)
  
  gi = fread(genes_ifn) %>% filter(gene_type %in% biotypes)
  
  gs = intersect(paste0(gene_pref, mat@genes), names(which(table(gi$gene_name) == 1)))
  
  gi_s = gi %>% filter(gene_name %in% gs)
  v_biotype = names(which(table(gi_s$gene_type) > 1))
  gi_s = gi_s %>% filter(gene_type %in% v_biotype) %>% 
    tibble::column_to_rownames("gene_name")
  
  m = mat@mat[gsub(gene_pref, "", gs), ]
  
  f_umis = t(t(tgstat::tgs_matrix_tapply(t(m), gi_s[paste0(gene_pref, rownames(m)), 'gene_type'], sum)) / colSums(m))
  
  invisible(f_umis)
}

###
## inner helper to plot composition of a frequency matrix with metadata and legends  
composition_matrix_barplot_with_metadata = function(x, ofn, val2col, md_fields_cols, samp_md_df) {  
  y = length(md_fields_cols)
  
  .plot_start(ofn, 1000 + 30 * y, max(100 + 20 * nrow(x), 400))
  layout(matrix(seq(1, y+3), nrow=1), widths=c(600, rep(30, y), 200, 200))
  par(mar=c(8,24, 4, 0.5))
  barplot(t(x), col=val2col[colnames(x)], las=2, horiz=T, cex.names=1.2)
  
  par(mar=c(8, 0.5, 4, 0.5))
  samps = rownames(x)
  for (md_nm in names(md_fields_cols)) {
    barplot(rep(1, nrow(x)), horiz=T, col=unlist(md_fields_cols[[md_nm]][samp_md_df[samps, md_nm]]), xaxt='n', las=2)  
    axis(1, 0.5, md_nm, las=2, tick = F)
  }
  
  plot.new()
  par(mar=c(1,2,4,1))
  plot.window(0:1, 0:1)
  legend("topleft", legend=colnames(x), fill=val2col[colnames(x)], bty='n', cex=1.5)
  
  concat_leg = c()
  for (md_nm in names(md_fields_cols))  { 
    concat_leg = c(concat_leg, setNames(c(NA, NA, md_fields_cols[[md_nm]]), c("", paste("--", md_nm, "--"), names(md_fields_cols[[md_nm]])))) 
  }
  
  plot.new()
  par(mar=c(1,2,4,1))
  plot.window(0:1, 0:1)
  legend('topleft', legend=names(concat_leg), fill=concat_leg, border=ifelse(is.na(concat_leg), NA, 'black'), bty='n', cex=1.5)
  
  dev.off()
}

####
# mc group composition plots by metadata field (e.g. patient). Supports meta-group (e.g. lymphocytes (T, B, NK..) level plots and meta-group specific plots.
#' Title
#'
#' @param mc_id 
#' @param mat_id 
#' @param sample_field name of field in cell_metadata to break cells by
#' @param md_fields_cols metadata fields color dictionary to show beside the sample composition
#' @param sort_by How to sort samples. clust (defauly, hclust compositions), mdCols (by metadata fields) sampleName (by sample name)boolean. If true (not default) md_fields_cols will be used to sort samples 
#' @param groups Used both for filering groups of interest and to order them 
#' @param group2meta grouping of mc groups to meta-groups (will gen a barplot of meta-groups and one per meta-group)
#' @param meta_grp2col color dictionary for meta_groups
#' @param min_sample_cells min number of cells in a sample
#' @param name add to output file names
#'
#' @return
#' @export
#'
#' @examples
mcell_mc_plot_group_composition = function(mc_id, mat_id, sample_field, 
                                           md_fields_cols=NULL, sort_by='clust', groups=NULL, group2meta=NULL, meta_grp2col=NULL, min_sample_cells = 1, selected_samples=NULL, name=NULL, cells=NULL)
{
  mc = scdb_mc(mc_id)
  mat = scdb_mat(mat_id)
  stopifnot(!is.null(mc) & !is.null(mat))
  
  if (is.null(cells)) {
    cells = names(mc@mc)
  }  else {
    stopifnot(all(cells %in% intersect(names(mc@mc), mat@cells)))
  }
  
  md = mat@cell_metadata[cells, ]
  stopifnot(all(c(sample_field, names(md_fields_cols)) %in% colnames(md)))
  mdu = unique(md[, c(sample_field, names(md_fields_cols))]) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var=sample_field)
  
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  
  samp_grp = table(md[, sample_field], col2grp[mc@colors[mc@mc[cells]]])
  
  if (!is.null(groups)) {
    samp_grp = samp_grp[, intersect(groups, colnames(samp_grp))]
  }
  
  if (!is.null(selected_samples)) {
    samp_grp = samp_grp[selected_samples, ]
    mdu = mdu[rownames(samp_grp), , drop=F]
  }
  
  samp_grp = samp_grp[rowSums(samp_grp) >= min_sample_cells, ]
  samp_grp_n = samp_grp / rowSums(samp_grp)
  
  samp_grp_ln = log2((samp_grp + 1) / rowSums(samp_grp + 1))
  
  # start with clustering samples composition correlation
  cc = cor(t(samp_grp_ln))
  
  .plot_start(scfigs_fn(mc_id, sprintf("%s_sample_cor", name), scfigs_dir(mc_id, "composition")), 400 + 15 * nrow(cc), 300 + 15 * nrow(cc))
  pheatmap(cc, breaks=seq(-1, 1, len=101), treeheight_col=10, treeheight_row=10, cellwidth=15, cellheight=15, annotation_col=mdu, annotation_color=md_fields_cols)
  dev.off()
  
  
  if (sort_by == 'clust') {
    hc = hclust(dist(cor(t(samp_grp_n))), method='ward.D2')
    samp_grp_n = samp_grp_n[hc$order, ]
    mdu = mdu[rownames(samp_grp_n), , drop=F]
  } else if (sort_by == 'mdCols') {
    mdu = mdu %>% tibble::rownames_to_column(var=sample_field) %>% arrange_at(c(names(md_fields_cols), sample_field)) %>% tibble::column_to_rownames(var=sample_field)
    samp_grp_n = samp_grp_n[intersect(rownames(mdu), rownames(samp_grp_n)), ]
  } else if (sort_by == 'sampleName') {
    samp_grp_n = samp_grp_n[sort(rownames(samp_grp_n)), ]
  } else if (sort_by == 'selectedSamples') {
    samp_grp_n = samp_grp_n[selected_samples, ]
  } else {
    stop("sort_by should be either: clust, mdCols, sampleName")
  }
  
  composition_matrix_barplot_with_metadata(x = samp_grp_n, 
                                           ofn = scfigs_fn(mc_id, sprintf('mc_group_comp%s_%s', ifelse(is.null(name), "", paste0("_", name)), sort_by), scfigs_dir(mc_id, "composition")), 
                                           val2col = grp2col, 
                                           md_fields_cols = md_fields_cols, 
                                           samp_md_df = mdu)
  
  # meta groups plots
  if (!is.null(group2meta)) {
    stopifnot(!is.null(meta_grp2col))
    stopifnot(all(names(group2meta) %in% col2grp))
    
    samp_mg = table(md[, sample_field], group2meta[col2grp[mc@colors[mc@mc]]])
    samp_mg_n = samp_mg / rowSums(samp_mg)
    
    if (sort_by == 'clust') {
      hc = hclust(dist(cor(t(samp_mg_n))), method='ward.D2')
      samp_mg_n = samp_mg_n[hc$order, ]
      mdu = mdu[rownames(samp_mg_n), ]
    } else {
      samp_mg_n = samp_mg_n[rownames(samp_grp_n), ]
    }
    
    composition_matrix_barplot_with_metadata(x = samp_mg_n, 
                                             ofn = scfigs_fn(mc_id, sprintf('mc_meta_group_comp%s_%s', ifelse(is.null(name), "", paste0("_", name)), sort_by), scfigs_dir(mc_id, "composition")), 
                                             val2col = meta_grp2col, 
                                             md_fields_cols = md_fields_cols, 
                                             samp_md_df = mdu)
    
    mg_s = split(names(group2meta), group2meta)
    
    for (mg in names(mg_s)) {
      cgrps = mg_s[[mg]]
      if (length(cgrps) > 1) {
        c_samp_mg = samp_grp[rownames(samp_grp_n), cgrps]
        c_samp_mg_n = c_samp_mg / rowSums(c_samp_mg)
        
        composition_matrix_barplot_with_metadata(x = c_samp_mg_n, 
                                                 ofn = scfigs_fn(mc_id, sprintf('mc_meta_group_%s_comp%s_%s', mg, ifelse(is.null(name), "", paste0("_", name)), sort_by), scfigs_dir(mc_id, "composition")), 
                                                 val2col = grp2col, 
                                                 md_fields_cols = md_fields_cols, 
                                                 samp_md_df = mdu)
        
      }      
    }
  }
}

####
mcell_mc_plot_group_enr_heatmaps = function(mc_id, groups=NULL, grp_ord=NULL, group2meta=NULL, meta_grp2col=NULL, name=NULL, zlim=3, n_glob_enr=5, n_grp_enr=5, n_grp_outliers=5, min_lfp_glob_enr=0.8, gene_width=15, mc_height=10, space_groups=T, export_plts_to_ppt=T, ppt_n_panels_in_row=4, show_marker_type=T)
{
  mc = scdb_mc(mc_id)
  stopifnot(!is.null(mc))
  lfp = log2(mc@mc_fp)
  lfp = lfp[apply(lfp, 1, max) >= min_lfp_glob_enr, ]
  
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  
  if (is.null(groups)) {
    groups = names(grp2col)
  }
  
  ## inner helper to plot group enrichment 
  enr_hm_helper = function(lfp, index, name, val2col, foc_groups, grp_ord=NULL) {
    stopifnot(is.null(grp_ord) || all(grp_ord %in% foc_groups))
    f = index %in% foc_groups
    fp_f = 2 ** lfp[, f]
    lfp_f = log2( fp_f / rowMeans(fp_f))
    
    x =  t(tgstat::tgs_matrix_tapply(lfp[,f], index[f], mean))
    glob_enr_gs = unique(as.vector(apply(x, 2, function(v) { v = v[v > 0]; if (length(v) > 0) { names(tail(sort(v), n_glob_enr)) }})))
    grp_outlier_gs =  unique(unlist(sapply(foc_groups, function(foc) { foc_mcs = which(index == foc); y = sort(apply(lfp[, foc_mcs, drop=F], 1, min) - apply(lfp[, setdiff(1:ncol(lfp), foc_mcs), drop=F], 1, max)); y = y[y > 0]; names(tail(y, n_grp_outliers)) } )))
    
    x_f =  t(tgstat::tgs_matrix_tapply(lfp_f, index[f], mean))
    grp_enr_gs = unique(as.vector(apply(x_f, 2, function(v) { v = v[v > 0]; if (length(v) > 0) { names(tail(sort(v), n_grp_enr)) }})))
    
    
    #gs_ann = table(c(glob_enr_gs, grp_outlier_gs, grp_enr_gs), c(rep('glob_enr', length(glob_enr_gs)), rep('grp_enr', length(grp_enr_gs)), rep('grp_outlier', length(grp_outlier_gs))))
    #gs_ann = ifelse(gs_ann > 0, 'in', 'out')
    
    foc_gs = unique(c(glob_enr_gs, grp_outlier_gs, grp_enr_gs))
    gs_ann = data.frame(row.names=foc_gs, glob_enr=ifelse(foc_gs %in% glob_enr_gs, 'yes', 'no'), grp_enr=ifelse(foc_gs %in% grp_enr_gs, 'yes', 'no'), grp_outlier=ifelse(foc_gs %in% grp_outlier_gs, 'yes', 'no'))
    x = x[foc_gs, ]
    
    gr_hc = hclust(dist(t(x)), method='ward.D2')
    grp_levels = colnames(x)[gr_hc$order]
    
    gs_hc = hclust(dist(cor(t(x))), method='ward.D2')
    gs_dd = as.dendrogram(gs_hc)
    gs_dd_reord = reorder(gs_dd, apply(x[, gr_hc$order], 1, which.max))
    gs_ord = as.hclust(gs_dd_reord)$order
    
    #gs_ord = order(apply(x[, gr_hc$order], 1, which.max) + 1e-3 * order(gs_hc$order))
    
    
    if (!is.null(grp_ord)) {
      x = x[, grp_ord]
    }
    
    val2col = val2col[colnames(x)]
    grp_ann = data.frame(row.names=names(val2col), grp=names(val2col))
    
    .plot_start(scfigs_fn(mc_id, sprintf("%s_grp_avg%s", name, ifelse(is.null(grp_ord), '_clust', '')), scfigs_dir(mc_id, "enr_heatmaps")), nrow(x) * gene_width + 400, ncol(x) * mc_height + 200)
    #hm = pheatmap(pmin(pmax(t(x), -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_rows=is.null(grp_ord), annotation_row=grp_ann, cellwidth=10, cellheight = 20, annotation_colors = list(grp=val2col[rownames(grp_ann)]), treeheight_row=10, treeheight_col=10, main=name)
    gs_cols = list()
    for (nm in colnames(gs_ann)) {
      gs_cols[[nm]] = c(yes='darkgray', no='white')
    }
    
    if (!show_marker_type) {
      gs_ann = NULL
    }  
    
    grp_levels = grp_ord
    if (is.null(grp_levels)) {
      grp_levels = colnames(x)[gr_hc$order]
    }
    
    gr_ord = colnames(x)[gr_hc$order]
    if (!is.null(grp_ord)) {
      gr_ord = grp_ord
    }
    
    hm = pheatmap(pmin(pmax(t(x[gs_ord, gr_ord]), -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_rows=F, cluster_cols=T, annotation_row=grp_ann, annotation_col=gs_ann, cellwidth=gene_width, cellheight = mc_height, annotation_colors = c(list(grp=val2col[rownames(grp_ann)]), gs_cols), treeheight_row=10, treeheight_col=10, main=name)
    dev.off()
    write.csv(x[gs_ord, gr_ord], scfigs_fn(mc_id, sprintf("%s_grp_avg%s", name, ifelse(is.null(grp_ord), '_clust', '')), scfigs_dir(mc_id, "enr_heatmaps"), ext='csv'), quote=F)
    
    
    lfp_d = lfp[foc_gs, f]
    hc = hclust(dist(t(lfp_d)), method='ward.D2')
    lfp_d_ord = order(as.numeric(factor(index[f], levels=grp_levels)) + 1e-3 * order(hc$order))
    lfp_d = lfp_d[, lfp_d_ord]
    
    grp_ann = data.frame(row.names=colnames(lfp_d), grp=index[as.numeric(colnames(lfp_d))])
    gaps_mcs = NULL
    if (space_groups) {
      gaps_mcs = cumsum(rle(grp_ann$grp)$lengths)
    }
    .plot_start(scfigs_fn(mc_id, sprintf("%s_mc%s", name, ifelse(is.null(grp_ord), '_clust', '')), scfigs_dir(mc_id, "enr_heatmaps")), nrow(lfp_d) * gene_width + 400, ncol(lfp_d) * mc_height + 200)
    hm = pheatmap(pmin(pmax(t(lfp_d), -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_rows=F, annotation_row=grp_ann, cellwidth=gene_width, cellheight = mc_height, annotation_colors = list(grp=val2col[grp_levels]), treeheight_col=10, main=name, gaps_row=gaps_mcs)
    dev.off()
    write.csv(lfp_d[hm$tree_col$order, ], scfigs_fn(mc_id, sprintf("%s_mc%s", name, ifelse(is.null(grp_ord), '_clust', '')), scfigs_dir(mc_id, "enr_heatmaps"), ext='csv'), quote=F)
    
    if (export_plts_to_ppt) {
      bw = 1.95
      sp = 0
      title_h = 0.52
      sub_title_h = 0.3
      pres = read_pptx()
      
      x =  t(tgstat::tgs_matrix_tapply(lfp[,f], index[f], mean))
      x_f =  t(tgstat::tgs_matrix_tapply(lfp_f, index[f], mean))
      
      ppt_odir = sprintf("%s/ppt_figs", scfigs_dir(mc_id, "enr_heatmaps"))
      dir.create(ppt_odir, showWarnings = F)
      
      legend_ofn = sprintf("%s/%s_legend.pdf", ppt_odir, name)
      leg_dim = mcell_mc_plot_legend(mc_id, ofn=legend_ofn, ncols=1, filter_names = foc_groups)
      
      for (nm in foc_groups) {
        
        g_df = NULL
        v = x[, nm] 
        v = v[v > 0]
        if (length(v) > 0) { 
          gs = names(tail(sort(v), ppt_n_panels_in_row * 2))
          g_df = data.frame(type=rep('glob', length(gs)), supmc=rep(nm, length(gs)), gene=rev(gs)) 
        }
        
        v = x_f[setdiff(rownames(x_f), g_df$gene), nm] 
        v = v[v > 0]
        if (length(v) > 0) { 
          gs = names(tail(sort(v), ppt_n_panels_in_row * 2))
          g_df = rbind(g_df, data.frame(type=rep('grp', length(gs)), supmc=rep(nm, length(gs)), gene=rev(gs)))
        }
        
        foc_mcs = which(index == nm)
        foc_gs = setdiff(names(which(apply(lfp, 1, max) > 0.5)), g_df$gene)
        y = sort(apply(lfp[foc_gs, foc_mcs, drop=F], 1, min) - apply(lfp[foc_gs, setdiff(1:ncol(lfp), foc_mcs), drop=F], 1, max))
        y = y[y > 0]
        if (length(y) > 0) {
          gs = names(tail(y, ppt_n_panels_in_row * 2))
          g_df = rbind(g_df, data.frame(type=rep('outlier', length(gs)), supmc=rep(nm, length(gs)), gene=rev(gs)))
        }
        
        message(paste("Adding slide for", nm))
        pres = add_slide(pres)
        pres = ph_with(pres, value=fpar(ftext(sprintf("%s (%s)", nm, name), fp_text(font.size = 26)), fp_p = fp_par(text.align = "center")), location=ph_location(sp, 0, 10, title_h))
        c_y = title_h
        leg_h = leg_dim[2] * 0.75
        leg_w = leg_dim[1] * 0.75
        pres = ph_with(pres, external_img(legend_ofn, leg_w, leg_h), location = ph_location(4 * (bw + sp), c_y, leg_w, leg_h))
        types = c('glob', 'grp', 'outlier')
        
        for (i in seq_along(types)) {
          pres = ph_with(pres, value=fpar(ftext(types[i], fp_text(font.size = 16))), location=ph_location(sp, c_y, bw , sub_title_h))
          #c_x = 0.3 + 2 * sp
          c_x = sp
          c_y = c_y + sub_title_h
          c_df = filter(g_df, type == types[i] & supmc == nm)
          if (nrow(c_df) > 0) {
            if (nrow(c_df) %% 2 == 1) {
              c_df = rbind(c_df, c_df[1,])
            }
            for (j in seq(1, nrow(c_df), by=2)) {
              c_ofn = sprintf("%s/%s_%s_%s_vs_%s.pdf", ppt_odir, nm, types[i], c_df[j, 'gene'], c_df[j+1, 'gene'])
              plt(c_df[j, 'gene'], c_df[j+1, 'gene'], lfp[,f], mc@colors[f], ofn=c_ofn, cex.lab=1.5)
              pres = ph_with(pres, external_img(c_ofn, bw, bw), location = ph_location(c_x, c_y, bw, bw))
              c_x = c_x + bw + sp
            }
          }
          c_y = c_y + bw
        }
      }
      print(pres, target=scfigs_fn(mc_id, sprintf("%s_gene_markers", name), scfigs_dir(mc_id, "enr_heatmaps"), ext='pptx'))
    }
    
  }
  
  enr_hm_helper(lfp, col2grp[mc@colors], sprintf('mc_group_enr%s', ifelse(is.null(name), "", paste0("_", name))), grp2col, foc_groups=groups, grp_ord=grp_ord)
  
  
  # meta groups plots
  if (!is.null(group2meta)) {
    stop("Not supported, requires update!")
    stopifnot(!is.null(meta_grp2col))
    stopifnot(all(names(group2meta) %in% col2grp))
    
    enr_hm_helper(lfp, group2meta[col2grp[mc@colors]], sprintf('mc_meta_group_enr%s', ifelse(is.null(name), "", paste0("_", name))), meta_grp2col, n_genes_from_each=16)
    
    mg_s = split(names(group2meta), group2meta)
    
    grps = col2grp[mc@colors]
    for (mg in names(mg_s)) {
      cgrps = mg_s[[mg]]
      ind = grps %in% cgrps
      if (length(cgrps) > 1) {
        enr_hm_helper(lfp[, ind], grps[ind], sprintf('mc_meta_group_%s_enr%s', mg, ifelse(is.null(name), "", paste0("_", name))), grp2col, n_genes_from_each=16)
      }      
    }
  }
}

####
mcell_mc2d_plot_group_subset = function(mc_id, groups, name, graph_id=mc_id, mc2d_id=NULL, show_mc_ids=T)
{
  mc = scdb_mc(mc_id)
  stopifnot(!is.null(mc))
  
  col2grp = get_mc_col2group(mc)
  mc_grps = col2grp[mc@colors]
  stopifnot(all(groups %in% mc_grps))
  
  if (is.null(mc2d_id)) {
    mc2d_id = paste(mc_id, name, sep="_")
  }
  sub_mc = mc_set_outlier_mc(mc, which(! mc_grps %in% groups))
  scdb_add_mc(mc2d_id, sub_mc)
  
  mcell_mc2d_force_knn(mc2d_id, mc2d_id, graph_id)
  mcell_mc2d_plot(mc2d_id, show_mc_ids=show_mc_ids)	
  
  mcell_mc_plot_legend(mc2d_id, name=name)
  
}

###
# Plot mc/supmc/meta-supmc composition by field (e.g. sample) with metadata fields annotation. Also: z-score on field + breakdown to md fields, Also: Inverse simpson score by (sample) field.
mcell_mc_plot_mc_composition_by_metadata_field = function(mc_id, break_by_field, md_fields_cols, mat_id=mc_id, grp_ord=NULL, group2meta=NULL, meta_grp2col=NULL, max_z=100, z_pval_breaks=c(1e-2, 1e-3, 1e-4), selected_samples=NULL)
{
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  stopifnot(!any(sapply(list(mat, mc), is.null)))
  stopifnot(is.null(md_fields_cols) || all(c(break_by_field, names(md_fields_cols)) %in% colnames(mat@cell_metadata)))
  
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  if (is.null(grp_ord)) {
    grp_ord = sort(unique(col2grp))
  }
  
  mc_ann = get_mc_pheatmap_ann(mc, 'group')
  
  cells = names(mc@mc)
  if (!is.null(selected_samples)) {
    samp2cells = split(cells, mat@cell_metadata[cells, break_by_field])
    cells = unlist(samp2cells[selected_samples])
  }
  
  ann_cols = list()
  md_ann = NULL
  if (!is.null(md_fields_cols)) {
    md_ann = unique(mat@cell_metadata[cells, c(break_by_field, names(md_fields_cols))]) %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var=break_by_field)
    ann_cols = md_fields_cols
  }
  ann_cols[['group']] = grp2col[grp_ord]
  
  mc_samp = table(mc@mc[cells], mat@cell_metadata[cells, break_by_field])
  if (ncol(mc_samp) > 1) {
    mc_samp_n = mc_samp / rowSums(mc_samp)
    
    mc_hc = hclust(dist(mc_samp_n), method='ward.D2')
    mc_ord = order(as.numeric(factor(col2grp[mc@colors[as.numeric(rownames(mc_samp_n))]], levels=grp_ord)) + 1e-3 * order(mc_hc$order))
    
    .plot_start(scfigs_fn(mc_id, sprintf("mc_%s_comp_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(mc_samp_n) * 4 + 800, ncol(mc_samp_n) * 20 + 800)
    pheatmap(t(mc_samp_n), cluster_rows=ncol(mc_samp_n) > 2, cluster_cols = nrow(mc_samp_n) > 2, color=colorRampPalette(c('white', 'darkred'))(100), annotation_col=mc_ann$ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=4, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    
    .plot_start(scfigs_fn(mc_id, sprintf("mc_%s_comp_heatmap_bySupmc", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(mc_samp_n) * 4 + 800, ncol(mc_samp_n) * 20 + 800)
    pheatmap(t(mc_samp_n[mc_ord, ]), cluster_cols=F, cluster_rows=ncol(mc_samp_n) > 2, color=colorRampPalette(c('white', 'darkred'))(100), annotation_col=mc_ann$ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=4, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    write.csv(mc_samp_n[mc_ord, ], scfigs_fn(mc_id, sprintf("mc_%s_comp_heatmap_bySupmc", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    
    e_mc_samp = (rowSums(mc_samp) %*% t(colSums(mc_samp))) / sum(mc_samp)
    z_mc_samp = (mc_samp - e_mc_samp) / sqrt(e_mc_samp + 1)
    mc_hc = hclust(dist(z_mc_samp), method='ward.D2')
    mc_ord = order(as.numeric(factor(col2grp[mc@colors[as.numeric(rownames(z_mc_samp))]], levels=grp_ord)) + 1e-3 * order(mc_hc$order))
    
    .plot_start(scfigs_fn(mc_id, sprintf("mc_%s_Z_heatmap_bySupmc", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_mc_samp) * 4 + 800, ncol(z_mc_samp) * 20 + 800)
    pheatmap(pmin(pmax(t(z_mc_samp[mc_ord, ]), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), cluster_cols=F, cluster_rows=ncol(z_mc_samp) > 2, annotation_col=mc_ann$ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=4, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    write.csv(z_mc_samp[mc_ord, ], scfigs_fn(mc_id, sprintf("mc_%s_Z_heatmap_bySupmc", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    
    md_gap = NULL
    if (!is.null(md_fields_cols)) {
      z_mc_md = NULL
      md_val_ann = NULL
      md_val_cols = list(group=mc_ann$ann_col, val=c())
      for (md_nm in names(md_fields_cols)) {
        c_mc_md = table(mc@mc[cells], mat@cell_metadata[cells, md_nm])
        colnames(c_mc_md) = paste(md_nm, colnames(c_mc_md), sep='-')
        
        e_mc_md = (rowSums(c_mc_md) %*% t(colSums(c_mc_md))) / sum(c_mc_md)
        z_mc_md = cbind(z_mc_md, (c_mc_md - e_mc_md) / sqrt(e_mc_md + 1))
        
        md_val_ann = rbind(md_val_ann, data.frame(row.names=colnames(c_mc_md), val=colnames(c_mc_md)))
        
        c_md_cols = md_fields_cols[[md_nm]]
        names(c_md_cols) = paste(md_nm, names(c_md_cols), sep='-')
        md_val_cols[['val']] = c(md_val_cols[['val']], c_md_cols)
        md_gap = c(md_gap, ncol(c_mc_md))
      }
      md_gap = cumsum(md_gap)
      
      mc_hc = hclust(dist(z_mc_md), method='ward.D2')
      mc_ord = order(as.numeric(factor(col2grp[mc@colors[as.numeric(rownames(z_mc_md))]], levels=grp_ord)) + 1e-3 * order(mc_hc$order))
      
      .plot_start(scfigs_fn(mc_id, "mc_md_Z_heatmap_bySupmc", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_mc_md) * 4 + 800, ncol(z_mc_md) * 20 + 800)
      pheatmap(pmin(pmax(t(z_mc_md[mc_ord, ]), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), cluster_cols=F, cluster_rows=F, annotation_col=mc_ann$ann, annotation_row=md_val_ann, annotation_colors=md_val_cols, cellwidth=4, cellheight=20, gaps_row=md_gap)
      dev.off()
      write.csv(z_mc_md[mc_ord, ], scfigs_fn(mc_id, "mc_md_Z_heatmap_bySupmc", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    }
    
    # supmc  
    smc_samp = table(col2grp[mc@colors[mc@mc[cells]]], mat@cell_metadata[cells, break_by_field])
    smc_samp_n = smc_samp / rowSums(smc_samp)
    
    grp_ann = data.frame(row.names=grp_ord, group=grp_ord)
    .plot_start(scfigs_fn(mc_id, sprintf("smc_%s_comp_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(smc_samp_n) * 20 + 800, ncol(smc_samp_n) * 20 + 800)
    pheatmap(t(smc_samp_n), color=colorRampPalette(c('white', 'darkred'))(100), annotation_col=grp_ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=20, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    
    .plot_start(scfigs_fn(mc_id, sprintf("smc_%s_comp_heatmap", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(smc_samp_n) * 20 + 800, ncol(smc_samp_n) * 20 + 800)
    pheatmap(t(smc_samp_n[grp_ord, ]), cluster_cols=F, color=colorRampPalette(c('white', 'darkred'))(100), annotation_col=grp_ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=20, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    write.csv(smc_samp_n[grp_ord, ], scfigs_fn(mc_id, sprintf("smc_%s_comp_heatmap", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    
    e_smc_samp = (rowSums(smc_samp) %*% t(colSums(smc_samp))) / sum(smc_samp)
    z_smc_samp = (smc_samp - e_smc_samp) / sqrt(e_smc_samp + 1)
    
    .plot_start(scfigs_fn(mc_id, sprintf("smc_%s_Z_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_smc_samp) * 20 + 800, ncol(z_smc_samp) * 20 + 800)
    pheatmap(pmin(pmax(t(z_smc_samp), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), annotation_col=grp_ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=20, cellheight=20, treeheight_col=10, treeheight_row=10)
    dev.off()
    write.csv(z_smc_samp[grp_ord, ], scfigs_fn(mc_id, sprintf("smc_%s_Z_heatmap", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    
    if (!is.null(md_fields_cols)) {
      z_smc_md = NULL
      md_val_cols[['group']] = grp2col
      for (md_nm in names(md_fields_cols)) {
        c_smc_md = table(col2grp[mc@colors[mc@mc[cells]]], mat@cell_metadata[cells, md_nm])
        colnames(c_smc_md) = paste(md_nm, colnames(c_smc_md), sep='-')
        
        e_smc_md = (rowSums(c_smc_md) %*% t(colSums(c_smc_md))) / sum(c_smc_md)
        z_smc_md = cbind(z_smc_md, (c_smc_md - e_smc_md) / sqrt(e_smc_md + 1))
      }
      
      p_vals = matrix(p.adjust(pnorm(abs(z_smc_md), lower.tail = F), method='fdr'), nrow=nrow(z_smc_md))
      p_vals_txt = F
      if (!is.null(z_pval_breaks)) {
        p_vals_txt = matrix("", nrow(p_vals), ncol(p_vals))
        for (i in seq_along(z_pval_breaks)) {
          p_vals_txt[p_vals <= z_pval_breaks[i]] = paste0(rep("*", i), collapse="")
        }
        p_vals_txt = t(p_vals_txt)
      }
      
      .plot_start(scfigs_fn(mc_id, "smc_md_Z_heatmap_clust", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_smc_md) * 20 + 800, ncol(z_smc_md) * 20 + 800)
      pheatmap(pmin(pmax(t(z_smc_md), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), cluster_rows=F, annotation_col=grp_ann, annotation_row=md_val_ann, annotation_colors=md_val_cols, cellwidth=20, cellheight=20, gaps_row=md_gap, display_numbers=p_vals_txt)
      dev.off()
      write.csv(z_smc_md, scfigs_fn(mc_id, "smc_md_Z_heatmap_clust", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
    }
    
    if (!is.null(group2meta)) {
      smc_samp = table(group2meta[col2grp[mc@colors[mc@mc[cells]]]], mat@cell_metadata[cells, break_by_field])
      smc_samp_n = smc_samp / rowSums(smc_samp)
      
      meta_ann = data.frame(row.names=names(meta_grp2col), group=names(meta_grp2col))
      ann_cols[['group']] = meta_grp2col
      .plot_start(scfigs_fn(mc_id, sprintf("meta_smc_%s_comp_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(smc_samp_n) * 20 + 800, ncol(smc_samp_n) * 20 + 800)
      pheatmap(t(smc_samp_n), color=colorRampPalette(c('white', 'darkred'))(100), annotation_col=meta_ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=20, cellheight=20, treeheight_col=10, treeheight_row=10)
      dev.off()
      write.csv(smc_samp_n, scfigs_fn(mc_id, sprintf("meta_smc_%s_comp_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F )
      
      e_smc_samp = (rowSums(smc_samp) %*% t(colSums(smc_samp))) / sum(smc_samp)
      z_smc_samp = (smc_samp - e_smc_samp) / sqrt(e_smc_samp + 1)
      
      .plot_start(scfigs_fn(mc_id, sprintf("meta_smc_%s_Z_heatmap_clust", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_smc_samp) * 20 + 800, ncol(z_smc_samp) * 20 + 800)
      pheatmap(pmin(pmax(t(z_smc_samp), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), annotation_col=meta_ann, annotation_row=md_ann, annotation_colors=ann_cols, cellwidth=20, cellheight=20, treeheight_col=10, treeheight_row=10)
      dev.off()
      write.csv(z_smc_samp, scfigs_fn(mc_id, sprintf("meta_smc_%s_Z_heatmap", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
      
      if (!is.null(md_fields_cols)) {
        z_smc_md = NULL
        md_val_cols[['group']] = meta_grp2col
        for (md_nm in names(md_fields_cols)) {
          c_smc_md = table(group2meta[col2grp[mc@colors[mc@mc[cells]]]], mat@cell_metadata[cells, md_nm])
          e_smc_md = (rowSums(c_smc_md) %*% t(colSums(c_smc_md))) / sum(c_smc_md)
          z_smc_md = cbind(z_smc_md, (c_smc_md - e_smc_md) / sqrt(e_smc_md + 1))
        }
        
        p_vals = matrix(p.adjust(pnorm(abs(z_smc_md), lower.tail = F), method='fdr'), nrow=nrow(z_smc_md))
        p_vals_txt = F
        if (!is.null(z_pval_breaks)) {
          p_vals_txt = matrix("", nrow(p_vals), ncol(p_vals))
          for (i in seq_along(z_pval_breaks)) {
            p_vals_txt[p_vals <= z_pval_breaks[i]] = paste0(rep("*", i), collapse="")
          }
          p_vals_txt = t(p_vals_txt)
        }
        
        .plot_start(scfigs_fn(mc_id, "meta_smc_md_Z_heatmap_clust", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), nrow(z_smc_md) * 20 + 800, ncol(z_smc_md) * 20 + 800)
        pheatmap(pmin(pmax(t(z_smc_md), -max_z), max_z), breaks=seq(-max_z, max_z, len=101), cluster_rows=F, annotation_col=meta_ann, annotation_row=md_val_ann, annotation_colors=md_val_cols, cellwidth=20, cellheight=20, gaps_row=md_gap, display_numbers=p_vals_txt)
        dev.off()
        write.csv(z_smc_md, scfigs_fn(mc_id, "meta_smc_md_Z_heatmap", dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)
      }
      
    } 
    
    # plot Simpson index over mc's
    simp = apply(mc_samp, 1, function(v) { 1 / sum( (v / sum(v))^2) })
    simp_df = data.frame(row.names=rownames(mc_samp), group=col2grp[mc@colors], inv_simp=simp)
    
    smc = factor(col2grp[mc@colors], levels=grp_ord)
    smc2 = factor(col2grp[mc@colors], levels=names(sort(tapply(simp, col2grp[mc@colors], median))))
    
    .plot_start(scfigs_fn(mc_id, sprintf("mc_Simpson_by_%s", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), max(400, 200 + 20 * length(col2grp)), 600)
    par(mar=c(20,4,4,1))
    boxplot(simp ~ smc, outline=F, main=sprintf("MCs by %s", break_by_field), las=2, col=grp2col[grp_ord], xlab='', ylab='Inverse Simpson')
    abline(h=c(1, ncol(mc_samp)), lty=2, col='darkgrey')
    stripchart(simp ~ smc, method='jitter', jitter=0.2, pch=19, vertical=T, add=T, cex=0.3)
    dev.off()
    
    .plot_start(scfigs_fn(mc_id, sprintf("mc_Simpson_by_%s_ord", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), max(400, 200 + 20 * length(col2grp)), 600)
    par(mar=c(20,4,4,1))
    boxplot(simp ~ smc2, outline=F, main=sprintf("MCs by %s", break_by_field), las=2, col=grp2col[levels(smc2)], xlab='', ylab='Inverse Simpson')
    abline(h=c(1, ncol(mc_samp)), lty=2, col='darkgrey')
    stripchart(simp ~ smc2, method='jitter', jitter=0.2, pch=19, vertical=T, add=T, cex=0.3)
    dev.off()
    
    if (!is.null(group2meta)) {
      
      smc2 = factor(group2meta[col2grp[mc@colors]], levels=names(sort(tapply(simp, group2meta[col2grp[mc@colors]], median))))
      
      .plot_start(scfigs_fn(mc_id, sprintf("meta_mc_Simpson_by_%s_ord", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field))), max(400, 200 + 20 * length(meta_grp2col)), 600)
      par(mar=c(20,4,4,1))
      boxplot(simp ~ smc2, outline=F, main=sprintf("MCs by %s", break_by_field), las=2, col=meta_grp2col[levels(smc2)], xlab='', ylab='Inverse Simpson')
      abline(h=c(1, ncol(mc_samp)), lty=2, col='darkgrey')
      stripchart(simp ~ smc2, method='jitter', jitter=0.2, pch=19, vertical=T, add=T, cex=0.3)
      dev.off()
      
      simp_df$meta_group = group2meta[col2grp[mc@colors]]
    } 
    
    write.csv(simp_df, scfigs_fn(mc_id, sprintf("mc_Simpson_by_%s", break_by_field), dir=scfigs_dir(mc_id, sprintf("composition_by_%s", break_by_field)), ext='csv'), quote=F)  
  }
}

###
# Common plots over an annotated mc object
mcell_common_plots = function(mc_id, sample_field, md_fields_cols, 
                              mat_id = mc_id, graph_id = mc_id, 
                              grp_ord = NULL, group2meta=NULL, meta_groups_cols=NULL, meta_groups_mc2d_show_mc_ids=F,
                              comp_sort_by='clust', comp_min_sample_cells = 1,
                              hm_n_glob_enr=5, hm_n_grp_enr=5, hm_n_grp_outliers=5, hm_show_marker_type = T, hm_zlim=4,
                              selected_samples=NULL,
                              n_legend_cols=1)
{
  
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  stopifnot(!any(sapply(list(mat, mc), is.null)))
  stopifnot(!is.null(md_fields_cols) && all(c(sample_field, names(md_fields_cols)) %in% colnames(mat@cell_metadata)))
  
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  
  if (is.null(grp_ord)) {
    grp_ord = sort(unique(col2grp))
  }
  
  # sample composition
  mcell_mc_plot_group_composition(mc_id, mat_id, sample_field, md_fields_cols=md_fields_cols, sort_by=comp_sort_by, groups=grp_ord, min_sample_cells = comp_min_sample_cells, name="all", selected_samples=selected_samples)
  mcell_mc_plot_legend(mc_id, ncols=n_legend_cols, names_ord=grp_ord)
  
  # mc and supmc composition by sample
  mcell_mc_plot_mc_composition_by_metadata_field(mc_id=mc_id, break_by_field=sample_field, md_fields_cols=md_fields_cols, mat_id=mat_id, grp_ord=grp_ord, group2meta = group2meta, meta_grp2col = meta_groups_cols, selected_samples=selected_samples)
  for (i in seq_along(md_fields_cols)) {
    mcell_mc_plot_mc_composition_by_metadata_field(mc_id=mc_id, break_by_field=names(md_fields_cols)[i], md_fields_cols=NULL, mat_id=mat_id, grp_ord=grp_ord)
  }  
  
  # Enriched gene heatmaps
  mcell_mc_plot_group_enr_heatmaps(mc_id, groups=grp_ord, grp_ord=grp_ord, name='all', n_glob_enr=hm_n_glob_enr, n_grp_enr=hm_n_grp_enr, n_grp_outliers=hm_n_grp_outliers, show_marker_type=hm_show_marker_type, zlim=hm_zlim )  
  
  if (!is.null(group2meta)) {
    for (gnm in unique(group2meta)) { 
      nms = names(which(group2meta == gnm))
      
      if (length(nms) > 1) {
        mcell_mc_plot_group_composition(mc_id, mat_id, sample_field, md_fields_cols=md_fields_cols, sort_by=comp_sort_by, groups=nms, min_sample_cells = comp_min_sample_cells, name=gnm, selected_samples=selected_samples)
        mcell_mc_plot_group_enr_heatmaps(mc_id, groups=nms, grp_ord=intersect(grp_ord, nms), name=gnm, n_glob_enr=hm_n_glob_enr, n_grp_enr=hm_n_grp_enr, n_grp_outliers=hm_n_grp_outliers, show_marker_type=hm_show_marker_type, zlim=hm_zlim)
        mcell_mc_plot_group_enr_heatmaps(mc_id, groups=nms, name=gnm, n_glob_enr=hm_n_glob_enr, n_grp_enr=hm_n_grp_enr, n_grp_outliers=hm_n_grp_outliers, show_marker_type=hm_show_marker_type, zlim=hm_zlim)  
        mcell_mc2d_plot_group_subset(mc_id, groups=nms, name=gnm, graph_id=graph_id, show_mc_ids=meta_groups_mc2d_show_mc_ids)
      }
    }
  }  
}


