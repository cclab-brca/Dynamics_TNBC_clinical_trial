# PAR1006 PDX scRNAseq analysis -----
require(readr, quietly = T)
require(dplyr, quietly = T)
require(ggplot2, quietly = T)
require(MASS, quietly = T)
require(RColorBrewer, quietly = T)
require(devtools, quietly = T)
require(Matrix, quietly = T)
require(pheatmap, quietly = T)
require(metacell, quietly = T)

source("PAR1006_scRNAseq_utils.R")

rebuild = F

# Set output directory (change if needed)
wd = "output"
dir.create(wd, showWarnings = F, recursive = T)

# Color and other definitions
hg_pref = "GRCh38"
models = 'PAR1006'

cdict = list(
  Treatment = c(Untreated="darkgray", Olaparib="red"),
  Condition = c(Un="darkgray", "T"="lightblue", PT="darkblue"),
  Experiment = c(ASHEA_PAR_EXP1_v2="#1B9E77",  ASHEA_PAR_EXP2b_v2="#D95F02", ASHEA_PAR_EXP3_v2="#7570B3"))

grp_ord = c('Epithelial', 'EM-hybrid', 'EM-hybrid-IER', 'Mesenchymal', 'Myoepithelial', 'Hypoxic', 'Immune-Act')

figs_dir = sprintf("%s/figs", wd)
scdb_dir = sprintf("%s/scdb", wd)

dir.create(figs_dir, showWarnings = F)
copyDirectory("data/scdb", scdb_dir)
scdb_init(scdb_dir, T)
scfigs_init(figs_dir)

tgconfig::override_params("config/PAR1006.yaml", package="metacell")

# build count matrix -----
sinfo = read_csv("data/PAR1006_metadata.csv")

exp_cols = setNames(brewer.pal(n=3, 'Set2'), unique(sinfo$Experiment))

hg_filt_id = "ASHEA_PAR_hg_filt"

if (rebuild | !scdb_obj_exists("mat", hg_filt_id)) {
  message("Building count matrix ...")
  
  hg_mat = NULL
  hg_md = NULL
  
  for (i in 1:nrow(sinfo)) {
    message(sprintf("Adding mat %d of %d from %s...", i, nrow(sinfo), sinfo[i, 'Sample_ID']))
    cmat = scdb_mat(sprintf("%s", sinfo[i, 'Submission_ID']))
    
    if (cmat@ncells > 0) {
      md = cbind(cmat@cell_metadata, sinfo[i,])
      
      hg_cells = rownames(md)[!is.na(md$soc_fine_assign) & md$soc_fine_assign == "singlet Human" & md$hg_tot_umis >= md$mm_tot_umis]
      
      hg_genes = grep(sprintf("^%s_", hg_pref), cmat@genes, v=T)
      
      if (i > 1) {
        stopifnot(all(rownames(hg_mat) == hg_genes))
      }
      
      if (length(hg_cells) > 0) {
        hgcmat = cmat@mat[hg_genes, hg_cells]
        hgcmd = md[hg_cells, ]
        hg_cells = paste(sinfo[i, 'Sample_ID'], hg_cells, sep="_")
        colnames(hgcmat) = hg_cells
        rownames(hgcmd) = hg_cells
        hg_mat = cbind(hg_mat, hgcmat)
        hg_md = rbind(hg_md, hgcmd)
      }
    }
  }
  
  rownames(hg_mat) = gsub(sprintf("%s_", hg_pref), "", rownames(hg_mat))

  scdb_add_mat(hg_filt_id, tgScMat(hg_mat, "umi", hg_md))
}

filt_id = hg_filt_id
filtD_id = paste0(filt_id, "_noD")

# remove doublets with Scrublet -----
if (rebuild || !scdb_obj_exists("mat", filtD_id)) {
  message("marking doublets with scrublet...")
  doub = mcell_mat_filter_doublets(filt_id, new_mat_id=filtD_id, break_by="Sample_ID", break_by_cutoffs=setNames(rep(0.3, nrow(sinfo)), sinfo$Sample_ID), plot_hist=T, update_mat_metadata=T)
}

mat = scdb_mat(filtD_id)
md = mat@cell_metadata[mat@cells, ]
sinfo = mat@cell_metadata %>% group_by(Sample_ID) %>% 
  summarise(f_doub = mean(doublet_by_Scrublet), n_cells = length(Sample_ID)) %>% 
  inner_join(x=sinfo)

.plot_start(scfigs_fn(filtD_id, "fDoub_by_nCells_per_samp"), 400, 400)
plot(log2(sinfo$n_cells), sinfo$f_doub, col=exp_cols[sinfo$Experiment], pch=19, cex=2, xlab='#cells (log2)', ylab="%doublet (estim)")
legend("topright", legend=names(exp_cols), fill=exp_cols, title='Exp')
dev.off()

# Generate laterel gene sets (blacklisted from feature genes) -----
meta_build_blacklist_gsets_by_gene_nms(filtD_id, "hg19")	
if (rebuild | !scdb_obj_exists("gset", "hg19_lateral")) {
  message("Lateral gene sets...")
  meta_generate_lateral_gene_sets(filtD_id, "hg19", selected_sets=c('cell_cycle', 'ifn', 'stress'))
}

mat = scdb_mat(filtD_id)
lat_gset = scdb_gset("hg19_mito")
lat_gset = gset_add_genes(lat_gset, names(scdb_gset("hg19_ncrna")@gene_set), 'ncRNA')
lat_gset = gset_add_genes(lat_gset, names(scdb_gset("hg19_lateral")@gene_set), 'cc_ifn_stress')

glob_lat_id = "hg19_mito_ncRNA_lat"
scdb_add_gset(glob_lat_id, lat_gset)

# Generate metacell object ----
if (rebuild || !scdb_obj_exists("mc", filtD_id)) {
  meta_mat2mc(filtD_id, lateral_gset_id = glob_lat_id, mc_K = 50)
}

# Remove suspected nuclei and Epithelial-Mesenchymal doublet MCs ----
filt2_id = paste0(filtD_id, "2")
if (rebuild || !scdb_obj_exists("mc", filt2_id)) {
  message("Removing low quality mcs ...")
  
  mc = scdb_mc(filtD_id)
  lfp = log2(mc@mc_fp)
  mc2d = scdb_mc2d(filtD_id)
  
  uc = colSums(mat@mat)
  mc_mu_umi = log2(tapply(uc[names(mc@mc)], mc@mc, mean))
  
  mc_f_mito = tapply(colSums(mat@mat[grep("^MT-", mat@genes), names(mc@mc)]) / uc[names(mc@mc)], mc@mc, mean)
  mc_f_ribo = tapply(colSums(mat@mat[grep("^RP[LS]", mat@genes), names(mc@mc)]) / uc[names(mc@mc)], mc@mc, mean)
  mc_f_malat_neat = tapply(colSums(mat@mat[c('MALAT1', 'NEAT1'), names(mc@mc)]) / colSums(mat@mat[, names(mc@mc)]), mc@mc, mean)
  
  f_biotype = mcell_mat_f_umi_by_gene_biotype(filtD_id, gene_pref="", genes_ifn="data/refdata-gex-GRCh38-2020-A_genes.tsv")
  f_biotype = f_biotype[, names(mc@mc)]
  
  mc_f_lnc = tapply(f_biotype['lncRNA', ], mc@mc, mean)
  mc_f_lnc_neto = mc_f_lnc - mc_f_malat_neat
  
  sus_h = which(mc_f_lnc > 0.2)
  sus_grps = rep('rest', ncol(lfp))
  sus_grps[sus_h] = 'high'
  sus_cdict = c(rest='darkgray', high='red')
  
  plt("MALAT1", "NEAT1", lfp, mc@colors, ofn=scfigs_fn(filtD_id, "NEAT1_vs_MALAT1"))
  plt("f_lncRNA", "f_ribo", lfp, sus_cdict[sus_grps], x=mc_f_lnc, y=mc_f_ribo, ofn=scfigs_fn(filtD_id, "fRibo_vs_flncRNA"))
  plt("f_lncRNA_neto", "f_ribo", lfp, sus_cdict[sus_grps], x=mc_f_lnc_neto, y=mc_f_ribo, ofn=scfigs_fn(filtD_id, "fRibo_vs_flncRNAneto"))
  plt("f_lncRNA", "mu_UMIs", lfp, sus_cdict[sus_grps], x=mc_f_lnc, y=log2(mc_mu_umi), ofn=scfigs_fn(filtD_id, "muUMIs_vs_flncRNA"))
  plt("f_mito", "mu_UMIs", lfp, sus_cdict[sus_grps], x=mc_f_mito, y=mc_mu_umi, ofn=scfigs_fn(filtD_id, "muUMIs_vs_fMito"))
  
  .plot_start(scfigs_fn(filtD_id, "mc2d_sus_nuclei"), 400, 400)
  plot(mc2d@mc_x, mc2d@mc_y, pch=21, bg=sus_cdict[sus_grps], cex=2)
  dev.off()
 
  sus_doub = which(lfp['KRT81', ] > 2 & lfp['VIM', ] > 1)
  sus_grps2 = rep('rest', ncol(lfp))
  sus_grps2[sus_doub] = 'Doublets'
  sus_cdict2 = c(rest='darkgray', Doublets='red')
  plt("KRT81", "VIM", lfp, sus_cdict2[sus_grps2], ofn=scfigs_fn(filtD_id, "VIM_vs_KRT81_sus_doub_MCs"))
  
  mc_samp = table(mc@mc, mat@cell_metadata[names(mc@mc), 'Sample_ID'])
  mc_samp_e = rowSums(mc_samp) %*% t(colSums(mc_samp)) / sum(mc_samp)
  mc_samp_z = (mc_samp - mc_samp_e) / sqrt(mc_samp_e + 1)
  .plot_start(scfigs_fn(filtD_id, "mc_sample_composition_z_score"), nrow(mc_samp) * 7 + 600, ncol(mc_samp) * 20 + 100)
  z_max = 25
  pheatmap(pmin(pmax(t(mc_samp_z), -z_max), z_max), breaks=seq(-z_max, z_max, len=101), treeheight_col=10, treeheight_row=10, cellwidth=7, cellheight=20, annotation_col=data.frame(row.names=1:ncol(lfp), Type=paste(sus_grps, sus_grps2)))
  dev.off()
  
  mc2 = mc_set_outlier_mc(mc, which(sus_grps != 'rest' | sus_grps2 != 'rest'))
  scdb_add_mc(filt2_id, mc2)
  for (nm in c('mat', 'gstat')) {
    file.symlink(sprintf("%s/%s.%s.Rda", .scdb_base, nm, filtD_id), sprintf("%s/%s.%s.Rda", .scdb_base, nm, filt2_id))
  }
}

# Annotate metacells and generare some basic figures ----
if (rebuild || !scdb_obj_exists("mc2d", filt2_id)) {
  supmc = meta_mc_plots(mc_id=filt2_id, graph_id=filtD_id, mat_id=filtD_id, lateral_gset_id = glob_lat_id, 
                        metadata_fields_to_export = c('Treatment', 'Condition', 'Experiment'), color_by_conf = T, max_gene_len=12)
              
}

# Paper plots ------
message("Plots...")
fig5b_genes = c('KRT81', 'CLDN3', 'SLPI', 'EGFR', 'VIM', 'ITGA6', 'FOS', 'JUNB', 'TAGLN', 'ACTA2', 'CA9', 'HILPDA', 'HLA-B', 'CD74', 'CCNB1', 'CDK1')

mc_plot_e_gc_barplots(mc_id = filt2_id, 'key_and_marker', genes=fig5b_genes,
                      grp_ord=grp_ord, intra_grp_ord_by=list("EM-hybrid"=c('VIM'), "EM-hybrid-IER"=c('FOS', 'JUNB'), "Mesenchymal"=c('VIM'),
                                                             "Myoepithelial"=c('VIM', 'TAGLN'), "Hypoxic"='CA9', "Immune-Act"='HLA-B', "Epithelial"=c('KRT81', 'CLDN3')))


mc = scdb_mc(filt2_id)
lfp = log2(mc@mc_fp)
mat = scdb_mat(filtD_id)

sig_sets = sapply(gsub("gset.", "", grep("gset.hg19_sig_.*_filt$", scdb_ls("gset"), v=T, perl=T)), function(gset_nm) names(scdb_gset(gset_nm)@gene_set) )
names(sig_sets) = gsub("hg19_sig_", "", gsub("_filt", "", names(sig_sets)))
sig_sets_filt = sapply(sig_sets, intersect, y=rownames(mc@e_gc))
sig_ann = do.call('rbind', lapply(names(sig_sets_filt), function(nm) data.frame(row.names=sig_sets_filt[[nm]], set=rep(nm, length(sig_sets_filt[[nm]])))))

uc = colSums(mat@mat)
f_mito = colSums(mat@mat[grep("^MT-", mat@genes), names(mc@mc)]) / uc[names(mc@mc)]
f_ribo = colSums(mat@mat[grep("^RP[LS]", mat@genes), names(mc@mc)]) / uc[names(mc@mc)]
mc_f_mito = tapply(f_mito, mc@mc, mean)
mc_f_ribo = tapply(f_ribo, mc@mc, mean)
mc_mu_umi = log2(tapply(uc[names(mc@mc)], mc@mc, mean))

f_biotype = mcell_mat_f_umi_by_gene_biotype(filtD_id, gene_pref="", genes_ifn="data/refdata-gex-GRCh38-2020-A_genes.tsv")
f_biotype = f_biotype[, names(mc@mc)]
mc_f_lnc = tapply(f_biotype['lncRNA', ], mc@mc, mean)

col2grp = get_mc_col2group(mc)

grp2col = get_mc_group2col(mc)
mc_ann = get_mc_pheatmap_ann(mc, 'grp')
mcs_grp = col2grp[mc@colors]

df = data.frame(mcs = colnames(lfp), type=mcs_grp, set='mean_UMIs', value_lab="mean UMIS (log2)", value=mc_mu_umi) 
df = rbind(df, data.frame(mcs = colnames(lfp), type=mcs_grp, set='f_lncRNA', value_lab="frac lncRNA UMIs", value=mc_f_lnc))
df = rbind(df, data.frame(mcs = colnames(lfp), type=mcs_grp, set='f_riboRNA', value_lab="frac Ribo UMIs", value=mc_f_ribo))
df = rbind(df, data.frame(mcs = colnames(lfp), type=mcs_grp, set='f_mitoRNA', value_lab="frac Mito UMIs", value=mc_f_mito))

df = df %>% filter(type %in% grp_ord) %>% mutate(type=factor(type, levels=grp_ord))
p = ggplot(df, aes(x=type, y=log2(value))) + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.9) +  geom_point(aes(fill=type), pch = 21, position = position_jitter(width=0.2)) + scale_fill_manual(values=grp2col[grp_ord]) + scale_x_discrete(guide=guide_axis(angle=30)) + labs(x="") + facet_wrap(~set, scales="free_y")
ggsave(scfigs_fn(filt2_id, "qcs_per_mc", ext='pdf'), p, width=7, height=4)  

# %cycling per supmc and condition
is_cycling_cutoff = 0.05
min_cells_in_supmc_for_cycling = 50

f_cc = get_gset_f_umis(filtD_id, "hg19_cell_cycle_filt")

.plot_start(scfigs_fn(filt2_id, "f_cc_hist"), 300, 300)
plot(density(f_cc), xlab='Cell Cycle module (%UMIs)', lwd=2, col='steelblue', main="")
abline(v=is_cycling_cutoff, lty=2)
dev.off()

df = data.frame(Cond=mat@cell_metadata[names(mc@mc), 'Condition'], supmc=col2grp[mc@colors[mc@mc]], f_cc=f_cc[names(mc@mc)]) %>% 
  mutate(Cond=ifelse(is.na(Cond), 'Un', Cond), cycling=f_cc >= is_cycling_cutoff) %>%
  mutate(Supmc = factor(supmc, levels=grp_ord), 
         Condition = factor(Cond, levels=c('Un', 'T', 'PT')))

cells_per_cond = table(df$Cond)

df_f = df %>% group_by(Condition, Supmc) %>% summarise(n=length(Supmc))  %>% filter(n >= min_cells_in_supmc_for_cycling) %>% inner_join(df)
df_g = df_f %>% group_by(Condition, Supmc, n) %>% summarise(f_cc = mean(cycling)) %>% mutate(f_cells=n / cells_per_cond[as.character(Condition)]) 

p = ggplot(df_g, aes(x=Condition, fill=Supmc, y=f_cc)) + geom_bar(stat="identity", position='dodge', color='black')  +  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25, size=1.5) + scale_fill_manual(values=grp2col[grp_ord]) + scale_x_discrete(guide=guide_axis(angle=30)) + labs(x="", y='% cycling')
ggsave(scfigs_fn(filt2_id, sprintf("f_cycling_by_supmc_cond_co%.2f", is_cycling_cutoff), ext='pdf'), p, width=6, height=3)

df_g = df %>% group_by(Condition, Supmc) %>% summarise(n=length(Supmc), f_cc = mean(cycling)) %>% mutate(f_cells = n / cells_per_cond[as.character(Condition)]) 

cond_levels = levels(df_g$Condition)
max_cc = max(df_g$f_cc)

.plot_start(scfigs_fn(filt2_id, sprintf("f_cycling_and_composition_co%.2f", is_cycling_cutoff)), length(cond_levels) * 180, 300)
par(mfrow=c(1, length(cond_levels)), mar=c(4,2.5,3,0.5))
for (cond in cond_levels) {
  df_c = filter(df_g, Condition == cond)
  
  cs = cumsum(df_c$f_cells)
  mids = (c(0, cs[-length(cs)]) + cs)/2
  
  plot.new()
  plot.window(xlim=c(0, max_cc * 1.1), ylim=c(0,1))
  
  rect(rep(0, length(mids)), c(0, cs[-length(cs)]), max_cc/10, cs, col=grp2col[as.character(df_c$Supmc)])
  segments(rep(max_cc/10, length(mids)), mids, max_cc/10 + df_c$f_cc, mids, lwd=2, col=grp2col[as.character(df_c$Supmc)])
  segments(max_cc/10 + df_c$f_cc, mids-0.03, max_cc/10 + df_c$f_cc, mids + 0.03, lwd=2, col=grp2col[as.character(df_c$Supmc)])
  axis(2)
  f_cc_at = seq(0, max_cc*1.1, by=0.04)
  axis(1, at=max_cc/10 + f_cc_at, labels = f_cc_at)
  title(main=cond, xlab = "% prolif cells")
}
dev.off()

supmc_cond_cc = dcast(df_g %>% select(Condition, Supmc, f_cc), Supmc ~ Condition) %>% tibble::column_to_rownames(var='Supmc')
.plot_start(scfigs_fn(filt2_id, sprintf("f_cycling_by_supmc_barplot_co%.2f", is_cycling_cutoff)), nrow(supmc_cond_cc) * 100 + 50, 300)
barplot(t(as.matrix(supmc_cond_cc)), beside=T, col=cdict$Condition, legend.text=names(cdict$Condition), ylim=c(-0.01, max_cc), ylab="% prolif cells")
rect(seq(1, by=4, length=length(mids)), -0.008, seq(4, by=4, length=length(mids)), -0.002, col=grp2col[rownames(supmc_cond_cc)])
dev.off()

# supmc freq by condition/treatment
samp_supmc = table(mat@cell_metadata[names(mc@mc), 'Sample_ID'], col2grp[mc@colors[mc@mc]])
samp_supmc_n = samp_supmc / rowSums(samp_supmc)

sann = sinfo %>% select(Treatment, Condition, Sample_ID, Experiment, n_cells) %>% 
  mutate(n_cells=log2(n_cells), Experiment=gsub("ASHEA_PAR_EXP", "", gsub("_v2", "", Experiment)), Condition=ifelse(is.na(Condition), "Un", Condition)) %>% 
  tibble::column_to_rownames('Sample_ID') %>% 
  arrange(Condition, Treatment, Experiment) %>% 
  as.data.frame

y = reshape2::melt(samp_supmc_n, varnames=c('Sample_ID', 'supmc'), value.name='frac')  %>% 
  left_join(sann %>% tibble::rownames_to_column(var="Sample_ID")) %>% 
  mutate(Condition=factor(Condition, levels=c('Un', 'T', 'PT')),
         supmc=factor(supmc, levels=grp_ord))

p = ggplot(y, aes(x=Condition, y=frac, fill=Treatment)) + geom_boxplot(outlier.shape = NA, aes(group=Condition, fill=supmc)) + geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.2), width=0.9, dotsize = 1.2) + facet_wrap(~supmc, nrow=1, scales="free_y") + scale_fill_manual(values=c(grp2col, cdict$Treatment)) + guides(fill="none") + labs(x="", title="%supmc per sample") + geom_vline(xintercept = 1.5, linetype='dashed', color='darkgray')
ggsave(scfigs_fn(filt2_id, "supmc_f_by_cond", ext='pdf'), p, width=ncol(samp_supmc)*1.75, height=3)

p = ggplot(y, aes(x=supmc, y=frac, fill=supmc)) + geom_boxplot(outlier.shape = NA, aes(group=supmc, fill=supmc)) + geom_jitter(position=position_jitter(0.2), show.legend=F, shape=21, colour='black', size=1.2) + facet_wrap(~Condition, nrow=1) + scale_fill_manual(values=c(grp2col, cdict$Treatment)) + guides(fill="none") + scale_x_discrete(guide=guide_axis(angle=30)) + labs(x="", title="%supmc per sample")
ggsave(scfigs_fn(filt2_id, "supmc_f_by_cond2", ext='pdf'), p, width=6, height=3)

# common plots
mcell_common_plots(filt2_id, mat_id=filtD_id, 'Sample_ID', cdict, grp_ord=grp_ord, hm_n_glob_enr=8, hm_n_grp_enr=8, hm_n_grp_outliers=5, hm_show_marker_type=F, hm_zlim=2, selected_samples=grep('AZD', unique(mat@cell_metadata[names(mc@mc), 'Sample_ID']), v=T, invert=T))

####
# gene cor over Mesen-* and Epithelial
emt_min_enr = 1.25
emt_min_enr_diff = 2.5
foc_mcs = grep('Mesenchymal|hybrid|Epithelial', mcs_grp)
gs_by_max =  names(which(apply(lfp[, foc_mcs], 1, max) > emt_min_enr))
gs_by_diff = names(which(apply(apply(lfp[, foc_mcs], 1, range), 2, diff) > emt_min_enr_diff))
foc_gs = intersect(gs_by_max, gs_by_diff)
cc = cor(t(lfp[foc_gs, foc_mcs]))
cw = 9; cutree_k = 4; treeheight = 15;

.plot_start(scfigs_fn(filt2_id, sprintf('gene_cor_over_Mesen_and_Epi_max%.2f_diff%.2f', emt_min_enr, emt_min_enr_diff)), nrow(cc) * cw + 200, ncol(cc) * cw + 100)
hm = pheatmap(cc, breaks=seq(-1, 1, len=101), treeheight_col=treeheight, treeheight_row=treeheight, cellwidth=cw, cellheight=cw, fontsize=7, cutree_row=cutree_k, cutree_col=cutree_k)
dev.off()

gmods = cutree(hm$tree_row, cutree_k)
mc_gmods = tgs_matrix_tapply(t(lfp[names(gmods), foc_mcs]), gmods, mean)
gmods_odir = scfigs_dir(filt2_id, 'Mesen_Epithelial_gmods')

.plot_start(scfigs_fn(filt2_id, "gmods_plts", dir=gmods_odir), 250 * cutree_k, 250 * cutree_k)
par(mfrow=rep(cutree_k, 2), mar=c(4,4,1,1))
for (i in seq(cutree_k)) {
  for (j in seq(cutree_k)) {
    if (i == j) {
      plot.new()
      plot.window(0:1, 0:1)
      text(0.1, 0.5, i, cex=3)
    } else {
      plt(j, i, mc_gmods, mc@colors[foc_mcs])
    }
  }
}
dev.off()  

for (i in seq(cutree_k)) {
  for (j in seq(cutree_k)) {
    if (i != j) {
      plt(j, i, mc_gmods, mc@colors[foc_mcs], ofn=scfigs_fn(filt2_id, sprintf("gmods_plts_%d_vs_%d", i, j), dir=gmods_odir))
    }
  }
}

# Signatures by cor to anchor genes - Epi: KRT81, Mesen: VIM, IER: IER2
foc_g = names(which(apply(lfp, 1, max) >= 0.5))
max_genes_in_sig = 50
min_cor_to_anchor = 0.5
sig_anchors = list(Epithelial='KRT81', Mesenchymal='VIM', IER='IER2')
sig_genes = sapply(names(sig_anchors), function(nm) rev(tail(names(which(sort(cor(lfp[sig_anchors[[nm]], foc_mcs], t(lfp[foc_g, foc_mcs]))[1,]) >= min_cor_to_anchor)), max_genes_in_sig)))
write.csv(sig_genes, scfigs_fn(filt2_id, "Epi_Mesen_IER_sig_genes", ext='csv'))

mc_sigs = t(apply(sig_genes, 2, function(v) colMeans(lfp[v, ])))
plt("Mesenchymal", "IER", mc_sigs[, foc_mcs], mc@colors[foc_mcs], ofn=scfigs_fn(filt2_id, "mu_sigs_IER_vs_Mesenchymal", dir=gmods_odir))
plt("Mesenchymal", "Epithelial", mc_sigs[, foc_mcs], mc@colors[foc_mcs], ofn=scfigs_fn(filt2_id, "mu_sigs_Epithelial_vs_Mesenchymal", dir=gmods_odir))

# cor with TFs
g_max_th = 0.2
n_tfs_per_sig = 50
tf_min_cor_to_sig = 0.5
g_max = apply(lfp, 1, max)
foc_gs = names(which(g_max >= g_max_th))

tfs = intersect(read.table("data/Lambert2018_TF_names_v_1.01.txt")$V1, foc_gs)

mc_sigs_wo_tfs = t(apply(sig_genes, 2, function(v) { del_tfs = intersect(v, tfs); message(sprintf("Removing %d TFs (%s)", length(del_tfs), paste0(del_tfs, collapse=", "))); colMeans(lfp[setdiff(v, tfs), ]) }))

gm_tf_cc = cor(t(mc_sigs_wo_tfs[, foc_mcs]), t(lfp[tfs, foc_mcs]))
foc_tfs = unique(unlist(lapply(rownames(gm_tf_cc), function(i) tail(names(which(sort(gm_tf_cc[i, ]) >= tf_min_cor_to_sig)), n_tfs_per_sig))))
tfs_ord = foc_tfs[order( apply(gm_tf_cc[, foc_tfs], 2, which.max) - 1e-3 * apply(gm_tf_cc[, foc_tfs], 2, max))]

g_ann = data.frame(row.names=foc_tfs, max_lfp=g_max[foc_tfs], dominant_sig=rownames(gm_tf_cc)[apply(gm_tf_cc[, foc_tfs], 2, which.max)])

.plot_start(scfigs_fn(filt2_id, sprintf('gmods_TFs_cor_min%.2f_%dperSig', g_max_th, n_tfs_per_sig), gmods_odir), length(tfs_ord)*12 + 300, nrow(gm_tf_cc)*40 + 200) 
pheatmap(gm_tf_cc[, tfs_ord], cluster_cols=F, breaks=seq(-1, 1, len=101), cluster_rows=F, annotation_col=g_ann, cellwidth=12, cellheight=40)
dev.off()  

# TFs in mcs ordered by sig strength
n_strats_l = list(EMT=5, IER=4)
foc_mcs_l  = list(EMT=which(mcs_grp %in% c('Epithelial', 'EM-hybrid', 'Mesenchymal')), 
                IER=which(mcs_grp %in% c('EM-hybrid', 'EM-hybrid-IER')))

tf_peaks_at = NULL
for (run_name in names(n_strats_l)) {
  message(sprintf("Processing %s...", run_name))
  n_strats = n_strats_l[[run_name]]
  foc_mcs = foc_mcs_l[[run_name]]
  mcs_grad = switch(run_name,
                    EMT=mc_sigs['Mesenchymal', foc_mcs] - mc_sigs['Epithelial', foc_mcs],
                    IER=mc_sigs['IER', foc_mcs])
  
  mc_strat = cut(mcs_grad, breaks=quantile(mcs_grad, seq(0, 1, len=n_strats+1)), labels=paste0(run_name, seq_len(n_strats)), include.lowest=T)
  tfs_by_strat = t(tgs_matrix_tapply(lfp[foc_tfs, foc_mcs], mc_strat, mean))
  tfs_by_strat_n = t(scale(t(tfs_by_strat)))
  max_strat = apply(tfs_by_strat_n, 1, which.max)
  tf_peaks_at[[run_name]] = data.frame(gene=names(max_strat), peak=levels(mc_strat)[max_strat])
                      
  tfs_ord = rownames(tfs_by_strat_n)[order(max_strat + 1e-3 * apply(tfs_by_strat_n, 1, max))]
  
  .plot_start(scfigs_fn(filt2_id, sprintf('TFs_on_%s_min%.2f_%dper_%dqs', run_name, g_max_th, n_tfs_per_sig, n_strats), gmods_odir), 600, length(tfs_ord) * 8 + 100) 
  pheatmap(tfs_by_strat_n[tfs_ord, ], cluster_cols=F, cluster_rows=F, annotation_row=g_ann, fontsize=7)
  dev.off()
  
  gmods_dir2 = sprintf("%s/TFs_by_%s_strats", gmods_odir, run_name)
  dir.create(gmods_dir2, showWarnings = F)
  for (i in seq_along(tfs_ord)) {
    nm = tfs_ord[i]
    .plot_start(scfigs_fn(filt2_id, sprintf("%s%d_%d_%s_by_%dqs", run_name, max_strat[nm], i, nm, n_strats), gmods_dir2), 400, 300)
    boxplot(lfp[nm, foc_mcs] ~ mc_strat, outline=F, col=NA, main=nm, xlab='', ylab='enr')
    stripchart(lfp[nm, foc_mcs] ~ mc_strat, method='jitter', jitter=0.2, vertical=T, add=T, pch=19)
    dev.off()
  }
  
  strat_comp = table(mc_strat, mcs_grp[foc_mcs])
  strat_comp_n = strat_comp / rowSums(strat_comp)
  strat_comp_n = strat_comp_n[, intersect(grp_ord, colnames(strat_comp_n))]
  .plot_start(scfigs_fn(filt2_id, sprintf('%s_%dqs_supmc_composition', run_name, n_strats), gmods_odir), 500, 300)
  barplot(t(strat_comp_n), col=grp2col[colnames(strat_comp_n)], ylab='% of mcs')
  dev.off()
}

tf_peaks_pairs = inner_join(tf_peaks_at[[1]], tf_peaks_at[[2]], by='gene')
peaks_t = table(tf_peaks_pairs$peak.x, tf_peaks_pairs$peak.y)
write.csv(peaks_t, scfigs_fn(filt2_id, 'gene_peak_strat', gmods_odir, ext='csv'))
