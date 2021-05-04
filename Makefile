### RECIPES FOR MKDOCS
### phumph.github.io/pleiotropy/

.PHONY: site serve

site:
	cp output/figures/*.png docs/img/
	mkdocs build docs

serve:
	mkdocs serve

###
### VARIABLES
###

hBFA_S_FILE=data/fitness_data/fitness_estimation/hBFA1_s_03_23_18_GC_cutoff_5.csv
dBFA_S_FILE=data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv
MUTATIONS_DIR=data/mutation_data
MUTATIONS_FILE=$(MUTATIONS_DIR)/master_mutation_calls.txt

hBFA_BASENAME=hBFA1_cutoff-5
hBFA_NEUTRAL_ENV=YPD_alpha
dBFA_BASENAME=dBFA2_cutoff-5
dBFA_NEUTRAL_ENV=Ancestor_YPD_2N
OUTDIR=data/fitness_data/fitness_calls
COMBINED=data/combined
OUTPUT=output
TABLEDIR=$(OUTPUT)/tables
FIGDIR=$(OUTPUT)/figures


.PHONY: WGS

WGS: data/mutation_data/mutations_by_bc.csv

data/mutation_data/mutations_by_bc.csv: scripts/combine_BCs_and_WGS.R
	Rscript $<

.PHONY: figs

# cluster-level plots
fig_deps := $(FIGDIR)/$(dBFA_BASENAME)_pleio_plot.pdf
fig_deps += $(FIGDIR)/$(dBFA_BASENAME)_clusts_plot.pdf
fig_deps += $(FIGDIR)/$(dBFA_BASENAME)_mu_plot.pdf
fig_deps += $(FIGDIR)/$(dBFA_BASENAME)_mu_v_sigma_plot.pdf
fig_deps += $(FIGDIR)/$(hBFA_BASENAME)_pleio_plot.pdf
fig_deps += $(FIGDIR)/$(hBFA_BASENAME)_clusts_plot.pdf
fig_deps += $(FIGDIR)/$(hBFA_BASENAME)_mu_plot.pdf
fig_deps += $(FIGDIR)/$(hBFA_BASENAME)_mu_v_sigma_plot.pdf

# barcode-level plots
fig_deps += $(FIGDIR)/heatmap_by_bc_plot.pdf
fig_deps += $(FIGDIR)/heatmap_with_muts.pdf
fig_deps += $(FIGDIR)/dendro_by_bc_ploidy.pdf
fig_deps += $(FIGDIR)/dendro_by_bc_source.pdf

figs: $(fig_deps)

$(FIGDIR)/heatmap_by_bc_plot.pdf $(FIGDIR)/heatmap_with_muts.pdf $(FIGDIR)/dendro_by_bc_ploidy.pdf $(FIGDIR)/dendro_by_bc_source.pdf:
	Rscript scripts/plot_by_bc.R \
		"data/combined/hBFA1_cutoff-5_compiled_data_by_barcode.csv data/combined/dBFA2_cutoff-5_compiled_data_by_barcode.csv" \
		"data/mutation_data/mutations_by_bc.csv" \
		-o $(FIGDIR)

$(FIGDIR)/$(dBFA_BASENAME)_mu_plot.pdf $(FIGDIR)/$(dBFA_BASENAME)_pleio_plot.pdf $(FIGDIR)/$(dBFA_BASENAME)_clusts_plot.pdf $(FIGDIR)/$(dBFA_BASENAME)_mu_v_sigma_plot.pdf: scripts/plot_by_cluster.R
	Rscript scripts/plot_by_cluster.R \
		--outdir=$(FIGDIR)\
		$(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_plot-data.csv \
		$(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_table.csv

$(FIGDIR)/$(hBFA_BASENAME)_mu_plot.pdf $(FIGDIR)/$(hBFA_BASENAME)_pleio_plot.pdf $(FIGDIR)/$(hBFA_BASENAME)_clusts_plot.pdf $(FIGDIR)/$(hBFA_BASENAME)_mu_v_sigma_plot.pdf: scripts/plot_by_cluster.R
	Rscript scripts/plot_by_cluster.R \
		--outdir=$(FIGDIR) \
		$(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_plot-data.csv \
		$(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_table.csv


.PHONY: summaries

summary_deps := $(COMBINED)/$(hBFA_BASENAME)_compiled_data_by_barcode.csv
summary_deps += $(COMBINED)/$(dBFA_BASENAME)_compiled_data_by_barcode.csv

summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_table.csv

summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_table.csv

summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_table.csv

summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_table.csv

summary_deps += $(MUTATIONS_DIR)/dBFA2_mutations_by_cluster.csv
summary_deps += $(MUTATIONS_DIR)/hBFA1_mutations_by_cluster.csv

summary_deps += $(FIGDIR)/$(hBFA_BASENAME)_source_summaries_plot.pdf
summary_deps += $(FIGDIR)/$(dBFA_BASENAME)_source_summaries_plot.pdf

summaries: $(summary_deps)

$(COMBINED)/$(hBFA_BASENAME)_compiled_data_by_barcode.csv:
	Rscript scripts/compile_data_by_barcode.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(COMBINED) \
		$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv \
		$(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv

$(COMBINED)/$(dBFA_BASENAME)_compiled_data_by_barcode.csv:
	Rscript scripts/compile_data_by_barcode.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(COMBINED) \
		$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv \
		$(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv

$(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_plot-data.csv $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_table.csv $(FIGDIR)/$(hBFA_BASENAME)_source_summaries_plot.pdf:
	Rscript scripts/summarise_sources.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTPUT) \
		$(COMBINED)/$(hBFA_BASENAME)_compiled_data_by_barcode.csv


$(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_plot-data.csv $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_table.csv $(FIGDIR)/$(dBFA_BASENAME)_source_summaries_plot.pdf:
	Rscript scripts/summarise_sources.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTPUT) \
		$(COMBINED)/$(dBFA_BASENAME)_compiled_data_by_barcode.csv


$(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_plot-data.csv $(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_table.csv $(MUTATIONS_DIR)/hBFA1_mutations_by_cluster.csv:
	Rscript scripts/summarise_clusters.R \
		--exclude=48Hr \
		--outdir=$(OUTPUT) \
		$(COMBINED)/$(hBFA_BASENAME)_compiled_data_by_barcode.csv \
		data/mutation_data/mutations_by_bc.csv


$(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_plot-data.csv $(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_table.csv $(MUTATIONS_DIR)/dBFA2_mutations_by_cluster.csv:
	Rscript scripts/summarise_clusters.R \
		--exclude=48Hr \
		--outdir=$(OUTPUT) \
		$(COMBINED)/$(dBFA_BASENAME)_compiled_data_by_barcode.csv \
		data/mutation_data/mutations_by_bc.csv


PHONY: hBFA

hBFA: $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv


$(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv:
	Rscript scripts/call_adapteds.R \
		--use_iva \
	  --exclude=CLM\|FLC4\|Stan\|48Hr \
	  --gens=8 \
	  --cutoff=0.05 \
	  --base_name=$(hBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(hBFA_S_FILE) \
		$(hBFA_NEUTRAL_ENV)


$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv:
	Rscript scripts/filter_autodiploids.R \
	  --use_iva \
	  --exclude=CLM\|FLC4\|Stan\|48Hr \
	  --gens=8 \
	  --cutoff=0.01 \
	  --base_name=$(hBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv \
	  autodiploids


$(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv:
	Rscript scripts/cluster_lineages.R \
		--use_iva \
		--exclude=CLM\|FLC4\|Stan\|48Hr \
		--gens=8 \
		--outdir=$(OUTDIR) \
		$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv


.PHONY: dBFA

dBFA: $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv


$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv :
	Rscript scripts/call_adapteds.R \
		--use_iva \
	  --exclude=CLM\|FLC4\|Stan\|48Hr \
	  --gens=8 \
	  --cutoff=0.05 \
	  --base_name=$(dBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(dBFA_S_FILE) \
		$(dBFA_NEUTRAL_ENV)


$(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv:
	Rscript scripts/cluster_lineages.R \
		--use_iva \
		--exclude=CLM\|FLC4\|Stan\|48Hr \
		--gens=8 \
		--outdir=$(OUTDIR) \
		$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv


.PHONY: clean

clean_all:

	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv
	rm -f $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv
	rm -f $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv
	rm -f $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv
	rm -f $(summary_deps)


.PHONY: clean_hBFA

clean_hBFA:

	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv
	rm -f $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv


.PHONY: clean_dBFA

clean_dBFA:

	#rm $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv
	rm -f $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv
	rm -f $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv


.PHONY: clean_summaries

clean_summaries:

	rm -f $(summary_deps)
