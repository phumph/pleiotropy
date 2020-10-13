###
### VARIABLES
###

hBFA_S_FILE=data/fitness_data/fitness_estimation/hBFA1_s_03_23_18_GC_cutoff_5.csv
dBFA_S_FILE=data/fitness_data/fitness_estimation/dBFA2_s_03_23_18_GC_cutoff_5.csv
MUTATIONS_FILE=data/mutation_data/master_mutation_calls.txt

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


data/mutation_data/mutations_by_bc.csv:
	Rscript scripts/combine_BCs_and_WGS.R


.PHONY: summaries

summary_deps := $(COMBINED)/$(hBFA_BASENAME)_compiled_data_by_barcode.csv
summary_deps += $(COMBINED)/$(dBFA_BASENAME)_compiled_data_by_barcode.csv
summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_table.csv
summary_deps += $(FIGDIR)/$(hBFA_BASENAME)_source_summaries_plot.pdf
summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_plot-data.csv
summary_deps += $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_table.csv
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
		$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv \
		$(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv

$(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_plot-data.csv $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_table.csv $(FIGDIR)/$(dBFA_BASENAME)_source_summaries_plot.pdf:
	Rscript scripts/summarise_sources.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTPUT) \
		$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv \
		$(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv

$(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_plot-data.csv $(TABLEDIR)/$(hBFA_BASENAME)_cluster_summaries_table.csv $(FIGDIR)/$(hBFA_BASENAME)_cluster_summaries_plot.pdf:
	Rscript scripts/summarise_clusters.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTPUT) \
		$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv \
		$(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv

$(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_plot-data.csv $(TABLEDIR)/$(dBFA_BASENAME)_cluster_summaries_table.csv $(FIGDIR)/$(dBFA_BASENAME)_cluster_summaries_plot.pdf:
	Rscript scripts/summarise_clusters.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTPUT) \
		$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv \
		$(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv \
		data/mutation_data/mutations_by_bc.csv


PHONY: hBFA

hBFA: $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv


$(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv:
	Rscript scripts/call_adapteds.R \
		--use_iva \
	  --exclude=X48Hr \
	  --gens=8 \
	  --cutoff=0.05 \
	  --base_name=$(hBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(hBFA_S_FILE) \
		$(hBFA_NEUTRAL_ENV)


$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv:
	Rscript scripts/filter_autodiploids.R \
	  --use_iva \
	  --exclude=X48Hr \
	  --gens=8 \
	  --cutoff=0 \
	  --base_name=$(hBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv \
	  autodiploids


$(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv:
	Rscript scripts/cluster_lineages.R \
		--use_iva \
		--exclude=X48Hr \
		--gens=8 \
		--outdir=$(OUTDIR) \
		$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv


.PHONY: dBFA

dBFA: $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv


$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv :
	Rscript scripts/call_adapteds.R \
		--use_iva \
	  --exclude=CLM\|FLC4\|Stan \
	  --gens=8 \
	  --cutoff=0.05 \
	  --base_name=$(dBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(dBFA_S_FILE) \
		$(dBFA_NEUTRAL_ENV)


$(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv:
	Rscript scripts/cluster_lineages.R \
		--use_iva \
		--exclude=CLM\|FLC4\|Stan \
		--gens=8 \
		--outdir=$(OUTDIR) \
		$(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv


.PHONY: clean

clean_all:

	rm $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv
	rm $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv
	rm $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv
	rm $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv

	rm $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv
	rm $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv
	rm $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv


.PHONY: clean_hBFA

clean_hBFA:

	#rm $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv
	#rm $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv
	rm $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv
	rm $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv


.PHONY: clean_dBFA

clean_dBFA:

	#rm $(OUTDIR)/$(dBFA_BASENAME)_adapteds.csv
	rm $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clusts.csv
	rm $(OUTDIR)/$(dBFA_BASENAME)_adapted_w_clust_means.csv


.PHONY: clean_summaries

clean_summaries:

	rm $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_plot-data.csv
	rm $(TABLEDIR)/$(hBFA_BASENAME)_source_summaries_table.csv
	rm $(FIGDIR)/$(hBFA_BASENAME)_source_summaries_plot.pdf
	rm $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_plot-data.csv
	rm $(TABLEDIR)/$(dBFA_BASENAME)_source_summaries_table.csv
	rm $(FIGDIR)/$(dBFA_BASENAME)_source_summaries_plot.pdf
