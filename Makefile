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

.PHONY: WGS

WGS: data/mutation_data/mutations_by_bc.csv


data/mutation_data/mutations_by_bc.csv :
	Rscript scripts/combine_BCs_and_WGS.R


.PHONY: hBFA

hBFA: $(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv $(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clusts.csv $(OUTDIR)/$(hBFA_BASENAME)_adapted_w_clust_means.csv


$(OUTDIR)/$(hBFA_BASENAME)_adapteds.csv :
	Rscript scripts/call_adapteds.R \
		--use_iva \
	  --exclude=X48Hr \
	  --gens=8 \
	  --cutoff=0.05 \
	  --base_name=$(hBFA_BASENAME) \
	  --outdir=$(OUTDIR) \
	  $(hBFA_S_FILE) \
		$(hBFA_NEUTRAL_ENV)


$(OUTDIR)/$(hBFA_BASENAME)_adapteds_autodips.csv :
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
