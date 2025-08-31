# Tillandsia phylogenomics and hybridization

Code to reproduce results from Yardeni et al., [2025](https://doi.org/10.1093/sysbio/syaf039)
This was used to generate whole-genome phylogenomics and infer hybridization in _Tillandsia_ subgenus _Tillandsia_.

* raw_to_vcf_GATK.md - processing BAM files to VCF calls and filtering.
* KING_analysis.sh - short script to infer kinship coefficients with KING
* ASTAL_gene_trees_wins.sh - script to generate VCF for genomic intervals and infer a species tree with ASTRAL
* example_plot_mod_180121.R - script to create *Twisst* plots, modified slightly from [original](https://github.com/simonhmartin/twisst)
* plot_twisst_mod.R - modified *Twisst* functions
* iqtree.sh - generic script to produce trees with iqtree2
* makewindows.sh - generic trees to infer ML trees on genomic windows from a VCF file.
