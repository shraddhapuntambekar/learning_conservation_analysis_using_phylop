#### Ref: This work is published in "Puntambekar, Shraddha, Rachel Newhouse, Jaime San-Miguel, Ruchi Chauhan, Grégoire Vernaz, Thomas Willis, Matthew T. Wayland, Yagnesh Umrania, Eric A. Miska, and Sudhakaran Prabakaran. "Evolutionary divergence of novel open reading frames in cichlids speciation." _Scientific reports 10_, no. 1 (2020): 1-18." (https://www.nature.com/articles/s41598-020-78555-0 )

# Conservation scores analysis using phyloP

Cichlid fishes are remarkable example of adaptive radiation as they exhibit large phenotypic diversity and rapid speciation. The objective of this study is to determine the genomic segments underlying phenotypic variation in cichlids. 

I used phyloP to identify the genomic regions in Oreochromis niloticus (ON, a cichlid fish) that have diverged at a rate different than that expected at a neutral drift. The five-way whole genome alignment of 5 cichlid fishes (*O. niloticus, P. nyererei, A. burtoni, M. zebra, N. brichardi*) constructed by Brawand et al was used to detect signatures of nonneutral evolution in aligned genomic sequences. 

Prior to running phyloP on the 5-way whole genome multiple alignment (maf); to produce genome-wide, per base, conservation scores; the alignments were first filtered to remove potential alignment and assembly errors and to minimize the effect of missing data.  Next few steps involve constructing a neutral substitution model by fitting  a time reversible substitution 'REV' model on the phylogeny obtained from (four fold degenerate) 4D sites, which are extracted based on ON's protein coding sequences. Then using a neutral model phyloP is run on the filtered, sorted mafs with an 'all branches' test to predict conservation-acceleration (CONACC) score for every site in the whole genome. All these steps are documented in the file **'1_run_all_steps_phylop.sh'**.

The next few steps are to map the predicted conservation scores on the ON's features like cds, exon, intron, intergenes and ancestral repeats. Conservation score for a particular feature (suppose say for exon 1 of a particular gene) is measured as score average over regions mapped contigously with scores; or the score is calculated as average of average calculate over 10bp window. Score mapping steps are documented in file **'2_mapping_scores_on_features_using_bedops.sh'**.

The final R script **'3_making_cdf_plots_ofmappedphylopscores.R'** contains steps to plot and to perform a Welch t-test (or Wilcoxon test) to compare the distribution of conservation scores for different features.
