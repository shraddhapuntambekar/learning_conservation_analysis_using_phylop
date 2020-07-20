#### Mapping the per base, genome-wide conservation scores predicted by phylop on ON's different annotation features like CDS, exons, introns, 5'UTR, 3'UTR, intergenes and ancestral repeats. 

# (1) Download the ON's annotations from Cambridge cichlid browser. The annotations are downloaded in bed format. 
# (i) For CDS, exons, introns and UTRs of ON -> Open the table browser - select 'Vertebrate', 'Orechromis niloticus', 'Gene and Gene Prediction', 'Broad Genes', 'Bed format' and select the region to download ((http://em-x1.gurdon.cam.ac.uk/cgi-bin/hgTables?hgsid=21982&clade=vertebrate&org=O.+niloticus&db=on11&hgta_group=genes&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=LG2%3A1959784-2269783&hgta_outputType=bed&hgta_outFileName=)).

# (ii) Intergenic regions are assumed to be regions that are not annotated in the whole genome, the regions complementary to the features annotated in gff (https://www.biostars.org/p/112251/).  
# Download information regarding the chromosome sizes. Make a two columned file having chromosome name in the first column and the respective chromosome size in the second column. 
	awk '{print $1"\t"$3}' chromSizes.bed > chromSizes_2cols.bed 
	sort -k1,1 -k2,2n chromSizes_2cols.bed > chromSizes_2cols_sorted.bed 

# Using bedtools to find regions complementary to the features in gff file.
	bedtools complement -i wholeGene_sorted.bed -g chromSizes_2cols_sorted.bed > ON_intergenic_complementof_wholegene_and_chrmSizes.bed 

# (iii) Ancestral repeats (ARs) are assumed to be repeat masked sequences from ON that are also conserved in teleosts. The AR regions were downloaded from cichlid genome browser by taking an intersection (having at least 80% overlap) between repeat masked regions from ON and 8-way cichlids multiple alignments (On_Mz_Pn_Ab_Nb_oryLat2_gasAcu1_danRer7_maf).

# (2) Replace refseq chr to ucsc chr

# (3) Using 'awk' just keep the first three columns of bed files 
	awk -F "\t" '{ print $1"\t"$2"\t"$3)' cds.bed > cds_ON.bed 

# (4) Add a unique id to each row as id1, id2..

# (5) Using bedops to sort the bed file 
	sort-bed introns_ON_id.bed > introns_ON_id_sorted.bed

# (6) Splitting each annotation bed files perbase, as is needed for bedops 
awk ' \
    { \
         regionChromosome = $1; \
         regionStart = $2; \
         regionStop = $3; \
         regionID = $4; \
         baseIdx = 0; \
         for (baseStart = regionStart; baseStart < regionStop; baseStart++) { \
             baseStop = baseStart + 1; \
             print regionChromosome"\t"baseStart"\t"baseStop"\t"regionID"-"baseIdx; \
             baseIdx++; \
         } \
    }' introns_ON_id_sorted.bed  > input_to_bedops/introns_ON_id_sorted_perbase.bed

# (7) Using bedops to map scores on the features:
### (i) averaging over contigous regions
cd input_to_bedops/perbase/
WF=phylopoutput/all
outF=out_map_roi/all 
for f in *_perbase.bed; do
	bedmap --skip-unmapped --echo --echo-map-score --delim '\t' $f $WF/ACC/sorted_phylop_ONmod_Allbranch_ACC.bed > $outF/ACC/bm_${f}
	bedops --merge $outF/ACC/bm_${f} > $outF/ACC/merged_${f}
	bedmap --echo --mean --stdev --min --max --median --sum --variance --delim '\t' $outF/ACC/merged_${f} $outF/ACC/bm_${f} > $outF/ACC/meanAndSd_${f}

	bedmap --skip-unmapped --echo --echo-map-score --delim '\t' $f $WF/CONACC/sorted_phylop_ONmod_Allbranch_CONACC.bed > $outF/CONACC/bm_${f}
	bedops --merge $outF/CONACC/bm_${f} > $outF/CONACC/merged_${f}
	bedmap --echo --mean --stdev --min --max --median --sum --variance --delim '\t' $outF/CONACC/merged_${f} $outF/CONACC/bm_${f} > $outF/CONACC/meanAndSd_${f}

	bedmap --skip-unmapped --echo --echo-map-score --delim '\t' $f $WF/GERP/sorted_phylop_ONmod_Allbranch_GERP.bed > $outF/GERP/bm_${f}
	bedops --merge $outF/GERP/bm_${f} > $outF/GERP/merged_${f}
	bedmap --echo --mean --stdev --min --max --median --sum --variance --delim '\t' $outF/GERP/merged_${f} $outF/GERP/bm_${f} > $outF/GERP/meanAndSd_${f}
done

### (i) averaging over 10bp window size
## All branch - ACC and CONACC and GERP
cd input_to_bedops/perbase/
outF=out_map_roi/all
outF10bp=out_map_roi/all_10bpsw
for f in *_perbase.bed; do
   bedops --chop 10 --stagger 1 $outF/ACC/merged_${f} > $outF10bp/ACC/chopped10_${f}
   bedmap --echo --mean --delim '\t' $outF10bp/ACC/chopped10_${f} $outF/ACC/bm_${f} > $outF10bp/ACC/allDetails_${f}
 bedmap --echo --mean --delim '\t' $outF/ACC/merged_${f} $outF10bp/ACC/allDetails_${f} > $outF10bp/ACC/avgOFavg_of10bpSW_${f}
 sudo rm $outF10bp/ACC/chopped10_${f}
 sudo rm $outF/ACC/merged_${f}
 sudo rm $outF/ACC/bm_${f}
        
   bedops --chop 10 --stagger 1 $outF/CONACC/merged_${f} > $outF10bp/CONACC/chopped10_${f}
   bedmap --echo --mean --delim '\t' $outF10bp/CONACC/chopped10_${f} $outF/CONACC/bm_${f} > $outF10bp/CONACC/allDetails_${f}
 bedmap --echo --mean --delim '\t' $outF/CONACC/merged_${f} $outF10bp/CONACC/allDetails_${f} > $outF10bp/CONACC/avgOFavg_of10bpSW_${f}
   sudo rm $outF10bp/CONACC/chopped10_${f}
 sudo rm $outF/CONACC/merged_${f}
 sudo rm $outF/CONACC/bm_${f}


  bedops --chop 10 --stagger 1 $outF/GERP/merged_${f} > $outF10bp/GERP/chopped10_${f}
  bedmap --echo --mean --delim '\t' $outF10bp/GERP/chopped10_${f} $outF/GERP/bm_${f} > $outF10bp/GERP/allDetails_${f}
 bedmap --echo --mean --delim '\t' $outF/GERP/merged_${f} $outF10bp/GERP/allDetails_${f} > $outF10bp/GERP/avgOFavg_of10bpSW_${f}
  sudo rm $outF10bp/GERP/chopped10_${f}
 sudo rm $outF/GERP/merged_${f}
 sudo rm $outF/GERP/bm_${f}

done


