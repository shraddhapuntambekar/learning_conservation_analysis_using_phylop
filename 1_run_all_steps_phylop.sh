### Following are the steps to calculate genome-wide, per base conservation scores using 'phylop'. The five-way whole genome alignment of ​ O. niloticus (ON)​ , ​ N. brichardi , ​ ​ A. burtoni , ​ ​ M. zebra and ​ P. nyererei (PN) genomes provided by Brawand et al and ON's annotation file downloaded from Cichlids genome browser is used in this analysis.

## Processing steps before running phylop: 
## (i) phylop requires chromosome-wise maf and gtf files for conservation score calculations, hence the need to chromosome-wise split the maf and gtf files.
## (ii) filtering the multiple alignments (maf file) using mafFilter to discard blocks which have sequences less than five and to remove gap only columns from maf blocks. The filtered mafs are then sorted using 'mafs-sort.sh' script from LAST (https://github.com/UCSantaCruzComputationalGenomicsLab/last.git). The alignment blocks are filtered to remove potential alignment and assembly errors and to minimize the effect of missing data. 

## Step 1 - Extracting names of chromosomes and scaffolds from Orenil1.1 assembly
sudo wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/235/GCF_000188235.2_Orenil1.1/GCF_000188235.2_Orenil1.1_assembly_report.txt

awk '{print $1}' <(awk 'NR > 34 { print }' < GCF_000188235.2_Orenil1.1_assembly_report.txt) > chrm_from_ON.txt

## Step 2: Splitting ON's gtf file (downloaded from cichlid genome browser) into multiple files; each having information per chromosome.
awk '{print "Niloticus."$0}' OreochromisNiloticus_fromCichlidBrowser_gtf.gtf > OreochromisNiloticus_fromCichlidBrowser_gtf_1.gtf

while read chr; do
  grep -P "^Niloticus.$chr\t" OreochromisNiloticus_fromCichlidBrowser_gtf_1.gtf > chrwise_gtf_ON/OreochromisNiloticus_fromCichlidBrowser_gtf_$chr.gtf
done < chrm_from_ON.txt

## Step 3: Editing the MAF files so that it is readable for maffilt
while read eachMaf; do
sudo sed -i 's/^s\tNiloticus:/s\tNiloticus./' $eachmaf.maf
sudo sed -i 's/^s\tBrichardi:/s\tBrichardi./' $eachmaf.maf 
sudo sed -i 's/^s\tNyererei:/s\tNyererei./' $eachmaf.maf
sudo sed -i 's/^s\tZebra:/s\tZebra./' $eachmaf.maf
sudo sed -i 's/^s\tBurtoni:/s\tBurtoni./' $eachmaf.maf
done < chrm_from_ON.txt

## Running maffilt for filtering the low quality regions from maf. 'MinBlockSize' of maffilter discards blocks with less sequences than a given threshold. 'XFullGap' removes any column containing gaps in the specified species
## Step 4: First create a text file named 'options_file' and save following lines in it. 

input.file=LG1.maf
input.format=Maf
 maf.filter=(MinBlockSize(min_size=5),XFullGap(species=(Niloticus,Nyererei,Brichardi,Zebra,Burtoni)),Output(file=processedMAFS/filteredMAFS/filt_LG1.maf))

## Step 5: then run the following loop
while read chr; do
 sed "s/LG1/$chr/g" options_file > all_optionfiles/options_file_$chr
 
MafFilter/maffilter param=all_optionfiles/options_file_$chr
done < chrm_from_ON.txt

## Step 6: sort the filtered mafs
cd processedMAFS/filteredMAFS/
for eachfile in filt_*.maf
do
  mafTools/last/scripts/maf-sort.sh $eachfile > ../sortedMAFS/s_$eachfile
done

### The next following steps involve in constructing a neutral substitution model using phyloFit (from PHAST) by fitting  a time reversible substitution 'REV' model on the phylogeny obtained from (four fold degenerate) 4D sites. msa_view (from PHAST) is used to extract 4D sites based on ON's protein coding sequences (excluding sites from mitochondria). 
 
## Step 7 - Extracting 4D sites  
while read chr; do
  	sudo msa_view processedMAFS/sortedMAFS/s_filt_$chr.maf --in-format MAF --4d --features chrwise_gtf_ON/OreochromisNiloticus_fromCichlidBrowser_gtf_$chr.gtf > msa_view_output_ON/codons/Oniloticus_4D_codons_$chr.ss

	sudo msa_view msa_view_output_ON/codons/Oniloticus_4D_codons_$chr.ss --in-format SS --out-format SS --tuple-size 1 > msa_view_output_ON/sites/${chr}_Oniloticus_4Dsites.ss
done < chrm_from_ON.txt

## Step 8 - Aggregating the 4D sites
cd msa_view_output_ON/sites/
awk '{if ($1!="0") print $2}'  <(wc -l *_Oniloticus_4Dsites.ss ) > chr_with_4D_sites_ON.txt

sed -i '823d' chr_with_4D_sites_ON.txt

Fnames=$(awk 'FNR!=1{print l}{l=$0};END{ORS="";print l}' ORS=' ' chr_with_4D_sites_ON.txt)

msa_view --unordered-ss --in-format SS --out-format SS --aggregate Niloticus,Brichardi,Nyererei,Zebra,Burtoni  $Fnames > ../aggregated_allchromosomes_ON_4Dsites.ss

## Step 9 - Constructing a neutral model 
phyloFit --tree "(Niloticus, (Brichardi, (Burtoni, (Nyererei, Zebra))))" --subst-mod REV --out-root neutralNiloticus_allChr_n_scaff_7Apr19 --msa-format SS aggregated_allchromosomes_ON_4Dsites.ss

## Running phyloP (chromosome-wise), on the filtered mafs using the likelihood ratio test (LRT) method and an 'all branches' test to predict conservation-acceleration (CONACC) score for every site in the whole genome. The output of phylop is stored in 'wig' format.

## Step 10 - Calculating ACC, CONACC and GERP scores using phylop
while read chr; do
	phyloP --wig-scores --method LRT --chrom $chr --mode ACC msa_view_output_ON/neutralNiloticus_allChr_n_scaff_7Apr19.mod --msa-format MAF processedMAFS/sortedMAFS/s_filt_$chr.maf > phylopoutput/all/ACC/${chr}_allbranch_ACC_ONmod_7Apr.wig

	phyloP --wig-scores --method LRT --chrom $chr --mode CONACC msa_view_output_ON/neutralNiloticus_allChr_n_scaff_7Apr19.mod --msa-format MAF processedMAFS/sortedMAFS/s_filt_$chr.maf > phylopoutput/all/CONACC/${chr}_allbranch_CONACC_ONmod_7Apr.wig

	phyloP --wig-scores --method LRT --chrom $chr --mode GERP msa_view_output_ON/neutralNiloticus_allChr_n_scaff_7Apr19.mod --msa-format MAF processedMAFS/sortedMAFS/s_filt_$chr.maf > phylopoutput/all/GERP/${chr}_allbranch_GERP_ONmod_7Apr.wig
done < chrm_from_ON.txt

## Step 11 - Converting wig files to bed
for fn in `ls phylopoutput/all/ACC/*.wig`;  
do 
	bedops/bin/wig2bed < $fn > ${fn%%.*}.bed;
done

for fn in `ls phylopoutput/all/CONACC/*.wig`;  
do 
	bedops/bin/wig2bed < $fn > ${fn%%.*}.bed;
done

for fn in `ls phylopoutput/all/GERP/*.wig`;  
do 
	bedops/bin/wig2bed < $fn > ${fn%%.*}.bed;
done

## Step 12 - Concantenate and sort the bed files
# All branch -  ACC
cd phylopoutput/all/ACC/ 
cat *.bed > sorted_phylop_ONmod_Allbranch_ACC_1.bed
sort -k1,1 -k2,2n sorted_phylop_ONmod_Allbranch_ACC_1.bed > sorted_phylop_ONmod_Allbranch_ACC.bed
rm sorted_phylop_ONmod_Allbranch_ACC_1.bed

# All branch - CONACC
cd phylopoutput/all/CONACC/  
cat *.bed > sorted_phylop_ONmod_Allbranch_CONACC_1.bed
sort -k1,1 -k2,2n sorted_phylop_ONmod_Allbranch_CONACC_1.bed > sorted_phylop_ONmod_Allbranch_CONACC.bed
rm sorted_phylop_ONmod_Allbranch_CONACC_1.bed

# All branch - GERP
cd phylopoutput/all/GERP/  
cat *.bed > sorted_phylop_ONmod_Allbranch_GERP_1.bed
sort -k1,1 -k2,2n sorted_phylop_ONmod_Allbranch_GERP_1.bed > sorted_phylop_ONmod_Allbranch_GERP.bed
rm sorted_phylop_ONmod_Allbranch_GERP_1.bed


 


