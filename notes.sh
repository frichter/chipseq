
cp epi_roadmap_heart/source_files/heart_fetal/H3K4me3/GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.bed.gz ../../chipseq/
cp epi_roadmap_heart/source_files/ ../../chipseq/
cp epi_roadmap_heart/source_files/heart_fetal/H3K4me3/GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.bed.gz ../../chipseq/
cp epi_roadmap_heart/source_files/heart_fetal/H3K4me3/GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.bed.gz ../../chipseq/
cp epi_roadmap_heart/source_files/heart_fetal/H3K4me3/GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.bed.gz ../../chipseq/

cp ../wgs/bed_annotations/epi_roadmap_heart/GSM906406_UCSD.Left_Ventricle.H3K4me3.STL003.merged.sorted.bed .
cp ../wgs/bed_annotations/epi_roadmap_heart/GSM772735_BI
cp ../wgs/bed_annotations/epi_roadmap_heart/GSM906406
cp ../wgs/bed_annotations/epi_roadmap_heart/GSM910580_UCSD.Left_Ventricle.H3K4me3.STL003.merged.sorted.bed .

fibroblast_primary_cell_line

# E013 mesoderm
# E004 Mesoendoderm
# download RA and RV
# figure out what the difference is between all iPS cells: run with 1 then take all 5
E020 iPS-20b

# E004, E013, E020, E083, E095, E105
# careful with right atrial E104 due to antibody issues
# E086 is fetal KIDNEY but seems to have permated into analysis

wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_sample/iPS-20b_cell_line/H3K4me3/GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.bed.gz

wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/roadmapepigenomics/by_sample/iPS-20b_cell_line/H3K4me3/GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.bed.gz


wget http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/

http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/regions_prom_E004.bed.gz
http://egg2.wustl.edu/roadmap/data/byDataType/dnase/BED_files_prom/regions_prom_E095.bed.gz

sort -V -k1,1 -k2,2

module load bedtools/2.21.0

sed 's/^/chr/' Chip_3_WT_input.bed | head
cd chdiTrios/Felix/chipseq/
bedtools jaccard -a Chip_3_WT_input.bed -b E086_15_coreMarks_1.bed

sed 1d Chip_3_WT_input.bed | head


Integrative analysis of 350 patients with sporadic congenital heart defects (CHD)
