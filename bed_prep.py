# module load bedtools/2.21.0
# module unload python
# module load python/2.7.6
# module load py_packages/2.7
# module unload gcc
# module load ucsctools/linux.x86_64
# python


module load bedtools/2.21.0
cd chdiTrios/Felix/chipseq/

# remove first line if necessary
sed -i 1d Chip_3_WT_input.bed
sed -i 's/^/chr/' Chip_3_WT_input.bed
FILE="Chip_3_WT_input"
FILE="Chip_4_WT_input"
FILE="Chip_7_MT_input"
FILE="Chip_9_MT_input"
sed 1d "$FILE.bed" | sed 's/^/chr/' > "$FILE.sorted.bed"
# add chr to start of every line if necessary
# sort and merge if necessary
sort -V -k1,1 -k2,2 Chip_3_WT_input.bed > Chip_3_WT_input.sorted.bed

gunzip GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.bed.gz
sort -V -k1,1 -k2,2 GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.bed > GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.sorted.bed
bedtools merge GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.sorted.bed > GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.merged.sorted.bed

# find jaccard index. a and b order does not matter
bedtools jaccard -a Chip_3_WT_input.sorted.bed -b E086_15_coreMarks_1.sorted.bed
# 0.213752
bedtools jaccard -a E086_15_coreMarks_1.sorted.bed -b Chip_3_WT_input.sorted.bed
# 0.213752
time bedtools jaccard -a Chip_3_WT_input.sorted.bed -b GSM906406_UCSD.Left_Ventricle.H3K4me3.STL001.sorted.bed.gz
# 30sec, same result with merged.sorted only 4sec
time bedtools jaccard -a Chip_3_WT_input.sorted.bed -b GSM906406_UCSD.Left_Ventricle.H3K4me3.STL001.merged.sorted.bed
# 0.0145867
time bedtools jaccard -a Chip_3_WT_input.sorted.bed -b GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.merged.sorted.bed
# 14715506

gunzip GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.bed.gz
sort -V -k1,1 -k2,2 GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.bed > GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.sorted.bed
bedtools merge GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.sorted.bed > GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353.sorted.merged.bed


def merge_bed(bed_name):
    """ MERGES a bed file
    """
    pybedtools.set_tempdir('/sc/orga/scratch/richtf01')
    bed_in = bed_name + '.sorted.bed'
    bed_out = bed_name + '.merged.sorted.bed'
    if not os.path.isfile(bed_out):
        bed = BedTool(bed_in)
        print "Merging " + bed_in + "..."
        bed_merged = bed.merge()
        bed_merged.saveas(bed_out)
        print bed_name + " done!"
    else:
        print bed_out + " already merged"

def calculate_length(bed_name):
    """read bedfile per line as running sum
    """
    print "Counting", bed_name
    bed_in = bed_name
    running_sum = 0
    with open(bed_in, 'r') as bed_iter:
        for bed_line in csv.reader(bed_iter, delimiter="\t"):
            running_sum += (float(bed_line[2]) - float(bed_line[1]))
    return running_sum


def intersect_bed(bed_new_name, bed_db_name):
    """Intersect new and known ChIP-Seq data
    """
    pybedtools.set_tempdir('/sc/orga/scratch/richtf01')
    if not os.path.isfile(bed_new_name + '.' + bed_db_name + '.bed'):
        bed_new = BedTool(bed_new_name + '.sorted.bed')
        bed_db = BedTool(bed_db_name + '.sorted.bed')
        print "Intersecting", bed_new_name, bed_db_name
        bed_overlap = bed_new.intersect(bed_db)
        bed_overlap.saveas(bed_new_name + '.' + bed_db_name + '.bed')
        print bed_new_name, "done for", bed_db_name
    else:
        print "Already intersected", bed_new_name, bed_db_name


def sort_bed(bed_name):
    """ SORTS a bed file
    """
    bed_in = bed_name + '.bed'
    bed_out = bed_name + '.sorted.bed'
    if not os.path.isfile(bed_out):
        print "Sorting " + bed_in + "... "
        sort_cmd = ("sort -V -k1,1 -k2,2 %s > %s"
            % (bed_in, bed_out))
        print sort_cmd
        subprocess.call(sort_cmd, shell = True)
        print bed_name + " sorted!"
    else:
        print bed_out + " already sorted"


import glob
import csv
from pybedtools import BedTool
import pybedtools
import os
import subprocess


bed_name = "GSM1127085_UCSF-UBC.Breast_Fibroblast_Primary_Cells.H3K4me3.RM071"
bed_name = "GSM772844_BI.iPS-20b.H3K4me3.DNA_Lib_353"
bed_name = "GSM906409_UCSD.Left_Ventricle.H3K9me3.STL001"
merge_bed(bed_name)

bed_name = "Chip_3_WT_input.sorted.bed"
bed_new_name = "Chip_3_WT_input"
bed_db_name = "GSM772735_BI.Fetal_Heart.H3K4me3.UW_H23914.merged"
bed_db_name = "GSM906406_UCSD.Left_Ventricle.H3K4me3.STL001.merged"
bed_db_name = "GSM910580_UCSD.Left_Ventricle.H3K4me3.STL003.merged"
bed_db_name = "GSM621423_BI.iPS-20b.H3K4me3.Lib_MC_20100119_05--ChIP_MC_20100113_05_hiPS-20b_H3K4Me3.merged"
bed_db_name = "E095_15_coreMarks_1"
bed_db_name = "regions_prom_E095"
bed_db_name = "GSM908951_UCSD.Left_Ventricle.H3K27ac.STL001.merged"
bed_db_name = "GSM906404_UCSD.Left_Ventricle.H3K4me1.STL001.merged"
bed_db_name = "GSM906396_UCSD.Left_Ventricle.H3K27ac.STL003.merged"
bed_pct_dict = {}
bed_db_name = "regions_enh_E095"
merge_bed(bed_db_name)

bed_new_iterable = glob.iglob('/hpc/users/richtf01/chipseq/H3k4me3_peaks/*.bed')
bed_new_name = bed_new_iterable.next()[:-4]

for bed_new_name in bed_new_iterable:
    sort_bed(bed_new_name[:-4])

bed_pct_dict = {}
bed_new_iterable = glob.iglob('input/*.sorted.bed')
for bed_new_name in bed_new_iterable:
    bed_new_name = bed_new_name[:-11]
    print bed_new_name[6:]
    bed_new_length = calculate_length(bed_new_name + '.sorted.bed')
    bed_iterable = glob.iglob('*merged.sorted.bed') #wgEncodeBroadHmm
    for bed_name in bed_iterable:
        bed_db_name = bed_name[:-11]
        print bed_db_name
        intersect_bed(bed_new_name, bed_db_name)
        bed_inter_length = calculate_length(bed_new_name + '.' + bed_db_name + '.bed')
        bed_pct_dict[bed_new_name + '___' + bed_db_name] = bed_inter_length/bed_new_length


file_name = "/hpc/users/richtf01/chipseq/pct_overlap.txt"
with open(file_name, 'w+') as f:
    w = csv.DictWriter(f, bed_pct_dict.keys() )
    w.writeheader()
    w.writerow(bed_pct_dict)

# 33445740.0
14715506.0/33445740.0
0.439981
10650494.0/33445740.0


del bed_pct_dict[bed_db_name]
