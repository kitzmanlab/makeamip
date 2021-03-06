# load the 'makeamip' conda environment that you previously set up
conda activate makeamip

# handy shortcut:
alias span='python -c "import numpy as np; import sys; x=np.array([ int(l.split(\"\t\")[2])-int(l.split(\"\t\")[1]) for l in sys.stdin]) ; print(x.sum())"'

echo "\
TP53
BRCA1
TTN
MSH2
MSH6
MLH1 
PMS2
XIST" > test_genelist.txt

# get coding exons for named genes

gene_list_to_bed_hg19.py \
    --inGeneList test_genelist.txt \
    --knownGeneTbl ${REFS}/annots/hg19_knownGene.txt.gz \
    --refGeneTbl ${REFS}/annots/hg19_refGene.txt.gz  \
    --gencodeTbl ${REFS}/annots/hg19_gencodeGeneV24lift37.txt.gz \
    --outBed test_genes_codingexons.bed \
    --inclExclusivelyNoncoding \
> get_gene_coords.log

# fix chromosome names to match hs37d5 reference
cat test_genes_codingexons.bed | perl -pi -e 's/chr//;s/^M\b/MT/'  > test_genes_codingexons.chromnamefix.bed

# include common SNP fingerprint panel 
cp ${REFS}/popsnps.target.bed ./

# merge 
cat popsnps.target.bed test_genes_codingexons.chromnamefix.bed | sort -k1,1 -k2,2n >  coding_exons_plus_popsnps.bed

# count total unique targeted bp 
cat coding_exons_plus_popsnps.bed | bedtools merge | span
# 155789

# pad out target regions by +/- 5 bp (min); some will end up having more padding dependent upon probe placement. 
cat coding_exons_plus_popsnps.bed | bedtools slop -b 5 -g ${REFS}/refs/hs37d5.contigLengths.txt  | bedtools merge > coding_exons_plus_popsnps.pad5bp.bed

cat coding_exons_plus_popsnps.pad5bp.bed | span 
# 161969

####

# split target regions into ~5kbp chunks to parallelize
mkdir target_split

split_target_bed.py \
    --inBed coding_exons_plus_popsnps.pad5bp.bed     --fxnChunkbedNumToFn "lambda i:'target_split/design_merge.%03d.bed'%i" \
    --outKey design_chunks.k --targetBpPerChunk 5000 --maxbpSingleInterval 2000


# need to have xa2multi.pl in path -- on my system this is included with bwa
export PATH=${PATH}:/opt/apps/rhel6/bwa/0.7.12/bin/
# add mrfast to path  on our system load module 
#export PATH=${PATH}:/mnt/userdata/annots_and_files/software/mrfast-2.6.1.0/
module load mrfast/2.6.1.0
# add jellyfish to path
export PATH=${PATH}:/nfs/kitzman2/lab_software/flux_rhel6/jellyfish-2.1.4/bin/
# point env var to picard - on our system load module
module load picard/2.9.0
# export PICARD_DIR=/mnt/userdata/tools 



##
################################################################
################################################################
# longer arms - 23-29

mkdir design

export CPUlo=24
export CPU=32

export RUN_KEY_IN=design_chunks.k
export RUN_KEY_OUT=design_chunks.pass1.k

export BASE_DIR=design

export TMPDIR=`readlink -f ${BASE_DIR}/temp/`
mkdir -p TMPDIR

export GENOME_PATH=${REFS}/refs/hs37d5.fa
export GENOME_DICT=${REFS}/refs/hs37d5.fa.dict

export JELLYFISH_PATH=${REFS}/indices/jellyfish/hs37d5_k15_canon.jf
export BWA_INDEX_PATH=${REFS}/indices/bwa/hs37d5.fa
export MRF_INDEX_PATH=${REFS}/indices/mrfast/hg19_rmskToN.fa
export MRF_CHROMNAME_TRANSLATE_TBL=${REFS}/refs/hg19_to_g1k.txt

export SNP_VCF=${REFS}/annots/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz

export COL_BASENAME=chunkidx
export COL_TARGET_BED=targetbed

export DO_BED_ALLARMS=yes

export ONLY_PICK=no

export TRAIN_PATH=${REFS}/parameters/trainsets/

# parameters that get passed into the candidate MIP arm and pair generation scripts. 

export GEN_CANDS_DESIGN_OPTS_STAGE1=" --ligArmLengthRange 25,32 --extArmLengthRange 25,32  --gapFillRange 100,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"

export GEN_PAIRS_DESIGN_OPTS_STAGE1=" --gapFillRange 100,120 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly --generateMipSeqFxn \"lambda picker:partial(picker.makeseq_dsmip_bbsi,adaptorLeft=new2019mipsPrep1f,adaptorRight=new2019mipsPrep1r_rc,tagLen=9,normalizeLength=169,maxLigArmLen=35 )\""

# this runs up until but not including picking.
pipeline_template.sh > ${BASE_DIR}/out.log 2>${BASE_DIR}/out.err

# now, given all probes in design database,
# select 10000 probes (--targetNumProbes) to make up this panel.
# 
# require each target base be covered by minimum of 1 probe's
# gap-fill (--coverage)
#
# 
greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.pad5bp.bed \
    --coverage 1 \
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_10kprobes_cov1.h5 \
    --targetNumProbes 10000 &

# make a bed file of probes. 
pairtbl_to_bed.py  \
    --table pairsPlus \
    --armStore design/arms_filtuniq.h5\
    --inStore design/pairs_picked_10kprobes_cov1.h5 \
    --bedOut design/pairs_picked_10kprobes_cov1.bed &

# picked_panel_check_coverage.py\
#     --inTargetBed coding_exons_plus_popsnps.pad5bp.bed \
#     --inMipBed design/pairs_picked_10kprobes_cov1.bed \
#     --out design/pairs_picked_10kprobes_cov1.stats.txt \
#     --outByTarg design/pairs_picked_10kprobes_cov1.stats_by_targ.txt \
#     --libname 10kprobes_cov1

evalulate_design.py \
    --inStorePairs design/pairs_all_joined.h5 \
    --inStorePairsPicked design/pairs_picked_10kprobes_cov1.h5  \
    --outPlotBase design/pairs_picked_10kprobes_cov1 \
    --considerUnpickedFracOvl 0.75 

# extract just the capture (gap-fill) regions from probe bed file
cut -f1,7,8 design/pairs_picked_10kprobes_cov1.bed > design/pairs_picked_10kprobes_cov1.gapfill.bed


# output final table w/ annotated probe pairs -- 
# this is the table require by the MIPS alignment pipeline. 

pair_table_out.py \
    --inStorePairs design/pairs_picked_10kprobes_cov1.h5 \
    --tablePairs pairsPlus\
    --armAlignHitsSuffix mrfAndBwa \
    --refFasta ${REFS}/refs/hs37d5.fa \
    --gcFlankingBp 100 --gcFlankingBp 1000 --gcFlankingBp 2500 --gcFlankingBp 5000 --gcFlankingBp 10000 --gcFlankingBp 20000 --gcFlankingBp 100000 --gcFlankingBp 500000  \
    --inStoreArms design/arms_filtuniq.h5\
    --mipTableOut design/FINAL_PROBE_TABLE_pairs_picked_10kprobes_cov1.txt &

# output final list of sequences to synthesize
pairtbl_to_oligo_list.py \
    --inStore design/pairs_picked_10kprobes_cov1.h5 \
    --table pairsPlus \
    > design/pairs_picked_10kprobes_cov1.oligos.txt
