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
cat coding_exons_plus_popsnps.bed | bedtools slop -b 5 -g /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.contigLengths.txt  | bedtools merge > coding_exons_plus_popsnps.pad5bp.bed

cat coding_exons_plus_popsnps.pad5bp.bed | span 
# 161969

####

# split target regions into ~5kbp chunks to parallelize
mkdir target_split

python $JP/proj/capdesign/split_target_bed.py \
    --inBed coding_exons_plus_popsnps.pad5bp.bed     --fxnChunkbedNumToFn "lambda i:'target_split/design_merge.%03d.bed'%i" \
    --outKey design_chunks.k --targetBpPerChunk 5000 --maxbpSingleInterval 2000

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

export GEN_CANDS_DESIGN_OPTS_STAGE1=" --ligArmLengthRange 23,29 --extArmLengthRange 23,30  --gapFillRange 100,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"

export GEN_PAIRS_DESIGN_OPTS_STAGE1=" --gapFillRange 100,120 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly --generateMipSeqFxn \"lambda picker:partial(picker.makeseq_dsmip_eari,adaptorLeft=yoonOlEvenFshort,adaptorRight=yoonOlEvenRrcshort,normalizeLength=150 )\""

# this runs up until but not including picking.
pipeline_template.sh > ${BASE_DIR}/out.log 2>${BASE_DIR}/out.err


#<- editing here need to update df.sort 



# more stringent with flakning 
python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.nomito.pad5bp.bed \
    --coverage 1 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_12kpbs_params2.h5 \
    --targetNumProbes 12400 &

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.nomito.pad5bp.bed \
    --coverage 1 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_18kpbscov1_params2.h5 \
    --targetNumProbes 18000 &

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.nomito.pad5bp.bed \
    --coverage 2 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_18kpbs_params2.h5 \
    --targetNumProbes 18000 &

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.nomito.pad5bp.bed \
    --coverage 3 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_45kpbs_params2.h5 \
    --targetNumProbes 45000 &

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed coding_exons_plus_popsnps.nomito.pad5bp.bed \
    --coverage 4 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design/pairs_all_joined.h5\
    --inArmStore  design/arms_filtuniq.h5  \
    --outPairStore design/pairs_picked_90kpbs_params2.h5 \
    --targetNumProbes 90000 &

for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs; do # 90kpbs; do
    python $JP/proj/capdesign/pairtbl_to_bed.py  \
        --table pairsPlus \
        --armStore design/arms_filtuniq.h5\
        --inStore design/pairs_picked_${arrayversion}_params2.h5 \
        --bedOut design/pairs_picked_${arrayversion}_params2.bed&
done

for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs; do # 90kpbs; do
    python $JP/proj/capdesign/picked_panel_check_coverage.py\
        --inTargetBed coding_exons_plus_popsnps.bed\
        --inMipBed design/pairs_picked_${arrayversion}_params2.bed\
        --out design/pairs_picked_${arrayversion}_params2.stats.txt\
        --outByTarg design/pairs_picked_${arrayversion}_params2.bytarg.stats.txt\
        --libname design_${arrayversion} &
done

for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs; do # 90kpbs; do
    python $JP/proj/capdesign/evalulate_design.py \
        --inStorePairs design/pairs_all_joined.h5 \
        --inStorePairsPicked design/pairs_picked_${arrayversion}_params2.h5  \
        --outPlotBase design/pairs_picked_${arrayversion}_params2 \
        --considerUnpickedFracOvl 0.75 &
done

for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs; do # 90kpbs; do
    cut -f1,7,8 design/pairs_picked_${arrayversion}_params2.bed > design/pairs_picked_${arrayversion}_params2.gapfill.bed;
done


for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs; do # 90kpbs; do
    python $JP/proj/capdesign/pair_table_out.py \
        --inStorePairs design/pairs_picked_${arrayversion}_params2.h5 \
        --tablePairs pairsPlus\
        --armAlignHitsSuffix mrfAndBwa \
        --refFasta /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
        --gcFlankingBp 100 --gcFlankingBp 1000 --gcFlankingBp 2500 --gcFlankingBp 5000 --gcFlankingBp 10000 --gcFlankingBp 20000 --gcFlankingBp 100000 --gcFlankingBp 500000  \
        --inStoreArms ${BASE_DIR}/arms_filtuniq.h5\
        --mipTableOut ${BASE_DIR}/pairs_picked_${arrayversion}_params2.plus.txt&
done

for arrayversion in 12kpbs 18kpbscov1 18kpbs 45kpbs 90kpbs; do
    bedtools intersect -a <( cat design/pairs_picked_${arrayversion}_params2.bed|awk '{print $1 "\t" $2 "\t" $7}' ) -b <(cat  design/pairs_picked_${arrayversion}_params2.bed |awk '{print $1 "\t" $8 "\t" $3}') | sort -k1,1 -k2,2n | awk '{if ($3 - $2 > 10) print $_}' | sort -k1,1 -k2,2n | bedtools merge  > design/pairs_picked_${arrayversion}_params2.daisychain.bed;
    
    bedtools subtract -a <(sort -k1,1 -k2,2n coding_exons_plus_popsnps.bed) -b <( bedtools merge -i <(sort -k1,1 -k2,2n design/pairs_picked_${arrayversion}_params2.gapfill.bed ) )  > design/not_covered_${arrayversion}_param2.bed;
done

cat design/*params?.stats.txt | column -t

"""
libname            perbase_average_coverage  perbase_average_coverage_distinct_probe  frac_gapfills_off_target  total_probes  total_probes_unique  total_probes_repeated
design_12kpbs      1.61184858995             1.61184858995                            0.122060907586            12648         12648                0
design_18kpbscov1  2.24315897858             2.22887168733                            0.139261812466            18000         17898                102
design_18kpbs      3.22228137202             2.86047895669                            0.150486462066            25115         22299                2816
design_45kpbs      5.57524596132             4.28807952887                            0.22364080417             45001         35038                9963
design_90kpbs      10.8639902319             7.27660740681                            0.386511021706            90003         61399                28604
"""


# by comparison, sampson v3 array:
"""
libname perbase_average_coverage        perbase_average_coverage_distinct_probe frac_gapfills_off_target        total_probes    total_probes_unique     total_probes_repeated
designv3        2.62352583917   2.60413813936   0.153386729367  12400   12318   82
"""  

# --> OK, we will go with 18kpbs cov2 (22k probes)


#######################################################
# DESIGN PROBES AGAINST ADDL GENES



# add:
#   APOA1 (was misnamed)
#   NPPB
#   CLCN6   
#   GUCY1B3
#   NPR1
#   ANGPTL8

# microsatellites
# /nfs/kitzman2/jacob/proj/mip_design_pipe/microsat_201508/microsats_slop.bed

echo "\
APOA1
NPPB
CLCN6  
GUCY1B3
NPR1
ANGPTL8" > scd_addl_genes.txt

python $JP/proj/capdesign/gene_list_to_bed_hg19.py \
--inGeneList scd_addl_genes.txt \
--knownGeneTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_knownGene.txt.gz \
--refGeneTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_refGene.txt.gz  \
--gencodeTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_gencodeGeneV24lift37.txt.gz \
--outBed scd_addl_genes_coding_exons.bed \
--inclExclusivelyNoncoding \
> get_addl_gene_coords.log

cat scd_addl_genes_coding_exons.bed | perl -pi -e 's/chr//;s/^M\b/MT/'  > scd_addl_genes_coding_exons.chromnamefix.bed;

cp /nfs/kitzman2/jacob/proj/mip_design_pipe/microsat_201508/microsats_slop.bed ./

cat scd_addl_genes_coding_exons.chromnamefix.bed microsats_slop.bed | sort -k1,1 -k2,2n >  addlgenes_plus_usats.bed

cat addlgenes_plus_usats.bed | bedtools merge | span
# 13401

cat addlgenes_plus_usats.bed | bedtools slop -b 5 -g /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.contigLengths.txt  | bedtools merge > addlgenes_plus_usats.pad5bp.bed

cat addlgenes_plus_usats.pad5bp.bed | span 
# 14791

mkdir target_split_addlgenes

python $JP/proj/capdesign/split_target_bed.py \
--inBed addlgenes_plus_usats.pad5bp.bed\
--fxnChunkbedNumToFn "lambda i:'target_split_addlgenes/design_merge.%03d.bed'%i" \
--outKey addlgenes_design_chunks.k --targetBpPerChunk 1000 --maxbpSingleInterval 2000


export PATH=$PATH:/nfs/kitzman2/lab_software/flux_rhel6/jellyfish-2.1.4/bin/
module load mrfast/2.6.1.0

##
################################################################
################################################################
# longer arms - 23-29

mkdir design_addlgenes

export TMPDIR=/scratch/kitzman-lab/tmp/

export CPUlo=24
export CPU=32

export RUN_KEY_IN=addlgenes_design_chunks.k
export RUN_KEY_OUT=addlgenes_design_chunks.pass1.k

export BASE_DIR=design_addlgenes
export GENOME_PATH=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa
export GENOME_DICT=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict

export JELLYFISH_PATH=/nfs/kitzman2/jacob/proj/mip_design_pipe/jellyfish/hs37d5_k15_canon.jf
export BWA_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/bwa-0.7.10/human/hs37d5/hs37d5.fa
export MRF_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/mrfast-2.6/ucsc_hg19/ws12/hg19_rmskToN.fa
export MRF_CHROMNAME_TRANSLATE_TBL=/nfs/kitzman2/lab_common/refs/human/hg19_to_g1k.txt

export SNP_VCF=/nfs/kitzman2/lab_common/1000genomes/downloaded_20141219/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz

export COL_BASENAME=chunkidx
export COL_TARGET_BED=targetbed

export DO_BED_ALLARMS=yes

export ONLY_PICK=no

export TRAIN_PATH=diaorch_trainsets/trainsets/

# 23,29 because script assumes (this is in common.py and not exposed) that the max ligation arm len is 30 but that counts the add'l stuffer base.
export GEN_CANDS_DESIGN_OPTS_STAGE1=" --ligArmLengthRange 23,29 --extArmLengthRange 23,30  --gapFillRange 100,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"

export GEN_PAIRS_DESIGN_OPTS_STAGE1=" --gapFillRange 100,120 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly --generateMipSeqFxn \"lambda picker:partial(picker.makeseq_dsmip_eari,adaptorLeft=yoonOlEvenFshort,adaptorRight=yoonOlEvenRrcshort,normalizeLength=150 )\""


# this runs up until but not including picking.
sh $JP/proj/capdesign/pipelines/pipeline_design_v2_adaboostmodel.sh > ${BASE_DIR}/out.log 2>${BASE_DIR}/out.err


# 14791/792386 * 18000 = 335 < 400

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed addlgenes_plus_usats.pad5bp.bed \
    --coverage 2 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design_addlgenes/pairs_all_joined.h5\
    --inArmStore  design_addlgenes/arms_filtuniq.h5  \
    --outPairStore design_addlgenes/pairs_picked_18kpbs_params2.h5 \
    --targetNumProbes 400 &

python $JP/proj/capdesign/pairtbl_to_bed.py  \
    --table pairsPlus \
    --armStore design_addlgenes/arms_filtuniq.h5 \
    --inStore design_addlgenes/pairs_picked_18kpbs_params2.h5 \
    --bedOut design_addlgenes/pairs_picked_18kpbs_params2.bed&

cut -f1,7,8 design_addlgenes/pairs_picked_18kpbs_params2.bed > design_addlgenes/pairs_picked_18kpbs_params2.gapfill.bed

python $JP/proj/capdesign/pair_table_out.py \
    --inStorePairs design_addlgenes/pairs_picked_18kpbs_params2.h5 \
    --tablePairs pairsPlus\
    --armAlignHitsSuffix mrfAndBwa \
    --refFasta /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
    --gcFlankingBp 100 --gcFlankingBp 1000 --gcFlankingBp 2500 --gcFlankingBp 5000 --gcFlankingBp 10000 --gcFlankingBp 20000 --gcFlankingBp 100000 --gcFlankingBp 500000  \
    --inStoreArms design_addlgenes/arms_filtuniq.h5 \
    --mipTableOut design_addlgenes/pairs_picked_18kpbs_params2.plus.txt


bedtools subtract -a <(sort -k1,1 -k2,2n addlgenes_plus_usats.bed) -b <( bedtools merge -i <(sort -k1,1 -k2,2n design_addlgenes/pairs_picked_18kpbs_params2.gapfill.bed ) )  > design_addlgenes/not_covered_18kpbs_params2.gapfill.bed;

#######################################################
# DESIGN MITOPROBES
# 
# will need to update GEN_PAIRS_DEISGN_OPTS_STAGE1

cat coding_exons_plus_popsnps.bed | grep MT  > mito_targets.bed

cat mito_targets.bed | bedtools merge | span
# 4048

cat mito_targets.bed | bedtools slop -b 5 -g /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.contigLengths.txt  | bedtools merge > mito_targets.pad5bp.bed

cat mito_targets.pad5bp.bed | span 
# 4104

mkdir target_split_mito

python $JP/proj/capdesign/split_target_bed.py \
--inBed mito_targets.pad5bp.bed\
--fxnChunkbedNumToFn "lambda i:'target_split_mito/design_merge.%03d.bed'%i" \
--outKey mito_design_chunks.k --targetBpPerChunk 1000 --maxbpSingleInterval 2000


##
################################################################
################################################################
# longer arms - 23-29

mkdir design_mito

export TMPDIR=/scratch/kitzman-lab/tmp/

export CPUlo=24
export CPU=32

export RUN_KEY_IN=mito_design_chunks.k
export RUN_KEY_OUT=mito_design_chunks.pass1.k

export BASE_DIR=design_mito
export GENOME_PATH=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa
export GENOME_DICT=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict

export JELLYFISH_PATH=/nfs/kitzman2/jacob/proj/mip_design_pipe/jellyfish/hs37d5_k15_canon.jf
export BWA_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/bwa-0.7.10/human/hs37d5/hs37d5.fa
export MRF_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/mrfast-2.6/ucsc_hg19/ws12/hg19_rmskToN.fa
export MRF_CHROMNAME_TRANSLATE_TBL=/nfs/kitzman2/lab_common/refs/human/hg19_to_g1k.txt

export SNP_VCF=/nfs/kitzman2/lab_common/1000genomes/downloaded_20141219/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz

export COL_BASENAME=chunkidx
export COL_TARGET_BED=targetbed

export DO_BED_ALLARMS=yes

export ONLY_PICK=no

export TRAIN_PATH=diaorch_trainsets/trainsets/

# 23,29 because script assumes (this is in common.py and not exposed) that the max ligation arm len is 30 but that counts the add'l stuffer base.
export GEN_CANDS_DESIGN_OPTS_STAGE1=" --ligArmLengthRange 23,29 --extArmLengthRange 23,30  --gapFillRange 100,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"

export GEN_PAIRS_DESIGN_OPTS_STAGE1=" --gapFillRange 100,120 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly --generateMipSeqFxn \"lambda picker:partial(picker.makeseq_dsmip_eari,adaptorLeft=yoonOlOddFshort,adaptorRight=yoonOlOddRrcshort,normalizeLength=150 )\""


# this runs up until but not including picking.
sh $JP/proj/capdesign/pipelines/pipeline_design_v2_adaboostmodel.sh > ${BASE_DIR}/out.log 2>${BASE_DIR}/out.err


# 4104/792386 * 18000 = 93 < 100

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed mito_targets.pad5bp.bed \
    --coverage 2 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design_mito/pairs_all_joined.h5\
    --inArmStore  design_mito/arms_filtuniq.h5  \
    --outPairStore design_mito/pairs_picked_18kpbs_params2.h5 \
    --targetNumProbes 100 

python $JP/proj/capdesign/pairtbl_to_bed.py  \
    --table pairsPlus \
    --armStore design_mito/arms_filtuniq.h5 \
    --inStore design_mito/pairs_picked_18kpbs_params2.h5 \
    --bedOut design_mito/pairs_picked_18kpbs_params2.bed

cut -f1,7,8 design_mito/pairs_picked_18kpbs_params2.bed > design_mito/pairs_picked_18kpbs_params2.gapfill.bed

python $JP/proj/capdesign/pair_table_out.py \
    --inStorePairs design_mito/pairs_picked_18kpbs_params2.h5 \
    --tablePairs pairsPlus\
    --armAlignHitsSuffix mrfAndBwa \
    --refFasta /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
    --gcFlankingBp 100 --gcFlankingBp 1000 --gcFlankingBp 2500 --gcFlankingBp 5000 --gcFlankingBp 10000 --gcFlankingBp 20000 --gcFlankingBp 100000 --gcFlankingBp 500000  \
    --inStoreArms design_mito/arms_filtuniq.h5 \
    --mipTableOut design_mito/pairs_picked_18kpbs_params2.plus.txt


bedtools subtract -a <(sort -k1,1 -k2,2n mito_targets.bed) -b <( bedtools merge -i <(sort -k1,1 -k2,2n design_mito/pairs_picked_18kpbs_params2.gapfill.bed ) )  > design_mito/not_covered_18kpbs_params2.gapfill.bed;


#######################################################
# DESIGN PROBES AGAINST FINAL LIST OF EXTRA GENES

echo "\
ITGA2B
LIPC
PLOD3
F12
F7" > scd_addl_genes_batch2.txt

python $JP/proj/capdesign/gene_list_to_bed_hg19.py \
--inGeneList scd_addl_genes_batch2.txt \
--knownGeneTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_knownGene.txt.gz \
--refGeneTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_refGene.txt.gz  \
--gencodeTbl /nfs/kitzman2/lab_common/annots/ucsc/hg19/hg19_gencodeGeneV24lift37.txt.gz \
--outBed scd_addl_genes_batch2_coding_exons.bed \
--inclExclusivelyNoncoding \
> get_addl_gene_batch2_coords.log

cat scd_addl_genes_batch2_coding_exons.bed | perl -pi -e 's/chr//;s/^M\b/MT/'  > scd_addl_genes_batch2_coding_exons.chromnamefix.bed;

# need to add extra rsids
# rs1077190
# rs3825807
# rs6601627
# rs118065347


echo "\
chr8    11799863    11799864    rs118065347 0   +
chr8    11778802    11778803    rs6601627   0   +
chr15   79071383    79071384    rs1077190   0   -
chr15   79089110    79089111    rs3825807   0   -" | perl -pi -e 's/ +/\t/gi' | cut -f1-3| perl -pi -e 's/chr//;s/^M\b/MT/'  > addl_willer_snps.bed

cat scd_addl_genes_batch2_coding_exons.chromnamefix.bed addl_willer_snps.bed | sort -k1,1 -k2,2n >  batch2_targets.bed

cat batch2_targets.bed | bedtools merge | span
# 10655

cat batch2_targets.bed | bedtools slop -b 5 -g /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.contigLengths.txt  | bedtools merge > batch2_targets.pad5bp.bed

cat batch2_targets.pad5bp.bed | span 
# 11525



mkdir target_split_batch2

python $JP/proj/capdesign/split_target_bed.py \
--inBed batch2_targets.pad5bp.bed\
--fxnChunkbedNumToFn "lambda i:'target_split_batch2/design_merge.%03d.bed'%i" \
--outKey addlgenes_design_chunks.k --targetBpPerChunk 1000 --maxbpSingleInterval 2000


export PATH=$PATH:/nfs/kitzman2/lab_software/flux_rhel6/jellyfish-2.1.4/bin/
module load mrfast/2.6.1.0

##
################################################################
################################################################
# longer arms - 23-29

mkdir design_batch2

export TMPDIR=/scratch/kitzman-lab/tmp/

export CPUlo=24
export CPU=32

export RUN_KEY_IN=addlgenes_design_chunks.k
export RUN_KEY_OUT=addlgenes_design_chunks.pass1.k

export BASE_DIR=design_batch2
export GENOME_PATH=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa
export GENOME_DICT=/nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict

export JELLYFISH_PATH=/nfs/kitzman2/jacob/proj/mip_design_pipe/jellyfish/hs37d5_k15_canon.jf
export BWA_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/bwa-0.7.10/human/hs37d5/hs37d5.fa
export MRF_INDEX_PATH=/nfs/kitzman2/lab_common/alignment_index/mrfast-2.6/ucsc_hg19/ws12/hg19_rmskToN.fa
export MRF_CHROMNAME_TRANSLATE_TBL=/nfs/kitzman2/lab_common/refs/human/hg19_to_g1k.txt

export SNP_VCF=/nfs/kitzman2/lab_common/1000genomes/downloaded_20141219/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz

export COL_BASENAME=chunkidx
export COL_TARGET_BED=targetbed

export DO_BED_ALLARMS=yes

export ONLY_PICK=no

export TRAIN_PATH=diaorch_trainsets/trainsets/

# 23,29 because script assumes (this is in common.py and not exposed) that the max ligation arm len is 30 but that counts the add'l stuffer base.
export GEN_CANDS_DESIGN_OPTS_STAGE1=" --ligArmLengthRange 23,29 --extArmLengthRange 23,30  --gapFillRange 100,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"

export GEN_PAIRS_DESIGN_OPTS_STAGE1=" --gapFillRange 100,120 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly --generateMipSeqFxn \"lambda picker:partial(picker.makeseq_dsmip_eari,adaptorLeft=yoonOlEvenFshort,adaptorRight=yoonOlEvenRrcshort,normalizeLength=150 )\""

# this runs up until but not including picking.
sh $JP/proj/capdesign/pipelines/pipeline_design_v2_adaboostmodel.sh > ${BASE_DIR}/out.log 2>${BASE_DIR}/out.err


# 11525/792386 * 18000 = 261 < 300

python $JP/proj/capdesign/greedy_pick_mip_set_v2.py \
    --targetBed batch2_targets.pad5bp.bed \
    --coverage 2 \
    --scorePenaltyForOverlapExisting 0.25 \
    --scorePenaltyForOverlapExistingPass2 0.70 \
    --scorePenaltyForNonOverlapTarget 0.65 \
    --maxBpDaisychain 7\
    --inPairStore design_batch2/pairs_all_joined.h5\
    --inArmStore  design_batch2/arms_filtuniq.h5  \
    --outPairStore design_batch2/pairs_picked_18kpbs_params2.h5 \
    --targetNumProbes 300 &

python $JP/proj/capdesign/pairtbl_to_bed.py  \
    --table pairsPlus \
    --armStore design_batch2/arms_filtuniq.h5 \
    --inStore design_batch2/pairs_picked_18kpbs_params2.h5 \
    --bedOut design_batch2/pairs_picked_18kpbs_params2.bed&

cut -f1,7,8 design_batch2/pairs_picked_18kpbs_params2.bed > design_batch2/pairs_picked_18kpbs_params2.gapfill.bed

python $JP/proj/capdesign/pair_table_out.py \
    --inStorePairs design_batch2/pairs_picked_18kpbs_params2.h5 \
    --tablePairs pairsPlus\
    --armAlignHitsSuffix mrfAndBwa \
    --refFasta /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
    --gcFlankingBp 100 --gcFlankingBp 1000 --gcFlankingBp 2500 --gcFlankingBp 5000 --gcFlankingBp 10000 --gcFlankingBp 20000 --gcFlankingBp 100000 --gcFlankingBp 500000  \
    --inStoreArms design_batch2/arms_filtuniq.h5 \
    --mipTableOut design_batch2/pairs_picked_18kpbs_params2.plus.txt

bedtools subtract -a <(sort -k1,1 -k2,2n batch2_targets.pad5bp.bed ) -b <( bedtools merge -i <(sort -k1,1 -k2,2n design_batch2/pairs_picked_18kpbs_params2.bed ) )  > design_batch2/not_covered_18kpbs_params2.gapfill.bed;


export JP=/nfs/kitzman1/jacob/dev/jkproj/


python $JP/proj/capdesign/pairtbl_to_oligo_list.py \
--inStore design/pairs_picked_18kpbs_params2.h5\
--table pairsPlus\
> design/pairs_picked_18kpbs_params2.oligos.txt

python $JP/proj/capdesign/pairtbl_to_oligo_list.py \
--inStore design_addlgenes/pairs_picked_18kpbs_params2.h5\
--table pairsPlus\
> design_addlgenes/pairs_picked_18kpbs_params2.oligos.txt

python $JP/proj/capdesign/pairtbl_to_oligo_list.py \
--inStore design_mito/pairs_picked_18kpbs_params2.h5\
--table pairsPlus\
> design_mito/pairs_picked_18kpbs_params2.oligos.txt

python $JP/proj/capdesign/pairtbl_to_oligo_list.py \
--inStore design_batch2/pairs_picked_18kpbs_params2.h5\
--table pairsPlus\
> design_batch2/pairs_picked_18kpbs_params2.oligos.txt




python $PP/proj/pals/libprep_v1/pad_with_junk.py \
--libIn <(echo "SEQ"; cut -f2 design/pairs_picked_18kpbs_params2.oligos.txt; cut -f2 design/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design/pairs_picked_18kpbs_params2.oligos.txt; cut -f2 design_addlgenes/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_addlgenes/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_addlgenes/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_mito/pairs_picked_18kpbs_params2.oligos.txt; cut -f2 design_mito/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_mito/pairs_picked_18kpbs_params2.oligos.txt; cut -f2 design_batch2/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_batch2/pairs_picked_18kpbs_params2.oligos.txt ; cut -f2 design_batch2/pairs_picked_18kpbs_params2.oligos.txt ) \
--libOut /dev/stdout \
--desiredLength 150 --randomPad --padding3p |grep -v SEQ> scd_v1_20170310.pad150.txt

#######################################################
# finally join_mutli_design_tables to make a final merged design table.

python /nfs/kitzman1/jacob/dev/mips_pipeline/novc/mimips/join_multi_design_tables.py  --mipTableIn design/pairs_picked_18kpbs_params2.plus.txt --inDesignName main --mipTableIn design_addlgenes/pairs_picked_18kpbs_params2.plus.txt --inDesignName addl --mipTableIn design_mito/pairs_picked_18kpbs_params2.plus.txt --inDesignName mito --inDesignName batch2 --mipTableIn  design_batch2/pairs_picked_18kpbs_params2.plus.txt --mipTableOut ./pairs_picked_joined.plus.txt


#######################################################
# make sure scoring model still works

cp /nfs/kitzman2/jacob/proj/mip_design_pipe/sampsonv3_20160726/designv3/arms_filtuniq.h5 $TMPDIR/sampsonv3_arms_thinned.h5

cp /nfs/kitzman2/jacob/proj/mip_design_pipe/sampsonv3_20160726/designv3/pairs_picked.h5 $TMPDIR/sampsonv3_pairs_picked.h5

python $JP/proj/capdesign/rescore_pairs.py \
--genome /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
--genomeDict /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict \
--inArmStore $TMPDIR/sampsonv3_arms_thinned.h5 \
--inPairStore $TMPDIR/sampsonv3_pairs_picked.h5 \
--outPairStore $TMPDIR/sampsonv3_pairs_picked.rescored.h5 \
--hitTblSufifx mrfAndBwa \
--snpTblSuffix 1kg \
--learningModel AdaBoost \
--trainingMode preTrained \
--trainSetPath diaorch_trainsets/trainsets/  \
--setPrefix trainset \
>rescore.log


cp /nfs/kitzman2/jacob/proj/mip_design_pipe/camperv1_20151130/design1_2_0p25_0p25_7/arms_filtuniq.h5  $TMPDIR/camperv1_arms_thinned.h5

cp /nfs/kitzman2/jacob/proj/mip_design_pipe/camperv1_20151130/design1_2_0p25_0p25_7/pairs_picked.h5 $TMPDIR/camperv1_pairs_picked.h5

python $JP/proj/capdesign/rescore_pairs.py \
--genome /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa \
--genomeDict /nfs/kitzman2/lab_common/refs/human/hs37d5/hs37d5.fa.dict \
--inArmStore $TMPDIR/camperv1_arms_thinned.h5 \
--inPairStore $TMPDIR/camperv1_pairs_picked.h5 \
--outPairStore $TMPDIR/camperv1_pairs_picked.rescored.h5 \
--hitTblSufifx mrfAndBwa \
--snpTblSuffix 1kg \
--learningModel AdaBoost \
--trainingMode preTrained \
--trainSetPath diaorch_trainsets/trainsets/  \
--setPrefix trainset \
>rescore.log
# !! recomputed scores are different than those computed at deisgn time!!

# model was trained on willer v2+v3.
# retrained

mkdir -p retrain/tables
cp diaorch_trainsets/trainsets/tables/trainset.tsv  retrain/tables/

python rdscript/scripts/pre-train_model.py \
--trainPath retrain/ \
--setPrefix trainset

# --> but weights are the same ... okay???

#
# see notebook 170301_check_stored_classifiers and 170301* (under sampson folder)
#
# tested the exisitng model vs sampsonv3, sampsonv1, camperv3 probe sets
# using  :
#  -design time stored scores for each
#  -manually recomputed scores for each
#  -commandline refit scores for each
# 
# --> all have AUC 0.68-0.83, excpet sampson v1 which has AUC 0.47
# 
# --> recomputed scores which have changed have nearly identical AUC to original scores -- must be getting transformed differently or something
