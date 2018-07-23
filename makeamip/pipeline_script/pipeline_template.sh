#!/bin/bash
set -e

# expected environment variables   
RUN_KEY_IN=${RUN_KEY_IN?:"required parameter RUN_KEY_IN"}
RUN_KEY_OUT=${RUN_KEY_OUT?:"required parameter RUN_KEY_OUT"}

BASE_DIR=${BASE_DIR?:"required parameter BASE_DIR, will create dirs under here"}

GENOME_PATH=${GENOME_PATH?:"required parameter GENOME_PATH - path to ref genome fasta"}
GENOME_DICT=${GENOME_DICT?:"required parameter GENOME_PATH - path to ref genome seqdict"}
JELLYFISH_PATH=${JELLYFISH_PATH?:"required parameter JELLYFISH_PATH - path to jellyfish histo file"}

BWA_INDEX_PATH=${BWA_INDEX_PATH?:"required parameter BWA_INDEX_PATH - path to bwa reference index"}

MRF_INDEX_PATH=${MRF_INDEX_PATH?:"required parameter MRF_INDEX_PATH - path to mrf reference index"}
MRF_CHROMNAME_TRANSLATE_TBL=${MRF_CHROMNAME_TRANSLATE_TBL?:"required parameter MRF_CHROMNAME_TRANSLATE_TBL - path to table of chromosome names to translate from mrfast-ese (eg hg19 -> grc37)"}

SNP_VCF=${SNP_VCF?:"required parameter SNP_VCF - path to common variant VCF"}



COL_BASENAME=${COL_BASENAME:?"required parameter COL_BASENAME = basename for this chunk"}

COL_TARGET_BED=${COL_TARGET_BED:?"required parameter COL_TARGET_BED = bed file of target intervals"}

GEN_CANDS_DESIGN_OPTS_STAGE1=${GEN_CANDS_DESIGN_OPTS_STAGE1:-" --ligArmLengthRange 23,30 --extArmLengthRange 23,30  --gapFillRange 80,120 --padBy 100 --maxHomopolLen 8 --gcRange 5,90 --enforceTopStrandOnly"}

GEN_PAIRS_DESIGN_OPTS_STAGE1=${GEN_PAIRS_DESIGN_OPTS_STAGE1:-" --gapFillRange 80,130 --tmDiffRange ' -20,20' --enforceGapFillOnTopOnly "}

PICK_PAIRS_DESIGN_OPTS=${PICK_PAIRS_DESIGN_OPTS:-" --coverage 2 --scorePenaltyForOverlapExisting 0.50 "}


ONLY_PICK=${ONLY_PICK:-"no"}  # yes|no

CPU=${CPU:-24}
CPUlo=${CPUlo:-8}

DO_BED_ALLARMS=${DO_BED_ALLARMS:-"no"}  # yes|no

#######################################

baseOfKeyOut=`basename $RUN_KEY_OUT | sed -e 's/\..*//'`;

# intermediate keys
# KEY_INT_1=`dirname $RUN_KEY_IN`/$( echo `basename $RUN_KEY_IN` | sed -e 's/\.[^\.]*$//' ).bwaPE.$( echo `basename $RUN_KEY_IN` | sed -e 's/.*\.//' )

#######################################

echo "mkdir   ${BASE_DIR}/temp/"
mkdir -p ${BASE_DIR}/temp/
mkdir -p ${BASE_DIR}/runscripts/
mkdir -p ${BASE_DIR}/bed/stage1_arms_hf/
mkdir -p ${BASE_DIR}/bed/pairs1/
mkdir -p ${BASE_DIR}/bed/pairs_picked/

mkdir -p ${BASE_DIR}/temp/uniq_screen/

echo "key     $RUN_KEY_IN --> $RUN_KEY_OUT"

t.py \
  -i $RUN_KEY_IN \
  -o $RUN_KEY_OUT \
  --fillInCols \
  armcands_stage1:"lambda l:'%s/%s.armsStage1.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  armcands_stage1_hf:"lambda l:'%s/%s.armsStage1.hardFilt.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  bed_extarm_armcands_stage1_hf:"lambda l:'${BASE_DIR}/bed/stage1_arms_hf/%s.extArm.armsStage1.hardFilt.bed'%(l[os.environ['COL_BASENAME']])"\
  bed_ligarm_armcands_stage1_hf:"lambda l:'${BASE_DIR}/bed/stage1_arms_hf/%s.ligArm.armsStage1.hardFilt.bed'%(l[os.environ['COL_BASENAME']])"\
  uniqcheck1_fq:"lambda l:'${BASE_DIR}/temp/uniq_screen/%s.arms.fq'%(l[os.environ['COL_BASENAME']])"\
  uniqcheck1_bwabam:"lambda l:'${BASE_DIR}/temp/uniq_screen/%s.arms.bwa.bam'%(l[os.environ['COL_BASENAME']])"\
  uniqcheck1_mrfbase:"lambda l:'${BASE_DIR}/temp/uniq_screen/%s.arms.mrf'%(l[os.environ['COL_BASENAME']])"\
  uniqcheck1_jointbam:"lambda l:'${BASE_DIR}/temp/uniq_screen/%s.arms.joint.bam'%(l[os.environ['COL_BASENAME']])"\
  \
  armcands_filtuniq:"lambda l:'%s/%s.armsStage1.filtUniq.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  armcands_thinned:"lambda l:'%s/%s.armsStage1.thinned.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  \
  paircands_stage1:"lambda l:'%s/%s.pairs1.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  \
  bed_pairs1_plus:"lambda l:'${BASE_DIR}/bed/pairs1/%s.plus.bed'%(l[os.environ['COL_BASENAME']])"\
  bed_pairs1_notcoverable:"lambda l:'${BASE_DIR}/bed/pairs1/%s.plus.notCoverable.bed'%(l[os.environ['COL_BASENAME']])"\
  \
  pairs_picked:"lambda l:'%s/%s.pairs_picked.h5'%(os.environ['BASE_DIR'],l[os.environ['COL_BASENAME']] )"\
  bed_pairs_picked:"lambda l:'${BASE_DIR}/bed/pairs_picked/%s.plus.bed'%(l[os.environ['COL_BASENAME']])"\
  bed_pairs_notcovered:"lambda l:'${BASE_DIR}/bed/pairs_picked/%s.plus.notCoverable.bed'%(l[os.environ['COL_BASENAME']])"


  # uniqcheck1_joint:"lambda l:'${BASE_DIR}/temp/uniq_screen/%s.arms.mrf'%(l[os.environ['COL_BASENAME']])"\

#####################################################################
# generate candidates


if [[ "$ONLY_PICK" == "yes" ]] ; then
    echo "candidates are picked; skipping"
else

echo "\
generate_candidates.py  \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --targetBed \${${COL_TARGET_BED}} \
    ${GEN_CANDS_DESIGN_OPTS_STAGE1} \
    --outTable \${armcands_stage1} \
    --dsmipLeftIndent 3 \
    --dsmipRightIndent 3 \
    --generateDsMIPs "|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_gencands_stage1.mf ;

echo "generate candidates  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_gencands_stage1.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_gencands_stage1.mf -j ${CPUlo} all 

#####################################################################
# add Tm

echo "\
cands_add_tm.py  \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --inStore \${armcands_stage1} \
    --table extArm --table ligArm --table extArmMinus --table ligArmMinus"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_tmcands_stage1.mf ;

echo "Tm candidates  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_tmcands_stage1.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_tmcands_stage1.mf -j ${CPU} all 

#####################################################################
# add freq info

echo "\
cands_add_nmerfreq.py \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --inStore \${armcands_stage1} \
    --table extArm --table ligArm --table extArmMinus --table ligArmMinus \
    -k 15 --jellyfishFile ${JELLYFISH_PATH}"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_jfcands_stage1.mf ;

echo "jellyfish candidates  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_jfcands_stage1.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_jfcands_stage1.mf -j ${CPU} all 

#####################################################################
# hard filter unsuitable arms

echo "\
filter_arm_cands.py \
    --inStore \${armcands_stage1} \
    --outStore \${armcands_stage1_hf} \
    --dsMIPs \
    --fxnPass hard_filter_2_ds"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_hardfilt_stage1.mf ;

echo "hard filt candidates  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_hardfilt_stage1.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_hardfilt_stage1.mf -j ${CPU} all 


#####################################################################
# make bed tracks from all arms
if [[ "$DO_BED_ALLARMS" == "yes" ]] ; then

    echo "armtbl_to_bed.py \
        --inStore \${armcands_stage1_hf} \
        --table extArm \
        --bedOut \${bed_extarm_armcands_stage1_hf} \
        --fxnDesc 'lambda c,h:\"%s:Tm:%.2f,freq:%.2f\"%(c.name,c.arm_tm,c.arm_mean_kmer_freq)'; \
    armtbl_to_bed.py \
        --inStore \${armcands_stage1_hf} \
        --table ligArm \
        --bedOut \${bed_ligarm_armcands_stage1_hf} \
        --fxnDesc 'lambda c,h:\"%s:Tm:%.2f,freq:%.2f\"%(c.name,c.arm_tm,c.arm_mean_kmer_freq)' "|\
    pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
        > ${BASE_DIR}/runscripts/${baseOfKeyOut}_armtbltobed_stage1hf.mf ;

    echo "write bed for candidate arms  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_armtbltobed_stage1hf.mf"
    gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_armtbltobed_stage1hf.mf -j ${CPU} all 

fi

#####################################################################
# arm uniqueness check 

echo "cands_to_fastq_for_align.py \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --inStore \${armcands_stage1_hf} \
    --table extArm --table ligArm --table extArmMinus --table ligArmMinus \
    --format fastq -o \${uniqcheck1_fq}"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_stage1tofq.mf ;

echo "candidates to fq  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_stage1tofq.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_stage1tofq.mf -j ${CPU} all 

echo "bwa aln -t 1 -l 15 -n 0.01 ${BWA_INDEX_PATH} \${uniqcheck1_fq}   > \${uniqcheck1_fq}.aln"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwa.mf ;

echo "bwa aln  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwa.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwa.mf -j ${CPU} all 

echo "bwa samse -n 500000 ${BWA_INDEX_PATH} \${uniqcheck1_fq}.aln \${uniqcheck1_fq}  | perl xa2multi.pl | samtools view -b -  > \${uniqcheck1_bwabam}"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwasamse.mf ;

echo "bwa samse  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwasamse.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_bwasamse.mf -j ${CPU} all 

echo "mrfast_multilength_runmaker.py \
    --workDir $TMPDIR \
    --inFastq \${uniqcheck1_fq} \
    --mrfastPath mrfast \
    --mrfastOptions \" --search ${MRF_INDEX_PATH} -e 4 \" \
    --outBam \${uniqcheck1_mrfbase}.bam \
    --outRunScript \${uniqcheck1_mrfbase}.mrfrun_initial.sh \
    --outJoinScript \${uniqcheck1_mrfbase}.mrfrun_join.sh \
    --outCleanupScript \${uniqcheck1_mrfbase}.mrfrun_cleanup.sh"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfrunmaker.mf ;

echo "mrf runmaker  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfrunmaker.mf"
gmake -f ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfrunmaker.mf -j 1 all  

echo "cat \${uniqcheck1_mrfbase}.mrfrun_initial.sh; echo ' '" | pipeline.py --rkIn $RUN_KEY_OUT --loopCmds --quiet | sh > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfinitial.sh

cat ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfinitial.sh | pipeline.py --taskListToMakefile > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfinitial.mf ;
echo "mrf align  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfinitial.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfinitial.mf -j 32

echo "sh \${uniqcheck1_mrfbase}.mrfrun_join.sh"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfjoin.mf ;
echo "mrf join  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfjoin.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfjoin.mf -j ${CPU}

echo "sh \${uniqcheck1_mrfbase}.mrfrun_cleanup.sh"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfcleanup.mf ;
echo "mrf cleanup  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfcleanup.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfcleanup.mf -j ${CPU}

echo "translate_chrom_names.py \
    --inSam \${uniqcheck1_mrfbase}.bam \
    --outSam /dev/stdout \
    --chromMap ${MRF_CHROMNAME_TRANSLATE_TBL} \
    --outGenomeDict ${GENOME_DICT} |\
samtools calmd -b - ${GENOME_PATH} > \
    \${uniqcheck1_mrfbase}.fixMDandChangeRef.bam 2>/dev/null"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfCalmd.mf ;
echo "mrf translate names  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfCalmd.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mrfCalmd.mf -j ${CPU}

echo "java -Xmx8g -Djava.io.tmpdir=${TMPDIR} -jar $PICARD_DIR/picard.jar MergeSamFiles I=\${uniqcheck1_mrfbase}.fixMDandChangeRef.bam I=\${uniqcheck1_bwabam} O=\${uniqcheck1_jointbam} SO=queryname  VALIDATION_STRINGENCY=LENIENT MSD=true MAX_RECORDS_IN_RAM=5000000"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_mergeSamfiles.mf ;
echo "merge mrf+bwa  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mergeSamfiles.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_mergeSamfiles.mf -j ${CPU}

echo "cands_add_align_hits.py \
    --inStore \${armcands_stage1_hf} \
    --table extArm --table ligArm \
    --table extArmMinus --table ligArmMinus \
    --maxedsCloseFrac 0.101 --newtblSuffix mrfAndBwa \
    --inBam \${uniqcheck1_jointbam} \
    --disregardDoubleHitsAtSameLocation"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_addAlignHits.mf ;
echo "count align hits  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_addAlignHits.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_addAlignHits.mf -j ${CPU}


#####################################################################
# add population SNPs

echo "cands_add_pop_snps.py \
    --inStore \${armcands_stage1_hf} \
    --inVcfList <( echo ${SNP_VCF} ) \
    --table extArm --table ligArm \
    --table extArmMinus --table ligArmMinus \
    --newtblSuffix 1kg"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_addPopSnps.mf ;
echo "popsnps  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_addPopSnps.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_addPopSnps.mf -j ${CPU}

#####################################################################
# filter by mapping count

echo "filter_arm_cands.py \
    --inStore \${armcands_stage1_hf} \
    --outStore \${armcands_filtuniq} \
    --hitTblSuffix mrfAndBwa\
    --snpTblSuffix 1kg\
    --dsMIPs \
    --fxnPass hard_filter_uniq_ds_5exactclose"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_filterUniq.mf ;
echo "filter uniq  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_filterUniq.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_filterUniq.mf -j ${CPU}

#####################################################################
# thin arm sets

# TODO export pickStep, pickNum?

echo "score_and_thin_arms.py \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --inArmStore \${armcands_filtuniq} \
    --hitTblSufifx mrfAndBwa \
    --snpTblSuffix 1kg \
    --outArmStore \${armcands_thinned}  \
    --scoringClass ScoringModelLearned \
    --learningModel AdaBoost \
    --trainingMode preTrained \
    --trainSetPath ${TRAIN_PATH} \
    --setPrefix trainset \
    --pickStep 5 --pickNbest 2"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_thinArms.mf ;
echo "thin arms  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_thinArms.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_thinArms.mf -j ${CPU}

#####################################################################
# generate pairs

echo "generate_arm_pairs.py \
    --genome ${GENOME_PATH} \
    --genomeDict ${GENOME_DICT} \
    --inArmStore \${armcands_thinned} \
    --outPairStore \${paircands_stage1} \
    --hitTblSufifx mrfAndBwa\
    --snpTblSuffix 1kg \
    ${GEN_PAIRS_DESIGN_OPTS_STAGE1} \
    --generateDsMIPs \
    --learningModel AdaBoost \
    --trainingMode preTrained \
    --trainSetPath ${TRAIN_PATH} \
    --setPrefix trainset"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_genPairs1.mf ;
echo "generate pairs  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_genPairs1.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_genPairs1.mf -j ${CPU}

echo "pairtbl_to_bed.py  --inStore \${paircands_stage1} --armStore \${armcands_thinned} --table pairsPlus --bedOut \${bed_pairs1_plus}; \
    bedtools subtract -a \${${COL_TARGET_BED}} -b <( cat \${bed_pairs1_plus} | cut -f1,7,8) > \${bed_pairs1_notcoverable}"|\
pipeline.py --rkIn $RUN_KEY_OUT --tasksAsMakefile --loopCmds \
    > ${BASE_DIR}/runscripts/${baseOfKeyOut}_pairs1ToBed.mf ;
echo "pairs to bed  :  ${BASE_DIR}/runscripts/${baseOfKeyOut}_pairs1ToBed.mf"
gmake -f  ${BASE_DIR}/runscripts/${baseOfKeyOut}_pairs1ToBed.mf -j ${CPU}

  echo "sort : SKIP";
fi  #if [[ "$ONLY_PICK" == "yes" ]] ; then

#####################################################################
# join generated pairs.

python join_chunked_arm_and_pair_tbls.py \
    --key ${RUN_KEY_OUT} \
    --colnameArmTbl armcands_filtuniq \
    --colnamePairTbl paircands_stage1 \
    --outArmStore ${BASE_DIR}/arms_filtuniq.h5 \
    --outPairStore ${BASE_DIR}/pairs_all_joined.h5 \
    --hitTblSuffix mrfAndBwa \
    --snpTblSuffix 1kg
