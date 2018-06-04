from setuptools import setup
#from distutils.core import setup
from distutils.extension import Extension
# from Cython.Build import cythonize
# from Cython.Distutils import build_ext

import pysam
import numpy 
import glob
import os.path as op

setup(
	name="makeamip",

    packages=['makeamip'],

    include_dirs = [numpy.get_include()]+pysam.get_include(),

    package_data = { 'makeamip': ['data/*' ] },

    scripts = [
        "makeamip/pipeline_script/pipeline_template.sh",
        "makeamip/annotate_miptbl_with_target.py",
        "makeamip/armtbl_to_bed.py",
        "makeamip/cands_add_align_hits.py",
        "makeamip/cands_add_nmerfreq.py",
        "makeamip/cands_add_pop_snps.py",
        "makeamip/cands_add_tm.py",
        "makeamip/cands_to_fastq_for_align.py",
        "makeamip/evalulate_design.py",
        "makeamip/filter_arm_cands.py",
        "makeamip/generate_arm_pairs.py",
        "makeamip/generate_candidates.py",
        "makeamip/gene_list_to_bed_hg19.py",
        "makeamip/greedy_pick_mip_set_v2.py",
        "makeamip/join_chunked_arm_and_pair_tbls.py",
        "makeamip/make_ado_insilico_ref.py",
        "makeamip/merged_tables_fix_flags.py",
        "makeamip/mkcands_from_yoon_tbl.py",
        "makeamip/mrfast_multilength_runmaker.py",
        "makeamip/pair_flanking_base_tbl.py",
        "makeamip/pair_metrics.py",
        "makeamip/pair_table_out.py",
        "makeamip/pairtbl_to_bed.py",
        "makeamip/pairtbl_to_oligo_list.py",
        "makeamip/pairtbl_to_partialread_mappable_mask.py",
        "makeamip/pair_to_mipgen_table.py",
        "makeamip/picked_panel_check_coverage.py",
        "makeamip/pipeline.py",
        "makeamip/regenerate_existing_pairs_with_new_adaptors.py",
        "makeamip/rescore_pairs_by_observed_eff.py",
        "makeamip/rescore_pairs.py",
        "makeamip/score_and_thin_arms.py",
        "makeamip/split_monolithic_armspairs_across_chunks.py",
        "makeamip/split_target_bed.py",
        "makeamip/suppress_bad_mrfast_alignments.py",
        "makeamip/suppress_duplicate_probes_in_joined_table.py",
        "makeamip/t.py",
        "makeamip/transfer_picked_info.py",
        "makeamip/translate_chrom_names.py",
        "makeamip/rescore_pairs.py",

    ],

    # entry_points = {
    # 	'console_scripts': [ 'rev_read_strip_mid = makeamip.rev_read_strip_mid:main',
    # 		 	 	         'merged_pe_strip_mid = makeamip.merged_pe_strip_mid:main',
    # 		 	 	         'flag_dup_mid_reads = makeamip.flag_dup_mid_reads:main',
    # 		 	 	         'exact_trim_mip_arms = makeamip.exact_trim_mip_arms_v2:main',
    # 		 	 	         'annotate_bam_by_mip_SE = makeamip.annotate_bam_by_mip_SE:main',
    # 		 	 	         'annotate_bam_by_mip = makeamip.annotate_bam_by_mip:main',
    # 		 	 	         'gather_target_coverage_at_thresholds = makeamip.uniformity.gather_target_coverage_at_thresholds:main',
    #                          'gather_perbp_uniformity = makeamip.uniformity.gather_perbp_uniformity:main',
    # 		 	 	         'probe_hit_count_from_bam = makeamip.uniformity.probe_hit_count_from_bam:main',
    # 		 	 	         'target_coverage_uniformity = makeamip.uniformity.target_coverage_uniformity:main',
    #                          'forcore_mip_sample_sheet = makeamip.forcore_mip_sample_sheet:main',
    #                          'join_mip_libsheet_core_iemsheet = makeamip.join_mip_libsheet_core_iemsheet:main'
    # 		 	 	          ]
    # }

)





