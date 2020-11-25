from gwf import Workflow, AnonymousTarget
import os, re, sys
import numpy as np
os.environ['NUMEXPR_MAX_THREADS'] = '8' # to suppress a warning
import pandas as pd

from templates import *

gwf = Workflow(defaults={'account': 'clues'})

#################################################################################
# Paths to data files
#################################################################################

chromosomes = list(map(str, range(1, 23))) + ['X']

# build dict of vcf files
vcf_files = dict()
vcf_dir = '/project/simons/faststorage/data/1000Genomes/'
for chrom in chromosomes:
    vcf_file_name = os.path.join(vcf_dir, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom))
    if chrom == 'X':
        vcf_file_name = os.path.join(vcf_dir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')
    vcf_files[chrom] = vcf_file_name

arg_sample_times_file = 'data/tennessen_times_fine.txt'
arg_sample_popsize_file = 'data/tennessen_popsize_fine.txt'
decode_recomb_map_file = 'data/decode_hg38_lifted_to_hg19.bed'

# produced by separate workflow_freqs.py:
freq_data_file = 'steps/freq_data/derived_pop_freqs.h5'

################################################################################
# Dummy run of ARGweaver steps to get the log file neaded for transition probs.
################################################################################

window_start, window_end = 29500000, 30500000
pop = 'CEU'
chrom = '3'

# # make sites file
# dummy_sites_task = gwf.target_from_template(
#     name=f'dummy_sites_{window_start}_{window_end}',
#     template=arg_sites_file(
#         start=window_start,
#         end=window_end,
#         sites_file=f'steps/dummy/dummy_{window_start}_{window_end}.sites',
#         fasta_files=west_eur_fasta_files
#     )
# )

# extract 1000 genomes vcf 
dummy_vcf_task = gwf.target_from_template(
    name=f'dummy_vcfwindow_{chrom}_{window_start}_{window_end}_{pop}',
    template=vcf_window(
        vcf_file=vcf_files[chrom],
        chrom=chrom, win_start=window_start, win_end=window_end, pop=pop,
        stepsdir='steps/recode_vcf/dummy'
    )
)

# make sites file for argweaver
dummy_sites_task = gwf.target_from_template(
    name=f'dummy_vcf2sites_{chrom}_{window_start}_{window_end}_{pop}',
    template=vcf2sites(
        vcffile=dummy_vcf_task.outputs['vcf_window_file'],
        chrom=chrom, win_start=window_start, win_end=window_end, pop=pop,
        stepsdir='steps/sites_files/dummy'
    )
)

# make a recombination file for argweaver
dummy_recomb_task = gwf.target_from_template(
    name=f'dummy_rec_{window_start}_{window_end}',
    template=arg_recomb_file(
        recomb_map=decode_recomb_map_file,
        chrom=chrom,
        start=window_start,
        end=window_end,
        stepsdir='steps/recomb_files/dummy'
    )
)

# run argsample && smc2bed-all
dummy_argsample_target = gwf.target_from_template(
    name=f'dummy_args_{window_start}_{window_end}',
    template=argsample(
        sites_file=dummy_sites_task.outputs['sites_file'], 
        times_file=arg_sample_times_file, 
        popsize_file=arg_sample_popsize_file, 
        recomb_file=dummy_recomb_task.outputs['recomb_file'], 
        stepsdir='steps/argsample/dummy'
    )
)

################################################################################
# Precompute transistion probabilities for CLUES
################################################################################

# make transition matrices for clues for a range of selection coeficients:
selection_coeficients = np.logspace(-5, np.log10(0.5), 50)
trans_task_list = list()
for i, sel_coef in enumerate(selection_coeficients):
    task = gwf.target_from_template(
        name=f'trans_mat_{i}', 
        template=transition_matrices(sel_coef, dummy_argsample_target.outputs['log_file']))
    trans_task_list.append(task)

# all transition matrix files:
trans_mat_files = [output for task in trans_task_list for output in task.outputs]

# make transition matrices conditional on snp sampling frequencies:
freqs = '0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98'
cond_trans_matrix_file = 'steps/trans/trans_tennesen_fine.hdf5'
cond_trans_matrix_file_no_suffix = modpath(cond_trans_matrix_file, suffix='')

gwf.target(name='cond_matrices',
     inputs=trans_mat_files, 
     outputs=[cond_trans_matrix_file], 
     walltime='4-00:00:00', 
     memory='36g') << f"""
mkdir -p {os.path.dirname(cond_trans_matrix_file)}
python ./clues/conditional_transition_matrices.py {dummy_argsample_target.outputs['log_file']} steps/trans/matrices/ \
    --listFreqs {freqs} -o {cond_trans_matrix_file_no_suffix} --debug
"""

# ################################################################################
# # Run ARGweaver and CLUES steps for each SNP
# ################################################################################


window_centers = [55000000, 110000000]
arg_win_size = 1500000
center_analyzed = 500000
flank = int((arg_win_size - center_analyzed) / 2)

clues_windows = [(int(pos - arg_win_size/2), int(pos + arg_win_size/2)) for pos in window_centers]

for window_start, window_end in clues_windows:

    # get recombination file
    recomb_task = gwf.target_from_template(
        name=f'rec_{chrom}_{window_start}_{window_end}_{pop}',
        template=arg_recomb_file(
            recomb_map=decode_recomb_map_file,
            chrom=chrom,
            start=window_start,
            end=window_end,
            stepsdir=f'steps/recomb_file'
        )
    )

    for pop in ['CEU', 'FIN']:


        # # make sites file
        # sites_task = gwf.target_from_template(
        #     name=f'sites_{window_start}_{window_end}',
        #     template=arg_sites_file(
        #         start=window_start,
        #         end=window_end,
        #         sites_file=f'steps/argsample/{window_start}_{window_end}/{window_start}_{window_end}.sites',
        #         fasta_files=west_eur_fasta_files
        #     )
        # )


        # extract vcf 
        vcf_task = gwf.target_from_template(
                    name=f'vcfwindow_{chrom}_{window_start}_{window_end}_{pop}',
                    template=vcf_window(
                        vcf_file=vcf_files[chrom],
                        chrom=chrom,
                        win_start=window_start,
                        win_end=window_end,
                        pop=pop,
                        stepsdir='steps/recode_vcf'
                    )
                )
    

        # make sites file
        sites_task = gwf.target_from_template(
            name=f'vcf2sites_{chrom}_{window_start}_{window_end}_{pop}',
            template=vcf2sites(
                vcffile=vcf_task.outputs['vcf_window_file'],
                chrom=chrom,
                win_start=window_start,
                win_end=window_end,
                pop=pop,
                stepsdir='steps/sites_files'

            )
        )

        min_freq = 0.25
        nr_snps = 5

        snps_start, snps_end = window_start + flank, window_end - flank
        snp_list = get_snps(freq_data_file, chrom, pop, snps_start, snps_end, min_freq, nr_snps)

        for chain in range(1, 3):

            argsample_target = gwf.target_from_template(
                name=f'args_{chrom}_{window_start}_{window_end}_{pop}_{chain}',
                template=argsample(
                    sites_file=sites_task.outputs['sites_file'], 
                    times_file=arg_sample_times_file, 
                    popsize_file=arg_sample_popsize_file, 
                    recomb_file=recomb_task.outputs['recomb_file'], 
                    stepsdir=f'steps/argsample/{chrom}_{window_start}_{window_end}_{pop}_{chain}'
                )
            )

            clues_task_list = list()
            stepsdir = f'steps/clues/{chrom}_{window_start}_{window_end}_{pop}_{chain}'
            for chrom, snp_pos, derived_allele, derived_freq in snp_list:
                clues_task = gwf.target_from_template(
                    name=f'clues_{chrom}_{window_start}_{window_end}_{pop}_{chain}_{snp_pos}',
                    template=clues(
                        bed_file=argsample_target.outputs['bed_file'], 
                        sites_file=sites_task.outputs['sites_file'], 
                        cond_trans_matrix_file=cond_trans_matrix_file, 
                        snp_pos=snp_pos, chrom=chrom, win_start=window_start, win_end=window_end, 
                        derived_allele=derived_allele, derived_freq=derived_freq,
                        chain=chain,
                        stepsdir=stepsdir
                    )
                )
                clues_task_list.append(clues_task)

            clues_files = [output for task in clues_task_list for output in task.outputs]


            # Extract info from all clues files for each window
            clues_csv_file_name = f'steps/csv/clues_{chrom}_{window_start}_{window_end}_{pop}_{chain}.csv'
            # steps_dir = os.path.dirname(clues_csv_file_name)

            clues_file_base_names = ' '.join([modpath(f, parent='') for f in clues_files])

            gwf.target(f'clues_{chrom}_{window_start}_{window_end}_{pop}_{chain}_csv', 
                inputs=clues_files, outputs=[clues_csv_file_name], walltime='1:00:00', memory='1g') << f"""

            mkdir -p {os.path.dirname(clues_csv_file_name)}
            python scripts/extract_clues_info.py {clues_csv_file_name} {stepsdir} {clues_file_base_names}
            """
