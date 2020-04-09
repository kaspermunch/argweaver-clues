from gwf import Workflow, AnonymousTarget
import os, re, sys
import numpy as np

gwf = Workflow()

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


################################################################################################
# Build data structure for derived variants that we can use for lookup:
################################################################################################

def extract_pop_frequencies(vcf_file_name, freq_file_name, males, females):

    inputs = [vcf_file_name]
    outputs = [freq_file_name + '.frq']
    options = {
        'cores': 1,
        'memory': '16g'
    }
    spec = '''
    mkdir -p steps/freq_data

    vcftools --gzvcf {} \
        --remove-indels --remove-filtered-all --max-alleles 2 \
        --non-ref-ac-any 1 \
        --keep {} --keep {} --freq --out {}
    '''.format(vcf_file_name, males, females, freq_file_name)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


data_dir = '/project/simons/faststorage/data/1000Genomes/'

# read population identifiers
pop_info_file_name = os.path.join(data_dir, 'metainfo/pop_names.tsv')
populations = []
with open(pop_info_file_name) as f:
    for line in f:
        pop, _ = line.split(maxsplit=1)
        populations.append(pop)

freq_task_list = list()
chromosomes = list(map(str, range(1, 23))) + ['X']
for chrom in chromosomes:

    # name of vcf file
    vcf_file_name = os.path.join(data_dir, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom))

    if chrom == 'X':
        vcf_file_name = os.path.join(data_dir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'.format(chrom))

    for pop in populations:

        # names of files for males and females from population
        males_file_name = os.path.join(data_dir, 'metainfo/{}_male.txt'.format(pop))
        females_file_name = os.path.join(data_dir, 'metainfo/{}_female.txt'.format(pop))

        freq_file_name = 'steps/freq_data/{}_{}'.format(pop, chrom)

        target = gwf.target_from_template("freqs_{}_{}".format(pop, chrom), extract_pop_frequencies(vcf_file_name, 
            freq_file_name, males_file_name, females_file_name))

        freq_task_list.append(target)

freq_files = [output for task in freq_task_list for output in task.outputs]

# Build data set:
hdf_file_name = 'steps/freq_data/derived_pop_freqs.h5'

gwf.target('compile_freq_data', 
    inputs=freq_files, 
    outputs=[hdf_file_name], 
    walltime='10:00:00', memory='36g') << f"""

python scripts/build_derived_freq_data.py {hdf_file_name}
"""

################################################################################################
# Example run of ARGweaver to check that it works:
# We also need the resulting log file to build transition matrices.
################################################################################################

chrom = '2'
window_start = 136000000
window_end = 137000000
sites_file = 'steps/sites/chr2.sites' # I just use the one Janne made

snp_pos = 136608646
argweaver_bed_file = 'steps/arg_sample/arg_sample.bed.gz'
argweaver_log_file = 'steps/arg_sample/arg_sample.log'
argweaver_trees_file = 'steps/arg_sample/arg_sample.trees'
argweaver_base_name = modpath(argweaver_log_file, suffix='')

gwf.target('argweaver', 
    inputs=[], 
    outputs=[argweaver_bed_file, argweaver_log_file, argweaver_trees_file],
    walltime='10:00:00', memory='10g') << f"""

mkdir -p steps/arg_sample

/home/kmt/anaconda3/envs/clues/bin/arg-sample -s {sites_file} --times-file data/tennessen_times_fine.txt \
    --popsize-file data/tennessen_popsize_fine.txt -r 1e-8 -m 1.2e-8 -c 25 -n 3000 --overwrite -o {argweaver_base_name}

../../software/argweaver/bin/smc2bed-all {argweaver_base_name}

~/anaconda3/envs/clues/bin/arg-summarize -a {argweaver_bed_file} -r {chrom}:{snp_pos}-{snp_pos} \
    -l {argweaver_log_file} -E > {argweaver_trees_file}
"""

#############################################################################################
# Build transition matrices using tennensen_fine demography:
################################################################################################

def make_transition_matrices_from_argweaver(selection_coef, arg_weaver_log_file):

    output_file = f'trans.s_{selection_coef}.h5'

    path = 'steps/trans/matrices'

    inputs=[arg_weaver_log_file]
    outputs=[os.path.join(path, output_file)]
    options = {'memory': '16g', 'walltime': '5-00:00:00'}

    spec = f'''
    mkdir -p {path}
    ORIGDIR=`pwd`
    cd {path}

    python $ORIGDIR/../../software/clues/make_transition_matrices_from_argweaver.py 10000 {selection_coef} \
        $ORIGDIR/{arg_weaver_log_file} {output_file} --breaks 0.95 0.025 --debug
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


trans_task_list = list()
for i, sel_coef in enumerate(np.logspace(-5, np.log10(0.5), 50)):
    task = gwf.target_from_template(f'trans_mat_{i}', make_transition_matrices_from_argweaver(sel_coef, argweaver_log_file))
    trans_task_list.append(task)

trans_mat_files = [output for task in trans_task_list for output in task.outputs]

#freqs = ' '.join([str(round(x, 2)) for x in np.linspace(0.05, 0.95, 19)])
freqs = '0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95'
cond_trans_matrix_file = 'steps/trans/trans_tennesen_fine.hdf5'
cond_trans_matrix_file_no_suffix = modpath(cond_trans_matrix_file, suffix='')
gwf.target('cond_matrices', inputs=trans_mat_files, outputs=[cond_trans_matrix_file], walltime='24:00:00', memory='36g') << f"""

python ../../software/clues/conditional_transition_matrices.py {argweaver_log_file} steps/trans/matrices/ \
    --listFreqs {freqs} -o {cond_trans_matrix_file_no_suffix}
"""

################################################################################################
# Run clues on one SNP:
################################################################################################

clues_output_file = 'steps/clues/chr2_136608646.h5'
clues_output_base_name = modpath(clues_output_file, suffix='')
derived_allele = 'G'
derived_freq = 75e-2
snp_pos = 136608646

gwf.target('clues', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""
mkdir -p steps/clues

python ../../software/clues/clues.py {argweaver_trees_file} {cond_trans_matrix_file} {sites_file} {derived_freq} \
    --posn {snp_pos} --derivedAllele {derived_allele} --noAncientHap \
    --thin 10 --burnin 100 --output {clues_output_base_name}
"""


# ################################################################################################
# # Run clues on many SNPs in a window:
# ################################################################################################

# freq_data_file = 'steps/freq_data/derived_pop_freqs.h5'
# chrom = '2'
# window_start = 136000000
# window_end = 137000000
# pop = 'CEU'
# min_freq = 0.25
# nr_snps = 5

# from subprocess import PIPE, Popen

# def execute(cmd, stdin=None):
#     process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
#     stdout, stderr = process.communicate(stdin)
#     return stdout, stderr

# def get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps):

# #    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} snps.txt --snppos {snp_pos}")
#     execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} snps.txt --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps}")
#     print(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} snps.txt --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps}")

#     snp_list = list()
#     with open('snps.txt', 'r') as snp_file:
#         for line in snp_file:
#             chrom, snp_pos, derived_allele, derived_freq = line.split()
#             snp_pos = int(snp_pos)
#             derived_freq = float(derived_freq)
#             snp_list.append((chrom, snp_pos, derived_allele, derived_freq))
#     return snp_list

# snp_list = get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps)

# for chrom, snp_pos, derived_allele, derived_freq in snp_list:

#     print(chrom, snp_pos, derived_allele, derived_freq)

#     clues_output_file = f'steps/clues/clues_{chrom}_{snp_pos}.h5'
#     clues_output_base_name = modpath(clues_output_file, suffix='')

#     gwf.target(f'clues_{chrom}_{snp_pos}', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""
#     mkdir -p steps/clues

#     python ../../software/clues/clues.py {argweaver_trees_file} {cond_trans_matrix_file} {sites_file} {derived_freq} \
#         --posn {snp_pos} --derivedAllele {derived_allele} --noAncientHap \
#         --thin 10 --burnin 100 --output {clues_output_base_name}
#     """













####################################################################


# clues_output_file = 'steps/clues/clues'
# clues_output_base_name = modpath(clues_output_file, suffix='')

# gwf.target('clues', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""

# mkdir -p steps/clues

# JOBDIR=/scratch/$GWF_JOBID

# JOBDIR=.

# python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} $JOBDIR/snps.txt \
#     --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps} 


# python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} $JOBDIR/snps.txt \
#     --snppos {snp_pos} 

# while IFS= read -r line
# do
#     POS=$(cut -f1 -d " " <<<$line)
#     BASE=$(cut -f2 -d " " <<<$line)
#     FREQ=$(cut -f3 -d " " <<<$line)

#     ~/anaconda3/envs/clues/bin/arg-summarize -a {argweaver_bed_file} -r {chrom}:$POS-$POS \
#         -l {argweaver_log_file} -E > $JOBDIR/$POS.trees

#     python ../../software/clues/clues.py $JOBDIR/$POS.trees {cond_trans_matrix_file} {sites_file} $FREQ \
#         --posn $POS --derivedAllele $BASE --noAncientHap \
#         --thin 10 --burnin 100 --output {clues_output_base_name}_{}_$POS

# done < "$JOBDIR/snps.txt"
# """ 
