



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
        vcf_file_name = os.path.join(data_dir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')

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
