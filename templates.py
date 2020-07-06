
import sys, os, re
from gwf import Workflow, AnonymousTarget
from subprocess import PIPE, Popen

def modpath(p, parent=None, base=None, suffix=None):
    """
    Modifies dir, basename, or suffix of path.
    """

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


# # def arg_sites_file(start, end, sites_file, fasta_files):
# def arg_sites_file(start, end, sites_file, fasta_files):

#     inputs = fasta_files

#     outputs = {'sites_file': sites_file}
#     options = {
#         'memory': '4g',
#         'walltime': '01:00:00'
#     }

#     spec = f'''
#     mkdir -p {os.path.dirname(sites_file)}
#     python scripts/argsample_sites_file.py X {start} {end} {sites_file} {" ".join(fasta_files)}
#     '''

#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def vcf_window(vcf_file, chrom, win_start, win_end, pop, stepsdir):
    """
    Extracts recoded VCF of a given window.
    """
    output_vcf_file = f'{stepsdir}/{chrom}_{win_start}_{win_end}_{pop}.recode.vcf'
    vcf_base_name = modpath(output_vcf_file, suffix=('.recode.vcf', ''))
    
    inputs = [vcf_file]
    outputs = {'vcf_window_file': output_vcf_file}
    options = {
        'memory': '8g',
        'walltime': '05:00:00'
    }
    
    spec = f'''
    
    mkdir -p steps/recode_vcf/{chrom}_{win_start}_{win_end}

    vcftools --gzvcf {vcf_file} --chr {chrom} --from-bp {win_start} --to-bp {win_end-1} \
        --keep ~/simons/faststorage/data/1000Genomes/metainfo/{pop}_male.txt \
        --keep ~/simons/faststorage/data/1000Genomes/metainfo/{pop}_female.txt \
        --remove-indels --remove-filtered-all --max-alleles 2 --recode --non-ref-ac-any 1 \
        --out {vcf_base_name}
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def vcf2sites(vcffile, chrom, win_start, win_end, pop, stepsdir):
    """
    Convert a VCF files to the sites format used by ARGweaver.
    """
    sitesfile=f'{stepsdir}/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.sites'

    inputs = [vcffile]
    outputs = {'sites_file': sitesfile}
    options = {
        'memory': '8g',
        'walltime': '01:00:00'
    }
    
    spec = f'''
    
    mkdir -p steps/sitesfiles/{chrom}_{win_start}_{win_end}

    python scripts/vcf2sites.py {vcffile} {sitesfile}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def arg_recomb_file(recomb_map, chrom, start, end, stepsdir):
    """
    Produce a recombination map for a window in bins with means for each 10kb.
    """
    recomb_file = f'{stepsdir}/{chrom}_{start}_{end}.rec'

    inputs = [recomb_map]
    outputs = {'recomb_file': recomb_file}
    options = {
        'memory': '4g',
        'walltime': '01:00:00'
    }

    spec = f'''
    mkdir -p {stepsdir}
    python scripts/argsample_rec_window.py {recomb_map} {chrom} {start} {end} {recomb_file}
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def argsample(sites_file, times_file, popsize_file, recomb_file, stepsdir):
    """
    Run ARGweaver sampling and concatetenate resulting files into one bed file.
    """
    bed_file = modpath(sites_file, parent=stepsdir, suffix='.bed.gz')

    arg_sample_base_name = modpath(bed_file, suffix='')
    log_file = modpath(arg_sample_base_name, suffix='.log')

    inputs = {'sites_file': sites_file, 'recomb_file': recomb_file}
    outputs = {'bed_file': bed_file, 'log_file': log_file}
    options = {
        'memory': '40g',
        'walltime': '14-00:00:00'
    }

    spec = f'''
    mkdir -p {stepsdir}
    arg-sample -s {sites_file} \
            --times-file {times_file} \
            --popsize-file {popsize_file} \
            --recombmap {recomb_file} \
            -m 1.247e-08 \
            -c 25 \
            -n 30000 \
            --overwrite \
            -o {arg_sample_base_name} \
    && \
    ./argweaver/bin/smc2bed-all {arg_sample_base_name}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def transition_matrices(selection_coef, arg_weaver_log_file):

    output_file = f'trans.s_{selection_coef}.h5'

    path = 'steps/trans/matrices'

    inputs=[arg_weaver_log_file]
    outputs=[os.path.join(path, output_file)]
    options = {'memory': '16g', 'walltime': '5-00:00:00'}

    spec = f'''
    mkdir -p {path}
    ORIGDIR=`pwd`
    cd {path}

    python $ORIGDIR/clues/make_transition_matrices_from_argweaver.py 10000 {selection_coef} \
        $ORIGDIR/{arg_weaver_log_file} {output_file} --breaks 0.95 0.025 --debug    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def clues(bed_file, sites_file, cond_trans_matrix_file, snp_pos, chrom, win_start, win_end, derived_allele, derived_freq, stepsdir):

    clues_output_file = modpath(bed_file, suffix=('.bed.gz', f'_{snp_pos}.h5'), parent=stepsdir)
    log_file = modpath(bed_file, suffix=('.bed.gz', '.log'))
    trees_file = modpath(bed_file, suffix=('.bed.gz', '.trees'))
    clues_output_base_name = modpath(clues_output_file, suffix=('.h5', ''))
    
    inputs = [bed_file]
    outputs = [clues_output_file]
    options = {
        'memory': '8g',
        'walltime': '1-00:00:00'
    }

    spec = f'''        
    mkdir -p {stepsdir}
    arg-summarize -a {bed_file} -r {chrom}:{snp_pos}-{snp_pos} -l {log_file} -E > {trees_file} \
    && \
    python ./clues/clues.py {trees_file} {cond_trans_matrix_file} {sites_file} {derived_freq} --posn {snp_pos} \
        --derivedAllele {derived_allele} --noAncientHap --approx 10000 --thin 10 --burnin 100 --output {clues_output_base_name}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def execute(cmd, stdin=None):
    process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(stdin)
    if not process.returncode == 0:
        print(cmd)
        print(stderr.decode())
    return stdout, stderr


def read_snp_info(snp_file):
    snp_list = list()
    with open('snps.txt', 'r') as snp_file:
        for line in snp_file:
            chrom, snp_pos, derived_allele, derived_freq = line.split()
            snp_pos = int(snp_pos)
            derived_freq = float(derived_freq)
            snp_list.append((chrom, snp_pos, derived_allele, derived_freq))
    return snp_list


def get_single_snp(freq_data_file, chrom, pop, snp_pos):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --snppos {snp_pos}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


def get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --start {window_start} --end {window_end} --minfreq {min_freq} --maxsnps {nr_snps}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


def clues2csv(clues_files, chrom, win_start, win_end, pop):
    clues_csv_file_name = f'steps/clues/csv/clues_{chrom}_{win_start}_{win_end}_{pop}.csv'
    clues_file_base_names = ' '.join([modpath(f, parent='', suffix='') for f in clues_files])

    inputs = clues_files
    outputs = [clues_csv_file_name]
    options = {
        'memory': '1g',
        'walltime': '1:00:00'
    }

    spec = f'''

    mkdir -p steps/clues/csv

    python scripts/extract_clues_info.py {clues_csv_file_name} steps/clues/clues_{chrom}_{win_start}_{win_end}_{pop} {clues_file_base_names}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
