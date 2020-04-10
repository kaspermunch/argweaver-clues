


cat steps/sim1.trees | perl -pi -e 's/n(\d+):/$1:/g' > ../kmt/sim1_fix.trees

python ../../software/clues/clues.py ../kmt/sim1_fix.trees steps/example.f_75.hdf5 steps/sim1.sites 75e-2 --thin 10 --burnin 100 -o ../kmt/sim1.clues


# Test

Here are the commands we run:

Sites file extracted 

She uses a sites file `chr2.sites` with CEU variation around the lactase gene, extracted from a 1000 genomes VCF file:

ARGweaver sampling:

    /home/kmt/anaconda3/envs/clues/bin/arg-sample -s chr2.sites --times-file tennessen_times_fine.txt --popsize-file tennessen_popsize_fine.txt -r 1e-8 -m 1.2e-8 -c 25 -n 3000 --overwrite -o arg_sample

    ../../../software/argweaver/bin/smc2bed-all arg_sample

Extracting SNPs for lactase SNP:

    ~/anaconda3/envs/clues/bin/arg-summarize -a arg_sample.bed.gz -r 2:136608646-136608646 -l arg_sample.log -E > arg_sample.trees

Add ancestral state to sites file:

    python ../scripts/add_ancestral_state.py chr2.sites ~/simons/faststorage/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/all_chromosomes.fa chr2_with_ancestral.sites

> When selecting SNPs to analyse, I should read sites file, select positions from that

> When running clues, you need to pick the right example.f_75.hdf5 and to round the derived allele frequency to the corresponding value.

Run clues:


    python ../../../software/clues/clues.py arg_sample.trees example.f_75.hdf5 chr2.sites 75e-2 -ancientHap ANCESTRAL --thin 10 --burnin 100 --output clues.hdf

    python ../../../software/clues/clues.py arg_sample.trees example.f_75.hdf5 chr2.sites 75e-2 -ancientHap ANCESTRAL --sitesFile chr2_with_ancestral.sites --thin 10 --burnin 100 --output clues.hdf

    python ../../../software/clues/clues.py arg_sample.trees example.f_75.hdf5 chr2.sites 75e-2 --noAncientHap --derivedAllele G --thin 10 --burnin 100 --output clues.hdf

    python ../../../software/clues/clues.py arg_sample.trees example.f_75.hdf5 chr2.sites 75e-2 --derivedAllele G --thin 10 --burnin 100 --output clues.hdf


This works:

    python ../../../software/clues/clues.py arg_sample.trees example.f_75.hdf5 chr2.sites 75e-2 --posn 136608646 --noAncientHap --derivedAllele G --thin 10 --burnin 100 --output chr2_136608646


conda install pyh5
conda install -c anaconda future
