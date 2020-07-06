
# Running CLUES (ARGweaver version)

This is a workflow for running CLUES on the computationally phased 1000 genomes data. There are two workflow files: `workflow_1000g_derived_info.py` and `workflow_clues.py`. The `workflow_1000g_derived_info.py` workflow needs to be run first to generate a HDF5 file that is used by the main workflow `workflow_clues.py` to look up information about derived variants. Both workflows hardcoded paths to data files in `/project/simons/faststorage/data/1000Genomes`. So membership of the `simons` project folder is required.

It is assumed that clues-v0 is installed (unpacked) in the base directory of the repository.

You need to make a conda envirommnet, E.g. `clues-v0`, with at least these packages:

    conda create -n <some name> -c anaconda -c conda-forge -c genomedk -c gwforg -c kaspermunch -c bioconda vcftools gwf pandas numpy h5py biopython future argweaver scikit-allel chromwindow

Activate environemnt and run like this:

    gwf -f workflow_1000g_derived_info.py run

Once it finishes you run:

    gwf -f workflow_clues.py

This workflow produces a table file with LR scores for all analysed SNPs.

The notebooks folder contains a notebook with hints to the interpretation of raw CLUES output.
