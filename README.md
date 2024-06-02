# Analyzing the fitness landscape of COV107-23, a SARS-CoV-2 spike antibody

Code to analyze deep sequencing files of COV107-23 combinatorial mutations, and ddG simulations of COV107-23 mutations

## Analysis of deep sequencing data following fluorescence-activated cell sorting of COV107-23 combinatorial mutants

### Dependencies

* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)

### Input files

Due to large size of fastq files, they cannot be uploaded to Github. Create a /fastq/ folder. Please download fastq files from NCBI, and transfer to fastq folder.
* [./Fasta/COV107_germline_ref.fa](./Fasta/COV107_germline_ref.fa): Amino acid sequence of COV107-23 germline
* [./fastq/Sample1_ATCACGAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801776): Input expression forward reads
* [./fastq/Sample1_ATCACGAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801777): Input expression reverse reads
* [./fastq/Sample2_CGATGTAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801778): Sorted expression forward reads (Replicate 1)
* [./fastq/Sample2_CGATGTAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801779): Sorted expression reverse reads (Replicate 1)
* [./fastq/Sample3_TTAGGCAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801780): Sorted expression forward reads (Replicate 2)
* [./fastq/Sample3_TTAGGCAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801781): Sorted expression reverse reads (Replicate 2)

### I. Calculate counts and fitness from fastq files

1. Calculate read counts and fitness from fastq files. Create a /fastq/ folder to store all downloaded fastq files.
``python scripts/COV107SHM_fq2fit.py``

- Input files
    - [./Fasta/COV107_germline_ref.fa](./Fasta/COV107_germline_ref.fa): Amino acid sequence of COV107-23 germline
    - [./fastq/Sample1_ATCACGAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801776): Input expression forward reads
    - [./fastq/Sample1_ATCACGAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801777): Input expression reverse reads
    - [./fastq/Sample2_CGATGTAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801778): Sorted expression forward reads (Replicate 1)
    - [./fastq/Sample2_CGATGTAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801779): Sorted expression reverse reads (Replicate 1)
    - [./fastq/Sample3_TTAGGCAT_L001_R1_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801780): Sorted expression forward reads (Replicate 2)
    - [./fastq/Sample3_TTAGGCAT_L001_R2_001.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRS9801781): Sorted expression reverse reads (Replicate 2)
- Output file
    - [./results/COV107_mutlib_fit.tsv](./results/COV107_mutlib_fit.tsv)

2. Filter results
``python scripts/COV107SHM_filter_result.py``

- Input file 
    - [./results/COV107_mutlib_fit.tsv](./results/COV107_mutlib_fit.tsv)
- Output file
    - [./results/COV107_mutlib_fit_exp.tsv](./results/COV107_mutlib_fit_exp.tsv)

3. Plot correlation of expression fitness between two independent experimental replicates
``Rscript scripts/COV107_filtered.R``

- Input file
    - [./results/COV107_mutlib_fit_exp.tsv](./results/COV107_mutlib_fit_exp.tsv)
- Output file
    - [./graph/cor_exp.pdf](./graph/cor_exp.pdf)


## Analysis of ddG simulations of COV107-23 mutations

### Dependencies

* [Python](https://www.python.org/) (version 3.9)
* [Rosetta Software Suite](https://www.rosettacommons.org/software/license-and-download)
    * Either an academic or commercial license is required. One can request a license in the link above.
* [PyMOL](https://pymol.org/2/)
* [pdb-tools](http://www.bonvinlab.org/pdb-tools/)

### Input files

* [Crystal structure of COV107-23](./structure/7lka.pdb)
* [Mutation files](./mut_files/)

### I. Renumber PDB file to prepare for mutagenesis

1. Using PyMOL and the pdb file (PDB: 7LKA), remove the solvent by running, in the PyMOL terminal
``remove solvent``

2. After removing the solvent, only select chain A (antibody heavy chain) and chain B (antibody light chain) by running, in the PyMOL terminal
``remove chain C+D+E+F+H+L``

3. Export the new molecule and retain atom IDs as [COV107.pdb](./structure/COV107.pdb).

4. Using pdb_reres.py in pdb-tools and the [COV107.pdb](./structure/COV107.pdb) file, run in the terminal
``python pdb_reres.py COV107.pdb > COV107_renum.pdb``

[COV107_renum.pdb](./structure/COV107_renum.pdb) is now ready to be used as input for ddG prediction using Rosetta.

### II. Predicting ddG using a modified high-resolution protocol of the ddG_monomer application in Rosetta
Link to ddG_monomer documentation: https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
Instead of 50 iterations, only 30 iterations were performed.

1. Pre-minimize the input structure [COV107_renum.pdb](./structure/COV107_renum.pdb)
``nohup /path/to/rosetta/main/source/bin/minimize_with_cst.static.linuxgccrelease -s /path/to/COV107_renum.pdb -in:file:fullatom -ignore_zero_occupancy false -ignore_unrecognized_res -fa_max_dis 9.0 -database /path/to/rosetta/main/database/ -ddg::harmonic_ca_tether 0.5 -score:weights /path/to/rosetta/main/database/scoring/weights/pre_talaris_2013_standard.wts -restore_pre_talaris_2013_behavior -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false -score:patch /path/to/rosetta/main/database/scoring/weights/score12.wts_patch > mincst.log 2>&1 </dev/null &``

- Input file
    - [COV107_renum.pdb](./structure/COV107_renum.pdb)
- Output file
    - [min_cst_0.5.COV107_renum_0001.pdb](./structure/min_cst_0.5.COV107_renum_0001.pdb)

2. Convert the .log file to a .cst file
``bash /path/to/rosetta/main/source/src/apps/public/ddg/convert_to_cst_file.sh ./mincst.log > ./Constraint.cst``

- Input file
    - [mincst.log](./data/mincst.log)
- Output file
    - [Constraint.cst](./data/Constraint.cst)

3. Perform ddG prediction in the background. Perform 3 independent replicates.
``nohup /path/to/rosetta/main/source/bin/ddg_monomer.static.linuxgccrelease -in:file:s /path/to/min_cst_0.5.COV107_renum_0001.pdb -ignore_zero_occupancy false -resfile F27I.resfile -ddg:weight_file soft_rep_design -ddg:minimization_scorefunction /path/to/rosetta/main/database/scoring/weights/pre_talaris_2013_standard.wts -restore_pre_talaris_2013_behavior -ddg::minimization_patch /path/to/rosetta/main/database/scoring/weights/score12.wts_patch -database /path/to/rosetta/main/database/ -fa_max_dis 9.0 -ddg::iterations 30 -ddg::dump_pdbs true -ignore_unrecognized_res -ddg::local_opt_only false -ddg::min_cst true -constraints::cst_file /path/to/Constraint.cst -ddg::suppress_checkpointing true -in::file::fullatom -ddg::mean false -ddg::min true -ddg::sc_min_only false -ddg::ramp_repulsive true -unmute core.optimization.LineMinimizer -ddg::output_silent false -out:path:all /path/to/F27I_rep1/ 2>&1 </dev/null &``

- Input file
    - [Point mutation resfiles](./mut_files/)
- Output file
    - ddg_predictions.out for each replicate

4. Compile total scores for all mutations.

- Input file
    - Scores from ddg_predictions.out
- Output file
    - [ddg_scores.csv](./data/ddg_scores.csv)
