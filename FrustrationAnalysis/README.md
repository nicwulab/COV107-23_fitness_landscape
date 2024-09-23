# Frustration Analysis 

Workflow and code to perform frustration analysis for COV107-23 WT and mutants

### Dependencies

* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [Rosetta](https://rosettacommons.org) (version 3.11)
* [PyMOL](https://www.pymol.org) (version 2.4.0)
* [spectrumany.py](./script/spectrumany.py) from https://pymolwiki.org/index.php/Spectrumany
* [data2bfactor.py](./script/data2bfactor.py) from https://gist.github.com/Croydon-Brixton/ca916769c285faed69740d384b5daf23 

### Input files

- [COV107_renum.pdb](./pdb/COV107_renum.pdb) - Identical to the one used for ddG calculation with Rosetta
- [Resfiles](./resfile/)

### I. Generate point mutations in COV107-23.
1. Point mutagenesis was run using `<rosetta_location>/main/source/bin/fixbb.static.macosclangrelease -s COV107_renum.pdb -resfile <Mutation>.resfile -nstruct 100`. One-hundred poses in total were generated.
    * Change `macosclangrelease` depending on your operating system.
2. Using the `score.sc` file, identify the pose with the lowest score. Use the lowest-scoring pose for downstream analysis.

### II. Submit frustration analysis to server.
1. Upload [wild-type COV107-23](./pdb/COV107_renum.pdb), or the lowest-scoring pose of mutants S53P, S35T, or F27L (found in [PDB][./pdb]) to the server http://frustratometer.qb.fcen.uba.ar. Select both chains A and B; do not use electrostatics. Run the frustration analysis through the server.
2. Download output files.

### III. Analyze frustration data.
1. Single residue frustration index.
    - Navigate to the FrustrationData folder of the server's output. Use the '_singleresidue' file, and add a '.txt' extension. These can be found in the [data](./data/) subfolder.
    - Run the R code `FrustrationAnalysis.R`. Note that the working directory should be set to the 'FrustrationAnalysis' folder.

2. Contact analysis.
    - Navigate to the VisualizationScripts folder of the server's output. Use the 'draw_links.py' file. These are located in the [script](./script) subfolder.
    - In the VisualizationScripts folder, open the '*_configurational.pml' file. Identify the contacts for each residue (site 27, 35, or 53; numbered 26, 34, and 52, respectively, in PDB files due to renumbering).
    - Run the following scripts in PyMOL: [Site53-WT.pml](./script/Site53-WT.pml), [Site53-Mut.pml](./script/Site53-Mut.pml), [Site35-WT.pml](./script/Site35-WT.pml), [Site35-Mut.pml](./script/Site35-Mut.pml), [Site27-WT.pml](./script/Site27-WT.pml), and [Site27-Mut.pml](./script/Site27-Mut.pml). Change directory to the 'FrustrationAnalysis' folder, and then run `run ./script/SiteX-WT/Mut.pml`
    
### IV. Structural analysis.
1. Run [Structure_Figs.pml](./script/Structure_Figs.pml) in PyMOL by changing directory to the 'FrustrationAnalysis' folder, and then running `run ./script/Structure_Figs.pml`. All image outputs should be in the [image](./image/) subfolder.
