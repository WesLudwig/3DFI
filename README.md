<p align="left"><img src="https://github.com/PombertLab/3DFI/blob/master/Misc/Logo.png" alt="3DFI - Three-dimensional function inference" width="800"></p>

The 3DFI pipeline predicts the 3D structure of proteins, and then searches for structural homology in the 3D space. Structures predicted in PDB format are searched against a local copy the [RSCB PDB](https://www.rcsb.org/) database with GESAMT (General Efficient Structural Alignment of Macromolecular Targets) from the [CCP4](https://www.ccp4.ac.uk/) package. Known PDB structures can also be searched against a set of predicted structures to identify potential structural homologs in predicted datasets.

## Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Howto](#howto)
  * [3D structure prediction](#3D-structure-prediction)
    * [RaptorX](#Raptorx---template-based-protein-structure-modeling)
    * [trRosetta](#trRosetta---deep-learning-based-protein-structure-modeling)
  * [Downloading PDB files from RCSB](#downloading-PDB-files-from-RCSB)
  * [Creating a list of PDB titles](#creating-a-list-of-PDB-titles)
  * [Creating/updating a GESAMT database](#creating/updating-a-GESAMT-database)
  * [Structural homology searches with GESAMT](#structural-homology-searches-with-GESAMT)
  * [Parsing the output of GESAMT searches](#Parsing-the-output-of-GESAMT-searches)
* [Miscellaneous](#miscellaneous)
* [Funding and acknowledgments](#Funding-and-acknowledgments)
* [References](#references)

### Introduction
###### About function inferences
Inferring the function of proteins using computational approaches usually involves performing some sort of homology search based on sequences or structures. In sequence-based searches, nucleotide or amino acid sequences are queried against known proteins or motifs using tools such as [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [DIAMOND](https://github.com/bbuchfink/diamond), [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi), or [Pfam](https://pfam.xfam.org/), but those searches may fail if the proteins investigated are highly divergent. In structure-based searches, proteins are searched instead at the 3D level for structural homologs.

###### Why structural homologs?
Because structure often confers function in biology, structural homologs often share similar functions, even if the building blocks are not the same (*i.e.* a wheel made of wood or steel is still a wheel regardless of its composition). Using this approach, we might be able to infer putative functions for proteins that share little to no similarity at the sequence level with known proteins, assuming that a structural match can be found.

###### What is needed for structure-based homology searches?
To perform structure-based predictions we need 3D structures —either determined experimentally or predicted computationally—that we can query against other structures, such as those from the [RCSB PDB](https://www.rcsb.org/). We also need tools that can search for homology at the structural level. Several tools are now available to predict protein structures, many of which are implemented as web servers for ease of use. A listing can be found at [CAMEO](https://www.cameo3d.org/), a website that evaluates their accuracy and reliability. Fewer tools are available to perform searches at the 3D levels (*e.g.* SSM and GESAMT). SSM is implemented in [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) while GESAMT is part of the [CCP4](https://www.ccp4.ac.uk/) package.

###### Why this pipeline?
Although predicting the structure of a protein and searching for structural homologs can be done online, for example by using [SWISS-MODEL](https://swissmodel.expasy.org/) and [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/), genomes often code for thousands of proteins and applying this approach on a genome scale would be time consuming and error prone. We implemented the 3DFI pipeline to enable the use of structure-based homology searches at a genome-wide level.

#### Requirements
1. [RaptorX](http://raptorx.uchicago.edu/)          ## Template-based predictions
2. [trRosetta](https://github.com/gjoni/trRosetta)  ## Deep-learning-based predictions 
3. [MODELLER](https://salilab.org/modeller/)        ## For RaptorX
4. [PyRosetta](http://www.pyrosetta.org/)           ## For trRosetta
5. [GESAMT](https://www.ccp4.ac.uk/)                ## From the CCP4 package
6. [Perl 5](https://www.perl.org/) modules - [File::Basename](https://perldoc.perl.org/File/Basename.html), [File::Find](https://perldoc.perl.org/File/Find.html), [Getopt::Long](https://perldoc.perl.org/Getopt/Long.html), [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)

### Howto
#### 3D structure prediction
##### RaptorX - template-based protein structure modeling
To perform 3D structure predictions locally with [RaptorX](http://raptorx.uchicago.edu/), the standalone programs should be [downloaded](http://raptorx.uchicago.edu/download/) and installed according to the authors’ instructions. Because of how RaptorX is implemented, 3D structure predictions should be performed from inside the RaptorX installation directory:
```Bash
cd RAPTORX_INSTALLATION_DIRECTORY/
```
Input files (in FASTA format) and the output location of the structures to be predicted can be specified with [raptorX.pl](https://github.com/PombertLab/3DFI/blob/master/raptorx.pl):
```Bash
raptorx.pl \
   -t 10 \
   -k 2 \
   -i ~/FASTA/ \
   -o ~/3D_predictions/
```
Options for raptorX.pl are:
```
-t (--threads)	Number of threads to use [Default: 10]
-i (--input)	Folder containing fasta files
-o (--output)	Output folder
-k (--TopK)	Number of top template(s) to use per protein for model building [Default: 1]
-m (--modeller)	MODELLER binary name [Default: mod10.1] ## Use absolute or relative path if not set in \$PATH
```

##### trRosetta - deep-learning-based protein structure modeling
To perform 3D structure predictions locally with [trRosetta](https://github.com/gjoni/trRosetta), [tensorflow](https://www.tensorflow.org/) version 1.15, [HH-suite3](https://github.com/soedinglab/hh-suite) and [PyRosetta](http://www.pyrosetta.org/) must be installed. A database for HHsuite3's hhblits, such as [Uniclust](https://uniclust.mmseqs.com/), should also be installed. For ease of use, [tensorflow](https://www.tensorflow.org/) 1.15 can be installed in a conda environment.

###### Tensorflow 1.15 in conda

To install tensorflow with GPU in conda:
```Bash
conda create -n tfgpu python=3.7
conda activate tfgpu
pip install tensorflow-gpu==1.15
```

To install tensorflow with CPU in conda: ## The files can eat through GPU VRAM very quickly...
```Bash
conda create -n tfcpu python=3.7
conda activate tfcpu
pip install tensorflow-cpu==1.15
```
###### Running trRosetta
Running [trRosetta](https://github.com/gjoni/trRosetta) involves 3 main steps: 1) searches with [HHsuite3](https://github.com/soedinglab/hh-suite)'s hhblits to generate alignments (.a3m); 2) prediction of protein inter-residue geometries (.npz) with [trRosetta](https://github.com/gjoni/trRosetta)'s predict.py; and 3) prediction of 3D structures (.pdb) with trRosetta.py and [PyRosetta](http://www.pyrosetta.org/). Performing these predictions on several proteins can be automated with 3DFI scripts.

1. Converting FASTA sequences to single string FASTA sequences with [fasta_oneliner.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/fasta_oneliner.pl):
```Bash
fasta_oneliner.pl \
   -f *.fasta \
   -o FASTA_OL
```

Options for [fasta_oneliner.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/fasta_oneliner.pl) are:
```
-f (--fasta)    FASTA files to convert
-o (--output)   Output folder
```

2. Running hhblits searches with [run_hhblits.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/run_hhblits.pl):
```Bash
## Running hhblits on multiple evalues independently
run_hhblits.pl \
   -t 10 \
   -f FASTA_OL/ \
   -o HHBLITS/ \
   -d /media/Data_3/Uniclust/UniRef30_2020_06 \
   -e 1e-40 1e-10 1e-03 1e+01

## Running hhblits on evalues sequentially, from stricter to more permissive
run_hhblits.pl \
   -t 10 \
   -f FASTA_OL/ \
   -o HHBLITS/ \
   -d /media/Data_3/Uniclust/UniRef30_2020_06 \
   -s \
   -se 1e-70 1e-50 1e-30 1e-10 1e-06 1e-04 1e+01
```
Options for [run_hhblits.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/run_hhblits.pl) are:
```
-t (--threads)	    Number of threads to use [Default: 10]
-f (--fasta)	    Folder containing fasta files
-o (--output)	    Output folder
-d (--database)     Uniclust database to query
-v (--verbosity)    hhblits verbosity; 0, 1 or 2 [Default: 2]

## E-value options
-e (--evalues)      Desired evalue(s) to query independently
-s (--seq_it)       Iterates sequentially through evalues
-se (--seq_ev)      Evalues to iterate through sequentially [Default:
                    1e-70 1e-60 1e-50 1e-40 1e-30 1e-20 1e-10 1e-08 1e-06 1e-04 1e+01 ]
-ne (--num_it)      # of hhblits iteration per evalue (-e) [Default: 3]
-ns (--num_sq)      # of hhblits iteration per sequential evalue (-s) [Default: 1] 
```

3. Create .npz files containing inter-residue geometries with [create_npz.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/create_npz.pl):
```Bash
create_npz.pl \
   -a HHBLITS/*.a3m \
   -o NPZ/ \
   -p /media/Data_3/opt/trRosetta/network/predict.py \
   -m /media/Data_3/opt/trRosetta/model2019_07
```

Options for [create_npz.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/create_npz.pl) are:
```
-a (--a3m)      .a3m files generated by hhblits
-o (--output)   Output folder
-p (--predict)  Path to predict.py from trRosetta
-m (--model)    Path to trRosetta model directory
```

4. Generate .pdb files containing 3D models from the .npz fil3s with [create_pdb.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/create_pdb.pl):
```Bash
create_pdb.pl \
   -c 10 \
   -n NPZ/*.npz \
   -o PDB/ \
   -f FASTA_OL/ \
   -t /media/Data_3/opt/trRosetta/pdb/trRosetta.py
```

Options for [create_pdb.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/create_pdb.pl) are:
```
-c (--cpu)	Number of cpu threads to use [Default: 10] ## i.e. runs n processes in parallel
-m (--memory)	Memory available (in Gb) required before thread launch [Default: 16] 
-n (--npz)	.npz files generated by hhblits
-o (--output)	Output folder [Default: ./]
-f (--fasta)	Folder containing the oneliner fasta files
-t (--trosetta)	Path to trRosetta.py from trRosetta
-p (--python)	Preferred Python interpreter [Default: python]
```

5. The .pdb files thus generated contain lines that are not standard and that can prevent applications such as [PDBeFOLD](https://www.ebi.ac.uk/msd-srv/ssm/) to run on the corresponding files. We can clean up the PDB files with [sanitize_pdb.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/sanitize_pdb.pl):
```Bash
sanitize_pdb.pl \
   -p PDB/*.pdb \
   -o PDB_clean
```

Options for [sanitize_pdb.pl](https://github.com/PombertLab/3DFI/blob/master/trRosetta/sanitize_pdb.pl) are:
```
-p (--pdb)      .pdb files generated by trRosetta
-o (--output)   Output folder
```

#### Downloading PDB files from RCSB
PDB files from the [Protein Data Bank](https://www.rcsb.org/) can be downloaded directly from its website. Detailed instructions are provided [here](https://www.wwpdb.org/ftp/pdb-ftp-sites). Because of the large size of this dataset, downloading it using [rsync](https://rsync.samba.org/) is recommended. This can be done as follows, wherein **/path/to/PDB/** should be replaced by the desired directory. PDB files (pdb*.ent.gz) will be located in subdirectories therein.

```Bash
rsync -rlpt -v -z --delete --port=33444
rsync.rcsb.org::ftp_data/structures/divided/pdb/ /path/to/PDB/
```
NOTE: For ease of use, [PDB_update.sh](https://github.com/PombertLab/3DFI/blob/master/PDB_update.sh) can also be modified, then used to rsync the PDB files.

#### Creating a list of PDB titles
To create a tab-delimited list of PDB entries and their titles from the downloaded PDB gzipped files (pdb*.ent.gz), we can use [PDB_headers.pl](https://github.com/PombertLab/3DFI/blob/master/PDB_headers.pl):
```Bash
PDB_headers.pl \
   -p /path/to/PDB/ \
   -o /path/to/PDB_titles.tsv
```
Options for [PDB_headers.pl](https://github.com/PombertLab/3DFI/blob/master/PDB_headers.pl) are:
```
-p (--pdb)	Directory containing PDB files downloaded from RCSB PDB/PDBe (gzipped)
-o (--output)	Output file in tsv format
```
The list created should look like this:
```
4hhg	CRYSTAL STRUCTURE OF THE PSEUDOMONAS AERUGINOSA AZURIN, RUH107NO YOH109
5xqj	CRYSTAL STRUCTURE OF A PL 26 EXO-RHAMNOGALACTURONAN LYASE FROM PENICILLIUM CHRYSOGENUM COMPLEXED WITH UNSATURATED GALACTURONOSYL RHAMNOSE SUBSTITUTED WITH GALACTOSE
1vqa	GENE V PROTEIN MUTANT WITH VAL 35 REPLACED BY ALA 35 AND ILE 47 REPLACED BY LEU 47 (V35A, I47L)
6rxa	EDDS LYASE VARIANT D290M/Y320M WITH BOUND FORMATE
```

#### Creating/updating a GESAMT database
Before performing structural homology searches with GESAMT, we should first create an archive to speed up the searches. We can also update the archive later as sequences are added (for example after the RCSB PDB files are updated with rsync). GESAMT archives can be created/updated with [run_GESAMT.pl](https://github.com/PombertLab/3DFI/blob/master/run_GESAMT.pl):
```Bash
## To create a GESAMT archive
run_GESAMT.pl \
   -cpu 10 \
   -make \
   -arch /path/to/GESAMT_ARCHIVE \
   -pdb /path/to/PDB/

## To update a GESAMT archive
run_GESAMT.pl \
   -cpu 10 \
   -update \
   -arch /path/to/GESAMT_ARCHIVE \
   -pdb /path/to/PDB/
```
Options for [run_GESAMT.pl](https://github.com/PombertLab/3DFI/blob/master/run_GESAMT.pl) are:
```
-c (--cpu)	CPU threads [Default: 10]
-a (--arch)	GESAMT archive location [Default: ./]
-m (--make)	Create a GESAMT archive
-u (--update)	Update existing archive
-p (--pdb)	Folder containing RCSB PDB files to archive
```

#### Structural homology searches with GESAMT
Structural homology searches with GESAMT can also be performed with [run_GESAMT.pl](https://github.com/PombertLab/3DFI/blob/master/run_GESAMT.pl):
```Bash
run_GESAMT.pl \
   -cpu 10 \
   -query \
   -arch /path/to/GESAMT_ARCHIVE \
   -input /path/to/*.pdb \
   -o /path/to/RESULTS_FOLDER \
   -mode normal
```
Options for [run_GESAMT.pl](https://github.com/PombertLab/3DFI/blob/master/run_GESAMT.pl) are:
```
-c (--cpu)	CPU threads [Default: 10]
-a (--arch)	GESAMT archive location [Default: ./]
-q (--query)	Query a GESAMT archive
-i (--input)	PDF files to query
-o (--outdir)	Output directory [Default: ./]
-d (--mode)	Query mode: normal of high [Default: normal]
```

Results of GESAMT homology searches will be found in the \*.gesamt files generated. Content of these should look like:
```
#  Hit   PDB  Chain  Q-score  r.m.s.d     Seq.  Nalign  nRes    File
#  No.   code   Id                         Id.                  name
     1   5HVM   A     0.8841   0.4410   0.4465    430    446   pdb5hvm.ent.gz
     2   5HVM   B     0.8639   0.5333   0.4419    430    452   pdb5hvm.ent.gz
     3   5HVO   D     0.8387   0.6602   0.4242    429    456   pdb5hvo.ent.gz
     4   5HVO   B     0.8356   0.6841   0.4252    428    454   pdb5hvo.ent.gz
```

#### Parsing the output of GESAMT searches
To add definitions/products to the PDB matches found with GESAMT, we can use the list generated by [PDB_headers.pl](https://github.com/PombertLab/3DFI/blob/master/PDB_headers.pl) together with [descriptive_GESAMT_matches.pl](https://github.com/PombertLab/3DFI/blob/master/descriptive_GESAMT_matches.pl):
```Bash
descriptive_GESAMT_matches.pl \
   -r /path/to/PDB_titles.tsv \
   -m *.gesamt \
   -q 0.3 \
   -o /path/to/GESAMT.matches
```
Options for [descriptive_GESAMT_matches.pl](https://github.com/PombertLab/3DFI/blob/master/descriptive_GESAMT_matches.pl) are:
```
-r (--rcsb)	Tab-delimited list of RCSB structures and their titles ## see PDB_headers.pl 
-p (--pfam)	Tab-delimeted list of PFAM stuctures and their titles (http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz)
-m (--matches)	Results from GESAMT searches ## see run_GESAMT.pl
-q (--qscore)	Q-score cut-off [Default: 0.3]
-b (--best)	Keep the best match(es) only (top X hits)
-o (--output)	Output name [Default: Gesamt.matches]
```

The concatenated list generated should look like:
```
### GPU96_01g00180-m1-3se6A
GPU96_01g00180-m1-3se6A       1   3SE6   B     0.8635   0.4638   0.2713    811    860   pdb3se6.ent.gz   CRYSTAL STRUCTURE OF THE HUMAN ENDOPLASMIC RETICULUM AMINOPEPTIDASE 2
GPU96_01g00180-m1-3se6A       2   3SE6   A     0.8611   0.4749   0.2699    815    870   pdb3se6.ent.gz   CRYSTAL STRUCTURE OF THE HUMAN ENDOPLASMIC RETICULUM AMINOPEPTIDASE 2
### GPU96_01g00200-m1-1j3iC
GPU96_01g00200-m1-1j3iC       1   5H39   A     0.8836   0.8169   0.6064    282    285   pdb5h39.ent.gz   STRUCTURAL ANALYSIS OF KSHV THYMIDYLATE SYNTHASE
GPU96_01g00200-m1-1j3iC       2   3IHI   B     0.8829   0.7784   0.6464    280    283   pdb3ihi.ent.gz   CRYSTAL STRUCTURE OF MOUSE THYMIDYLATE SYNTHASE
### GPU96_01g00840-m1-5hvmA
GPU96_01g00840-m1-5hvmA       1   5HVM   A     0.8841   0.4410   0.4465    430    446   pdb5hvm.ent.gz   STRUCTURE OF ASPERGILLUS FUMIGATUS TREHALOSE-6-PHOSPHATE SYNTHASE A INCOMPLEX WITH UDP AND VALIDOXYLAMINE A
GPU96_01g00840-m1-5hvmA       2   5HVM   B     0.8639   0.5333   0.4419    430    452   pdb5hvm.ent.gz   STRUCTURE OF ASPERGILLUS FUMIGATUS TREHALOSE-6-PHOSPHATE SYNTHASE A INCOMPLEX WITH UDP AND VALIDOXYLAMINE A

```

#### Miscellaneous 
###### Splitting multifasta files
Single fasta files for structure prediction with raptorx.pl can be created with [split_Fasta.pl](https://github.com/PombertLab/3DFI/blob/master/split_Fasta.pl):
```
split_Fasta.pl \
   -f file.fasta \
   -o output_folder \
   -e fasta
```

Options for [split_Fasta.pl](https://github.com/PombertLab/3DFI/blob/master/split_Fasta.pl) are:
```
-f (--fasta)	FASTA input file (supports gzipped files)
-o (--output)	Output directory; defaults to file name prefix
-e (--ext)	Desired file extension [Default: fasta]
```

###### Splitting PDB files
RCSB PDB files can be split per chain with [split_PDB.pl](https://github.com/PombertLab/3DFI/blob/master/split_PDB.pl):
```
split_PDB.pl \
   -p files.pdb \
   -o output_folder \
   -e pdb
```

Options for [split_PDB.pl](https://github.com/PombertLab/3DFI/blob/master/split_PDB.pl) are:
```
-p (--pdb)	PDB input file (supports gzipped files)
-o (--output)	Output directory. If blank, will create one folder per PDB file based on file prefix
-e (--ext)	Desired file extension [Default: pdb]
```

###### Renaming files
Files can be renamed using regular expressions with [rename_files.pl](https://github.com/PombertLab/3DFI/blob/master/rename_files.pl):
```
rename_files.pl \
   -o 'i{0,1}-t26_1-p1' \
   -n '' \
   -f *.fasta
```

Options for [rename_files.pl](https://github.com/PombertLab/3DFI/blob/master/rename_files.pl) are:
```
-o (--old)	Old pattern/regular expression to replace with new pattern
-n (--new)	New pattern to replace with; defaults to blank [Default: '']
-f (--files)	Files to rename
```

## Funding and acknowledgments
This work was supported by the National Institute of Allergy and Infectious Diseases of the National Institutes of Health (award number R15AI128627) to Jean-Francois Pombert. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

##### REFERENCES
1) [RCSB Protein Data Bank: Sustaining a living digital data resource that enables breakthroughs in scientific research and biomedical education](https://pubmed.ncbi.nlm.nih.gov/29067736/). **Burley SK, et al.** Protein Sci. 2018 Jan;27(1):316-330. PMID: 29067736 PMCID: PMC5734314 DOI: 10.1002/pro.3331

2) [RaptorX: exploiting structure information for protein alignment by statistical inference](https://pubmed.ncbi.nlm.nih.gov/21987485/). **Peng J, Xu J.** Proteins. 2011;79 Suppl 10:161-71. PMID: 21987485 PMCID: PMC3226909 DOI: 10.1002/prot.23175

3) [Template-based protein structure modeling using the RaptorX web server](https://pubmed.ncbi.nlm.nih.gov/22814390/). **Källberg M, et al.** Nat Protoc. 2012 Jul 19;7(8):1511-22. PMID: 22814390 PMCID: PMC4730388 DOI: 10.1038/nprot.2012.085

4) [Enhanced fold recognition using efficient short fragment clustering](https://pubmed.ncbi.nlm.nih.gov/27882309/). **Krissinel E.** J Mol Biochem. 2012;1(2):76-85. PMID: 27882309 PMCID: PMC5117261

5) [Overview of the CCP4 suite and current developments](https://pubmed.ncbi.nlm.nih.gov/21460441/). **Winn MD et al.** Acta Crystallogr D Biol Crystallogr. 2011 Apr;67(Pt 4):235-42. PMID: 21460441 PMCID: PMC3069738 DOI: 10.1107/S0907444910045749

6) [Improved protein structure prediction using predicted interresidue orientations](https://pubmed.ncbi.nlm.nih.gov/31896580/). **Yang J, Anishchenko I, Park H, Peng Z, Ovchinnikov S, Baker D.** Proc Natl Acad Sci USA. 2020 Jan 21;117(3):1496-1503.PMID: 31896580 PMCID: PMC6983395 DOI: 10.1073/pnas.1914677117
