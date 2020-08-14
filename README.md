<p align="left"><img src="https://github.com/PombertLab/3DFI/blob/master/Misc/Logo.png" alt="3DFI - Three-dimensional function inference" width="600"></p>

The 3DFI pipeline predicts the 3D structure of proteins, and then searches for structural homology in the 3D space. Structures predicted in PDB format are searched against a local copy the [RSCB PDB](https://www.rcsb.org/) database with GESAMT (General Efficient Structural Alignment of Macromolecular Targets) from the [CCP4](https://www.ccp4.ac.uk/) package. Known PDB structures can also be searched against a set of predicted structures to identify potential structural homologs in predicted datasets.

## Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Howto](#howto)
  * [3D structure prediction](#3D-structure-prediction)
  * [Downloading PDB files from RCSB](#downloading-PDB-files-from-RCSB)
  * [Creating a list of PDB titles](#creating-a-list-of-PDB-titles)
  * [Creating/updating a GESAMT database](#creating/updating-a-GESAMT-database)
  * [Structural homology searches with GESAMT](#structural-homology-searches-with-GESAMT)
* [Miscellaneous](#miscellaneous)
* [References](#references)

### Introduction
###### About function inferences
Inferring the function of proteins using computational approaches usually involves performing some sort of homology search based on sequences or structures. In sequence-based searches, nucleotide or amino acid sequences are queried against known proteins or motifs using tools such as [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [DIAMOND](https://github.com/bbuchfink/diamond), [CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi), or [Pfam](https://pfam.xfam.org/), but those searches may fail if the proteins investigated are highly divergent. In structure-based searches, proteins are searched instead at the 3D level for structural homologs.

###### Why structural homologs?
Because structure often confers function in biology, structural homologs often share similar functions, even if the building blocks are not the same (*i.e.* a wheel made of wood or steel is still a wheel regardless of its composition). Using this approach, we might be able to infer putative functions for proteins that share little to no similarity at the sequence level with known proteins, assuming that a structural match can be made.

###### What is needed for structure-based homology searches?
To perform structure-based predictions we need 3D structures —either determined experimentally or predicted computationally—that we can query against other structures, such as those from the [RCSB PDB](https://www.rcsb.org/). We also need tools that can search for homology at the structural level. Several tools are now available to predict protein structures, many of which are implemented as web servers for ease of use. A listing can be found at [CAMEO](https://www.cameo3d.org/), a website that evaluates their accuracy and reliability. Fewer tools are available to perform searches at the 3D levels (*e.g.* SSM and GESAMT). SSM is implemented in [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) while GESAMT is part of the [CCP4](https://www.ccp4.ac.uk/) package.

###### Why this pipeline?
Although predicting the structure of a protein and searching for structural homologs can be done online, for example by using [SWISS-MODEL](https://swissmodel.expasy.org/) and [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/), genomes often code for thousands of proteins and applying this approach on a genome scale would be time consuming and error prone. We implemented the 3DFI pipeline to enable the use of structure-based homology searches at a genome-wide level.

#### Requirements
1. RaptorX - http://raptorx.uchicago.edu/
2. MODELLER - https://salilab.org/modeller/
3. GESAMT -  https://www.ccp4.ac.uk/
4. Perl modules - [File::Basename](https://perldoc.perl.org/File/Basename.html), [File::Find](https://perldoc.perl.org/File/Find.html), [Getopt::Long](https://perldoc.perl.org/Getopt/Long.html), [PerlIO::gzip](https://metacpan.org/pod/PerlIO::gzip)

### Howto
#### 3D structure prediction
To perform 3D structure predictions locally with [RaptorX](http://raptorx.uchicago.edu/), the standalone programs should be [downloaded](http://raptorx.uchicago.edu/download/) and installed according to the authors’ instructions. Because of how RaptorX is implemented, 3D structure predictions should be performed from inside the RaptorX installation directory:
```
cd RAPTORX_INSTALLATION_DIRECTORY/
```
Input files (in FASTA format) and the output location of the structures to be predicted can be specified with [raptorX.pl](https://github.com/PombertLab/3DFI/blob/master/raptorx.pl):
```
raptorx.pl -t 10 -k 2 -i ~/FASTA/ -o ~/3D_predictions/
```
Options for raptorX.pl are:
```
-t (--threads)	Number of threads to use [Default: 10]
-i (--input)	Folder containing fasta files
-o (--output)	Output folder
-k (--TopK)	Number of top template(s) to use per protein for model building [Default: 1]
-m (--modeller)	MODELLER binary name [Default: mod9.23] ## Use absolute or relative path if not set in \$PATH
```

#### Downloading PDB files from RCSB
PDB files from the [Protein Data Bank](https://www.rcsb.org/) can be downloaded directly from its website. Detailed instructions are provided [here](https://www.wwpdb.org/ftp/pdb-ftp-sites). Because of the large size of this dataset, downloading it using [rsync](https://rsync.samba.org/) is recommended. This can be done as follows, wherein **/path/to/PDB/** should be replaced by the desired directory. PDB files (pdb*.ent.gz) will be located in subdirectories therein.

```
rsync -rlpt -v -z --delete --port=33444
rsync.rcsb.org::ftp_data/structures/divided/pdb/ /path/to/PDB/
```
NOTE: For ease of use, [PDB_update.sh](https://github.com/PombertLab/3DFI/blob/master/PDB_update.sh) can also be modified, then used to rsync the PDB files.

#### Creating a list of PDB titles
To create a tab-delimited list of PDB entries and their titles from the downloaded PDB gzipped files (pdb*.ent.gz), we can use [PDB_headers.pl](https://github.com/PombertLab/3DFI/blob/master/PDB_headers.pl):
```
PDB_headers.pl -p /path/to/PDB/ -o /path/to/PDB_titles.tsv
```
Options for PDB_headers.pl are:
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
Insert text...

```
run_GESAMT.pl -cpu 10 -make -arch /path/to/GESAMT_ARCHIVE -pdb /path/to/PDB/	## CREATE DB
```
```
run_GESAMT.pl -cpu 10 -update -arch /path/to/GESAMT_ARCHIVE -pdb /path/to/PDB/	## UPDATE DB
```

#### Structural homology searches with GESAMT
Insert text...

```
run_GESAMT.pl -cpu 10 -query -arch /path/to/GESAMT_ARCHIVE -input /path/to/*.pdb -o /path/to/RESULTS_FOLDER -mode normal
```

6) Parse the GESAMT output to add definitions/products to the PDB matches found
```
descriptive_GESAMT_matches.pl -t /path/to/PDB_titles.tsv -m *.gesamt -q 0.3 -o /path/to/GESAMT.matches
```

#### Miscellaneous 
1. Single fasta files for structure prediction with raptorx.pl can be created with split_Fasta.pl
2. Files can be renamed using regular expressions with rename_files.pl
3. RCSB PDB files can be split per chain with split_PDB.pl


##### REFERENCES
1) [RCSB Protein Data Bank: Sustaining a living digital data resource that enables breakthroughs in scientific research and biomedical education
Burley SK, et al.](https://pubmed.ncbi.nlm.nih.gov/29067736/) Protein Sci. 2018 Jan;27(1):316-330. PMID: 29067736 PMCID: PMC5734314 DOI: 10.1002/pro.3331

2) [RaptorX: exploiting structure information for protein alignment by statistical inference
Peng J, Xu J.](https://pubmed.ncbi.nlm.nih.gov/21987485/) Proteins. 2011;79 Suppl 10:161-71. PMID: 21987485 PMCID: PMC3226909 DOI: 10.1002/prot.23175

3) [Template-based protein structure modeling using the RaptorX web server
Källberg M, et al.](https://pubmed.ncbi.nlm.nih.gov/22814390/) Nat Protoc. 2012 Jul 19;7(8):1511-22. PMID: 22814390 PMCID: PMC4730388 DOI: 10.1038/nprot.2012.085

4) [Enhanced fold recognition using efficient short fragment clustering
Krissinel E.](https://pubmed.ncbi.nlm.nih.gov/27882309/) J Mol Biochem. 2012;1(2):76-85. PMID: 27882309 PMCID: PMC5117261

5) [Overview of the CCP4 suite and current developments
Winn MD et al.](https://pubmed.ncbi.nlm.nih.gov/21460441/) Acta Crystallogr D Biol Crystallogr. 2011 Apr;67(Pt 4):235-42. PMID: 21460441 PMCID: PMC3069738 DOI: 10.1107/S0907444910045749
