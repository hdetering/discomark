# DiscoMark

---

Discovery of phylogenetic markers from orthologous sequences.

---

## Citation

If you use DiscoMark in you research, please cite:

Detering H, Rutschmann S, Simon S, Fredslund J, Monaghan MT (2016) DiscoMark: Nuclear marker discovery from orthologous sequences using low coverage genome data. *BioRxiv* **047282** doi: 10.1101/047282


## Requirements

DiscoMark is a python script and depends on several programs. The below instructions should work well for Unix-based (Linux, Apple OS X, etc.) operation systems. It might also be possible to run the program under Windows, however this has not been tested.

DiscoMark, uses the following programs:
* [Python](https://www.python.org) (>= 2.7)
* [NCBI Blast+](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST) (>= 2.2.29)
* [MAFFT](http://mafft.cbrc.jp/alignment/software) (>= 6.864b)
* [TrimAl](https://github.com/scapella/trimal) (>= 1.4)

For Linux and Mac, we include binary files for BLASTn, MAFFT and TrimAl in the respective `bin/` folders. See the [Installation](#installation) section on how to use them.

Additionally, you'll need the following Python packages:
* [Biopython](http://biopython.org/) (>= 1.62)
* [SqlAlchemy](http://www.sqlalchemy.org/) (>= 0.9)

To facilitate the installation of these packages, we suggest to use the python module manager [pip](https://pypi.python.org/pypi/pip) (which normally comes with python3). To check which version is available on your computer type:
```
pip -V
```

To install the two python packages [Biopython](www.biopython.org/) and [SQLAlchemy](www.sqlalchemy.org/) you simply type the following:
```
sudo pip install biopython 
sudo pip install sqlalchemy
```

## Installation

Download DiscoMark (when setting the `PATH` environment variable, replace `linux` with `mac` if you're on a Mac):
```
git clone https://github.com/hdetering/discomark.git
export PATH=$PATH:$PWD/discomark/bin/linux
```

If you have the programs BLAST+, MAFFT and TrimAl installed on your system you don't need the `export PATH=[...]` statement. In case you want to use the provided binaries, make sure to execute the `export PATH=[...]` statement each time you login or make the change persistent in your user profile (works slightly different for [linux](http://superuser.com/questions/324832/how-can-i-permanently-add-a-path-to-my-bash-profile) and [mac](http://hathaway.cc/post/69201163472/how-to-edit-your-path-environment-variables-on-mac)).

## How to run DiscoMark

Make sure your [input data](https://github.com/hdetering/discomark/wiki#input-data) is formatted as FASTA format.   

Let´s say you want to discover markers for two species (i.e. species1 and species2) and use a reference (i.e. reference.fasta), you will call DiscoMark like this:

```
cd discomark
python run_project.py -i example/hamstr/species1 -i example/hamstr/species2 -r example/reference/reference.fasta -d output
```

If you want to add an annotation file for the input markers, you will call DiscoMark like this:

```
cd discomark
python run_project.py -i example/hamstr/species1 -i example/hamstr/species2 -r example/reference/reference.fasta -a input/co2go.ixosc.csv -d output
```

Under the defualt settings you will perform on online BLAST search for you primer pairs. Let´s say you do not have an internet connection or do not want to use this option you can call DiscoMark like this:

```
cd discomark
python run_project.py -i example/hamstr/species1 -i example/hamstr/species2 -r example/reference/reference.fasta -a input/co2go.ixosc.csv -d output --no-primer-blast
```

Please see the wiki for the complete information on the [command line options](https://github.com/hdetering/discomark/wiki/Command-Line-Options).


## Results

An interactive HTML file in the output directory at `7_report/discomark_results.html` contains all important information about the designed primer pairs.  
(In the example above that would be: `output/7_report/discomark_results.html`.)

A table with all primer pairs is also avialable at CSV format in `7_report/primers.xls`.



Notice: This program may contain errors. Please inspect results carefully.   
This program comes with ABSOLUTELY NO WARRANTY.   
This is free software, and you are welcome to redistribute it under certain conditions


## Experiencing Problems?

Have a look at the [FAQ](https://github.com/hdetering/discomark/wiki/Frequently-Asked-Questions). 

If the FAQ is not helpful, you can [open an issue](https://github.com/hdetering/discomark/issues/new) to inform us about your problem.
