# DiscoMark

---

Discovery of phylogenetic markers from orthologous sequences.

---

## Citation

Detering H, Rutschmann S, Simon S, Fredslund J, Monaghan MT. 2016. DiscoMark: Nuclear marker discovery from orthologous sequences using low coverage genome data. bioRxiv doi: 10.1101/047282


## Requirements

DiscoMark is a python script and depends on several programs. The below instructions should work well for Unix-based (Linux, Apple OS X, etc.) operation systems. It might also be possible to run the program under Windows, however this has not been tested.

DiscoMark, uses the following programs:
* [Python](https://www.python.org) (>= 2.7)
* [NCBI Blast+](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST) (>= 2.2.29)
* [MAFFT](http://mafft.cbrc.jp/alignment/software) (>= 6.864b)
* [TrimAl](https://github.com/scapella/trimal) (>= 1.2)

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
export PATH=$PATH:/path/to/discomark/bin/linux
```

If you have the programs BLAST+, MAFFT and TrimAl installed on your system you don't need the `export PATH=[...]` statement.

## How to run DiscoMark

Make sure your [input data](https://github.com/hdetering/discomark/wiki#input-data) is formatted as FASTA format. 
Let´s say you want to discover markers for two species (i.e. species1 and species2) and use a reference (i.e. reference.fasta), you will call DiscoMark like this:

```
cd discomark
python run_project.py -i example/hamstr/species1 -i example/hamstr/species2 -r example/reference/reference.fasta -d output
```

See the wiki for more info on the [command line options](https://github.com/hdetering/discomark/wiki/Command-Line-Options).


## Results

Inspect the marker report in the output directory at `7_report/discomark_results.html`.  
(In the example above that would be: `output/7_report/discomark_results.html`.)

The primer table can also be found in CSV format in `7_report/primers.xls`.
