### different functions concerning clustalw alignments.

from __future__ import division, print_function
import os
from os.path import exists
import subprocess
import sys
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import AlignInfo
from Bio import AlignIO
import config as cf



def columnsummary( col ):
    """returns a summary tupel of the given column - used to delimit the
    primer regions, the regions within which primers can be located."""

    # if we see an intron symbol, we must return (1, 0) - a dummy value
    # since there cannot be 0 nucleotides in total with 1 different nucs.
    # we return this value so we can score intron columns separately.
    # see the score matrix in config.py

    # in general mode (non-intron mode), X's are simply jokers, but we still
    # don't want them in our primers (nor do we want N's):

    if 'X' in col or 'x' in col or 'N' in col or 'n' in col:
        return (1, 0)


    d={}# holds different nucleotides
    n=0 # number of nucleotides in total

        # don't count gap symbols:
    for c in col:
        if c != '-':
            d[c] = 1
            n+=1
    return (len(d), n) # return (#different nucleotides, total #nucleotides)




def printslice( allseq, s, e ):
    """prints a slice of the given alignment (given as a list of strings), including consensus *'s"""
    # here's how the list of sequences might be retrieved:
    # allseq = alignment.get_all_seqs()
    for i in allseq:
        print(i.seq.data[s:e])

    stars = ''

    for i in range(s, e):
        col = ''
        for j in allseq[0:]:
            col = ''.join( [col, j.seq.data[i]] )
        a,b = columnsummary( col )
        if a == 1 and b>1:
            stars = ''.join( [stars, '*'] )
        else:
            stars = ''.join( [stars, ' '] )
    print(stars)




def create_alignment( filename, verbose=1, outfile = 'test.aligned.fasta' ):
    """This function performs a MAFFT alignment"""

    jfilename = filename
    if not exists( jfilename ):
        sys.exit( "No such file: %s"%jfilename )

    # let's see if the given file isn't already aligned
    try:

        alignment = AlignIO.read(filename, 'fasta')
        jfilename = filename.rsplit('.', 1)[0]
        if verbose:
            print("File '%s' is already aligned. Skipping alignment step." % filename)

    except:
        # deprecated:
        #cline = MultipleAlignCL( jfilename, cf.ClustalPath )
        # set output filename:
        #cline.set_output( jfilename+'.aln' )
        #cline.set_dna_matrix( cf.DNAMatrixPath )
        # ang. nedenstaaende: nuvaerende matrix er god i alle testede tilfaelde.
        ## current matrix is good in all tested cases.
        #cline = ClustalwCommandline( "clustalw", infile=jfilename, outfile=jfilename+'.aln', dnamatrix=cf.DNAMatrixPath )
        cline = ['mafft', '--localpair', '--maxiterate', '16', '--inputorder', '--preservecase', '--quiet', jfilename   ]

        #cline.gap_open_pen = 0.001
        # dette sammen med multalinDNAmatrix.clustal giver bedre alignment
        # hvis sekvenserne er meget ens (se Version2/test.fasta)
        ##  this along with multalinDNAmatrix.clustal provides better alignment
        ##  if the sequences are very similar (see Version2/test.fasta)

        # hvis linien udkommenteres og mymatrix_identity5.clustal bruges i
        # stedet, giver det bedre alignment hvis sekvenserne ikke ligner hinanden
        # helt saa meget (?) (f.eks. Paper/NAR/revieweralignment.fasta)
        ##  if the line is commented out and mymatrix_identity5.clustal used
        ##  instead, it provides better alignment if the sequences do not
        ##  resemble each other quite so much (?)
        ##  (eg. Paper / NAR / revieweralignment.fasta)

        #cline.gap_ext_pen = 0.00 # 0.01 seems to give other regions, better?

        # normalt var gappen ikke sat.
        ## normally gaps were not set.

        # cline.max_div = 100

        if verbose:
            print('-- running this command:', cline)

	    # deprecated:
        #alignment = Clustalw.do_alignment(cline)
        stdout = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=None).communicate()[0]
        with open(jfilename+'.mafft', 'wb') as handle:
            handle.write(stdout)
        alignment = AlignIO.read(jfilename+'.mafft', 'fasta')

    # make sure sequences are upper case
    for rec in alignment:
        rec.seq = rec.seq.upper()

    #
    # indices from 0 to len-1
    #
    # things you can do with the returned alignment object:
    #
    # deprecated:
    #allseq = alignment.get_all_seqs()
    allseq = list(alignment)
    summary = AlignInfo.SummaryInfo(alignment)
    l = alignment.get_alignment_length()

    return (allseq, summary, l)
