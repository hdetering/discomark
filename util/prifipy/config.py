from os.path import exists
from sys import modules
thismodule = modules[__name__] # a handle to this module

##
## global configuration parameters:
##



#Basedir = '/users/chili/BiRC/MolBio/Data/PrimerDesign/'
Basedir = '/home/harrem/dev/discomark/util/prifi/upgraded/'

# do our alignments have the XX intron markers and do we use that intron info
# in the scoring?
INTRONS = 'yes'

# how many bp two primers must overlap to be considered almost the same:
MatchOverlap = 10

PrimerPairSuggestions = 4 # how many non-overlapping primer pairs should be suggested

TailLength = 8 # length of tail in 3'-end checked for AT content

MaxPrimerLength = 35

MinPrimerLength = 15

SuggestedMaxTm = 77.0

MinTm = 55.0

CriticalTm = 58.0 # if both primer Tm's below this number, penalize the pair

MinTmWithMismatchesAllowed = 59.0

Min3endPerfectMatches = 2

MaxMismatches = 4

MinLengthWithThreeAmbiguities = 25 # shorter than 25: can have at most 2 amb.

MinLengthWithTwoAmbiguities = 21 # shorter than 21: can have at most 1 amb.

# if you have four, the first and last must be more than this number apart:
WindowWithFourMismatches = 17

# Min3endPerfectMatchesIfPrimerHasMaxMismatches = 6

# penalize ambiguities closer than this distance to 3'-end:
CriticalAmbiguityDistanceTo3End = 5

# MaxPrimerPrimerDistance = 3000

# MinLengthOfLongestInterPrimerExon = 50 # ID-exon

#MaxPrimerIDexonDistance = 300  # max distance from either primer to ID-exon

#IDexonBonus = 5 # reward for having a reachable ID exon

MaxDiversityPerColumn = 4

MaxMismatchesWithMaxDiversity = 1

# MinNumberOfNucleotidesPerAmbiguity = 8.0 # e.g. length 24 is too short for 3 amb.

MinNucleotidesPerColumn = 2

PrimerConcentrationNmolar = 250 # nanomolaer

SaltConcentrationMolar = .05 # [K+] i molaer

#MismatchTmDecrease = 2  # approximate presumed Tm decrease per mismatch

MaxPrimerPairTmDifference = 15 # max Tm difference between forw. and rev. primers

# optimal primer length interval, both inclusive:
OptimalPrimerLength = [ 25, 35 ] #MaxPrimerLength )

# optimal product length interval (optimum between two middle numbers).
# linear increase/decrease to/from optimum (i.e. linear from 800 to 1500,
# straight from 1500 to 1800, steep linear decrease from 1800 to 1900):
OptimalProductLength = [ 600, 800, 1800, 2200 ]

OptimalProductLengthReward = 60

MinProductLength = 450

MaxProductLength = 3000

# it's still optimal if we go 2 below the optimal minimum length, if we have
# no mismatches:
OptimalPrimerLengthDispensationWithNoMismatches = 2

# ambiguities that are flanked by conserved regions of at least
# this size are not penalized further:
GoodConservedRegionLength = 2

MinPreIntronLength = 50 #one primer must have at least this distance to nearest intron

OptimalPreIntronLength = [ 70, 150 ]

MinDistanceToIntron = 5  # we're not entirely sure of intron start/ends

#used to assign column values to alignment parts w/ 2 seq and for rewarding conservation
ConservationWindow = 80

NonConservationPenalty = -5 # penalty for not residing in a conserved region

# at least this high conservation-% is necessary in alignment parts with only 2 seqs
# and where a primer might be contained entirely:
MinConservationPercent = 90.0

MafftPath = '/usr/bin/mafft'
ClustalPath = '/users/chili/BiRC/MolBio/Clustal/ClustalW/clustalw1.83/clustalw'
DNAMatrixPath = '%smultalinDNA.clustal'%Basedir
#DNAMatrixPath = '%smymatrix_identity5.clustal'%Basedir
WorkingDir = '/users/chili/BiRC/MolBio/Data/Comparisons5'
ParameterFile = '%s/parameters.prifi'%Basedir # the file in which parameter are saved

DocFilename = '%s/documentation_PriFi.txt'%Basedir


from .meltingtemperature import Tm
TM = Tm( PrimerConcentrationNmolar, SaltConcentrationMolar )




def trueLengthOfIntronWithMarkerLength( markerlength ):
    # Assumed intron code:
    #
    # XXX : <= 200 nucleotides
    # XXXX : 201 - 500 nucleotides
    # XXXXX : 501 - 1000 nucleotides
    # XXXXXX : > 1000 nucleotides
    if markerlength == 3:
        return 200
    elif markerlength == 4:
        return 500
    elif markerlength == 5:
        return 1000
    elif markerlength == 6:
        return 1500 # 2000 is too precautious
    else:

        return markerlength #raise ValueError, "Unexpected intron marker of length%d; set the 'Introns in sequences' parameter to 'no' to analyze this alignment."%markerlength





# ---------------------

# score matrix for identifying candidate primer regions
#

INF = -99999999
if MaxMismatches > 0:
    #
    # the penalties below are probably too restrictive in that they hinder
    # any region which starts with e.g. ****-, where * is perfect match and
    # - is a mismatch, no matter what comes after. That is not good since e.g.
    # ****-******************** might be better than just the last part starting
    # after the mismatch.
    #
    # p1 = MaxPrimerLength / MaxMismatches
    # p2 = ( MaxPrimerLength / MaxMismatchesWithMaxDiversity ) - 2
    #
    # so use these instead for region identification (note that bad primers
    # are reeped out later anyway):
    #
    p1 = 2  # at most every second column can be a mismatch with 2 diff. nucs.
    p2 = 2
    p3 = 2
    #
else:
    p1 = INF
    p2 = INF
    p3 = INF

# columns with more than one sequence represented and only one nucleotide
# (perfect alignment) are rewarded; others are penalized. Some with
# -infinity if this column should never be part of any primer, some with
# lesser numbers if they might participate in a primer site.

#total #nucs:     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23      in alignment column.
scorematrix = ( [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF], # 0 different nucleotides
                [INF, INF,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1], # 1 different nucleotides
                [INF, INF, -p3, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1, -p1], # 2 different nucleotides
                [INF, INF, INF, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2], # 3 different nucleotides
                [INF, INF, INF, INF, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2, -p2], # 4 different nucleotides
                [INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF]) # 5 different nucleotides

# scorematrix = ( [INF, INF, INF, INF, INF, INF], # 0 different nucleotides
#                 [INF, INF,   1,   1,   1,   1], # 1 different nucleotides
#                 [INF, INF, -p3, -p1, -p1, -p1], # 2 different nucleotides
#                 [INF, INF, INF, -p2, -p2, -p2], # 3 different nucleotides
#                 [INF, INF, INF, INF, -p2, -p2], # 4 different nucleotides
#                 [INF, INF, INF, INF, INF, INF]) # 5 different nucleotides

# (NB: the following note no longer holds, the entry has been corrected from
# INF to -p3)
# Note that it is disallowed to have 2 nucleotides in total that are different

# NB: entry (1, 0) - 0 nucleotides in total, 1 different, which of course
# is impossible - is used to hold the penalty for a column which has an X
# or N in it, i.e. is an intron in one of the sequences. Such a column may not
# be part of a primer region, so we need to be able to score it with -infinity.



# No alignment columns with too few nucleotides may be included, change any entries
# in the scorematrix that otherwise might allow that to happen:
for tuple in scorematrix:
    for j in range( 1, MinNucleotidesPerColumn ):# entry 0 assumed to be INF a priori
        tuple[j] = INF





#
#
# here's a list of the above names in a format that the program can use:
#
#

# each parameter is defined by a short explanatory text, the corresponding variable
# name, a long explanation, a minimal value and a maximal value (for singular values),
# and a default value:

Parameters = [ ('Similarity overlap', 'MatchOverlap', """If two single primers overlap with at least this many nucleotides, they are considered similar and regarded, for all practical purposes, as the same primer. This parameter therefore has a profound impact on which primers are reported. Consider two primer sets P and Q consisting of P.fw/P.rv and Q.fw/Q.rv. If both P.fw and Q.fw are similar and P.rv and Q.rv are similar, and Q has a higher score than P, then P is not reported (since it is considered to be more or less the same set as Q). Furthermore, if some primer p is similar to two other primers already reported, p (i.e., any set that includes p) is not reported (note that such sets would have lower scores than the sets already reported). Thus, the primer sets eventually reported are the highest scoring sets under these restrictions. This way suggesting several almost identical pairs is avoided.\nIf the Similarity Overlap parameter is increased, the reported primer sets most likely will be highly overlapping. If the parameter is decreased, the program is forced to suggest primers that don't overlap as much and hence are located in different regions. That leads to higher variation among the suggested primer sets but may end up including sets which have a much lower score compared to sets not reported. Care is needed when changing this parameter. """, 1, 100, 10),
               ('Tail length', 'TailLength', """The length of the 3'-end tail in which the number of A's and T's are counted.""", 0, MinPrimerLength, 8),
               #('Primer pair suggestions', 'PrimerPairSuggestions', """The maximum number (at most 4) of non-overlapping primer pairs to be suggested. "Non-overlapping" is defined so that two single parameters are considered to overlap if they overlap with at least Match Overlap nucleotides; if both forward and reverse primers in pair P overlap with the forward and reverse primers of a pair Q with higher score, P is not reported. Further, at most two overlapping single primers are reported. Thus, the suggested pairs are the highest scoring pairs under these restrictions. This way suggesting several almost identical pairs is avoided.""", 1, 4, 4),
               ('Maximum primer length', 'MaxPrimerLength', """Maximum primer length allowed.""", 10, 80, 35),
               ('Minimum primer length', 'MinPrimerLength', """Minimum primer length allowed.""", 10, 80, 18),
               ('Maximum melting temperature', 'SuggestedMaxTm', """Maximum primer melting temperature allowed.""", 20, 1000, 77.0),
               ('Minimum melting temperature', 'MinTm', """Minimum primer melting temperature allowed.""", 20, SuggestedMaxTm, 55.0),
               ('Critical melting temperature', 'CriticalTm', """If both primer melting temperatures are below this value, penalize the pair.""", 20, 1000, 58.0),
               ('Minimum melting temperature with ambiguity positions', 'MinTmWithMismatchesAllowed', """If a primer melting temperature is below this value, the primer can have no ambiguity positions.""", 20, 1000, 58.0),
               ("Minimum number of 3'-end matches", 'Min3endPerfectMatches', """A primer has to have at least this many perfect matches in its 3'-end.""", 0, MinPrimerLength, 2),
               ("Maximum number of ambiguity positions", 'MaxMismatches', """A primer may have at most this many ambiguity positions.""", 0, 100, 4),
               ("Minimum length with three or more ambiguity positions", 'MinLengthWithThreeAmbiguities', """A primer which is shorter than this value may have at most two ambiguity positions.""", 0, 100, 25),
               ( "Minimum length with two ambiguity positions", 'MinLengthWithTwoAmbiguities', """A primer which is shorter than this value may have at most one ambiguity position.""", 0, 100, 21),
               ( "Critical ambiguity position distance from 3'-end", 'CriticalAmbiguityDistanceTo3End', """Penalize ambiguity positions closer than this distance in nucleotides to the 3'-end.""", 0, 100, 5),
               ( "Maximum diversity per nucleotide position", 'MaxDiversityPerColumn', """A primer cannot cover a nucleotide position in the alignment with more than this many different nucleotides (gaps are ignored).""", 1, 4, 4 ),
               ( "Maximum number of ambiguity positions with maximum diversity", 'MaxMismatchesWithMaxDiversity', """A primer may cover at most this many nucleotide positions of the alignment with the maximum diversity allowed.""", 0, 100, 1 ),
               ( "Minimum number of nucleotides per nucleotide position", 'MinNucleotidesPerColumn', """In all nucleotide positions in the alignment covered by a primer, at least this many sequences must have a nucleotide. Thus, with the value 2 no primer may cover a part of the alignment in which only 1 sequence is represented (has a nucleotide).""", 1, 100, 2 ),
               ( "Primer concentration in nanomolar", 'PrimerConcentrationNmolar', """The concentration of the primers in the PCR reaction.""", 0, 10000000000000, 250 ),
               ( "Salt concentration in molar", 'SaltConcentrationMolar', """The salt concentration in the 1 x PCR buffer.""", 0, 10000000000000, 0.05 ),
               ( "Maximum melting temperature difference", 'MaxPrimerPairTmDifference', """Maximum difference allowed between the melting temperatures of the two primers in a pair.""", 0, 100, 15 ),
               ( "Optimal primer length interval", 'OptimalPrimerLength', """A primer length between these two values (both inclusive) is considered optimal.""", -1, -1, [25, 35] ),
               ( "Optimal PCR product length interval", 'OptimalProductLength', """A PCR product length between the middle two of these four values (both inclusive) is considered optimal. A length between the first two values is considered acceptable but less than optimal, and likewise with a length between the last two values. A length outside the full range is penalized. The product length is the length of the two primers plus the distance between them, measured in nucleotides.""", -1, -1, [600, 800, 1800, 2200] ),
               ( "Minimum PCR product length", 'MinProductLength', """Minimum length of the PCR product allowed.""", 0, 100000, 450 ),
               ( "Maximum PCR product length", 'MaxProductLength', """Maximum length of the PCR product allowed.""", 0, 100000, 3000 ),
               ( "Optimal primer length dispensation with no ambiguity positions", 'OptimalPrimerLengthDispensationWithNoMismatches', """If a primer has no ambiguity positions, its length can be this many nucleotides shorter than the otherwise smallest optimal length value and still be considered optimal.""", 0, 100, 2 ),
               ( "Good conserved region length", 'GoodConservedRegionLength', """Ambiguity positions in a primer preferably are surrounded by conserved regions without ambiguity positions. A region of at least this length is considered an acceptable conserved region; if an ambiguity positions is flanked by a conserved region of less than this length, it is penalized.""", 0, 100, 2 ),
               #( "XXXud med mig 1", 'MinPreIntronLength', """.""" ),
               #( "XXXud med mig 2", 'OptimalPreIntronLength', """.""" ),
               #( "XXXud med mig 3", 'MinDistanceToIntron', """.""" ),
               ( "Conservation window length", 'ConservationWindow', """Primers (partly) based on only two sequences must lie in a highly conserved region. The window length used as a basis for the conservation calculation has this length. See Minimum conservation %.""", 1, 10000, 80 ),
               ( "Minimum conservation-% in conservation window", 'MinConservationPercent', """When calculating the overall conservation in the conservation window, there must be at least this percentage of perfect matches (gaps are allowed) in the window.""", 0, 100, 90.0) ,
               ( "Path to local clustalw executable", 'ClustalPath', """If you have a local clustalw application, PriFi can do the alignment of your fasta file for you if you give the path to the executable file.""", -1, -1, "/users/chili/BiRC/MolBio/Clustal/ClustalW/clustalw1.83/clustalw" ),
               #( "XXXud med mig 4", 'DNAMatrixPath', """.""" ),
               ( "Path to working directory", 'WorkingDir', """The path to the default directory when loading alignment or Fasta files.""", -1, -1, "/users/chili/BiRC/MolBio/Data/Comparisons5" ),
               ( "Introns in sequences", 'INTRONS', """If this option is set to 'yes', X'es in the sequences are interpreted as special intron markers following this translation scheme:\nXXX : intron, <= 200 bp\nXXXX : 201 - 500 bp intron\nXXXXX : 501 - 1000 bp\nXXXXXX : > 1000 bp\nFurther, the primers in a primer set *must* have an intron between them to increase the probability of finding polymorphisms in the PCR product, and some other evaluation criteria apply when scoring the found primer sets.\nIf the option is set to 'no', X's are regular wild card symbols, and primer sets are not required to have an intron between them.""", -1, -1, "yes" )
               ]




def assertTextMatchesTupel( txt, tupel ):
    """Auxiliary function that splits the txt argument, converts the elements into floats and asserts that each float is within min and max values as defined in the tuple of pairs."""

    values = txt.strip()[1:-1].split(',')
    if len(values) != len(tupel):
        return 0
    prev=-999999
    fvalues = []
    try:
        for i in xrange(len(values)):
            v =  float(values[i])
            if not '.' in values[i]:
                v = int(v)

            fvalues.append(v)
            if v < tupel[i][0] or v > tupel[i][1] or v < prev:
                return 0
            prev = v # values must be increasing!
    except ValueError:
        return 0 # one of the values was not a number

    return fvalues # return new values



def parseStringAndAssignToParameter( v, m1, m2, txt ):
    """v is the parameter name for which txt is a value - i.e., txt is a string
    which needs to be parsed to get the value which is then assigned to the
    parameter with the name v. m1 and m2 are minimum and maximum allowed
    for the value (only applied if the value turns out to be a number)."""

    global WorkingDir, ClustalPath, OptimalProductLength, OptimalPrimerLength, INTRONS

    # three cases:

    # the tuple values:
    if v=='OptimalPrimerLength':
        w = assertTextMatchesTupel( txt, [(1, 150), (1, 150)] )
        if not w:
            return 0
        OptimalPrimerLength = w
    elif v=='OptimalProductLength':
        w = assertTextMatchesTupel( txt, [(20, 10000), (20, 10000), (20, 10000), (20, 10000) ] )
        if not w:
            return 0
        OptimalProductLength = w

      #elif v=='OptimalPreIntronLength':
      #   if not self.assertTextMatchesTupel( txt, [(1, 150), (1, 150)] ):
      #      break

    # text values:

    elif v=='INTRONS':
        txt = txt.lower()
        if txt != 'no' and txt != 'yes':
            return 0
        INTRONS = txt


    elif v=='ClustalPath':
        if not exists( txt ):
            return 0
        ClustalPath = txt

    elif v=='WorkingDir':
        if not exists( txt ):
            return 0
        WorkingDir = txt


    # single values:
    else:
        try:
            value = float(txt)
        except ValueError:
            return 0
        if value < m1 or value > m2:
            return 0
        if not '.' in txt:
            value = int(value)

        # set the variable to this value

        #print 'before/after parsing: ',thismodule.__dict__[v],
        thismodule.__dict__[v] = value
        #print "setting something:",thismodule.__dict__[v]


    return 1 # it went well
