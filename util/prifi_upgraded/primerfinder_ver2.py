#! /usr/bin/env python


#
# PRIFI
#




## May/June 2005: This criterion discarded:
#      if one primer is based on 2 sequences only, the other has
#      to be based on at least 3.
## Apr 21: introduced the explain argument in scoreprimerpair(), saves memory = time!
## Apr 22: keep left/right info in primer tupel
## May 5th: introduced a primer class instead of the horribly long primer tuples
## May 6th: introduced pruning step in the primer candidate finder; speedup factor = 500-1000!
##

from __future__ import division, print_function
import time
from alignment import *
import config as cf
import sys
from math import log10, atan
from reversecomplement import reverse_and_complement


# needed for support of Python versions 2 _and_ 3
try:
    xrange
except NameError:
    xrange = range

# ------ auxiliary definitions --------------------------------------


class primer:
    """represents an individual primer"""

    def __init__( self, s, e, i, j, part, Tm, mm, nearinl, nearinr, sia,ATr, ATl,ts):
        # score used to be index 0
        self.regionstart = s
        self.regionend = e
        self.start = i
        self.end = j
        self.seq = part                    # sequence of primer (copy of one of the sequences)
        self.tm = Tm                       # melting temperature
        # tmr used to be index 7
        self.ambiguities = mm              # list of primer alignment mismatch positions
        self.nearestLeftIntron = nearinl   # distance to nearest intron on either side
        self.nearestRightIntron = nearinr
        self.avSeqInAlignment = sia        # average number of sequences in primer's alignment
        self.ATsInLeftTail = ATl           # number of As and Ts in 3'-end tail
        self.ATsInRightTail = ATr
        self.twoseqs = ts # number of bases in the primer based on only two seqs.

        # different contributors to the total score, used for explanation:
        self.siascore = 0
        self.tmscore = 0
        self.gcscore = 0
        self.score = 0
        self.lengthscore = 0
        self.lengthexplanation = None  # a text explanation for the contributed length score
        self.highdivscore = 0
        self.mmclusterscore = 0
        self.mmscore = 0
        self.scscore = 0 # self complementarity
        self.consscore = 0

        # the following are not contributors to the total score and have to be added separately
        # (using the helper sum variables below):
        self.fwpreintronscore = 0 # preintronlength score if p serves as forward primer
        self.rvpreintronscore = 0 # preintronlength score if p serves as reverse primer
        self.fwATscore = 0       # penalty for too many As and Ts in 3'-end tail if p fw primer
        self.rvATscore = 0
        self.fw3GCscore = 0
        self.rv3GCscore = 0
        self.fw3endmmscore = 0
        self.rv3endmmscore = 0
        # - the above 8 scores are added (fw scores with fw scores, rv scores with rv score)
        # and kept in these extra scores:
        self.fwscore = 0
        self.rvscore = 0






def doAlignment( file, verbose=0 ):
    """Perform clustalw alignment from the given file. The verbose parameter determines whether the command line given to clustalw is printed."""


    allseq, summary, l = create_alignment( file, verbose )

    return allseq, summary, l




def print_N_best_primers( primers, lastX ):
    print('\nBest %d primers overall:'%lastX)
    qq=len(primers) - lastX
    for score, s, e, i, j, primer, tmf, tmr, mm in primers[-lastX:]:
        print("\n(%d):\n[%d, %d[: (%d, %d). Score %d, Tmf %d, Tmr %d, length %d, #mism %d"%(len(primers)-qq, s, e, i, j, score, tmf, tmr, j-i, mm))
        printslice( allseq, i, j )
        qq+=1




def print_best_primer_in_each_region( regionprimerlist ):
    for rpl in xrange( len(regionprimerlists) ):
        if len( regionprimerlists[rpl])>0:
            print('\nBest primer in region %d:'%rpl)

        qq=0
        for score, s, e, i, j, primer, tmf, tmr, mm in regionprimerlists[rpl][-1:]:
            print("\n[%d, %d[: (%d, %d). Score %d, Tmf %d, Tmr %d, length %d, #mism %d"%(s, e, i, j, score, tmf, tmr, j-i, mm))
            printslice( allseq, i, j )
            qq+=1


def findPrimerOverlap( fw1, fw2 ):
    """find pairwise overlap in basepairs (or 0) between the two given primers (parameter
    names not clever, they don't have to be forward primers)."""

    fwoverlap = fw1.start < fw2.end-1 and fw2.start < fw1.end-1

    if fwoverlap:
        fwoverlap = min( fw1.end, fw2.end ) - max( fw1.start, fw2.start )

    return fwoverlap



def findOverlap(a, b, c, d, e, f, g, h):
    """find pairwise overlap in basepairs (or 0) for fw and rv primers of the two pairs given (first pair is indices a, b and c, d; second pair is e, f and g, h. (Ending indices are not inclusive). I.e., return overlap of forward primers, and overlap of reverse primers."""

    fwoverlap = a < f-1 and e < b-1
    rvoverlap = c < h-1 and g < d-1

    if fwoverlap:
        fwoverlap = min( b, f ) - max( a, e )

    if rvoverlap:
        rvoverlap = min( d, h ) - max( c, g )


    return fwoverlap, rvoverlap



def factorial( n ):
    if n == 0:
        return 1

    f = 1
    for i in xrange(1, 1+n):
        f *= i

    return f



def ambiguityCode( col ):
    """col is a string holding a mismatch column of the alignment. Based on the nucleotides in this column, return a one-letter ambiguity code."""
    a = 'A' in col
    c = 'C' in col
    g = 'G' in col
    t = 'T' in col
    x = 'X' in col
    n = 'N' in col

    if a and c and not g and not t:
        return 'M'
    if a and g and not c and not t:
        return 'R'
    if a and t and not c and not g:
        return 'W'
    if c and g and not a and not t:
        return 'S'
    if c and t and not a and not g:
        return 'Y'
    if g and t and not c and not a:
        return 'K'
    if a and c and g and not t:
        return 'V'
    if a and c and t and not g:
        return 'H'
    if a and g and t and not c:
        return 'D'
    if c and g and t and not a:
        return 'B'
    if x:  # intron column with gaps in all other sequences
        return 'intron-in-primer-error'
    if n or (a and c and g and t): # 'N' is also used if all 4 nucleotides are observed
        return 'N'

    return 'no-mismatch-error' # no mismatch in this column OR all gaps




def insertAmbiguities( sq, s, e, dir, sm, colsum ):
    """returns the sequence sq with inserted ambiguities (looked up in sm and colsum - i, j are the start, end indices of the string in the alignment), in rev.compl. form if dir==1"""

    seq = sq
    for i in range( s, e ):
        if colsum[i][0] > 1:
            # here's a mismatch
            # deprecated:
            #seq = ambiguityCode(sm.get_column(i)).join( [seq[:i-s], seq[i+1-s:]] )
            seq = ambiguityCode(sm.alignment[:,i]).join( [seq[:i-s], seq[i+1-s:]] )

    if dir:
        return reverse_and_complement(seq)
    else:
        return seq


def primer2string( p, dir, sm, colsum ):
    """p is a primer, dir is the direction, 0 is forward, 1 is reverse. sm is a
    summary of the alignment, colsum is a list of column summaries.
    We return the primer in string format WITH AMBIGUITY CODES INSERTED.
    A primer is a tuple on this form:
        (score, s, e, i, j, part, Tmf, Tmr, mm )
    where s, e are the regional start, end indices in the alignment, i, j are
    the start, end indices of the primer itself (e and j excluded), part is the actual
    string, Tmf is melting temperature in the forward direction, Tmr is the
    melting temp. in the reverse direction, mm is a list of ambiguity indices."""

    # seq = p.seq

    # now substitute ambiguity codes:
    # for i in range( p.start, p.end ):
    #    if colsum[i][0] > 1:
    #        # here's a mismatch
    #        seq = ambiguityCode(sm.get_column(i)).join( [seq[:i-p.start], seq[i+1-p.start:]] )


    seq = insertAmbiguities( p.seq, p.start, p.end, dir, sm, colsum)

    if dir:
        return 'rv: [%4d,%4d[, Tm=%.1f, #amb=%d: %s (not r.c.)'%(p.start, p.end, p.tm, len(p.ambiguities), reverse_and_complement(seq))+'\n                                  %s (r.c.)'%(seq)
    else:
        return 'fw: [%4d,%4d[, Tm=%.1f, #amb=%d: %s'%(p.start, p.end, p.tm, len(p.ambiguities), seq)

# ----------------------------------------------------------------------





def IntronsBetweenRegions( intronindices, a, b ):
    """intronindices is a list of introns (start,end indices), a and b are
    start indices of two primer regions."""

    intronsbetweenprimers = [] # list of all inter-region introns
    x6 = 0 # #XXXXXX intron markers (indicating length>1000)

    for ii in intronindices:
        if ii[0] > b:
            break  # all remaining introns are to the right of a and b
        if ii[0] > a and ii[0] < b:
            intronsbetweenprimers.append( ii )
            if ii[1] - ii[0] + 1 == 6:
                x6 += 1

    return (intronsbetweenprimers, x6)






def countATGCs( seq ):
    a=c=g=t=0
    for i in seq:
        if i == 'A':
            a += 1
        elif i == 'C':
            c += 1
        elif i == 'G':
            g += 1
        elif i == 'T':
            t += 1
    return a, t, g, c






def scoreprimerpair( p1, p2, realindices, intronsbetweenprimers, explain=0 ):
    """p1 is the forward primer, p2 is the reverse primer (needs to be reverse complemented). Thus p1 is assumed to reside to the left of p2 in the alignment.
    Return 3-tupel of (score, explanation for rewards, explanation for penalties) where the last two are strings and the first is a number. If the explain argument is 0, the explanations are empty strings."""

    # A primer is a tuple on this form:
    #   (score, s, e, i, j, part, Tmf, Tmr, mm )
    # where s, e are the regional start, end indices in the alignment (e excluded),
    # i, j are the start, end indices of the primer itself, part is the actual
    # string, Tmf is melting temperature in forward direction, Tmr in reverse
    # direction; mm is the number of mismatches.

    mmf = len(p1.ambiguities) # mismatches in forward primer
    mmr = len(p2.ambiguities) # mismatches in reverse primer



    # NB: this rule-of-thumb has been cut out:
    # cut: ------>

    # Eg.: Tm = 64 with three mismatches not okay but okay with 2.
    #if p1.tm - MismatchTmDecrease * mmf < MinTm or \
    #       p2.tm - MismatchTmDecrease * mmr < MinTm:
    #    return (-9, '', '')

    # <-------








    # cut: ------->


    # create a list of index markers with true values (substituting the
    # intron markers XXX etc. with their actual lengths)
    #totalintronl = 0 # total intron length so far

    #exonstart =  realindices[p1.end] # here starts the first exon
    #exonindices = [ ]  # real start, end (end not included) indices of exons
    #for ii in intronsbetweenprimers:
    #    exonindices.append( [exonstart, realindices[ii[0]] ] )

    #    # find real intron length and add to total:
    #    totalintronl += trueLengthOfIntronWithMarkerLength( ii[1]+1-ii[0] )
    #    # subtract marker length:
    #    totalintronl -= ( ii[1]+1 - ii[0] )

    #    exonstart = realindices[ii[1]] + 1
    # don't forget last exon:
    #exonindices.append( [ exonstart, realindices[p2.start] ] )

    # <--------



    #productlength = totalintronl + p2.end - p1.start

    # results in an error:
    #productlength = realindices[p2.end] - realindices[p1.start]
    productlength = realindices[p2.end-1] - realindices[p1.start]


    # product length okay?
    if productlength < cf.MinProductLength or productlength > cf.MaxProductLength:
        return (-12, '', '')





    # FORSLAG: pil de par fra, som har samlet preintrondist < 100 OG som ikke
    # kan naa en ID-exon.
    # Ikke implementeret lige nu.


    # cut: ------>

    # search for long enough exon between the introns:



    # if there is a contiguous exon of at least a certain length between primers
    # and within a certain distance from the closest primer, reward it:
    #foundIDexon = 0
    #IDexoncloseenough = 0
    # don't check first and last exon, the ID-exon must be another one:
    #if intronsbetweenprimers[0][0] - p1.end >= MinLengthOfLongestInterPrimerExon or \
    #       p2.start-intronsbetweenprimers[-1][1]-1 >=MinLengthOfLongestInterPrimerExon:
    #    foundIDexon = 1
    #    IDexoncloseenough = 1



    #distIDexonstring = "None found!"

    #for estart, eend in exonindices:
    #    if eend - estart >= MinLengthOfLongestInterPrimerExon:
    #        foundIDexon = 1
    #        # has to be 'within reach' of primers too (NB: first exon doesn't
    #        # count if it is only reachable from the fw primer, last exon
    #        # doesn't count if it is only reachable from the rv primer):

    #        lll = estart - realindices[p1.end]
    #        rrr = realindices[p2.start] - eend

    #        if (0 < lll <= MaxPrimerIDexonDistance) or (0 < rrr <= MaxPrimerIDexonDistance):
    #            IDexoncloseenough = 1
    #            #if 0 < lll and ( lll <= rrr or rrr == 0 ):
    #            #    distIDexonstring = "%d (fw)"%lll
    #            #else:
    #            #    distIDexonstring = "%d (rv)"%rrr
    #            break



    # <------------- cut






    # --------------------------------------------------------------------

    # OKAY, this primer pair should be considered, now rate it.




    pro = []      # explanation for rewards
    contra = []   # explanation for penalties
    score = 1000.0 # must stay above 0 for all primer pairs reaching this point
    score = 0.0 # must stay above 0 for all primer pairs reaching this point




    # add scores of individual primers:

    score += p1.score
    score += p2.score

    score += p1.fwscore + p2.rvscore

    preintronscore = p1.fwpreintronscore + p2.rvpreintronscore
    # score += preintronscore

    ATscore = p1.fwATscore + p2.rvATscore
    # score += ATscore

    # score += p1.fw3GCscore + p2.rv3GCscore

    endmmscore = p1.fw3endmmscore + p2.rv3endmmscore
    # score += endmmscore




    p1len = p1.end - p1.start
    p2len = p2.end - p2.start


    if explain:
        pro.append( 'Tm = %4.1f / %4.1f'%(p1.tm, p2.tm) )
        pro.append( 'Primer lengths: %d / %d'%(p1len, p2len) )
        pro.append( 'Avg. #sequences in primer alignments: %.1f / %.1f'%(p1.avSeqInAlignment, p2.avSeqInAlignment))
        if cf.INTRONS == 'yes':
            pro.append( 'Estimated product length: %d'%productlength)
            pro.append( 'Primer/intron distances: %d / %d'%(p1.nearestRightIntron,p2.nearestLeftIntron))
        else:
            pro.append( 'Product length: %d'%productlength)
        # pro.append( 'Estimated distance to closest ID-exon: %s'%distIDexonstring )
        pro.append("A/T's among last %d bp of 3'-end: %d / %d"%(cf.TailLength,p1.ATsInRightTail,p2.ATsInLeftTail))
        pro.append( 'Ambiguities: %d / %d'%(mmf, mmr) )

        if mmf > 0:
            ambs = ["Fw ambiguity positions:"]
            for q in p1.ambiguities:
                ambs.append('%d'%(q-p1.start+1))
            pro.append( ' '.join(ambs) )
        if mmr > 0:
            ambs = ["Rv ambiguity positions:"]
            for i in range(len(p2.ambiguities)-1, -1, -1):
                q = p2.ambiguities[i]
                ambs.append('%d'%(p2.end-q))
            pro.append( ' '.join(ambs) )


    pro.append('')




    # add explanations for individual primers:
    if explain:
        P = p1.tmscore + p2.tmscore
        if P!=0:
            pro.append('%5.1f: High-Tm bonus'%P)

        if p1.lengthscore > 0:
            pro.append(p1.lengthexplanation%'Fw')
        if p2.lengthscore > 0:
            pro.append(p2.lengthexplanation%'Rv')

        P = p1.siascore + p2.siascore
        if P!=0:
            pro.append( "%5.1f: bonus for #sequences in primer alignments"%P)

        if p1.gcscore != 0:
            contra.append("%5.1f: Fw GC content >= 75%"%p1.gcscore )

        if p2.gcscore != 0:
            contra.append("%5.1f: Rv GC content >= 75%"%p2.gcscore )

        if p1.scscore != 0:
            contra.append("%5.1f: Fw self-complementarity"%p1.scscore )

        if p2.scscore != 0:
            contra.append("%5.1f: Rv self-complementarity"%p2.scscore )

        if p1.consscore != 0:
            contra.append( "%5.1f: Fw in unconserved region or based mostly on 2 seqs"%p1.consscore )
        if p2.consscore != 0:
            contra.append( "%5.1f: Rv in unconserved region or based mostly on 2 seqs"%p2.consscore )


        if cf.INTRONS == 'yes':
            if preintronscore != 0:
                contra.append("%5.1f: Primer/intron distance(s) outside %d-%d bp"%(preintronscore, cf.OptimalPreIntronLength[0], cf.OptimalPreIntronLength[1]))

        if ATscore != 0:
            contra.append("%5.1f: Too high AT content in 3'-ends"%ATscore)

        if p1.lengthscore < 0:
            contra.append(p1.lengthexplanation%'Fw')

        if p2.lengthscore < 0:
            contra.append(p2.lengthexplanation%'Rv')


        if p1.fw3GCscore > 0:
            pro.append("%5.1f: Fw has G/C terminal in 3'-end"%p1.fw3GCscore)
        if p2.rv3GCscore > 0:
            pro.append("%5.1f: Rv has G/C terminal in 3'-end"%p2.rv3GCscore)


        if p1.fw3endmmscore != 0:
            contra.append("%5.1f: Fw has ambiguities near 3'-end"%p1.fw3endmmscore)

        if p2.rv3endmmscore != 0:
            contra.append("%5.1f: Rv has ambiguities near 3'-end"%p2.rv3endmmscore)

        if p1.highdivscore != 0:
            contra.append("%5.1f: Fw has high diversity in ambiguities"%p1.highdivscore)
        if p2.highdivscore != 0:
            contra.append("%5.1f: Rv has high diversity in ambiguities"%p2.highdivscore)

        if p1.mmclusterscore != 0:
            contra.append("%5.1f: Fw ambiguities flanked by short conserved regions"%p1.mmclusterscore)

        if p2.mmclusterscore != 0:
            contra.append("%5.1f: Rv ambiguities flanked by short conserved regions"%p2.mmclusterscore)



    # cut:------------->

    # Is there a reachable ID exon?

    #if foundIDexon and IDexoncloseenough:
    #    score += IDexonBonus
    #    pro.append("%5.1f: reachable ID exon"%IDexonBonus)

    # <-------------- cut




    # penalize complementarity between fw and rv primers (so far,
    # only self-complementarity between a primer and itself has been
    # tested for):

    # ** currently not checked **





    # both primers can't have too low Tm's:

    if p1.tm < cf.CriticalTm and p2.tm < cf.CriticalTm:
        pen = -( 9 * ( 2*cf.CriticalTm - p1.tm - p2.tm ) )
        if mmf == mmr == 0:
            pen /= 2        # if no mismatches, penalize less
        score += (pen)
        if explain:
            contra.append("%5.1f: Both primers have low Tm"%pen)



    # penalize too short and too long products:

    if productlength < cf.OptimalProductLength[0]:
        pen = (-(cf.OptimalProductLength[0] - productlength)/7) # var /2
        if explain:
            contra.append("%5.1f: Very short estimated product"%pen)
        score += pen

    elif productlength > cf.OptimalProductLength[3]:
        pen = (-(productlength - cf.OptimalProductLength[3])/10)
        if explain:
            contra.append("%5.1f: Very long estimated product"%pen)
        score += pen

    else:         # reward longer products within the limits:
        pen = cf.OptimalProductLengthReward

        if productlength <= cf.OptimalProductLength[1]:
            diff = float(cf.OptimalProductLength[1] - cf.OptimalProductLength[0])
            pen *= ((productlength-cf.OptimalProductLength[0])/diff)

        elif productlength >= cf.OptimalProductLength[2]:
            diff = float(cf.OptimalProductLength[3] - cf.OptimalProductLength[2])
            pen *= ((cf.OptimalProductLength[3]-productlength)/diff)

        if explain and pen != 0:
            pro.append("%5.1f: Good product length"%pen)
        score += pen



    if explain:
        mmscore = p1.mmscore + p2.mmscore
        if mmscore != 0:
            contra.append("%5.1f: %d ambiguities"%(mmscore, len(p1.ambiguities)+len(p2.ambiguities)) )


    # finalize report strings:

    if len(contra) > 0:
        contra.append('')
        contrastring = '\n'.join( contra )
    else:
        contrastring = ""

    if len(pro) > 0:
        pro.append('')
        prostring = '\n'.join( pro )
    else:
        prostring = ''



    return ( ( score ), prostring, contrastring )














def scoreIndividualPrimer( p, colsum, conservation, explain=0 ):
    """p is a primer, give it its individual score which is used when scoring primer pairs.
    If explain==1, add a textual explanation for the score to the primer."""


    #global cf.MinTm
    #global cf.OptimalPrimerLength
    #global cf.OptimalPreIntronLength
    #global cf.MinPreIntronLength

    # how many sequences in the primer alignments?
    # p1.avSeqInAlignment is the avg. number of seqs in primer alignment)

    p.score = 0
    plen = (p.end - p.start)


    # reward primers based on 3 or 4 sequences; don't reward alignments with only 2:


    pen = atan((p.avSeqInAlignment-2.9)*10)*11+16.06246 # sigmoid function
    #pen = (p.avSeqInAlignment - 2) * 10 # was 30
    p.score += pen

    if explain and abs(pen) > .1:
        p.siascore = pen




    # high melting temperatures better than low:
    pen = (p.tm - cf.MinTm) * 4
    p.score += pen

    if explain and pen != 0:
        p.tmscore = pen



    # GC content shouldn't be above 75%:

    w, x, y, z = countATGCs( p.seq )
    GC = (y+z)/plen

    # penalize primers with more than 75% GC and which are short also:
    # - BUT DON'T REWARD THEM OTHERWISE

    #Pp = 10
    if GC < .75:
        #if explain:
        #    pro.append("%5.1f: Fw GC content < 75%s" %(Pp, '%') )
        #score += pen
        pass
    else:
        pen = -abs( plen - cf.OptimalPrimerLength[1] )
        #pen = -abs( plen - OptimalPrimerLength[1] )
        p.score += pen
        if explain and pen != 0:
            p.gcscore = pen



    # complexity:

    # prod = factorial(w) * factorial(x) * factorial(y) * factorial(z)
    # K = (1/L) * log10( factorial(L) / prod ) / log10(L)
    # score += K*3


    # GC content:
    # outside40_60percent = ( abs( (y+z)/L - 0.5 ) - 0.1 ) # 0.1-0.4
    # if outside40_60percent > .1:
    #   score -= outside40_60percent * 30






    # score the pre-intron length in both directions:

    if cf.INTRONS == 'yes':
        if p.nearestRightIntron < cf.MinPreIntronLength:
            p.fwpreintronscore = -(cf.OptimalPreIntronLength[0] - p.nearestRightIntron) / 3.0

        elif p.nearestRightIntron < cf.OptimalPreIntronLength[0]:
            p.fwpreintronscore = -(cf.OptimalPreIntronLength[0] - p.nearestRightIntron) / 6.0

        elif p.nearestRightIntron > cf.OptimalPreIntronLength[1]:
            p.fwpreintronscore = -(p.nearestRightIntron - cf.OptimalPreIntronLength[1]) / 10.0


        if p.nearestLeftIntron < cf.MinPreIntronLength:
            p.rvpreintronscore = -(cf.OptimalPreIntronLength[0] - p.nearestLeftIntron) / 3.0

        elif p.nearestLeftIntron < cf.OptimalPreIntronLength[0]:
            p.rvpreintronscore = -(cf.OptimalPreIntronLength[0] - p.nearestLeftIntron) / 6.0

        elif p.nearestLeftIntron > cf.OptimalPreIntronLength[1]:
            p.rvpreintronscore = -(p.nearestLeftIntron - cf.OptimalPreIntronLength[1]) / 10.0

        # don't add to the score, only one of them is used depending on whether this
        # primer will be the fw or rv primer in a pair.




    # penalize AT content in 3'-end:

    t = cf.TailLength/2.0
    t1 = p.ATsInRightTail - t # only penalize if AT fraction is more than 50% of 3' tail.
    t2 = p.ATsInLeftTail - t

    if t1>0:
        p.fwATscore = -3 * t1*t1  # at most factor*9

    if t2>0:
        p.rvATscore = -3 * t2 * t2




    lenmm = len(p.ambiguities)

    # penalize primers with less than optimal length:

    optmax = cf.OptimalPrimerLength[1]
    optmin = cf.OptimalPrimerLength[0]
    if lenmm == 0:
        optmin -= cf.OptimalPrimerLengthDispensationWithNoMismatches

    HigherThanOptimum = 1
    LowerThanOptimum = 1

    if plen < optmin:
        pen = -( optmin - plen ) * LowerThanOptimum
        p.score += pen
        p.lengthscore = pen
        if explain and pen != 0:
            p.lengthexplanation = '%5.1f:'%pen + ' %s shorter than ' + str(optmin)

    elif plen > optmax:
        pen = -( plen - optmax ) * HigherThanOptimum
        p.score += pen
        p.lengthscore = pen
        if explain and pen != 0:
            p.lengthexplanation = '%5.1f:'%pen + ' %s longer than ' + str(optmax)

    else:  # reward longer primers a little (within limits):
        pen = (plen - optmin)/2.0
        if explain and pen != 0:
            p.lengthexplanation = '%5.1f:'%pen + ' %s primer length'
        p.lengthscore = pen
        p.score += pen






    # score G/C terminal (both orientations)

    P = 3.0
    if p.seq[-1] == 'G' or p.seq[-1] == 'C':
        p.fw3GCscore = P

    if p.seq[0] == 'C' or p.seq[0] == 'G':
        p.rv3GCscore = P






    # penalize ambiguities near 3'-end (both orientations):

    highdiversitypenalty = 0
    # find distance from each mismatch to the 3'-end of the primer:
    p.fw3endmmscore = p.rv3endmmscore = 0

    for i in p.ambiguities:
        # don't diversity-penalize mismatch columns with only two different
        # nucleotides:
        highdiversitypenalty -= (colsum[i][0]-2)*10

        # don't penalize ambiguities far from 3'-end:
        if p.end - i <= cf.CriticalAmbiguityDistanceTo3End:
            # p.fw3endmmscore += -( 30.0 / ( p.end - i ) + 5 )
            p.fw3endmmscore += -5

        if i - p.start <= cf.CriticalAmbiguityDistanceTo3End:
            # p.rv3endmmscore += -( 30.0 / ( i - p.start + 1 ) + 5 )
            p.rv3endmmscore += -5



    # penalize diversity in mismatch columns:
    if highdiversitypenalty != 0:
        p.score += highdiversitypenalty
        p.highdivscore = highdiversitypenalty





    # penalize 'clustering' of mismatches or mismatches near the ends
    # (should give 0 if mismatches are spread out as nicely as possible):

    # ideal number of nucleotides between mismatches (and ends):
    #     ideal = min( 7, (p1len - mmf) / (mmf + 1) )
    pp = p.start   # initialization (pp points to position of last ambiguity plus one)
    pen = 0
    factor = 6.0
    for i in p.ambiguities:  # go through list of ambiguity positions
        pen -=max( 0, ( cf.GoodConservedRegionLength - ( i - pp ) ) * factor )
        pp = i+1
    pen -= max( 0, ( cf.GoodConservedRegionLength - ( p.end - pp ) ) * factor )
    if pen != 0:
        p.score += pen
        p.mmclusterscore = pen



    # penalize primers not in highly conserved areas OR IN highly conserved
    # ares but based mostly on 2 seqs:

    if 0 in conservation[p.start:p.end] or p.twoseqs >= .67*plen:
        p.score += cf.NonConservationPenalty
        p.consscore = cf.NonConservationPenalty





    # penalize ambiguities (Max counts double):
    # (perhaps the penalty should depend on the Tm?)
    pen = -11 * lenmm

    if lenmm == cf.MaxMismatches:
        pen -= 11



    if pen < 0:
        p.score += pen
        p.mmscore = pen




    # don't check for self-complementarity, it seems to have only minor impact:

    # what we should do here is penalize IF:
    #  - # bindings is high (at least 8-10)
    #  - # not too many non-matches of the overlapping dimer
    # This would ensure that the erroneous product would survive.


    # penalize self-complementarity:
    #pen = 0
    #for f in xrange( plen ):
    #    cm = 0 # count complementary matches
    #    for g in xrange( 0, plen-f ):
    #        c1 = p.seq[f+g].lower()
    #        c2 = p.seq[plen-1-g].lower()
    #        #print '\n%s %s'%(c1, c2),
    #        if c1 == 'a' and c2 == 't' or \
    #           c1 == 't' and c2 == 'a' or \
    #           c1 == 'c' and c2 == 'g' or \
    #           c1 == 'g' and c2 == 'c':
    #            # self-pairing here
    #            cm += 1
    #            #print '#',
    #    #print '\nsc: %d\n'%(cm)
    #    if cm > pen: # remember the highest number of matching columns for any dimer
    #        pen = cm

    #pen *= (-0.5) # cf.SelfComplementarityFactor?
    #p.score += pen
    #p.scscore = pen


    # now add up the two extra scores for use if this primer will be a fw or
    # a rv primer, respectively:

    p.fwscore = p.fwpreintronscore + p.fwATscore + p.fw3GCscore + p.fw3endmmscore
    p.rvscore = p.rvpreintronscore + p.rvATscore + p.rv3GCscore + p.rv3endmmscore










# def scoreprimer( primer, colsum, i, j, mm ):
#     """score this primer (alignment columns [i, j[) using the column summaries. The higher the score, the better.


#     NB: not used at the moment!


#     """

#     L = j-i
#     score = 100.0

#     # length   ~ around 14 - 35:
#     score += L

#     # mismatches ~ 0-Max
#     score -= mm*12


#     # how many sequences in the alignment?
#     t = 0.0
#     for q in xrange(i, j):
#         t += colsum[q][1]*colsum[q][1]*3

#     # avg. #nucleotides per alignment column:
#     score += (t/L )


#     # complexity:
#     w, x, y, z = countATGCs( primer )

#     prod = factorial(w) * factorial(x) * factorial(y) * factorial(z)
#     K = (1/L) * log10( factorial(L) / prod ) / log10(L)

#     score += K*3


#     # GC content:
#     outside40_60percent = ( abs( (y+z)/L - 0.5 ) - 0.1 ) # 0.1-0.4
#     if outside40_60percent > .1:
#         score -= outside40_60percent * 30




#     return score




def has_x_perfect_end_matches( x, end, colsum, i, j):
    """each primer has to have at least Min3endPerfectMatches perfect alignment matches at the 3'-end. This function checks if there are at least x perfect matches. Argument 'end' determines which end is checked: 0 is left, 1 is right. Return 0 if that end has enough contiguous perfect matches, 0 otherwise."""

    # colsum[i][0] is the number of different nucleotides in
    # column i of the alignment.

    if end:
        r = 1 # assume all is well
        for q in xrange( x ):
            if colsum[j-1-q][0] > 1:
                r = 0
                break
        return r
    else:
        l = 1
        for q in xrange( x ):
            if colsum[i+q][0] > 1:
                l = 0
                break
        return l







def too_high_diversity_in_mismatches( colsum, i, j ):
    """do the mismatches have too high a diversity?"""

    #global cf.MaxDiversityPerColumn, cf.MaxMismatchesWithMaxDiversity

    colswithmaxdiversity = 0
    colswithmt2n = 0 # more than 2, less than max
    mm = 0

    for q in xrange(i, j):

        if colsum[q][0] > 1:
            # mismatch!

            mm += 1
            if colsum[q][0] > cf.MaxDiversityPerColumn:
                return 1
            elif colsum[q][0] == cf.MaxDiversityPerColumn:
                colswithmaxdiversity += 1
            elif colsum[q][0] > 2:
                colswithmt2n += 1

    if colswithmaxdiversity > cf.MaxMismatchesWithMaxDiversity:
        return 1

    # if there is the maximally allowed number of columns with max diversity,
    # there can only be one other ambiguity, and the other can only have two
    # different nucleotides:
    if colswithmaxdiversity == cf.MaxMismatchesWithMaxDiversity:
        if mm > 2 or colswithmt2n > 0:
            return 1

    return 0



def too_many_mismatches( colsum, i, j, tm ):
    """too many mismatches in this primer? if melting temperature is close to the minimal value, we need to have 0 mismatches.
    Return (a, b) where a is 0 or 1 (false or true as to the question 'too many mismatches..') and b is the number of mismatches."""

    #global cf.MaxMismatches
    #global cf.MinTmWithMismatchesAllowed

    mm = []       # list of mismatch column indices
    # twoseqs = 1   # yes, we think we have only 2 seqs in this alignment part

    for q in xrange(i, j):

        if colsum[q][0] > 1:
            #if q>i:
            #    if colsum[q-1][0] > 1:
            #        return (1, []) # two consecutive mismatches
            mm.append( q )

        # if colsum[q][1] > 2:
        #    twoseqs = 0


    if len(mm) > cf.MaxMismatches:
        return (1, mm)

    if tm < cf.MinTmWithMismatchesAllowed and (len(mm) > 0): # or twoseqs):

        return (1, mm)

    return (0, mm)

    # note that the mismatch list is sorted, smallest index first





def find_primer_regions_old( a ):
    """a is a list of columnscores"""

    # with column scores 1, -2, -infinity, this function delimits regions
    # when it no longer 'pays off' to include a column in the current
    # region. I use the Maximum Partial Sum algorithm, sort of: add column
    # score, if sum falls below 0, start new region and reinitialize sum to 0.
    # Unfortunately, the same list reversed might give different results since
    # we sum from the left. See following function.
    maxendinghere = 0
    # indices for primer regions
    indices = []
    start = 0
    end = 0
    for i in xrange(len(a)):
        maxendinghere = maxendinghere+a[i]
        #print maxendinghere,
        if maxendinghere < 0:
            maxendinghere = 0
            if i-start >= cf.MinPrimerLength:
                indices.append( (start, i) )
            start = i+1 # begin new block

    # remember last bit if necessary:
    if i-start >= cf.MinPrimerLength:
        indices.append( (start, i) )

    # return a list of tuples: (start, end) of interval, where
    # end is NOT included.
    return indices



def find_primer_regions( a ):
   """a is a list of columnscores. Go through a, locate negative entries (not
   -infinity), check whether a window of (MinPrimerLength) can be placed
   around each which holds in total at most (MaxMismatches) mismatches. If not,
   this entry/position cannot be part of any primer region."""


   # The list 'a' holds these values: cf.INF, -cf.p1, -cf.p2, -cf.p3 and 1.

   i = 0
   la = len(a)
   pr = [True]*la
   for ii in range(la):
       if a[ii] == cf.INF:
           pr[ii] = False
   # pr is a list of True's with False's where the 'a' list has INF values.

   a.append(cf.INF) # add 'stop value' to end of list
   #printslice(allseq, 69, 74)

   # locate mismatch positions, check if each resides in window
   # of at least MinPrimerLength with no INFs and has at most MaxMismatches :
   while 1:

       while i < la and (a[i] == cf.INF or a[i] == 1):
           i += 1
       if i >= la:
           break


       # now index i is a mismatch column
       for start in range(i, i - cf.MinPrimerLength, -1):
           if a[start] == cf.INF:
               start += 1
               break

       # now start is the start index of a window ending in position i which has at
       # most the length MinPrimerLength (and shorter if an INF was met).

       for end in range(i, i+cf.MinPrimerLength):
           if a[end] == cf.INF:
               break


       # now end is the end index of a window starting in position i which has at
       # most the length MinPrimerLength (and shorter if an INF was met).

       if end-start+1 < cf.MinPrimerLength:
           # impossible to find a sufficiently long window around position i
           pr[i] = False

       else:
           ok = False
           for ws in range(start, end-cf.MinPrimerLength+2):
               if a[ws:ws+cf.MinPrimerLength].count(-cf.p1) <= cf.MaxMismatches:
                   ok = True
                   break

           if not ok:
               pr[i] = False
       i += 1

   indices = []

   start = end = 0
   for i in range(len(pr)):
       if not pr[i]:
           if end - start >= cf.MinPrimerLength:
               indices.append( (start, end) )
           start = end = i+1
       else:
           end = i+1
   if end - start >= cf.MinPrimerLength:
       indices.append( (start, end) )

   #printslice(allseq, 0, la)
   return indices




# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------


def findprimers( verbose, allseq, summary, l ):
    """verbose is 1 if we want comments printed, 0 otherwise"""






    # calculate summary information for all columns and
    # create vector of columnscores from score matrix:
    mpsvector = []
    colsum = []
    # list of intron start, end indices (BOTH inclusive!), initialized
    # with dummy value:
    intronindices = [ (-2, -2) ]

    # create a 'real indices' vector which in entry i holds the estimated
    # real index (estimated because of the unspecific intron markers) of
    # position i in the alignment:
    realindices = []
    totalintronlength = 0
    lastintronlength = 0
    lenallseq = len(allseq)

    for i in xrange(l):
		# deprecated:
        #column = summary.get_column(i)
        column = summary.alignment[:,i]

        if cf.INTRONS == 'yes':

            # remember last index only of each intron:
            if ('X' in column or 'x' in column):

                realindices.append(i+totalintronlength) # a dummy value really
                lastintronlength += 1

                if i>intronindices[-1][1] + 1:
                    intronindices.append( [i, i] ) # new intron
                else:
                    intronindices[-1][1] = i     # same intron, update last index
            else:
                if lastintronlength > 0:
                    # intron ended.

                    # IF several sequences in the alignment have intron markers,
                    # and they are misaligned - like this:
                    #
                    # AAGGTTXXXAATTCC
                    # AAGGTTCCXXXXGCAT
                    #
                    # - there might be more than 6 consecutive columns holding
                    # X'es. In that case, find the longest series of X'es in the
                    # same sequence

                    lil = 0
                    for tse in range(lenallseq):
                        # deprecated:
                        #xs = allseq[tse].seq.data[i-lastintronlength:i].count('X')
                        xs = allseq[tse].seq[i-lastintronlength:i].count('X')
                        if xs > lil:
                            lil = xs

                    if lil > 6:
                        raise(ValueError, "Weird marker length at index %d: %d"%(lil, i-lastintronlength))

                    lastintronlength = lil


                    trull = (cf.trueLengthOfIntronWithMarkerLength(lastintronlength))
                    totalintronlength += trull-lastintronlength # don't count marker length
                    lastintronlength = 0
                realindices.append(i+totalintronlength)
        else:
            # we don't use special intron markers
            realindices.append(i) # same as regular index

        colsum.append( columnsummary( column ) )

        if colsum[-1][1]>23: # this is out of scope of the score matrix
            mpsvector.append( cf.scorematrix[colsum[-1][0]][23] )
        else:
            mpsvector.append( cf.scorematrix[colsum[-1][0]][colsum[-1][1]] )




    # delete dummy value:
    del intronindices[0]





    #printslice(allseq, 1895, 1905)


    # calculate conservation vector (if index i is 1 it means that column i
    # sits in a ConservationWindow with a conservation-% of at least
    # MinConservationPercent, otherwise it has value 0):

    colsumlen = len(colsum)
    conservation = [0]*colsumlen

    if colsumlen >= cf.ConservationWindow:
        # init:
        mism = 0
        for i in xrange(cf.ConservationWindow):
            if colsum[i][0] > 1:
                mism += 1

        # now mism holds number of mismatches in first ConservationWindow columns
        ws = 0
        we = cf.ConservationWindow # not inclusive
        lastone = -1 #index of last 1 in conservation list

        while 1:
            if (100.0*(cf.ConservationWindow - mism))/cf.ConservationWindow >= cf.MinConservationPercent:
                # here's a good window!
                for i in xrange(max(lastone+1, ws), we):
                    conservation[i] = 1
                lastone = we-1

            if we >= colsumlen:
                break

            # Invariant: mism holds #ambiguities in [ws, we[
            if colsum[ws][0] > 1:
                mism -= 1
            if colsum[we][0] > 1:
                mism += 1
            ws += 1
            we += 1

    # now conservation is a list of 0's and 1's: a 1 in index i means that column
    # i resides in a window of sufficient conservation, a 0 that it doesn't.





    # now make a 'safety zone' around introns such that no primers are
    # allowed too close to an intron:
    for ii in intronindices:
        i1 = max(0, ii[0]-cf.MinDistanceToIntron)
        i2 = min(len(mpsvector)-1, ii[1]+cf.MinDistanceToIntron+1 )
        for jj in range( i1, i2):
            mpsvector[jj] = cf.INF




    # now mpsvector is a list of column scores
    # and colsum is a list of tuples, one per column, with the information
    #  (#different nucleotides in this column, #nucleotides in total in this clmn)


    # find the primer regions:

    primerregions = find_primer_regions( mpsvector )
    if verbose:
        if cf.INTRONS == 'yes':
            print(' found introns:', intronindices)
        print(' found regions:', primerregions)


    # Now we have a list of primer regions *with no introns in them*.




    # primers = []

    tmlow = manymism = highdiv = mmend = 0


    # list of lists of primers; each inner list is associated to a primer region:
    regionprimerlists = [ ]


    for s,e in primerregions:

        regionprimerlists.append( [] )
        #print

        # initialize counting variable:
        inthisregion = 0

        for i in xrange(s, e):

            # only consider primers with length inside min and max requirements:
            endpoint = min( e+1, 1+i+cf.MaxPrimerLength )

            for j in xrange(i + cf.MinPrimerLength, endpoint):


                # i, j are start, end indices (end NOT inclusive) of
                # alignment which we must now test as a primer candidate.

                # look at the chosen part of the first sequence with
                # no gaps and use that as primer (candidate):


                ppp = []
                for q in xrange(i, j):
                    # deprecated:
                    #if allseq[0].seq.data[q] != '-':
                    if allseq[0].seq[q] != '-':
                        # deprecated:
                        #ppp.append(allseq[0].seq.data[q])
                        ppp.append(allseq[0].seq[q])
                    else:
                        # pick symbol from one of the other seqs:
                        for rl in xrange(1, lenallseq):
                            # deprecated:
                            #if allseq[rl].seq.data[q] != '-':
                            if allseq[rl].seq[q] != '-':
                                # deprecated:
                                #ppp.append(allseq[rl].seq.data[q])
                                ppp.append(allseq[rl].seq[q])
                                break # found non-gap symbol
                        else:
                            # all-gap column! shouldn't ever happen,
                            # but to avoid crash insert random symbol:
                            ppp.append('C')

                part = "".join(ppp)


                ##
                ## old method of defining the 'part':
                ##
                ##
                # part = '--'
#                 q=0
#                 while '-' in part and q < len(allseq):
#                     part = allseq[q].seq.data[i:j]
#                     q+=1

#                 for ci in xrange(len(part)):
#                     if part[ci] == '-':
#                 if '-' in part:
#                     print i, j, part
#                     sys.exit("Please edit the alignment to remove all-gap columns.")





                # part is now the primer candidate sequence.

                # NOTE that at this point we don't have any ambiguity symbols
                # in the primer, we have simply copied a piece from one of the
                # sequences. Thus the melting temperatures calculated below are
                # only correct for one version of this primer.


                #
                # is this primer candidate a keeper?
                #
                # if one of these conditions is true, don't keep the primer:
                #
                #

                if colsum[i][0] > 1 or colsum[j-1][0] > 1:
                    mmend += 1
                    continue  # we can't end in a mismatch


                # melting temp. le Novere/sugimoto (see config.py):


                Tm = cf.TM.tm( part )



                # check if Tm is too low or high:

                if Tm < cf.MinTm:
                    tmlow += 1
                    continue

                if Tm > cf.SuggestedMaxTm:
                    continue




                (tmm, mm) = too_many_mismatches( colsum, i, j, Tm )
                if tmm:
                    # print "too many mismatches (%d, %d), tm=%f:"%(i, j, Tm)
                    # printslice(allseq, i, j)
                    manymism += 1
                    continue


                # note that mm is a *list* of mismatch positions
                lenmm = len(mm)
                partlen = j-i



                # mismatches/length ratio:

                # can't have more than 1 ambiguity if length is below 21
                if lenmm > 1 and partlen < cf.MinLengthWithTwoAmbiguities:
                    continue

                # can't have more than 2 ambiguities if length is below 25
                if lenmm > 2 and partlen < cf.MinLengthWithThreeAmbiguities:
                    continue


                if too_high_diversity_in_mismatches( colsum, i, j ):
                    highdiv += 1
                    continue


                # if four mismatches, they must be spread out over a window of
                # at least WindowWithFourMismatches:
                if lenmm > 3 and abs( mm[-1] - mm[0] ) < cf.WindowWithFourMismatches:
                    continue



                # find average number of sequences in the primer's alignment part:
                ii = 0.0
                twoseqs = 0 # number of columns with only two sequences represented
                mmt = 0     # number of mismatches in cols with only two seqs
                for q in xrange(i, j):
                    ii += colsum[q][1]
                    if colsum[q][1] == 2:
                        twoseqs += 1
                        if colsum[q][0] > 1:
                            mmt += 1

                sia = ii/partlen # avg. number of seqs in alignment per column



                # if it is entirely based on two seqs, it can have at most
                # two mismatches:
                if twoseqs == partlen and lenmm > 2:
                    continue





                # if at least 2/3 of the primer is based on only two sequences,
                # it must be in a highly conserved region, and the 2-seq part
                # can have at most 2 mismatches ( and there may be only 3 in total):
                if twoseqs >= .67*partlen and (0 in conservation[i:j] or mmt > 2 or lenmm > 3):
                    continue

                score = None # scoreprimer( part, colsum, i, j, len(mm) )


                # find AT content in last X nucleotides (in both ends):
                ATl = part[:cf.TailLength].count('A') + part[:cf.TailLength].count('T')
                ATr = part[-cf.TailLength:].count('A') + part[-cf.TailLength:].count('T')

                if ATl == cf.TailLength and ATr == cf.TailLength:
                    continue # can't have all AT 3'-end



                if cf.INTRONS == 'yes':

                    # find nearest intron in both directions:
                    nearestintronl = nearestintronr = -1 # dummy values

                    for ii in range(len(intronindices)):
                        if intronindices[ii][0] >= j:
                            nearestintronr = intronindices[ii][0] - j
                            break
                        nearestintronl = i - intronindices[ii][1] - 1
                else:
                    nearestintronl = nearestintronr = cf.MinPreIntronLength+1
                    # these values won't affect anything










                # create primer object:
                pp = primer(s, e, i, j, part, Tm, mm, nearestintronl, nearestintronr, sia, ATr, ATl, twoseqs )
                scoreIndividualPrimer( pp, colsum, conservation )
                #print i, j, pp.score

                # if this candidate overlaps with another candidate already found
                # and obeys to the same restrictions as the other (first set of
                # level three criteria), keep only the better one.

                foundbettersimilar = 0
                jf = 0

                # NB: to skip the pruning step here (which shaves perhaps a factor
                # of 1000 off the time..), just put a '>' instead of the '<' in
                # the while condition..:
                # (debug:)
                while jf < len(regionprimerlists[-1]):

                    other = regionprimerlists[-1][jf]
                    if findPrimerOverlap( pp, other ) >= cf.MatchOverlap:


                       # yes, this other one overlaps. Can we throw one of them out?

                       # ( overlap is important since it ensures that the product
                       # length
                       # of a pair with one primer won't differ significantly from
                       # the product length if the other primer were used)


                       # The other is better than pp if:
                       #  if pp's distance to nearest left (right) intron is good,
                       #     then so must other's distance to nearest left (right)
                       #     intron be;
                       #  if pp has at least one G or C in its left (right) tail,then
                       #     so must other.
                       #  if pp has a perfect match in its left (right) terminal,then
                       #     so must other.
                       #  other must have a better score than pp.

                       # all the conditions must be fulfilled so we won't throw out a
                       # candidate with a perhaps slightly worse score but with
                       # better chances of being part of a succesful PAIR of primers
                       # later on. Thus, if pp has a chance of being part of a pair,
                       # then other's chance must be at least as good, and vice versa

                       if ( pp.nearestLeftIntron >= cf.MinPreIntronLength ) <= \
                          ( other.nearestLeftIntron >= cf.MinPreIntronLength ) and \
                          ( pp.nearestRightIntron >= cf.MinPreIntronLength ) <= \
                          ( other.nearestRightIntron >= cf.MinPreIntronLength ) and \
                          ( pp.ATsInLeftTail < cf.TailLength ) <= \
                          ( other.ATsInLeftTail < cf.TailLength ) and \
                          ( pp.ATsInRightTail < cf.TailLength ) <= \
                          ( other.ATsInRightTail < cf.TailLength ) and \
                          ( colsum[pp.start][0] == 1 ) <= ( colsum[other.start][0] == 1 ) and \
                          ( colsum[pp.end-1][0] == 1 ) <= ( colsum[other.end-1][0] == 1 ) and \
                          other.fwscore >= pp.fwscore and \
                          other.rvscore >= pp.rvscore and \
                          other.score > pp.score:

                           # yes, other is better, don't keep pp
                           foundbettersimilar = 1
                           break
                       elif ( other.nearestLeftIntron >= cf.MinPreIntronLength ) <= \
                            ( pp.nearestLeftIntron >= cf.MinPreIntronLength ) and \
                            ( other.nearestRightIntron >= cf.MinPreIntronLength ) <= \
                            ( pp.nearestRightIntron >= cf.MinPreIntronLength ) and \
                            ( other.ATsInLeftTail < cf.TailLength ) <= \
                            ( pp.ATsInLeftTail < cf.TailLength ) and \
                            ( other.ATsInRightTail < cf.TailLength ) <= \
                            ( pp.ATsInRightTail < cf.TailLength ) and \
                            ( colsum[other.start][0] == 1)<= ( colsum[pp.start][0] == 1 ) and \
                            ( colsum[other.end-1][0] == 1)<= ( colsum[pp.end-1][0] == 1 ) and \
                            pp.fwscore >= other.fwscore and \
                            pp.rvscore >= other.rvscore and \
                            pp.score > other.score:

                           # pp is better, remove other:
                           del regionprimerlists[-1][jf]
                       else: # they overlap but neither is better than the other
                           jf += 1


                    else:
                        jf += 1


                if not foundbettersimilar:
                    regionprimerlists[-1].append( pp )
                    # primers.append( pp )
                    inthisregion +=1

        if len( regionprimerlists[-1] ) > 0:

            # now prune the region and keep only the best of overlapping candidates:

            ## we do it each time we add a candidate instead, see above.

            if verbose:
                print('  region [%4d, %4d[ : %d candidates found'%(s, e,inthisregion))
        else:
            del regionprimerlists[-1] # don't keep an empty region list




    if 0 and verbose: # taken out for now
        if tmlow>0:
            print(' Too low Tm in both directions: %d'%tmlow)
        if mmend>0:
            print(' Terminal ambiguity: %d'%mmend)
        if manymism>0:
            print(' Too many ambiguities: %d'%manymism)
        if highdiv>0:
            print(' Too high diversity in mismatch: %d'%highdiv)
        print(' %d single primer candidates found in total'%len(primers))



    if regionprimerlists == []:
        #    if len(primers) == 0:
        return None



    # we compare all with all regardless of individual primer scores, so
    # skip the sorting:

    # sort overall primer list:
    # primers.sort()

    # sort each regional list of primers:
    # for primerlist in regionprimerlists:
    #     primerlist.sort()



    # these two lists have references to the same tuples (primers), one is
    # just a list of references to the tuples, the other is a list of lists
    # of references to the tuples. Each inner list represents a primer region.





    # print_N_best_primers( primers, 3 )



    # print_best_primer_in_each_region( regionprimerlists )




    # Now try all primer pairs



    if 0 and verbose:
        print(' ---------------------------\n checking all primer pairs..')
    primerpairs = []

    # stats variables:
    tmth = sameregion= difftm = nointrons = x6p2 = noID = IDtoofar = ne3m = tsmt = tmm = spi = pl = oot = 0

    pairsscored = 0
    tid1 = time.time()





    # now pair all with all:



    for a in range( len( regionprimerlists ) ):

        for p1 in regionprimerlists[a]:

            # primer temperature can't be too small or large:
            # CHECKED ABOVE
            #if p1.tm < cf.MinTm or p1.tm > cf.SuggestedMaxTm:
            #    continue

            # We have to have at least X perfect matches in the 3'):

            # if mmf == MaxMismatches: # has the maximum allowed number of mismatches
            #    x = Min3endPerfectMatchesIfPrimerHasMaxMismatches
            # else:
            x = cf.Min3endPerfectMatches

            if not has_x_perfect_end_matches( x, 1, colsum, p1.start, p1.end ):
                continue


            if p1.ATsInRightTail == cf.TailLength:
                continue # AT tail

            # p1 is okay as forward primer, find a partner from another region:

            # if we're looking at sequences with no intron symbols, we have
            # to pair primers from the same region, otherwise we don't:
            if cf.INTRONS == 'no':
                sreg = a
            else:
                sreg = a+1

            for b in range( sreg, len( regionprimerlists ) ):

                intronsbetweenregions, x6 = IntronsBetweenRegions( intronindices, regionprimerlists[a][0].regionstart, regionprimerlists[b][0].regionstart )

                # (these are set to [] and 0 if cf.INTRONS == 'yes' is 0)



                if cf.INTRONS == 'yes':
                    # there has to be at least one intron between the regions:
                    if len( intronsbetweenregions ) == 0:
                        # no introns, don't consider the a region together with
                        # the b region.
                        continue

                # primers can't be too far from each other:
                # if there is an XXXXXX intron between the primers,
                # there can be at most
                # two introns in total between them.

                if x6 > 0 and len( intronsbetweenregions ) > 2:
                    continue




                # okay, go on to score all pairs from these regions:

                for p2 in regionprimerlists[b]:

                    # p1 is fw primer, p2 is rv primer


                    # We have to have at least X perfect matches in the 3' end
                    # (i.e. left end of right primer):

                    # if mmf == MaxMismatches:
                    #    x = Min3endPerfectMatchesIfPrimerHasMaxMismatches
                    # else:
                    x = cf.Min3endPerfectMatches

                    if not has_x_perfect_end_matches(x, 0, colsum, p2.start, p2.end):
                        continue





                    # AT tail in 3'-end:
                    if p2.ATsInLeftTail == cf.TailLength:
                        continue


                    # at least one primer must have a certain distance to
                    # the nearest intron:
                    if p1.nearestRightIntron < cf.MinPreIntronLength and p2.nearestLeftIntron < cf.MinPreIntronLength:
                        continue
                    # both are set right if cf.INTRONS == 'yes' is 0



                    # CHECKED ABOVE:
                    # primer temperature can't be too small:
                    #if p2.tm < cf.MinTm or p2.tm > cf.SuggestedMaxTm:
                    #    print '\nSHOULD NEVER HAPPEN I THINK'
                    #    # jeg tror nok Tm er ens i begge retninger, og dem
                    #    # der er for lave i begge retninger er renset ud tidligere
                    #    continue


                    # primer melting temperatures have to be similar:
                    if abs( p1.tm - p2.tm ) > cf.MaxPrimerPairTmDifference:
                        continue



                    # find score of primer pair (obtain no textual report)
                    score, pro, contra = scoreprimerpair( p1,p2, realindices, intronsbetweenregions)

                    pairsscored += 1

                    # debugging:
                    # for checking a particular primer pair:
                    if 0:



                        if p1.start==456 and p1.end == 487 and \
                               p2.start==703 and p2.end==727:


                        #checkstring1 =insertAmbiguities( p1.seq,p1.start,p1.end,0,summary, colsum)
                        #checkstring2 =insertAmbiguities( p2.seq,p2.start,p2.end,1,summary, colsum)
                        #if checkstring1 == 'GCACYATAATTTATTTGCTTGGGCARCAGCT'\
                        #       and checkstring2 == 'TGGAGCACCACTTTGCTTAATAGC':

                            # get explanations:
                            scoreIndividualPrimer( p1, colsum, conservation, 1 )
                            scoreIndividualPrimer( p2, colsum, conservation, 1 )
                            score1,pro1,contra1=scoreprimerpair( p1,p2, realindices, intronsbetweenregions, 1)

                            print('\nCHECKED PRIMER PAIR:')
                            #print p1
                            #print p2


                            print('score',score1)
                            print(primer2string( p1, 0, summary, colsum ))
                            print(primer2string( p2, 1, summary, colsum ))
                            print(pro1, contra1)

                            # printslice(allseq, 1447, 1478)





                    if score > 0:
                        primerpairs.append( [score, p1, p2, pro, contra] )

                        # elif score == -1:
                        #     sameregion += 1
                    elif score == -2:
                        difftm += 1
                    elif score == -3:
                        nointrons += 1
                    elif score == -4:
                        x6p2 += 1
                    elif score == -5:
                        noID += 1
                    elif score == -6:
                        IDtoofar += 1
                    elif score == -7:
                        ne3m += 1
                    elif score == -8:
                        tsmt += 1
                    elif score == -9:
                        tmm += 1
                    elif score == -10:
                        tmth += 1
                    elif score == -11:
                        spi += 1
                    elif score == -12:
                        pl += 1
                    elif score == -13:
                        oot += 1

                        # if score <-1:
                        #    print primer2string(p1, 0)
                        #    print primer2string(p2, 1),' score %d\n'%score

    tid2 = time.time()
    if verbose:
        print('loop time: %3.3f s.'%((tid2-tid1)))

    if 0 and verbose:
        print('done')

    if 0 and verbose: # taken out for now
        # print ' Pairs from same region     :%5d'%sameregion
        if difftm>0:
            print(" Tm's too different          : %6d"%difftm)
        if tmth>0:
            print(' Tm too high                 : %6d'%tmth)
        if nointrons>0:
            print(' No introns between primers  : %6d'%nointrons)
        if x6p2>0:
            print(' XXXXXX intron plus two      : %6d'%x6p2)
        if noID>0:
            print(' No valid ID-exon            : %6d'%noID)
        if IDtoofar>0:
            print(' ID-exon unreachable         : %6d'%IDtoofar)
        if ne3m>0:
            print(" Too few 3'-end matches      : %6d"%ne3m)
        if tsmt>0:
            print(" Too small Tm                : %6d"%tsmt)
        if tmm>0:
            print(" Many ambiguities + small Tm : %6d"%tmm)
        if spi>0:
            print(" Introns too close to primers: %6d"%spi)
        if pl>0:
            print(" Invalid est. product length : %6d"%pl)
        if oot>0:
            print(" ID exons too short & too far: %6d"%pl)

        print(" Possibly valid pairs        : %6d"%len(primerpairs))
        print('                              ------')

    if verbose:
        print(' Total primer pairs considered   : %6d'%(pairsscored))



    if len(primerpairs) == 0:
        return None

    primerpairs.sort() # the score has to be the first entry of the primer pair tupel




    # report PrimerPairSuggestions best pairs that don't overlap:

    if verbose:
        print("\nBest %d non-overlapping primer pairs:\n"%cf.PrimerPairSuggestions)
    rep = 0
    reported = [] # remember indices of pairs to report
    count = -1 # overall best pair is last in the list
    while 1:
        if rep == cf.PrimerPairSuggestions: # found enough
            break
        s, p1, p2, pro, contra = primerpairs[count]
        if -count == len(primerpairs):
            break # didn't find PrimerPairSuggestions different pairs..

        seenalready = 0
        sawleft = sawright = 0
        for i in reported:
            #p33 = primerpairs[i][1].start
            #p34 = primerpairs[i][1].end
            #p43 = primerpairs[i][2].start
            #p44 = primerpairs[i][2].end


            # fwoverlap, rvoverlap = findOverlap( p1.start, p1.end, p2.start, p2.end, p33, p34, p43, p44 )
            fwoverlap = findPrimerOverlap( p1, primerpairs[i][1] )
            rvoverlap = findPrimerOverlap( p2, primerpairs[i][2] )

            #heq = p1.start==1069 and p1.end==1100 and p2.start==1316 and p2.end==1345
            #if heq: print '\nfwoverlap: %s  rvoverlap: %s\n'%(fwoverlap,rvoverlap)

            # if fw primers overlap and rv primers overlap
            # with a pair already reported, don't report this one:
            if fwoverlap >= cf.MatchOverlap and rvoverlap >= cf.MatchOverlap:
                # both this pair's primers overlap with the current primers:
                seenalready = 1
                break

            elif fwoverlap >= cf.MatchOverlap:
                # fw primers overlap
                if sawleft >= 1: # note, overlap means we already have one..
                    seenalready = 1
                    break # report max 2 pairs with overlapping fw primers
                sawleft += 1

            elif rvoverlap >= cf.MatchOverlap:
                # rv primers overlap
                if sawright >= 1:
                    seenalready = 1
                    break # report max 2 pairs with overlapping rv primers
                sawright += 1


        if not seenalready:
            reported.append( count )
            # print primerpairs[count][1].start,primerpairs[count][1].end, primerpairs[count][2].start,primerpairs[count][2].end

            rep += 1
        count -= 1


    #reported.append(-2) # also report 2nd and 3rd best overall
    #reported.append(-3)

    returnprimers = [] # list of primer pairs to return
    for count in reported:
        score, p1, p2, pro, contra = primerpairs[count]

        # now score the pairs again but obtain the explanations:
        intronsbetweenprimers, x6 =IntronsBetweenRegions(intronindices, p1.regionstart, p2.regionstart)
        scoreIndividualPrimer( p1, colsum, conservation, 1 )
        scoreIndividualPrimer( p2, colsum, conservation, 1 )
        s,pro,contra=scoreprimerpair( p1,p2,realindices,intronsbetweenprimers, 1)

        if s != score:
            print('something is WRONG %f %f\n\n'%(s, score))
            sys.exit(1)




        q1 = insertAmbiguities( p1.seq, p1.start, p1.end, 0, summary, colsum)
        q2 = insertAmbiguities( p2.seq, p2.start, p2.end, 1, summary, colsum)


        # return the explanations and indices too (indices are needed by the gui):
        txt = "Fw 5'-%s\nRv 5'-%s\n\n%s%s\nScore: %d"%(q1, q2, pro, contra, s)
        #txt = txt + pro + contra+"\nScore: %d"%s


        returnprimers.append((q1, q2, (txt, p1.start, p1.end, p2.start, p2.end, p1.tm, p2.tm)))



        if verbose:
            print('-'*50)
            print('score %.1f (rank %d)'%(s, -count))
            print(primer2string( p1, 0, summary, colsum ))
            # print p2.start, p2.end, reverse_and_complement(p2.seq), p2.tm, p2.tm
            print(primer2string( p2, 1, summary, colsum ))

            print(pro, contra)


            #         print '\n\nbedste sidst:\n\n'
            #         for s, p1, p2, pro, contra in primerpairs[reported[-1]:]:

            #             print '-'*50
            #             print 'score %d'%(s)
            #             print primer2string( p1, 0, summary, colsum )
            #
            #             print primer2string( p2, 1, summary, colsum  )

            #             print pro, contra

            # print 'Best pair (score %d):\n'%primerpairs[-1][0],primerpairs[-1][3], primerpairs[-1][4]

    # s, p1, p2, pro, contra = primerpairs[-1]
    # q1 = primer2string( p1, 0, summary, colsum ).split()[-1]
    # q2 = primer2string( p2, 1, summary, colsum ).split()[-2]
    # use a more convenient function:
    # q1 = insertAmbiguities( p1.seq, p1.start, p1.end, 0, summary, colsum)
    # q2 = insertAmbiguities( p2.seq, p2.start, p2.end, 1, summary, colsum)
    # return (q1,q2)
    return returnprimers




def test( tfil, start = '0'):
    """This test function takes a spreadsheet saved as a text file as argument. This file has a special format with .. well, found primers should appear in columns 7 and 9. Further, the start argument specifies the lowest number of primer pair (legume ID) to test agains. It is a string on the form 'N' or 'N-'; in the first case, only legN is tested."""

    try:
        fin = open(tfil)
    except IOError:
        sys.exit( "Couldn't open test file "+tfil )

    fin.readline()

    stands = 0
    if start[-1] == '-':
        start = int(start[:-1])
    else:
        start = int(start)
        if start > 0:
            stands = 1 # stop after this one

    while 1:
        line = fin.readline()
        if line=='':
            break
        f = line.replace('"','').split(',')


        if f[0] == '' or int(f[0][3:])<start:
            continue

        if stands and int(f[0][3:]) > start:
            break

        i = f[1] # filename

        print("\n[31;2m%s (%s):[0m"%(i, f[0])),
        o1 = f[6].strip()
        o2 = f[8].strip()


        allseq, summary, l = doAlignment( i, 0 )


        primerpairs = findprimers( 0, allseq, summary, l ) # get list of primer pair suggestions

        nomatch = 'Disagreement'
        nomatchstring = []
        if not primerpairs:
            nomatch = 'Disagreement, no primer pairs found'
            nomatchstring.append('None  None\n')
        else:
            ruth = 0
            for p1, p2,qz in primerpairs:
                ruth += 1
                if p1 != o1 or p2 != o2:
                    nomatchstring.append('%s  %s\n'%(p1, p2))
                    if (p1 != None and p2 != None):


                        # if both fw/fw and rv/rv primers overlap with at least
                        # MatchOverlap, consider the pairs almost identical:
                        fwoverlap = rvoverlap = 0
                        for startpoint in range(len(p1)-cf.MatchOverlap+1):
                            if p1[startpoint:startpoint+cf.MatchOverlap] in o1:
                                fwoverlap = 1
                                break

                        for startpoint in range(len(p2)-cf.MatchOverlap+1):
                            if p2[startpoint:startpoint+cf.MatchOverlap] in o2:
                                rvoverlap = 1
                                break

                        if fwoverlap and rvoverlap:
                            nomatch = 'almost identical (to %d. suggestion)'%ruth
                            break

                else:
                    print('identical')
                    nomatch = 0
                    break

        if nomatch:
            print(nomatch)
            if nomatch[0] == 'D':
                print('lene / prg:')
                print(o1,' ', o2)
                print()
                print(''.join(nomatchstring))


                # save Lenes fw and rv primers in two fasta files so that they
                # can easily be aligned with the others with clustalx:
                f1 = i+'.LHM.fw'
                f2 = i+'.LHM.rv'
                try:
                    fw = open( f1, 'w' )
                    rv = open( f2, 'w' )
                except IOException:
                    sys.exit("Couldn't write Lene's primer files")

                print(">LHM_fw_%s"%i, file=fw)
                print(o1, file=fw)
                print(">LHM_rv_%s"%i, file=rv)
                print(reverse_and_complement(o2), file=rv)

                fw.close()
                rv.close()
                print('wrote primer files %s and %s (for aligning)'%(f1, f2))
                print('(note that the rv primer is not revcomplemented in this file)')



    fin.close()


def writePrimersToFiles( path, primerpairs, verbose = 0 ):
    """Write report file and primers in fasta files. The method may throw exception"""


    ll = 60 # length of lines in report file
    f3 = path+'.rep'
    rep = open ( f3, 'w' )
    print("PriFi report. Suggested primers after analysis of this file:\n%s\n%s\n"%(path, "="*ll), file=rep)
    for i in range(len(primerpairs)):

        # save forward and reverse primers in two fasta files so that they
        # can easily be aligned with the others with clustalx:

        p1, p2, qz = primerpairs[i]

        f1 = path+'.fw%d'%(1+i)
        f2 = path+'.rv%d'%(1+i)

        fw = open( f1, 'w' )
        rv = open( f2, 'w' )


        print(">fw%d_%s (Tm = %4.1f) 5'-3'"%(1+i, path, qz[5]), file=fw)
        print(p1, file=fw)
        print(">rv%d_%s (Tm = %4.1f) 5'-3'"%(1+i, path,qz[6]), file=rv)
        print(p2, file=rv) #reverse_and_complement(p2) # p2 is already rev.complemented
        print("Primer set %d\t\t\t\t(%d-%d / %d-%d)\n"%(i+1, qz[1], qz[2], qz[3], qz[4]), file=rep)
        print("%s\n\n%s\n"%(qz[0], "-"*ll), file=rep)
        fw.close()
        rv.close()

        if verbose:
            print('%s\n%s'%(f1, f2))

    rep.close()
    if verbose:
        print(f3)



def hasIntrons( fname ):
    from Bio import SeqIO
    for rec in SeqIO.parse( open(fname, 'rt'), "fasta" ):
        if rec.seq.count( 'X' ) > 0:
            return True
    return False



if __name__ == '__main__':
    """given argument '-t' means run test function, otherwise needs multiple
    fasta file. Possible command lines:
    primerfinder.py At5g38900.1
    primerfinder.py -t
    primerfinder.py -t 66
    primerfinder.py -t 66-
    """

    if len(sys.argv) < 2:
        sys.exit( "Please provide a fasta file as argument" )

    if sys.argv[1] == '-t':
        if len(sys.argv) > 2:
            start = sys.argv[2]
        else:
            start = '0'
        test( "/users/chili/BiRC/MolBio/Data/ReadyAndFormatted/leg60-174_primers.csv", start)
        sys.exit(0)

    # check input for intron hints
    cf.INTRONS = 'yes' if hasIntrons( sys.argv[1] ) else 'no'

    allseq, summary, l = doAlignment( sys.argv[1], 1 )

    primerpairs = findprimers( 1, allseq, summary, l )



    if not primerpairs:
        print("No valid primer pair found")
        sys.exit( 0 )

    print('Found %d primer pair suggestions. Writing primer files:'%len(primerpairs))
    try:
        writePrimersToFiles( sys.argv[1], primerpairs, 1 )
    except:
        print("Couldn't write primer files")
