from math import log10


class Tm:
    """class for calculating melting temperature of a sequence"""
    
    def __init__( self, primerconcNmolar = 250, saltconcMolar = 0.05 ):
        # -----------------------------------------------
        # tables found here:
        # http://nar.oupjournals.org/cgi/content/full/24/22/4501
        # (Sugimoto)
        # kcal/molar
        self.dH = { 'AA': -8.0, 'TT': -8.0, 'AT': -5.6, 'TA': -6.6,
                    'CA': -8.2, 'TG': -8.2, 'GT': -9.4, 'AC': -9.4,
                    'CT': -6.6, 'AG': -6.6, 'GA': -8.8, 'TC': -8.8,
                    'CG':-11.8, 'GC':-10.5, 'GG':-10.9, 'CC':-10.9,
                    'NG': -9.4, 'NC': -9.9, 'NA': -7.9, 'NT': -7.4,
                    'CN': -9.4, 'GN': -9.9, 'AN': -7.4, 'TN': -7.9,
                    'NN': -8.6, 'initGC': 0.6, 'initAT':0.6 }
        
        # cal/(K*molar)
        self.dS = { 'AA':-21.9, 'TT':-21.9, 'AT':-15.2, 'TA':-18.4,
                    'CA':-21.0, 'TG':-21.0, 'GT':-25.5, 'AC':-25.5,
                    'CT':-16.4, 'AG':-16.4, 'GA':-23.5, 'TC':-23.5,
                    'CG':-29.0, 'GC':-26.4, 'GG':-28.4, 'CC':-28.4,
                    'NG':-23.7, 'NC':-25.9, 'NA':-21.2, 'NT':-19.8,
                    'CN':-23.7, 'GN':-25.9, 'AN':-19.8, 'TN':-21.2,
                    'NN':-22.6, 'initGC': -9.0, 'initAT':-9.0 }
        
        GasConstant = 1.987   # cal / (K*mol)
        
        # constants for calculating Tm - le Novere:
        Ct =  primerconcNmolar / 4000000000.0  # primer conc. in molar / 4

        # Ct burde maaske aendres for degenerate primers hvor koncentrationen
        # er mindre af hver primer, eller hvad..
        # se http://www.changbioscience.com/primo/pconc.html, som jo dog
        # er et kommercielt produkt.
        
        self.konstant1 = GasConstant * log10(Ct)
        logsalt = log10(saltconcMolar)
        self.logsaltf = logsalt * .368
        self.konstant2 = 16.6 * logsalt - 273.15



    def tm( self, seq ):
        """calculate the le Novere melting temperature, formula found in:
   http://www.pasteur.fr/recherche/unites/neubiomol/SOFTWARES/melting/melting.html"""

        dh = ds = 0
        for j in range(len(seq)-1):

            bid = seq[j:j+2]

            # handle ambiguity code characters
            # ( see http://mbcr.bcm.tmc.edu/Guide/Sequences/iupac.html):
            if 'U' in bid:
                bid = bid.replace('U', 'T')

            char = [seq[j], seq[j+1]]
            
            for i in [0,1]:
                if bid[i] == 'R':
                    char[i] = 'GA'
                elif bid[i] == 'Y':
                    char[i] = 'TC'
                elif bid[i] == 'M':
                    char[i] = 'AC'
                elif bid[i] == 'K':
                    char[i] = 'TG'
                elif bid[i] == 'S':
                    char[i] = 'GC'
                elif bid[i] == 'W':
                    char[i] = 'TA'
                elif bid[i] == 'H':
                    char[i] = 'TAC'
                elif bid[i] == 'B':
                    char[i] = 'TCG'
                elif bid[i] == 'V':
                    char[i] = 'GCA'
                elif bid[i] == 'D':
                    char[i] = 'TGA'

            if char != [0,0]:
                comb = (len(char[0])*len(char[1]))
                sumH = sumS = 0.0
                for i in char[0]:
                    for j in char[1]:
                        sumH += self.dH[ i+j ]
                        sumS += self.dS[ i+j ]

                sdh = sumH / comb
                sds = sumS / comb
            else:
                sdh = self.dH[ seq[j:j+2] ]
                sds = self.dS[ seq[j:j+2] ]
                
            dh += sdh*1000.0
            ds += sds


        # initiation values (for sugimoto the value is the same regardless
        # of terminal base pair):
        
        dh += self.dH['initGC']
        ds += self.dS['initGC']
        
        return dh / (ds + self.logsaltf*(len(seq)-1) + self.konstant1)+self.konstant2
        
