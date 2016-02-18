from Tkinter import *
from Bio import Clustalw
from Bio.Align import AlignInfo
from clustalalignment import columnsummary

class Alignmentviewer( Frame ):
   """Class for viewing an alignment (in .aln format, e.g. from Clustal) in a Frame.
   Has methods for adding sequences, e.g. primers. The linecolors argument is a list
   of pairs of colors"""
   
   def __init__( self, parent, bgcolor="LightBlue", slidercolor="#ca4c4a", linecolors = None ):

      self.W = 15
      self.H = 4
            
      Frame.__init__( self, parent )
      self.config( bg=bgcolor )
      self.pack( pady=10, padx=10,expand = YES, fill = BOTH )
      #self.master.title( "alignment viewer" )
      #self.master.geometry( "525x100" )  # width x length

      
      # create scrollbar:
      self.scrollbar = Scrollbar( self, orient=HORIZONTAL  )
      self.scrollbar.pack( fill=X, expand = YES )
      self.scrollbar.config( bg=slidercolor )
      
      self.frame = Frame(self)
      self.frame.config( bg=bgcolor )
      # create text field for sequence names:
      self.namefield = Text( self.frame, wrap=NONE, height = self.H, width = self.W )
      #self.namefield.insert( INSERT, """navn1\nnavn2""" )
      self.namefield.config( font = "Courier 15", bg="White", fg="Black", state = DISABLED )
      self.namefield.pack(side=LEFT)

      # define tags:
      self.namefield.tag_config( "p0", foreground="Red" )
      self.namefield.tag_config( "p1", foreground="Blue" )
      self.namefield.tag_config( "p2", foreground="#fe9a00" )
      self.namefield.tag_config( "p3", foreground="DarkGreen" )
      
      
      # create alignmentdata field:
      self.text1 = Text(self.frame, xscrollcommand=self.scrollbar.set, wrap=NONE, height = self.H )
      self.text1.config( font = "Courier 15", bg="White", fg="Black", state = DISABLED )
      self.text1.pack( fill=X, expand = YES, side = LEFT )
      
      # define tags:
      ii = 0
      if linecolors:
         for c1, c2 in linecolors:
            self.text1.tag_config( "p%d"%ii, foreground=c1, background=c2 )
            ii += 1
      self.text1.tag_config( "match", background="#d9d984" )
      self.text1.tag_config( "intron", background="#aaaeb2" )
      self.text1.tag_config( "A", background="#d97e84" )
      self.text1.tag_config( "C", background="#d97e34" )
      self.text1.tag_config( "G", background="#b4b89a" )
      self.text1.tag_config( "T", background="#56789a" )
      
      self.scrollbar.config(command=self.text1.xview)
      self.frame.pack( fill=X, expand = YES   )

      self.allseq = self.summary = self.l = None # alignment variables




   def getAlignment( self ):
      if self.allseq and self.summary and self.l:
         return (self.allseq, self.summary, self.l )
      else:
         return None


   def loadAlignment( self, path ):
      """ path is a path to an alignment file in .aln format"""
      alignment = Clustalw.parse_file( path )
      self.allseq = alignment.get_all_seqs()
      self.summary = AlignInfo.SummaryInfo(alignment)
      self.l = alignment.get_alignment_length()
      self.insertLoadedBioAlignment()
      


   def setHeight( self, n ):
      self.namefield.config( height = n )
      self.text1.config( height = n )
      self.H = n
      
   

      



   def insertLoadedBioAlignment( self ):
      """ aln is a list of alignment sequences, loaded via the Biopython module Clustalw"""

      names = []
      seqdata = []

      for seq in self.allseq:
         names.append( seq.description )
         seqdata.append( seq.seq.data )

      self.text1.config( state = NORMAL )
      self.namefield.config( state = NORMAL )
      self.setSequenceData( "\n".join(seqdata) )
      self.setNameData( "\n".join(names) )
      self.setHeight( len(self.allseq) )
      self.text1.config( state = DISABLED )
      self.namefield.config( state = DISABLED )
      


   def setSequenceData( self, txt ):


      # matching columns in different color (or each letter its own color)
      highlightmatches = 1 
      self.insertText( self.text1, "", clear = 1 )

      if highlightmatches:
         colsum = []
         for i in range(self.l):
            column = self.summary.get_column(i)
            colsum.append( columnsummary( column ) )

      highlightcolumns = [] # list of column indices to highlight
      introncolumns = []

      for i in xrange( self.l ):
         if colsum[i][0] == 1 and colsum[i][1] > 1:
            highlightcolumns.append(i)

            
      lines = txt.split('\n')
      for l in range(len(lines)):
         for i in range(len(lines[l])):

            if not highlightmatches:
               self.insertText( self.text1, lines[l][i], tag = lines[l][i] )
            else:
               if colsum[i][0] == 1 and colsum[i][1] > 1:
                  self.insertText( self.text1, lines[l][i], tag = "match" )
                  
               elif lines[l][i] == 'X' or lines[l][i] == 'x':
                  self.insertText( self.text1, lines[l][i], tag = 'intron' )
               else:
                  self.insertText( self.text1, lines[l][i] )

               
         if l<len(lines)-1:
            self.insertText( self.text1, "\n" )
         
      #self.insertText( self.text1, txt, clear = 1 )

      # if this view has an associated lineviewer, the lineviewer will need this list:
      if self.lineview:
         self.lineview.highlightcolumns = highlightcolumns 
         #self.lineview.introncolumns = introncolumns

   def setNameData( self, txt ):
      self.insertText( self.namefield, txt, clear = 1 ) 
      


   def addSequence( self, name, seqdata, color ):
      """insert an extra sequence in the alignment"""

      self.insertText( self.namefield, "\n" ) 
      self.insertText( self.text1, "\n" ) 
      self.insertText( self.namefield, name, tag = "p%d"%color ) 
      self.insertText( self.text1, seqdata, tag = "p%d"%color )
         
      self.setHeight( self.H+1 )




   def addSequenceData( self, seqdata, color ):
      """insert some sequencedata in the alignment"""

      self.insertText( self.text1, seqdata, tag = "p%d"%color ) # using tag for color
               
      if '\n' in seqdata:
         self.setHeight( self.H+1 )


   def scrollTo( self, offset ):
      """set scrollbar such that the given offset (between 0.0 and 1.0) of the
      window is visible in the middle."""


      self.text1.xview( MOVETO, offset - 10.0/self.l)
      return
      
      
      

   


   def insertText( self, widget, txt, clear=0, tag="p-1" ):
      widget.config( state = NORMAL )
      if clear:
         widget.delete(1.0, END)
      widget.insert( INSERT, txt, tag )
      
      widget.config( state = DISABLED )

      

def main():
   Alignmentviewer().mainloop()

if __name__ == "__main__":
   main()
