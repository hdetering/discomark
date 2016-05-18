from Tkinter import *
from Bio import Clustalw
from Bio.Align import AlignInfo
from clustalalignment import columnsummary
from config import INTRONS

LINEDISTANCE = 8  # distance between lines
HIGHLIGHTCOLOR = "#d9d984"
INTRONCOLOR = "#b2b2e8"

class AlignmentviewerLines( Frame ):
   """Class for viewing an alignment (in .aln format, e.g. from Clustal) in a Frame.
   Each sequences is represented by a line only."""
   
   def __init__( self, parent, mymaster, bgcolor="LightBlue", slidercolor ="LightBlue", linecolors=None ):
      """master is the application object that creates this frame and which supposedly
      has a method called handleMouseClick. The linecolors argument is a list of pairs
      of colors."""
      
      #self.W = 15
      self.H = 70



      # the Frame class already has a variable called 'master' which is set below
      # when the constructor is called, so this one needs another name in order not
      # to be overwritten:
      self.__master = mymaster 
      
      # this list is calculated by the regular alignment viewer if one is associated
      # with the application:
      self.highlightcolumns = None

      Frame.__init__( self, parent )

      self.config( bg=bgcolor )
      self.pack( pady=10, padx=10,expand = YES, fill = BOTH )
      #self.__master.title( "alignment viewer" )
      #self.__master.geometry( "525x100" )  # width x length

      
      # create scrollbar:
      self.scrollbar = Scrollbar( self, orient=HORIZONTAL  )
      self.scrollbar.pack( fill=X, expand = YES )
      self.scrollbar.config( bg=slidercolor )
      
      self.frame = Frame(self)
      self.frame.config( bg=bgcolor )

      # create text field for sequence names:
      #self.namefield = Text( self.frame, wrap=NONE, height = self.H, width = self.W )
      #self.namefield.insert( INSERT, """navn1\nnavn2""" )
      #self.namefield.config( font = "Courier 15", bg="White", fg="Black", state = DISABLED )
      #self.namefield.pack(side=LEFT)

      
      # create alignmentdata field:
      self.canvas = Canvas(self.frame, xscrollcommand=self.scrollbar.set,height =self.H)
      self.canvas.config( bg="White" )
      self.canvas.pack( fill=X, expand = YES, side = LEFT )

            
      # set colors for added lines:
      self.colors = linecolors

      
      self.scrollbar.config(command=self.canvas.xview)
      self.frame.pack( fill=X, expand = YES   )

      self.allseq = self.summary = self.l = None # alignment variables
      self.lines = [] # line objects

      

      self.canvas.bind( "<1>", self.__master.handleMouseClick )
      
      

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
      self.insertLoadedBioAlignment(self.allseq)
      

   def setAlignment( self, a, s, le ):

      self.allseq = a
      self.summary = s
      self.l = le
      


   #def setHeight( self, n ):
   #   self.namefield.config( height = n )
   #   self.text1.config( height = n )
   #   self.H = n
      
   

   def clear( self ):
      allitems = self.canvas.find_all()
      for item in allitems:
         self.canvas.delete(item)

      



   def insertLoadedBioAlignment( self ):
      """ aln is a list of alignment sequences, loaded via the Biopython module Clustalw"""
      allitems = self.canvas.find_all()
      for item in allitems:
         self.canvas.delete(item)

      # draw interval lines every 100 bp:
      for i in range(0, self.l, 100):
         self.canvas.create_line( i, 0, i, self.H, fill = "#e6e6e6" )

      self.canvas.create_text( 120, self.H-7, fill = "#a6a6a6", font = "Arial 10", text = '100 bp')
         
      # draw loaded alignment as lines:
      self.lines = []
      self.h = 10     # object variable that keeps track of height of last drawn line
      for seq in self.allseq:
         txt = seq.seq.data
         start = 0
         while txt[start] == '-':
            start+=1
         end = len(txt)-1
         while txt[end] == '-':
            end -= 1
         
         newline =self.canvas.create_line(start, self.h, end, self.h, width=3 )
         self.lines.append(newline)
         self.h += LINEDISTANCE

      # set the region in which the image can be scrolled:
      self.canvas.config( scrollregion = (-10, 0, 10+self.l, self.H ) )
      self.highlightColumns()

      self.canvas.update_idletasks()
      #self.printImage()




      


   def addPrimerPair( self, start1, end1, start2, end2, color ):
      """insert an extra sequence in the alignment"""
      
      self.lines.append( self.canvas.create_line( start1, self.h, end1, self.h,
                                                  fill=self.colors[color][0], width=3 ))
      self.lines.append( self.canvas.create_line( end1, self.h, start2, self.h,
                                                  fill=self.colors[color][0], width=1 ))
      self.lines.append( self.canvas.create_line( start2, self.h, end2, self.h,
                                                  fill=self.colors[color][0], width=3 ))
      self.h += LINEDISTANCE
      #self.printImage()

      


   def printImage( self ):
      
      self.canvas.postscript( file = "lines.ps" )
      print '\nwrote postscript file lines.ps'

      

   


   def highlightColumns( self ):
      """Highlight columns with perfect matches"""
      
      if INTRONS == 'yes': # if X'es are special intron markers, highlight them
         for i in xrange(self.l):
            column = self.summary.get_column(i)

            if ('X' in column or 'x' in column):
               newline = self.canvas.create_line( i, 9, i, 4+LINEDISTANCE*(len(self.allseq)), fill =INTRONCOLOR )
               self.canvas.lower( newline ) # move the highlighting under the lines

      # this list is calculated by the regular alignment viewer if one is associated
      # with the application:
      if self.highlightcolumns:
         
         for col in self.highlightcolumns:
            newline = self.canvas.create_line( col, 9, col, 4+LINEDISTANCE*(len(self.allseq)), fill =HIGHLIGHTCOLOR )
            self.canvas.lower( newline ) # move the highlighting under the lines

      
            

   def setScrollBar( self, offset ):
      """set scrollbar such that the given offset (between 0.0 and 1.0) of the
      window is visible in the middle."""

      lo, hi = self.scrollbar.get()
      
      l = hi - lo # bar length
      o1, o2 = offset - l/2.0, offset + l/2.0
      if o1<0:
         o1, o2 = 0, l
      elif o2>1:
         o1, o2 = 1-l, 1
      self.scrollbar.set( o1, o2 )
      
      

      

def main():
   Alignmentviewer().mainloop()

if __name__ == "__main__":
   main()
