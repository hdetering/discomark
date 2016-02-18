#!/usr/bin/env python


from Tkinter import *
from tkMessageBox import *   # dialogs, messages, ..
from alignmentviewer import Alignmentviewer
from alignmentviewerLines import AlignmentviewerLines
import tkFileDialog
import primerfinder_ver2 as primerfinder
import config as cf
#import tkSimpleDialog
from popup_configure import ConfigureDialog
from os.path import exists
from reversecomplement import reverse_and_complement
import re
from os import listdir
from os.path import join
import ScrolledText

fgcolor = "#9ad0ca"
bgcolor = "Black"#"#74a8ac"
buttoncolor = "#74a8ac"
W = 1600 # pixels
H = 1000



class PrimerConfigureDialog( ConfigureDialog ):

   
      

   def body(self, master): 
      #self.config( bg = self.bgcolor )
      
      
      from config import Parameters
      self.parameters = Parameters

      
      self.entries = {}
      for i in range(len(Parameters)):
         nametext, variable, description, m1, m2, default = Parameters[i]
         self.entries[nametext] =self.createEntryInGrid( master, nametext, getattr( cf, variable), i, description)





   def saveParameters( self ):

      from config import ParameterFile

      try:
         if not self.validate():
            raise Exception
         fout = open( ParameterFile, 'w' )
         for p in self.parameters:
            fout.write("%s\n"%self.entries[ p[0] ].get())
         fout.close()
         
         
         showinfo( "Saved values", "Parameter values were saved. They will be loaded the next time you run PriFi", parent=self )
         # self.loadedparameters = 1
      except:
         showinfo( "Save error", "Couldn't save parameter values", parent=self )


   



   def reset(self):
      """overridden method. Reset all parameters to their original values"""
      for p in self.parameters:
         self.entries[ p[0] ].delete(0, END)
         self.entries[ p[0] ].insert( INSERT, p[5] )
         








   def validate(self): # (this method is called before window is closed)
      
      from config import parseStringAndAssignToParameter

      for n, v, d, m1, m2, default in self.parameters:
         
         txt = self.entries[n].get()

         
         itwentokay = parseStringAndAssignToParameter( v, m1, m2, txt )

         if not itwentokay:
            showinfo( v, 'Invalid value for "%s". Explanation for this parameter:\n%s'%(n, self.ee[n.lower()]), parent=self )
            self.initial_focus = self.entries[n]
            return 0
                     
      return 1
       
      




class prifiGUI( Frame ):
   

   


   def __init__( self, aln = None ):
      """Constructor for the gui. aln is a .aln alignment file which is loaded initially if given."""


   

      Frame.__init__( self )
      self.pack( expand = YES, fill = BOTH )
      self.config( bg=bgcolor )
      self.master.title( "  --  PriFi  --  " )
      self.master.geometry( "%dx%d"%(W, H) )  # width x length

      primercolors = ( ( "Red", "Pink" ),
                       ( "Blue", "#CCCCFF" ),
                       ( "#fe9a00", "#fefe86" ),
                       ( "DarkGreen", "LightGreen" ) )

      self.filepath = None # path to current (alignment) file
      
      self.docframe = None # will hold the documentation/about/help window
      # frame 2 holds buttons and info window:
      self.frame2 = Frame(self, bg=bgcolor )  
      self.frame2.pack( expand=NO, fill=NONE, pady=10 )

      # holds two button frames:
      self.frameo = Frame(self.frame2, bg=bgcolor ) 
      self.frameo.pack( side=TOP, expand=NO, fill=NONE, padx=10 )

      # button frame holds buttons:
      self.buttonframe = Frame(self.frameo, bg=bgcolor ) 
      self.buttonframe.pack( side=LEFT, expand=NO, fill=NONE, padx=10 )
      loadbutton = self.createButton( self.buttonframe, "Load alignment", None, self.loadAlignment, LEFT)
      alignbutton = self.createButton( self.buttonframe, "Align Fasta file", None, self.doAlignment, LEFT )
      batchbutton = self.createButton( self.buttonframe, "Batch run (Fasta files)", None, self.fastaBatch, LEFT )
      prifibutton = self.createButton( self.buttonframe, "Show PriFi analysis", None, self.loadPrifiFile, LEFT )

      #holds 3 buttons
      self.buttonframe2 = Frame(self.frameo, bg=bgcolor )
      self.buttonframe2.pack( side=LEFT, expand=NO, fill=NONE, padx=10 )
      
      primerbutton = self.createButton( self.buttonframe2,"Find primers", None, self.findprimers )
      configbutton = self.createButton( self.buttonframe2, "Configure", None, self.configure )
      aboutbutton = self.createButton( self.buttonframe2, "About PriFi", None, self.about )
      quitbutton = self.createButton( self.buttonframe2, "Quit", "noname", self.quit )
      
      self.infow = StringVar()
      self.infolabel = Label(self.frame2, font = "Arial 15", bg=bgcolor, fg = fgcolor, textvariable=self.infow )
      self.infolabel.pack( side=TOP )
      
      
      # frame L holds alignment line view:
      self.frameL = Frame( self )
      self.frameL.config( bg = bgcolor)
      
      self.alviewL = AlignmentviewerLines( self.frameL, self, bgcolor, buttoncolor, primercolors )
      self.frameL.pack( fill=X, expand=YES )



      # frame 1 holds alignment view:
      self.frame1 = Frame( self )
      self.frame1.config( bg = bgcolor)
      self.alview = Alignmentviewer( self.frame1, bgcolor, buttoncolor, primercolors  )
      self.frame1.pack( fill=X, expand=YES )
      # the ordinary viewer needs a reference to the line viewer:
      self.alview.lineview = self.alviewL 
      
      self.shownprimers = [] # list of primers currently shown
       
      # frame 3 holds primer reports:
      self.frame3 = Frame(self, bg=bgcolor ) # for primer reports
      self.frame3.pack( expand=NO)
         
      # create N primer report windows:

      self.report = []
      ii = 0
      for i in range(cf.PrimerPairSuggestions):
         frame = Frame(self.frame3, bg=bgcolor )
         label = Label(frame, text="Primer set %d"%(i+1), fg=fgcolor, bg=bgcolor, font = "Arial 19"   )
         label.pack()
         window = Text(frame, width = W / 33, height = 28, bg="Black", font = "Arial 13", state=DISABLED, fg = primercolors[ii][1]  )
         window.pack(pady=5)

         #button = self.createButton(frame, "Show", "show%d"%i, None, None )
         # bind this way to that the button name can be retrieved:
         #button.bind( "<1>", self.showprimers ) 
         #button.pack()
         frame.pack( expand = NO, side=LEFT, pady=10 )
         self.report.append((window))
         self.primerpairs = []
         ii += 1


      # load documentation file:
      try:
         fin = open ( cf.DocFilename )
         self.doctext = fin.read()
         fin.close()

      except:
         self.doctext = None

      # if an alignment file was given, load it:
      if aln and aln[-4:].lower()=='.aln':
         self.loadAlignmentFile( aln )
         self.filepath = aln


   
      if exists( cf.ParameterFile ):
         # parameters have been saved previously, load these values
         # (but only load them once)
         
         self.loadParameters()

      # we only use the working directory once, after first load we just load
      # from where we last loaded. 
      self.firstload = 1 
      
      self.batchpath = None # path to dir of last batch run

      
   
   def handleMouseClick( self, event ):
      canvas = event.widget
      x = float(canvas.canvasx(event.x))
      
      self.alview.scrollTo( x / self.alviewL.l )
      

   

   def showprimers( self, n ): #event ):
      """n is the suggestion number"""
      
      if self.primerpairs == []:
         self.infow.set( "No primers to show!" )
         return
      
      self.infow.set( "" )
      #n = int(event.widget.winfo_name()[-1]) # last character of the name is an integer
      
      if n >= len(self.primerpairs) or n in self.shownprimers:

         return
      
      #fw = " "*self.primerpairs[n][1] + self.primerpairs[n][0]# + " "*self.primerpairs[n][2]
      #rv = " "*self.primerpairs[n][4] + self.primerpairs[n][3]# + " "*self.primerpairs[n][5]

      self.shownprimers.append(n)
      self.alview.addSequence( "FW%d / RV%d"%(n+1,n+1), "", n )



      # add the primer in the line view (get start and end indices of both primers):
      self.alviewL.addPrimerPair( self.primerpairs[n][1], self.primerpairs[n][2], self.primerpairs[n][4], self.primerpairs[n][5], n )

      self.alview.addSequenceData( "%s"%(" "*self.primerpairs[n][1]), -1 ) # no bg color
      self.alview.addSequenceData( self.primerpairs[n][0], n ) # first primer
      spaces = self.primerpairs[n][4] - self.primerpairs[n][1] - len(self.primerpairs[n][0])
      self.alview.addSequenceData( "%s"%(" "*spaces), n ) # with bg color

      # second primer:
      self.alview.addSequenceData( reverse_and_complement( self.primerpairs[n][3] ), n )
      



   
   def createButton( self, parent, txt, Name, cmd=None, SIDE=LEFT, color=buttoncolor ):
      """create and return a button"""
      button = Button( parent, text=txt, bg=color, command=cmd, font = "Arial 15", name=Name )
      if SIDE == None:
         button.pack( padx=15, pady=10 )
      else:
         button.pack( padx=15, pady=10,side=SIDE )
      return button








   def showPrimerReports( self, primerpairs ):
         
         self.infow.set( "" )

         if self.shownprimers:
            # redraw alignment window so old primers don't remain:
            self.shownprimers = []
            self.alview.insertLoadedBioAlignment( )
            self.alviewL.insertLoadedBioAlignment( )
         c = 0
         while c<len(self.report) and c<len(primerpairs):
            for p1, p2, expltuple in primerpairs:

               #primer1, primer1 start and end indices, primer2, primer2 start and
               #stop indices:
               self.primerpairs.append( (p1, expltuple[1], expltuple[2],
                                         p2, expltuple[3], expltuple[4]) )
               self.report[c].config( state=NORMAL )
               self.report[c].delete(1.0, END)
               self.report[c].insert( INSERT, expltuple[0] )

               self.report[c].config( state=DISABLED )
               self.showprimers(c)
               c += 1

         while c<len(self.report):
            self.report[c].config( state=NORMAL )
            self.report[c].delete(1.0, END)
            self.report[c].config( state=DISABLED )
            c+=1

         if len(primerpairs) > 0:
            # set scrollbars so that first suggestion's FW primer is centered:
            self.alviewL.setScrollBar(float(self.primerpairs[0][1]) / self.alview.l )
            self.alview.scrollTo( float(self.primerpairs[0][1]) / self.alview.l )





   def findprimers( self, showreports = 1 ):
      """By default, the found primers are reported in the GUI. Giving 0 as an
      argument turns the option off."""

      if showreports:
         self.infow.set("working ...")
         self.infolabel.update()#
      
      
      

      self.primerpairs = [] # the current primer pairs shown
      t = self.alview.getAlignment()
      if t:
         
         try:
            primerpairs = primerfinder.findprimers( 0, t[0], t[1], t[2] )
         except ValueError, e:
            if cf.INTRONS == 'yes':
               self.infow.set( e )
            else:
               self.infow.set( "Error, unexpected marker in the alignment" )
            return 0 # error
         
         if not primerpairs:
            if showreports:
               self.showPrimerReports( [] )
               self.infow.set( "No valid primers found!" )
            return 0 # error


         # build summary tuple representing the best primer pair
         # (recall that primerpairs is a list of primer tuples, where
         # each tuple has this form: p1, p2, (explanation tupel), and
         # the explanation tupel's first entry is a string whose LAST
         # field is the primer pair's score :)
         #
         # put score first so that the list in which this tupel is put
         # can be sorted by it:
         fields = primerpairs[0][2][0].split()
         
         if cf.INTRONS == 'yes':
         
            stup = [ float(fields[-1]), # score
                     primerpairs[0][0],  # fw primer
                     primerpairs[0][1],  # rv primer
                     "%4.1f"%primerpairs[0][2][5], # fw Tm
                     "%4.1f"%primerpairs[0][2][6], # rv Tm
                     primerpairs[0][2][2] - primerpairs[0][2][1], # fw length
                     primerpairs[0][2][4] - primerpairs[0][2][3], # rv length
                     int(fields[25]), # product length
                     int(fields[28]), # fw pre-intronlength
                     int(fields[30]), # rv pre-intronlength
                     self.filepath+'.rep' ] # analysis file
         else:
            stup = [ float(fields[-1]), # score
                     primerpairs[0][0],  # fw primer
                     primerpairs[0][1],  # rv primer
                     "%4.1f"%primerpairs[0][2][5], # fw Tm
                     "%4.1f"%primerpairs[0][2][6], # rv Tm
                     primerpairs[0][2][2] - primerpairs[0][2][1], # fw length
                     primerpairs[0][2][4] - primerpairs[0][2][3], # rv length
                     int(fields[24]), # product length
                     self.filepath+'.rep' ] # analysis file
            
                  
         if showreports:
            self.showPrimerReports( primerpairs )

         
         try:
            primerfinder.writePrimersToFiles( self.filepath, primerpairs )

            
            if not showreports:
               return stup

            
            filedir = ""
            file = self.filepath
            uds = 'Files written'
            sindex = self.filepath.rfind('/')
            if sindex != -1:
               filedir = self.filepath[:sindex+1]
               file = self.filepath[sindex+1:]
               uds += ' in directory %s'%filedir
            uds += ': report file %s.rep and primer files'%file
               
            uds += '\n%s.fw1, %s.rv1'%(file, file )
            if len(primerpairs) > 1:
               uds += '\n%s.fw2, %s.rv2'%(file, file )
               if len(primerpairs) > 2:
                  uds += '\netc.'
            self.infow.set(uds)
         except Exception, e:
            self.infow.set("Couldn't write primers to files %s"%e)
         
         return stup
      
      else:
         # no alignment loaded
         self.infow.set( "Load or create an alignment first" )





   def about( self ):
      """Pop up a window which explains how PriFi works."""

      if not self.doctext:
         self.infow.set("Couldn't load documentation file" )
         return
      
      if self.docframe:
         self.docframe.destroy()

         
      self.docframe = Toplevel( self, bg=bgcolor )  
      #docframe.pack( expand=NO, fill=NONE, pady=10 )
      
      doc = ScrolledText.ScrolledText(self.docframe, width=80, height=30, bg="White", font="Arial 15", state=DISABLED, fg="Black" )
      doc.pack(padx=15, pady=15)
      doc.config(state=NORMAL)
      doc.insert(INSERT, self.doctext)
      doc.config(state=DISABLED)
      okb = self.createButton( self.docframe, "OK", None, self.removeDocWindow,LEFT )
            

   def removeDocWindow( self ):
      self.docframe.destroy()

   




   def loadPrifiFile( self ):
      """Load primer sets of PriFi analysis previously completed."""

      path = tkFileDialog.askopenfilename( initialdir = self.batchpath, title="Load PriFi analysis file", filetypes=[(".rep files", "*.rep")])

      if path:

         if path[-4:] == '.rep':

            primerpairs = []
            try:
               fin = open( path )
               fin.readline()
               aln = fin.readline().strip() # get alignment filename from report file

               # deduce the exact name of the alignment file:

               # NEW VERSION: .aln saettes bare efter navnet
               alf = aln+'.aln'
               if not exists(alf):
                  ri = aln.rfind('.')
                  if ri != -1:
                     alf = aln[:ri]+'.aln'

               if not exists(alf):
                  alf = '[unknown]'
            

               self.filepath = alf #path to ALIGNMENT file, not a fasta file
         
                              
               # parse report file:
               indexpattern = re.compile( '(\d+)-(\d+) / (\d+)-(\d+)' )
               
               while 1:
                  line = fin.readline()
                  if line == '':
                     break

                  if line.startswith( 'Primer set' ):
                     primer = []
                     indexmatch = indexpattern.search( line )
                     
                     p1start = int(indexmatch.group(1))
                     p1end = int(indexmatch.group(2))
                     p2start = int(indexmatch.group(3))
                     p2end = int(indexmatch.group(4))
                     fin.readline()

                     txt = []
                     line = fin.readline()
                     while not line.startswith( '---' ):
                        txt.append(line)
                        line = fin.readline()
                        
                     # now txt is a list of lines holding the report; we get Tm's and
                     # primer strings from this list:
                     
                     fw = txt[0].strip()[6:]
                     rv = reverse_and_complement( txt[1].strip()[6:] )

                     fields = txt[3][5:].split('/') # this line holds the Tm's
                     
                     tm1 = float(fields[0])
                     tm2 = float(fields[1])

                     # staple report :) :
                     txt = "".join( txt )

                     tupel = (fw, rv, (txt, p1start, p1end, p2start, p2end, tm1, tm2 ) )

                     

                     primerpairs.append( tupel  )
                     
               fin.close()
               self.loadAlignmentFile( self.filepath )
               self.primerpairs = []
               self.showPrimerReports( primerpairs )
               self.infow.set("Done.")
               self.infolabel.update()#
               
            except Exception, exc:
               self.infow.set( "Couldn't show PriFi analysis; %s"%(exc) )
               return
            
            
      
         else:
            # not a .rep file
            showinfo( ".rep files", "Filename has to have .rep extension", parent=self )



   def doAlignment( self, path = None, drawalignment = 1 ):
      """Do alignment from file chosen in dialog window. The path may also
      be given as an argument. By default, the alignment is srawn/written in the
      GUI; this option can be turned off for batch runs."""

      if not path:
         path = tkFileDialog.askopenfilename( title="Load fasta file", filetypes=[(".fasta files", "*.fasta"), ("all files", "*")])

      if path:

         if drawalignment:
            self.infow.set("working ...")
            self.infolabel.update()#
      
         try:
            self.alview.allseq,self.alview.summary,self.alview.l=primerfinder.doAlignment(path)
            if drawalignment:
               self.alviewL.setAlignment( self.alview.allseq, self.alview.summary, self.alview.l )
            
               
         except Exception, exc:
            if drawalignment:
               self.infow.set( "Couldn't do clustalw alignment on file %s"%path )
            return 1 # error

         # deduce the exact name of the alignment file:

         # NEW VERSION: .aln saettes nu bare efter filnavnet.
         alf = path+'.aln'

         #if not exists(alf):
         #  ri = path.rfind('.')
         #  if ri != -1:
         #    alf = path[:ri]+'.aln'
               
         if not exists(alf):
            #alf = '[unknown]'
            alf += '(?)'

         self.filepath = alf #path to ALIGNMENT file, not a fasta file
         
         if not drawalignment:
            return 0

         
         for w in self.report:
            self.alview.insertText( w, "", clear=1 )

         self.infow.set( "Alignment succeeded; alignment saved in file %s"%alf )
         self.shownprimers = []
         self.master.title( "  --  PriFi  --  "+ path[path.rfind('/')+1:])
         self.alview.insertLoadedBioAlignment()
         self.alviewL.insertLoadedBioAlignment()
         
         
      





   def fastaBatch( self ):
      """Retrieve a directory name from the user and do alignment and primer
      finding on all fasta files (.fasta) found in it."""

      #path = tkFileDialog.askopenfilenames( initialdir = cf.WorkingDir, title="Choose one or several .fasta files", filetypes=[(".aln files", "*.aln")])
      
      path = tkFileDialog.askdirectory( initialdir = cf.WorkingDir, title="Choose a directory holding Fasta files")

      
      if path:
         self.batchpath = path
         
         self.infow.set("working (skipping files for which primers were found already ...")
         self.infolabel.update()#
      
         for w in self.report:
            self.alview.insertText( w, "", clear=1 )

         self.shownprimers = []
         self.master.title( "  --  PriFi  --  ")

         self.alview.insertText( self.alview.text1, "", clear = 1 )
         #setSequenceData( "" )
         self.alviewL.clear()
         
         allallfiles = listdir( path )

         # DO redo already done alignments but SKIP files for
         # which primers were found:
         allfiless = []
         primerreports = []
         for af in allallfiles:
            if '.aln' in af or '.dnd' in af or 'prifi' in af:
               if af.endswith('.aln.rep'):
                  primerreports.append(af) #new version: .aln just added to filename
               continue
            allfiless.append(af)

         # allfiless inderholder nu alle fasta-filer
         # primerreports indeholder alle rep-navne uden .aln.rep
         
         allfiles = []
         for af in allfiless:
            afrepnavn = af+'.aln.rep' # from old version: [:af.find('.')]
            if afrepnavn in primerreports:
               primerreports.remove(afrepnavn)
               print 'skipping', af
               continue
            allfiles.append(af)
            
         
         missedalignments = 0
         missedprimers = 0
         pos = len(allfiles)
         pp = 0
         excel = [] # list of summary tuples, one per found primer set 
         for f in allfiles:

            if f.endswith( '.dnd' ) or f.endswith( '.aln' ):
               missedalignments += 1
               continue
            
            self.filepath = join( path, f)
            if self.doAlignment( self.filepath, 0): # 0: don't write alignments
               missedalignments += 1
            else:
               # alignment was successful, now find primers:
               self.infow.set("analyzing file %s"%f)
               self.infolabel.update()#

               res = self.findprimers( 0 ) # 0:don't write primer reports in the gui
               if not res:
                  missedprimers += 1
               else:
                  # res is a summary tupel of the best primer set
                  excel.append(res)
                  pp += 1

         # now sort all the best primer pairs by score (first entry in each tupel):
         excel.sort()
         excel.reverse()


         # write summary tupels to .xl file:
         xlfil = path+"/prifibatch.xl"
         if exists(xlfil):
            xlexists = True
         else:
            xlexists = False
            
         stuf = open ( xlfil, 'a' )

         if not xlexists:
            print >> stuf, "Score\tFW\tRV\tTm(FW)\tTm(RV)\tlen(FW)\tlen(RV)\tProduct length\tFW/intron distance\tRV/intron distance\tPriFi analysis file"
         for tup in excel:
            for en in tup[:-1]:
               stuf.write( "%s\t"%en )
            stuf.write( "%s\n"%tup[-1] )
         stuf.close()
         
            
         self.infow.set("Found %d files, performed %d successful alignments, found %d primer sets. See sheet of best primer sets in \n%s"%(pos, pos-missedalignments, pp, (path+"/prifibatch.xl")))
         self.infolabel.update()
      
         

      

   def loadAlignment( self, path = None ):
      """Retrieve a path to an .aln file from the user and load this alignment"""

      if not path:
         print cf.WorkingDir
         if self.firstload == 1:
            path = tkFileDialog.askopenfilename( initialdir = cf.WorkingDir, title="Load .aln file", filetypes=[(".aln files", "*.aln"), ("all files", "*")])
            self.firstload = 0
         else:
             path = tkFileDialog.askopenfilename( title="Load .aln file", filetypes=[(".aln files", "*.aln")])
      
      if path:

         self.infow.set("working ...")
         self.infolabel.update()#
      
         for w in self.report:
            self.alview.insertText( w, "", clear=1 )
         self.loadAlignmentFile( path )
         self.filepath = path
      #else:
      #   showinfo( "Missing alignment file", "Couldn't find alignment file %s"%path, parent=self )



   def loadAlignmentFile( self, path ):
      self.infow.set( "" )
      self.shownprimers = []
      self.alview.loadAlignment( path )
      self.alviewL.setAlignment( self.alview.allseq, self.alview.summary, self.alview.l)
      self.alviewL.insertLoadedBioAlignment()
      self.master.title( "  --  PriFi  --  "+ path[path.rfind('/')+1:])
      

   def quit( self ):
      sys.exit(0)
   

   def configure( self ):
      PrimerConfigureDialog( self, "Configure", "#7896a0" )

         
   def loadParameters( self ):
      """load previously saved configuration parameters"""
      
      try:
         fin = open( cf.ParameterFile )
         lines = fin.readlines()
         fin.close()
         i = 0
         for n, v, d, m1, m2, default in cf.Parameters:
            # the entries in this list are in the same order as the lines
            # in the file.

            txt = lines[i].strip()
            i += 1
            if not cf.parseStringAndAssignToParameter( v, m1, m2, txt ):
               raise Exception
            
         self.infow.set( "Parameter values loaded" )

      except Exception, e:
         self.infow.set( "Parameter values couldn't be loaded, using default values" )

         














def main():
   if len(sys.argv) == 2:
      gui = prifiGUI( sys.argv[1] )
   else:      
      gui = prifiGUI()

   gui.mainloop()
   

if __name__ == "__main__":
   main()
