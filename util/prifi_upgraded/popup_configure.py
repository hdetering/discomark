from Tkinter import *
import os
from tkMessageBox import *

class ConfigureDialog(Toplevel):

    def __init__(self, parent, title = None, bgcolor="LightBlue"):

        Toplevel.__init__(self, parent)
        self.transient(parent)
        
        if title:
            self.title(title)

        # entry explanations displayed when the parameter name button is pressed:
        self.ee = {}
        
        self.bgcolor = bgcolor
        self.config( bg = bgcolor )
        
        self.parent = parent

        self.result = None

        body = Frame(self, bg=bgcolor)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks


    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        
        box = Frame(self, bg=self.bgcolor)

        w = Button(box, bg=self.bgcolor, text="OK", width=10, command=self.ok)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, bg=self.bgcolor, text="Reset", width=10, command=self.reset)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, bg=self.bgcolor, text="Save", width=10, command=self.saveParameters)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, bg=self.bgcolor, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        self.bind("&lt;Return>", self.ok)
        self.bind("&lt;Escape>", self.cancel)

        box.pack()


    def reset( self ):    # should be overridden
        """reset all parameters"""
        pass


    
    def saveParameters( self ):    # should be overridden
        """save all parameters"""
        pass

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()




    def createEntryInGrid( self, master, txt, value, rowno, expl ):
        """ each label is here a button so that if you press it you can get information
        about the parameter (from a popup window arranged in the buttonpress method)"""

        nm = txt.lower() # window names can't start with uppercase letter
        b = Button(master, name=nm, bg = self.bgcolor, font = "Arial 15", text=txt)
        b.grid(row=rowno)
        b.bind( "<1>", self.buttonpress )
        self.ee[nm] = expl # save explanation for the entry under the button's name
      
        e1 = Entry(master, font = "Arial 15", bg = "White" )
        e1.insert( INSERT, value )
        e1.grid(row=rowno, column=1)
        return e1

   
    def buttonpress( self, event ):
        name = event.widget.winfo_name() # button's name
        
        # names are lower case versions of button text:
        txt = self.ee[name]
        
        showinfo( name.capitalize(), txt, parent=self )
      

    def body(self, master): # override
        pass


    def validate(self): # override (this method is called before window is closed)
        return 1 # values ok  


    def apply(self): # override
              pass

    
