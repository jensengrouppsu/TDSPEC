#!/usr/bin/env python

from __future__ import print_function

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
from modules import dressed_spectroscopy as spec
from modules.euler_rotation import rotation_matrix_3D as rotate
import pickle
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from chem.constants import HART2WAVENUM, NM2WAVENUM
import os

__version = '2.0'

class tdspec_input():
    '''TDSpec input class.
    Holds all options for any TDSpec calcultion.'''

    def __init__(self):
        self.runtypes =     {   'Absorbance': 'Abs',
                                'Raman': 'Raman',
                                'Anti-Stokes Raman': 'ASRaman',
                                'Two photon absorbance': 'TwoAbs',
                                'HyperRaman': 'HyperRaman',
                                'Anti-Stokes HyperRaman': 'ASHyperRaman',
                                'Sum-Frequency Generation': 'SFG',
                                'Fluorescence': 'Fluor',
                                'Circular Dichroism': 'CD',
                                'Raman Optical Activity': 'RVROA',
                                'Second HyperRaman': 'SecondHyperRaman'    }
        self.Runtype            = 'Raman'
        self.lfreq              = 532.0
        self.gam2               = 0.0
        self.sshift             = 0.0
        self.Datafile           = '_TDSPEC.data'
        self.Width              = 7.0
        self.lnormalize         = False
        self.lherztell          = False
        self.latermoff          = False
        self.lhtints            = False
        self.lhtpref            = False
        self.lvoverlap          = False
        self.lrhrsb2            = False
        self.lhpolavg           = False
        self.lunitsphere        = False
        self.lnoRes             = False
        self.lSFGFull           = False
        self.IRUV               = False
        self.lerror             = False
        self.lruntype           = False
        self.printlevel         = 'None'
        self.tstep              = 500
        self.tstop              = 0.300
        self.npoints            = 2000
        self.abspts             = 10000
        self.cbguess            = 9000
        self.FluorState         = 1
        self.freq_screen        = 8000.0
        self.fc_screen          = 0.0001
        self.solvmod            = 'Simple'
        self.temperature        = 300.
        self.excnum             = 1
        self.theta              = 0.
        self.psi                = 0.
        self.orientation        = 1
        self.lintensity         = False
        self.rrs_pref           = None

    def update(self, entry):
        '''Updates values based on input options.'''
        self.Runtype = entry['calctype'].value
        self.Width   = entry['width'].value
        self.tstep   = entry['tstep'].value
        self.tstop   = entry['tstop'].value
    
        if self.Runtype == 'Raman':
            self.lfreq  = entry['exfreq'].value
        elif self.Runtype == 'Absorbance':
            self.lfreq  = entry['frange'].value[1]
            self.nabpts = int(1. + (NM2WAVENUM(entry['frange'].value[0]) - NM2WAVENUM(self.lfreq)) / 20.)
        elif self.Runtype == 'Fluorescence':
#            self.lfreq    = entry['frange'].value[0]
            self.lfreq  = entry['exfreq'].value
            self.FluorState = entry['fstate'].value
            self.nabpts = int((NM2WAVENUM(self.lfreq) - NM2WAVENUM(entry['frange'].value[1])) / 2.)
        else:
            print ('WARNING: Only "Raman", "Absorption", and "Fluorescence" calculations')
            print ('         are currently implemented in the GUI.')
            print ('         NO OTHER METHODS IMPLIMENTED YET.')
            return

        self.sshift = entry['sshift'].value
        self.gam2   = entry['gamma2'].value

        if entry['smodel'].value == 'I. Simple':
            self.solvmod = "Simple"
        elif entry['smodel'].value == 'II. High temperature':
            self.solvmod = "HighTemp"
        elif entry['smodel'].value == 'III. General temperature':
            self.solvmod = "GenTemp"
        else:
            print ('Unrecognized solvent model')
        if entry['ht'].value == 'I. A term':
            self.latermoff = False
            self.lherztell = False
        elif entry['ht'].value == 'II. A+B terms':
            self.latermoff = False
            self.lherztell = True
        elif entry['ht'].value == 'III. B term':
            self.latermoff = True
            self.lherztell = True

        printlevel = entry['ppol'].value
        if printlevel == 'alpha':
            self.printlevel = 'Pol'
        elif printlevel == 'alpha + A + C':
            self.printlevel = 'Atensor'
        elif printlevel == 'alpha + G':
            self.printlevel = 'Gtensor'
        elif printlevel == 'alpha + A + G + C + D':
            self.printlevel = 'AllTensors'
        else:
            self.printlevel = 'None'

def refresh_options(ctype):
    '''Changes the options displayed based on selected calctype.'''

    for key in entry.entry.keys():
        if entry[key].copt == 'all' or ctype in entry[key].copt:
            entry[key].grid()
        else:
            entry[key].ungrid()

# The "Calculate" button funciton
def execute():
    '''Generates a TDSpec input file from the given data and runs
    the TDSpec calculation.'''

    # Update TDSpec input object with entered options
    ob = tdspec_input()
    ob.update(entry)

    # Return if method not implemented yet
    if ob.Runtype!='Raman' and ob.Runtype!='Absorbance' and ob.Runtype!='Fluorescence':
        popup_message('Only "Raman", "Absorption", and "Fluorescence" calculations'
                      'are currently implemented in the GUI.\n'
                      'NO OTHER METHODS IMPLEMENTED YET.', title='BETA WARNING')
        return

    # Copy the datafile to this directory
    if len(entry['datafile'].value) < 5:
        popup_message("Plese select a valid datafile.")
        return
    string = "cp {0} {1}".format(entry['datafile'].value, ob.Datafile)
    os.system(string)

    # Open input file for writing
    infile = open('_TDSPEC.inp', 'w')

    print ("Runtype {0}".format(ob.runtypes[ob.Runtype]),   file=infile)
    print ("Freq {0:7.2f}".format(ob.lfreq),                file=infile)
    print ("Shift {0:12.8f}".format(ob.sshift),             file=infile)
    print ("Gamma2 {0:10.4f}".format(ob.gam2),              file=infile)
    print ("Datafile {0}".format(ob.Datafile),              file=infile)
    if ob.lherztell: print ("HerzTell Yes",                 file=infile)
    if ob.latermoff: print ("ATermOff",                     file=infile)
    if ob.Runtype == 'Raman':
        print ("Width {0:8.3f}".format(ob.Width),           file=infile)

    print ("Tstop {0:6.3f}".format(ob.tstop),               file=infile)
    print ("Tstep {0}".format(ob.tstep),                    file=infile)
    if ob.Runtype == 'Absorbance' or ob.Runtype == 'Fluorescence':
        print ("NABSpts {0}".format(ob.nabpts),             file=infile)
    if ob.Runtype == 'Fluorescence':
        print ("FluorState {0}".format(ob.FluorState),      file=infile)

    print ("SolvMod {0}".format(ob.solvmod),                file=infile)
    print ("Temp {0}".format(entry['temperature'].value),   file=infile)
    print ("ExcNum {0}".format(entry['excnum'].value),      file=infile)

    if ob.printlevel != 'None':
        print ("Print{0}".format(ob.printlevel),            file=infile)
    if entry['intensity'].value:
        print ("Int",                                       file=infile)

    # Close input file
    print ("End",                                           file=infile)
    infile.close()

    # Run TDSpec calculation
    string = "rm _TDSPEC.pol 2> /dev/null"
    os.system (string)
    string = "tdspec < _TDSPEC.inp > _TDSPEC.out"
    os.system (string)

    # Update Plot with data
    try:
        graph.deiconify()
    except (Tk.TclError,NameError):
        generate_graph()
        graph.deiconify()
    x, y = np.loadtxt("_TDSPEC.out", unpack=True)
    if entry['normalize'].var.get():
        m = y.max()
        y /= m
    line.set_ydata(y)
    line.set_xdata(x)
    if ob.Runtype == 'Raman':
        ax.set_xlabel("wave number / cm$^{-1}$", fontsize=16)
    else:
        ax.set_xlabel("wavelength / nm", fontsize=16)
    ax.set_xlim((x[-1],x[0]))
    if ob.Runtype == 'Fluorescence':
        ax.set_xlim((entry['frange'].value[0],entry['frange'].value[1]))
    ax.set_ylim((0.,y.max()*1.1))
    canvas.draw()

def getdatafile():
    '''Calculates Spectrum based on input.'''
    import tkFileDialog
    filetypes = [("TDSpec datafile","*.data"), ("All","*")]
    entry['datafile'].filename = tkFileDialog.askopenfilename(filetypes=filetypes)
    entry['datafile'].val1.delete(0, Tk.END)
    entry['datafile'].val1.insert(0, entry['datafile'].filename)

def saveoutputfile():
    '''Saves calculated data to disk.'''
    import tkFileDialog
    filetypes = [("TDSpec output","*.out")]
    outputfilename = tkFileDialog.asksaveasfilename(filetypes=filetypes)
    if outputfilename[-4:] != ".out":
        pass
    else:
        # Check if calculation has been done first
        if os.path.isfile("_TDSPEC.out"):
            print ("Output file saved as: {0}".format(outputfilename))
            string = "mv _TDSPEC.out {0}".format(outputfilename)
            os.system(string)
            polfile = outputfilename[:-4]+".pol"
            string = "mv _TDSPEC.pol {0} 2> /dev/null".format(polfile)
            os.system(string)
        else:
            print ("Run calculation first.")

def close_window():
    '''Closes the windows.'''
    try:
        root.destroy()
    except Tk.TclError:
        pass
    try:
        graph.destroy()
    except Tk.TclError:
        pass
    except NameError:
        pass

class options():
    '''A class for defining the different types of entries in GUI.'''

    def entry(self, root, pos, justify='right', **kwargs):
        '''Textbox for words or numbers.'''
        self.val1 = Tk.Entry(root, width=41, justify=justify)
        if 'default' in kwargs:
            self.val1.insert(0, kwargs['default'])
        self.val1.grid(row=pos,column=2,sticky=Tk.E+Tk.W,columnspan=2)

    def double_entry(self, root, pos, justify='right', **kwargs):
        '''Two textboxes, used for number ranges.'''
        self.val1 = Tk.Entry(root, justify=justify)
        if 'default1' in kwargs:
            self.val1.insert(0, kwargs['default1'])
        self.val1.grid(row=pos,column=2,sticky=Tk.W)
        self.val2 = Tk.Entry(root, justify=justify)
        if 'default2' in kwargs:
            self.val2.insert(0, kwargs['default2'])
        self.val2.grid(row=pos,column=3,sticky=Tk.E)

    def optionmenu(self, root, pos, default, lst, cmd):
        '''A dropdown list of options.'''
        self.vals = Tk.StringVar(root)
        self.vals.set(default)
        self.val1 = apply(Tk.OptionMenu, (root, self.vals) + tuple(lst), {'command': cmd})
        self.val1.grid(row=pos,column=2,sticky=Tk.E+Tk.W,columnspan=2)

    def datafilebox(self, root, pos, **kwargs):
        '''Box to read in datafile.'''
        self.filename   = ''
        self.val1       = Tk.Entry(root)
        self.val2       = Tk.Button(root, text='Browse', command=getdatafile, width=10)
        self.val1.grid(row=pos,column=2,sticky=Tk.W+Tk.E)
        self.val2.grid(row=pos,column=3,sticky=Tk.E)

    def checkbox(self, root, pos, **kwargs):
        '''Checkbox for True/False options.'''
        self.var        = Tk.IntVar()
        self.val1       = Tk.Checkbutton(root, variable=self.var)
        self.val1.grid(row=pos,column=2)

    @property
    def value(self):
        '''Get the value of this options object.'''
        if self.type == 'entry' or self.type == 'double_entry':
            val = self.val1.get()
            try:
                val = int(val)
            except ValueError:
                try:
                    val = float(val)
                except ValueError:
                    pass
            if self.type == 'double_entry':
                val2 = self.val2.get()
                try:
                    val2 = int(val2)
                except ValueError:
                    try:
                        val2 = float(val2)
                    except ValueError:
                        pass
                return [val, val2]
            else:
                return val
        elif self.type == 'datafilebox':
            return self.filename
        elif self.type == 'optionmenu':
            return self.vals.get()
        elif self.type == 'checkbox':
            if self.var.get():
                return True
            else:
                return False

    def grid(self):
        '''Show option.'''
        self.label.grid()
        if self.type == 'entry':
            self.val1.grid()
        elif self.type == 'double_entry':
            self.val1.grid()
            self.val2.grid()
        elif self.type == 'optionmenu':
            self.val1.grid()

    def ungrid(self):
        '''Remove option.'''
        self.label.grid_remove()
        if self.type == 'entry':
            self.val1.grid_remove()
        elif self.type == 'double_entry':
            self.val1.grid_remove()
            self.val2.grid_remove()
        elif self.type == 'optionmenu':
            self.val1.grid_remove()
        elif self.type == 'checkbox':
            self.val1.grid_remove()

    def __init__(self, root, label, type, pos, copt='all', **kwargs):
        self.type   = type
        self.labels = Tk.StringVar()
        self.labels.set(label)
        self.label  = Tk.Label(root, textvariable=self.labels, width=20, anchor=Tk.W)
        self.label.grid(row=pos,column=1)
        self.type   = type
        self.copt   = copt
        if self.type == 'entry': self.entry(root, pos, **kwargs)
        if self.type == 'double_entry': self.double_entry(root, pos, **kwargs)
        if self.type == 'optionmenu': self.optionmenu(root, pos, **kwargs)
        if self.type == 'datafilebox': self.datafilebox(root, pos, **kwargs)
        if self.type == 'checkbox': self.checkbox(root, pos, **kwargs)
        if 'grid' in kwargs and not kwargs['grid']: self.ungrid()

class entries():
    '''A class for storing all entries for the TDSpec GUI.'''
    def __init__(self):
        self.entry = {}
    def add(self, name, root, label, type, pos, **kwargs):
        self.entry.update({name: options(root, label, type, pos, **kwargs)})
    def __getitem__(self,name):
        return self.entry[name]

def generate_graph():
    '''Creates or refresh the output window.'''
    global graph, fig, ax, line, canvas

    graph       = Tk.Tk()
    graph.wm_title("TDSpec Output")
    fig         = Figure()
    ax          = fig.add_subplot(111)
    line,       = ax.plot(range(10),'b-',lw=2)
    canvas      = FigureCanvasTkAgg(fig,master=graph)
    canvas.show()
    canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
    toolbar = NavigationToolbar2TkAgg(canvas, graph)
    toolbar.update()
    graph.withdraw()

def popup_message(message, title='Error'):
    '''Creates a popup window with an error message.'''

    popup = Tk.Toplevel()
    popup.title(title)

    text = Tk.Message(popup, text=message, width=800, justify='center')
    text.pack()

    button = Tk.Button(popup, text='Close', command=popup.destroy)
    button.pack()

# ############
# MAIN PROGRAM
##############
def main():

    global root, entry

    # TCL main window
    root = Tk.Tk()
    root.wm_title("TDSpec v"+__version)
    
    # Set all entries
    entry = entries()

    # Calctype and datafile
    entry.add('datafile', root, 'TDspec datafile', 'datafilebox', 1)
    entry.add('calctype', root, 'Calculation Type', 'optionmenu', 2, default='Raman', lst=['Raman',
              'Absorbance', 'Fluorescence', 'Anti-Stokes Raman', 'Two-Photon Absorbance',
              'HyperRaman', 'Anti-Stokes HyperRaman', 'Sum-Frequency Generation',
              'Circular Dichroism', 'Second HyperRaman'], cmd=refresh_options)

    # Frequencies (range and excitation)
    entry.add('frange', root, 'Frequency range (nm)', 'double_entry', 3, default1=300.0,
              default2=800.0, grid=False, copt=['Absorbance','Fluorescence'])
    entry.add('exfreq', root, 'Excitation frequency (nm)', 'entry', 4, default=532.0,
              copt=['Raman', 'Fluorescence'])

    # Solvent and other broadening parameters
    entry.add('width', root, 'Lorentzian width (cm-1)', 'entry', 5, default=7.0, copt='Raman')
    entry.add('gamma2', root, 'Gamma2 (cm-1)', 'entry', 6, default=0.0)
    entry.add('sshift', root, 'Solvent shift (a.u.)', 'entry', 7, default=0.0)
    entry.add('smodel', root, 'Solvent model', 'optionmenu', 8, default='I. Simple', cmd=None,
              lst=['I. Simple', 'II. High temperature', 'III. General temperature'])
    entry.add('temperature', root, 'Temperature (K)', 'entry', 9, default=300.0)

    # FC and HT options
    entry.add('ht', root, 'Hertzberg-Teller', 'optionmenu', 10, default='I. A term', lst=['I. A term',
              'II. A+B terms', 'III. B term', 'IV. B1 term (only RHRS)'], cmd=None)
    entry.add('fstate', root, 'Fluorescence state', 'entry', 11, default=1, grid=False, copt='Fluorescence')
    entry.add('excnum', root, 'Overtones excitation nr.', 'entry', 12, default=1)

    # Integration paramters
    entry.add('tstep', root, 'Time Step', 'entry', 20, default=500, grid=False,
              copt=['Absorbance','Fluorescence'])
    entry.add('tstop', root, 'Time Stop', 'entry', 21, default=0.300, grid=False,
              copt=['Absorbance','Fluorescence'])

    # Print statements
    entry.add('ppol', root, 'Print level', 'optionmenu', 90, default='None', cmd=None, copt='Raman',
              lst=['None', 'alpha', 'alpha + A + C', 'alpha + G', 'alpha + A + G + C + D'])
    entry.add('intensity', root, 'Sticks spectrum', 'checkbox', 30, copt='Raman')
    entry.add('normalize', root, 'Normalize spectrum', 'checkbox', 31, copt='Raman')
    
    # Buttons:
    # Calculate
    b_execute   = Tk.Button(root, text='Calculate', command=execute, width=15)
    b_execute.grid(row=100,column=1)
    # Save
    b_save      = Tk.Button(root, text='Save', command=saveoutputfile, width=15)
    b_save.grid(row=100,column=2)
    # Close
    b_exit      = Tk.Button(root, text='Exit', command=close_window, width=15)
    b_exit.grid(row=100,column=3)
    
    # Cycle until windows are closed
    Tk.mainloop()
    close_window()
    
    # Remove files no longer needed
    string = "rm _TDSPEC.inp _TDSPEC.out _TDSPEC.data _TDSPEC.pol 2> /dev/null"
    os.system (string)

if __name__ == '__main__':
    main()
