import tkinter as tk
import copy
from tkinter import ttk, filedialog
import customtkinter
import csv
from plotly.graph_objs.layout import YAxis, XAxis, Margin
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from PIL import Image, ImageTk
import numpy as np
import matplotlib.patheffects as path_effects

# Window Params
W = 1900
H = 960

# Global Variables
Etypes = ['Free Energy', 'Relative Free Energy', 'Potential Energy', 'Relative Potential Energy', 'Enthalpy', 'Gibbs Free Energy']
# DEBUG FUNCTIONS
def combobox_callback(choice):
        print("combobox dropdown clicked:", choice)

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # DEFINE PARAMETERES:
        self.ngraph = 0
        self.nmechs = 1
        self.energylist = [0, 0]
        self.refs = [0, 0]
        self.S = [0, 0]
        self.RCoord = [0, 0]
        self.colorlist = []
        self.mechs = [0, 0]
        self.mechstypes = {}
        self.Units = [0, 0] # 0 for eV, 1 for kcal/mol, 2 for kJ/mol
        self.Titles = [0, 0]
        self.linedict = {}
        self.bardict = {}
        self.normenlist = [] #for normalised energies and plotting them

        #Frames
        self.LeftFrame = customtkinter.CTkFrame(self, width=W/3, height=2*H/3, border_color='red')
        self.LeftFrame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        self.RightFrame = customtkinter.CTkFrame(self, width=2*W/3, height=2*H/3, border_color='red')
        self.RightFrame.grid(row=0, column=1, padx=20, pady=20, sticky="nsew")
        self.TableFrame = customtkinter.CTkFrame(self, width=W, height=H/3)
        self.TableFrame.grid(row=1, column= 0, padx=20, pady=20, sticky='nsew', columnspan=2)


        self.title("m-PES")
        self.geometry('{}x{}'.format(W, H))
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.load_button = customtkinter.CTkButton(self.LeftFrame, text="Load file", command=lambda:self.readinput(self.ngraph, self.energylist, self.refs, self.RCoord, self.Titles, self.mechs))
        self.plot_button = customtkinter.CTkButton(self.LeftFrame, text="Plot PES", command=lambda:self.pltPES(self.graph_frame, self.ngraph, self.originalRC, self.Titles, self.originalEL, self.mechs, self.Units))
        self.load_button.grid(row=0, column=0,  padx=20, pady=20, sticky="w")
        self.plot_button.grid(row=3, column=0,  padx=20, pady=20, sticky="w")
        
        self.bg_image = Image.open("bg.jpg")  # Replace with your image path
        self.bg_image = self.bg_image.resize((int(W/3), int(2*H/3)))  # Resize if needed
        self.bg_image_tk = ImageTk.PhotoImage(self.bg_image)
        self.graph_frame = customtkinter.CTkFrame(self.RightFrame, width=W/3, height=2*H/3, border_color='blue')
        self.graph_frame.grid(row=0, column=1, padx=10, pady=10)
        self.gp_canvas = customtkinter.CTkCanvas(self.graph_frame, width=W/3, height=2*H/3)
        self.gp_canvas.pack(fill="both", expand=True)
        self.gp_canvas.create_image(0, 0, anchor="nw", image=self.bg_image_tk)





    def readinput(self, ngraph, energylist, refs, RCoord, Titles, mechs): # Read info from CSV
        #global energylist, refs, RCoord, Titles
        energydict = {}
        ref = []

        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("mPES files", "*.mpes")])
        print('FP', file_path)
        ma = {}
        if file_path:
            try:
                with open(file_path, 'r') as file:
                    reader = csv.reader(file)
                    T = next(reader)[0]
                    print('TITLE = ', T)
                    headers = next(reader)  
                    line = next(reader)  # Read gas reaction coords
                    for i in range(len(line)):
                        RC = [coord for coord in line[1:]]
                    U = line[0]
                    headers = next(reader) 
                    nPES = next(reader)  # Number of PES diagrams
                    nPES = int(nPES[0])
                    print(nPES)
                    for i in range(nPES):
                        energy = next(reader)
                        fenergy = [float(num) for num in energy[1:]]
                        energydict[energy[0]] = fenergy
                        ref.append(energy[0])
                        print(i)
                    refs[ngraph] = ref
                    Titles[ngraph] = T
                    RCoord[ngraph] = RC
                    energylist[ngraph] = energydict
                    self.Units[ngraph] = U.strip()
                    self.RCoord[ngraph] = RC
                    if file_path.endswith('.mpes'):
                        #mechs[ngraph] = dict.fromkeys((range(nPES)))
                        
                        print('READING MECHS')
                        headers = next(reader)
                        for i in range(nPES):
                            print('i', i)
                            me = next(reader)
                            print('me', me)
                            print(me[1].strip())
                            ma[me[0]] = Mech(color=me[1].strip(), linestyle=me[2].strip(), barstyle=me[3].strip())
                            print('ma', ma)
                        mechs[ngraph] = ma
                        print('ULO')
                    if not ma:
                        print("m is empty before calling createenergysets.")
                        self.createenergysets(ngraph, Titles[ngraph], refs[ngraph], mechs, energylist, self.Units, RCoord)
                    else:
                        self.createenergysets(ngraph, Titles[ngraph], refs[ngraph], mechs, energylist, self.Units, RCoord)
            except StopIteration as e:
                print(f"Reached the end of the CSV unexpectedly: {e}")
            except ValueError as e:
                print(f"Value error when converting data types: {e}")
            except Exception as e:
                print(f"Error loading file: {e}")

        print('Titles', Titles)
        print('Energylist', self.originalEL)
        print('Refs', refs)
        print('RCoord', self.originalRC)
        print('Units :', self.Units)
        print('MECHS', mechs)

    def normalize(self, ngraph, energylist, sel_ref):
        nlist = copy.deepcopy(energylist)
        if sel_ref.get() == 'no ref':
            pass
        elif sel_ref.get() == 'all zero':
            print('BEFORE el:', nlist)
            for key in nlist[ngraph]:
                nlist[ngraph][key] = [energy - nlist[ngraph][key][0] for energy in nlist[ngraph][key]]
            print('AFTER el:', nlist)
        else:
            print('BEFORE el:', nlist)
            S = nlist[ngraph][sel_ref.get()][0]
            print('S=', S)
            for key in nlist[ngraph]:
                nlist[ngraph][key] = [energy - S for energy in nlist[ngraph][key]]
            print('AFTER el:', nlist)
        return nlist

    def conversion(self, ngraph, energylist, uts, sel_conv):
        nlist = copy.deepcopy(energylist)
        u1 = Units_list.index(uts[ngraph])
        u2 = Units_list.index(sel_conv.get())
        print('nlist', nlist)
        print('u1 and u2', u1, u2)
        if u1 == u2:
            pass
        else:
            for key in nlist[ngraph]:
                nlist[ngraph][key] = [energy*Conv_mat[u1][u2] for energy in nlist[ngraph][key]]
        return nlist

    def createdatatable(self, energylist, RC):
        #print("ENERGY LIST:", energylist)
        rows = list(energylist[0].keys())
        #print('OOOOOOOOOO', self.cols)
        #print('RCRCRC', RC)
        TSs = {}

        h_scrollbar = ttk.Scrollbar(self.TableFrame, orient='horizontal')
        h_scrollbar.grid(row=1, column=0, sticky='ew')

        for i in rows:
            TSs[i], index = (self.calcacten(energylist[0][i]))
        #print('TSSSSS', TSs)
        
        cols = RC[0]

        endata = energylist[0]

        if TSs:
            # Get the length of the largest list
            largest_list_length = max(len(value) for value in TSs.values())
            for i in range(int(largest_list_length/2)):
                head1 = 'Ea dir ' + str(i+1)
                head2 = 'Ea inv ' + str(i+1)
                cols.append(head1)
                cols.append(head2)
            for i in rows:
                [endata[i].append(x) for x in TSs[i]]
        # Round decimal places
        for ref, values in endata.items():
            endata[ref] = [round(value, 6) for value in values]


        #############################
        #print('SELF_ENDATA', endata)
        #print('ENERGYLIST', energylist)
        #print('COLS', cols)
        treeview = ttk.Treeview(self.TableFrame, columns=(cols), xscrollcommand=h_scrollbar.set)
        h_scrollbar.config(command=treeview.xview)
        treeview.heading("#0", text="Mechanism")
        for i in cols:
            treeview.heading(i, text=i)
            treeview.column(i, width=150, anchor='center')
        for i in rows:
            treeview.insert(
                "", 
                tk.END,
                text=i, 
                values=endata[i]
            )

        treeview.grid(row=0, column=0, sticky="nsew")

    def calcacten(self, energies):
        TSlist = []
        TSindex = []
        loop = len(energies)-2
        for i in range(1, loop):
            a = energies[i] - energies[i-1]
            b = energies[i] - energies[i+1]
            if (a > 0 ) and (b > 0 ):
                TSlist.append(a)
                TSlist.append(b)
                TSindex.append(i+1)


        #print('TSs', TSlist)
        return TSlist, TSindex


    def createenergysets(self, ngraph, title, refs, mechs, energylist, units, RCoord):
        self.originalRC = copy.deepcopy(RCoord)
        self.originalEL = copy.deepcopy(energylist)
        print('hoa')
        self.energysets = EnergySettings(self.LeftFrame, ngraph, title, refs, mechs, self.originalEL, units, self.originalRC)
        print('hola')
        self.energysets.grid(row=1, column=0, padx=10, pady=(10, 0), sticky="nsw", columnspan=2)
        
        #self.originalRCA = []
        #self.originalRCA.extend(RCoord)

        #print('ORIGINALELantes', self.originalEL)
        print('ORIGINALRCantes', self.originalRC)
        #print('ELantes', energylist)
        print('RCantes', RCoord)
        self.createdatatable(energylist, RCoord)
        #print('ELdespues', energylist)
        print('RCdespues', RCoord)
        #print('ORIGINALELdespues', self.originalEL)
        print('ORIGINALELdespues', self.originalRC)

    def save_PES_png(self, fig):
        """Function to save the current figure as a PNG file."""
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
        if file_path:
            fig.savefig(file_path)
            print(f"Graph saved as {file_path}")

    def pltPES(self, frame, ngraph, RCoord, Titles, energylist, mechs, Us):
        self.normenlist = []
        self.convlist = []
        #energylist2 = copy.deepcopy(energylist)
        RC2 = copy.deepcopy(RCoord)
        self.mechs = mechs
        reaction_coordinates = [(i+1)*2 for i in range(len(RC2[0]))]
        for key in mechs[ngraph]:
            print(key, 'line =', self.mechs[ngraph][key].linestyle)
            print(key, 'color =', self.mechs[ngraph][key].color_name)
            print(key, 'SP =', self.mechs[ngraph][key].SPdict)

        for key in energylist[ngraph]:
            print('bartype', self.mechs[ngraph][key].barstyle)
            print('line', self.mechs[ngraph][key].linestyle)
            print('color', self.mechs[ngraph][key].color_name)
        # Normalization:

        if hasattr(self.energysets, 'sel_norm') and self.energysets.sel_norm.get() == 'on':
            if hasattr(self.energysets, 'refsubbox') and hasattr(self.energysets.refsubbox, 'sel_ref'):
                print('NORM?', self.energysets.sel_norm.get())
                print('SEL REF', self.energysets.refsubbox.sel_ref.get())
                self.normenlist = self.normalize(ngraph, energylist, self.energysets.refsubbox.sel_ref) # Hay que averiguar bien como normalizar esto.
                

        # Unit Conversion:
        if hasattr(self.energysets, 'sel_conv'):
            print('Actual Units: ', Us[ngraph])
            print('Obj Units:', self.energysets.sel_conv.get())
            if self.normenlist == []:
                convlist = self.conversion(ngraph, energylist, Us, self.energysets.sel_conv)
                self.convlist = convlist
            else:
                convlist = self.conversion(ngraph, self.normenlist, Us, self.energysets.sel_conv)
                self.convlist = convlist

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_title(f'{self.energysets.sel_title.get()}')
        ax.set_xlabel('Reaction Coordinate')
        ax.set_ylabel(f'{self.energysets.sel_Etype.get()} ({self.energysets.sel_conv.get()})')

        # Plot the energetic structures:
        for key in self.convlist[ngraph]:
            i = 0
            for i, value in enumerate(self.convlist[ngraph][key]):
                label = f"{key} ({self.mechs[ngraph][key].color_name}, {self.mechs[ngraph][key].barstyle})"
                if self.mechs[ngraph][key].SPdict == None:
                    ax.hlines(y=value, xmin=reaction_coordinates[i]-0.5, xmax=reaction_coordinates[i]+0.5, color=str(self.mechs[ngraph][key].color_code), linewidth=2, linestyles=self.mechs[ngraph][key].barstyle, label=label)
                elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'False':
                    ax.hlines(y=value, xmin=reaction_coordinates[i]-0.5, xmax=reaction_coordinates[i]+0.5, color=str(self.mechs[ngraph][key].color_code), linewidth=2, linestyles=self.mechs[ngraph][key].barstyle, label=label)
                elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'True':
                    ax.plot(reaction_coordinates[i], value, color=str(self.mechs[ngraph][key].color_code), marker='o', label=label)

                if i < len(self.convlist[ngraph][key]) - 1:
                    #plot lines connecting horizontal lines
                    if self.mechs[ngraph][key].SPdict == None and self.mechs[ngraph][key].SPdict == None:
                        x_values = [reaction_coordinates[i] + 0.5, reaction_coordinates[i+1] - 0.5]
                        y_values = [self.convlist[ngraph][key][i], self.convlist[ngraph][key][i+1]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle=self.mechs[ngraph][key].linestyle)
                    elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'False' and self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i+1]] == 'False':
                        x_values = [reaction_coordinates[i] + 0.5, reaction_coordinates[i+1] - 0.5]
                        y_values = [self.convlist[ngraph][key][i], self.convlist[ngraph][key][i+1]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle=self.mechs[ngraph][key].linestyle)
                    elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'True' and self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i+1]] == 'False':
                        x_values = [reaction_coordinates[i], reaction_coordinates[i+1] - 0.5]
                        y_values = [self.convlist[ngraph][key][i], self.convlist[ngraph][key][i+1]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle=self.mechs[ngraph][key].linestyle)
                    elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'False' and self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i+1]] == 'True':
                        x_values = [reaction_coordinates[i] + 0.5, reaction_coordinates[i+1]]
                        y_values = [self.convlist[ngraph][key][i], self.convlist[ngraph][key][i+1]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle=self.mechs[ngraph][key].linestyle)
                    elif self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i]] == 'True' and self.mechs[ngraph][key].SPdict[RC2[self.ngraph][i+1]] == 'True':
                        x_values = [reaction_coordinates[i], reaction_coordinates[i+1]]
                        y_values = [self.convlist[ngraph][key][i], self.convlist[ngraph][key][i+1]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle=self.mechs[ngraph][key].linestyle)
                i += 1
            if hasattr(self.energysets, 'sel_AE') and self.energysets.sel_AE.get() == 'on':
                #pass
                # Graph the Act En value:
                TSen, TSin = self.calcacten(self.convlist[ngraph][key])
                print('TSEN', TSen)
                print('TSINDEX', TSin)
                txtcoords = []
                for l, ind in enumerate(TSin):
                    if self.mechs[ngraph][key].SPdict[RC2[self.ngraph][ind-1]] == 'False':
                        # Vertical Line RIGHT:
                        x_values = [2*ind + 0.5, 2*ind + 0.5]
                        y_values = [self.convlist[ngraph][key][ind-1], self.convlist[ngraph][key][ind]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle='dotted')
                        # Vertical Line LEFT:
                        x_values = [2*ind - 0.5, 2*ind - 0.5]
                        y_values = [self.convlist[ngraph][key][ind-1], self.convlist[ngraph][key][ind-2]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle='dotted')
                        # Horizontal Line RIGHT:
                        x_values = [2*ind + 0.5, 2*(ind+1) - 0.5]
                        y_values = [self.convlist[ngraph][key][ind], self.convlist[ngraph][key][ind]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle='dotted')
                        # Horizonatal Line LEFT:
                        x_values = [2*ind - 0.5, 2*(ind-1) + 0.5]
                        y_values = [self.convlist[ngraph][key][ind-2], self.convlist[ngraph][key][ind-2]]
                        ax.plot(x_values, y_values, color=str(self.mechs[ngraph][key].color_code), linewidth=0.5, linestyle='dotted')
                        #print(self.convlist[ngraph][key][ind-1])
                        # Print the Activation Energy value in the graph
                        txt1 = plt.text(float(ind*2-1.5), np.mean([self.convlist[ngraph][key][ind-1], self.convlist[ngraph][key][ind-2]]), str(round(TSen[2*l],2)), color = str(self.mechs[ngraph][key].color_code), weight='bold')
                        txt2 = plt.text(float(ind*2+0.5), np.mean([self.convlist[ngraph][key][ind-1], self.convlist[ngraph][key][ind]]), str(round(TSen[2*l+1],2)), color = str(self.mechs[ngraph][key].color_code), weight='bold')
                        txt1.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
                        txt2.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
                        txt1_pos = txt1.get_position()
                        txt2_pos = txt2.get_position()
                        # If Ea text is overlapping, dispace them in y axis
                        evdiff = 1.0
                        for c in txtcoords:
                            if abs(txt1.get_position()[1] - c[1]) < evdiff and txt1.get_position()[1] - c[1] > 0:
                                txt1_pos = (txt1_pos[0], txt1_pos[1] + (evdiff-(txt1.get_position()[1] - c[1])))
                                txt1.set_position(txt1_pos)
                            elif abs(txt1.get_position()[1] - c[1]) < evdiff and txt1.get_position()[1] - c[1] < 0:
                                txt1_pos = (txt1_pos[0], txt1_pos[1] - abs(txt1.get_position()[1] - c[1]))
                                txt1.set_position(txt1_pos)

                            if abs(txt2.get_position()[1] - c[1]) < evdiff and txt2.get_position()[1] - c[1] > 0:
                                txt2_pos = (txt2_pos[0], txt2_pos[1] + (evdiff-(txt2.get_position()[1] - c[1])))
                                txt2.set_position(txt2_pos)
                            elif abs(txt2.get_position()[1] - c[1]) < evdiff and txt2.get_position()[1] - c[1] < 0:
                                txt2_pos = (txt2_pos[0], txt2_pos[1] - abs(txt2.get_position()[1] - c[1]))
                                txt2.set_position(txt2_pos)

                        txtcoords.append(txt1.get_position())
                        txtcoords.append(txt2.get_position())

                        print("DB1", txtcoords)

            #plt.text(5, 1, 'HOLA')
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

        # Embed Matplotlib figure into the tkinter frame
        for widget in frame.winfo_children():
            widget.destroy()  # Clear the frame

        # Create a canvas for the Matplotlib plot
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Add save button
        save_button = customtkinter.CTkButton(frame, text="Save as PNG", command=lambda: self.save_PES_png(fig))
        save_button.grid(row=1, column=0, pady=10)  # Place the button below the graph


        ## Set the frame to expand to fit the content
        #frame.grid_rowconfigure(0, weight=1)
        #frame.grid_columnconfigure(0, weight=1)

        plt.close(fig)
        self.createdatatable(convlist, RC2) # RC2 is needed so that the activation energies don't get inserted as new reaction coordinates.

class EnergySettings(customtkinter.CTkFrame):
    def __init__(self, master, ngraph, title, refs, mechs, energylist, Units, RC):
        super().__init__(master, width=W/3, height=3*H/8, border_color='green', border_width=4)
        self.grid_propagate(False)
        self.Title = title
        self.units = Units
        self.refs = refs
        self.ngraph = ngraph
        self.mechs = mechs
        self.RC = RC
        self.energylist = energylist
        self.refsubbox = None
        self.graphicalsets = None
        self.SPsettings = None
        print('HOLA')
        #DEFAULT SELECTIONS
        self.sel_ref = 'no ref'

        self.grid_columnconfigure(0, weight=1)
        self.title = customtkinter.CTkLabel(self, text=f'Energy settings for {self.Title}', font=('Helvetica', 18, 'bold'), fg_color="gray30", corner_radius=6)
        self.title.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="ew", columnspan=30)

        #DEFAULT SETTINGS
        self.values = []
        print(self.refs)
        for i in range(len(self.RC[self.ngraph])):
            print(i)
            print(self.values)
            self.values.append('False')
        print('TOT')
        if self.mechs[ngraph] == 0:
            m = {}
            j = 0
            for ref in self.refs:
                m[ref] = Mech(list(colors_html.keys())[j], 'dashed', 'solid', SPdict=dict(zip(self.RC[self.ngraph], self.values), S='no ref'))
                j += 1
            print('jelo', m[self.refs[0]].SPdict)
            self.mechs[self.ngraph] = m

# Conversion of units
        self.conv_label = customtkinter.CTkLabel(self, text='Unit Conversion:')
        self.conv_label.grid(row=1, column=1, padx=10, pady=(10, 0), sticky="w")
        self.sel_conv = tk.StringVar(value=Units[ngraph])
        self.conversion = customtkinter.CTkComboBox(self, values=['eV', 'kcal/mol', 'kJ/mol'], variable=self.sel_conv)
        self.conversion.grid(row=2, column=1, padx=10, pady=(10, 0), sticky="w")
# Change Title
        self.title_label = customtkinter.CTkLabel(self, text='Change Title:')
        self.title_label.grid(row=4, column=2, padx=10, pady=(10, 0), sticky="w")
        self.sel_title = tk.StringVar(value=title)
        self.chg_title = customtkinter.CTkEntry(self, placeholder_text=f'{title}', textvariable=self.sel_title)
        self.chg_title.grid(row=5, column=2, padx=10, pady=(10, 0), sticky="w")

# Energy Type
        self.Etype = customtkinter.CTkLabel(self, text='Energy Type:')
        self.Etype.grid(row=1, column=2, padx=10, pady=(10, 0), sticky="w")
        self.sel_Etype = tk.StringVar(value='Potential Energy')
        self.conversion = customtkinter.CTkComboBox(self, values=Etypes, variable=self.sel_Etype)
        self.conversion.grid(row=2, column=2, padx=10, pady=(10, 0), sticky="w", ipadx=5, ipady=5)

# Normalization
        self.sel_norm = tk.StringVar(value='off')
        self.normcheck = customtkinter.CTkCheckBox(self, text='Normalize', variable=self.sel_norm, onvalue='on', offvalue='off', command=lambda: self.createrefsubbox(self.sel_norm, self.refs)) #Comprobar si está pulsado cuando se genere el gráfico
        self.normcheck.grid(row=1, column=0, padx=10, pady=(10, 0), sticky="nsew")
# Create Graphical Settings Window
        self.graphsets = customtkinter.CTkButton(self, text="Graphical Settings", command=lambda:self.opengraphsets(self.Title, self.ngraph, self.refs, self.mechs, self.energylist, self.RC))
        self.graphsets.grid(row=3, column=1, padx=10, pady=(10, 0), sticky="w")
# Save as mpes button
        self.mpesbut = customtkinter.CTkButton(self, text="Save as .mpes", command=lambda:self.save_as_mpes(self.ngraph ,self.Title, self.units, self.refs, self.RC, self.energylist, self.mechs))
        self.mpesbut.grid(row=4, column=1, padx=10, pady=(10, 0), sticky="w")
# Show Activation Energy button
        self.sel_AE = tk.StringVar(value='off')
        self.AEbox = customtkinter.CTkCheckBox(self, text="Show Activation Energy", variable=self.sel_AE, onvalue='on', offvalue='off')
        self.AEbox.grid(row=3, column=0, padx=10, pady=(10, 0), sticky="nsew")

    def opengraphsets(self, title, ngraph, refs, mechs, energylist, RC):
        if not self.graphicalsets:
            print('Creating Graphical Settings')
            # Create Canvas
            self.newWin = customtkinter.CTkToplevel(self, fg_color="black")
            self.newWin.title('Graphical Settings')
            self.graphicalsets = GraphicalSettings(self.newWin, title, ngraph, refs, mechs, energylist, RC)
            self.graphicalsets.grid( padx=10, pady=(10, 0))
        else:
            if self.graphicalsets:
                self.graphicalsets.destroy()
                self.graphicalsets = None
            if self.newWin:
                self.newWin.destroy()
                self.newWin = None

    def createrefsubbox(self, sel_norm, refs):
        if sel_norm.get() == 'on':
            if not self.refsubbox:
                print('Creating RefsSubBox')
                self.refsubbox = RefSubBox(self, refs)
                self.refsubbox.grid(row=2, column=0, padx=10, pady=(10, 0), sticky="nsew")
        else:
            if self.refsubbox:
                self.refsubbox.destroy()
                self.refsubbox = None

    def save_as_mpes(self, ngraph, title, units, refs, RC, en, mechs):
        print('TITLE', title)
        print('UNTIS', units[ngraph])
        print('REFS', refs)
        print('RC', RC[ngraph])
        print('EN', en[ngraph])
        print('MECHS', mechs)
    
        file_path = filedialog.asksaveasfilename(defaultextension=".mpes", filetypes=[("MPES files", "*.mpes")])
        if file_path:
            with open(file_path, 'w') as f:
                f.write(title + '\n')
                f.write("Reaction Coordinates:"+ '\n')
                f.write(units[ngraph])
                for coor in RC[ngraph]:
                    f.write(', ' + (coor))
                f.write('\n')
                f.write("Energy:" + '\n')
                f.write(str(len(refs)) + '\n')
                for key, values in en[ngraph].items():
                    f.write(key)
                    for value in values:
                        f.write(', ' + str(value))
                    f.write('\n')
                f.write('Graph info' + '\n')
                for key, value in mechs[ngraph].items():
                    f.write(key)
                    f.write(', ' + value.color_name)
                    f.write(', ' + value.linestyle)
                    f.write(', ' + value.barstyle)
                    f.write('\n')

    
class RefSubBox(customtkinter.CTkFrame):
    def __init__(self, master, refs):
        super().__init__(master)
        # Setting Refference
        self.refbox = ['no ref', 'all zero'] + refs
        self.reflabel = customtkinter.CTkLabel(self, text='PES reference:')
        self.reflabel.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="nsew")
        self.sel_ref = tk.StringVar(self)
        self.sel_ref.set(self.refbox[0])
        self.refcombox = customtkinter.CTkOptionMenu(self, variable=self.sel_ref, values=self.refbox)
        self.refcombox.grid(row=0, column=1, padx=10, pady=(10, 0), sticky="nsew")

class GraphicalSettings(customtkinter.CTkFrame):
    def __init__(self, master, title, ngraph, refs, mechs, energylist, RC):
        super().__init__(master)
        self.ngraph = ngraph
        self.RC = RC
        self.title = title
        self.refs = refs # refs in a list
        self.mechs = mechs # list containing a dictionary with key=refs and values = Mech types containing graphical info. 
        self.m = mechs[ngraph]
        self.colordict = {}
        self.linedict = {}
        self.bardict = {}
        self.SPsubboxes = {}  # Dictionary to store SPsubboxes for each key
        print('MECHS', self.m)
        print('REFS', refs)
        # DEFAULT SETTINGS
        print('DEFAULTSETTING:', self.mechs[ngraph])
        print('RCin GRAPHSETTS:', RC)
        #CREATE SP grid:
        self.SPframe = customtkinter.CTkFrame(self)
        self.SPframe.grid(row=4+len(energylist[ngraph]) ,column=0, padx=10, pady=(10, 0), sticky='w', columnspan=len(self.refs))

        self.SP_vars = {}
        self.title = customtkinter.CTkLabel(self, text=f'Graphical Settings for {self.title}')
        self.title.grid(row=0 ,column=0, padx=10, pady=(10, 0), sticky='we', columnspan=3)
        self.subtitles = [f'Mechanism {title}', 'Color', 'Linestyle', 'Barstyle', 'Single Points']
        for i in range(len(self.subtitles)):
            stit = customtkinter.CTkLabel(self, text=self.subtitles[i])
            stit.grid(row=1 ,column=i, padx=10, pady=(10, 0), sticky='w')
        
        #print('len(refs[ngraph])', len(refs[ngraph]))

        for i in range(len(self.refs)):
            print('refs', self.refs)
            print('i', i)
            self.colordict[self.refs[i]] = 0
            print(self.refs[i])

        for index, (key, en) in enumerate(energylist[ngraph].items()):
            label = customtkinter.CTkLabel(self, text=key)
            label.grid(row = index+2, column = 0, padx=10, pady=5)
            # Colors
            color_var = tk.StringVar(self)
            color_var.set(self.m[key].color_name)
            self.colordict[key] = color_var
            color_menu = customtkinter.CTkOptionMenu(self,variable=color_var, values=list(colors_html.keys()), command=combobox_callback(color_var))
            color_menu.grid(row = index+2, column = 1, padx=10, pady=5)
            # Linestyle
            line_var = tk.StringVar(self)
            line_var.set(self.m[key].linestyle)
            self.linedict[key] = line_var
            line_menu =  customtkinter.CTkOptionMenu(self,variable=line_var, values=linestyles)
            line_menu.grid(row = index+2, column = 2, padx=10, pady=5)
            # Barstyle
            bar_var = tk.StringVar(self)
            bar_var.set(self.m[key].barstyle)
            self.bardict[key] = bar_var
            bar_menu =  customtkinter.CTkOptionMenu(self,variable=bar_var, values=linestyles)
            bar_menu.grid(row = index+2, column = 3, padx=10, pady=5)
            # SP option
            self.SP_vars[key] = tk.StringVar(value='off')
            SP_check = customtkinter.CTkCheckBox(self, text='',variable=self.SP_vars[key],  onvalue='on', offvalue='off', command=lambda k=key, i=index: self.createSPboxes(self.ngraph, self.SP_vars[k], energylist, k, i))
            SP_check.grid(row = index+2, column = 4, padx=10, pady=5, sticky='nsew')

        self.save_button = customtkinter.CTkButton(self, text='Save Settings', command=lambda:self.savesetts(ngraph, self.colordict, self.linedict, self.bardict, self.mechs, self.SPsubboxes))
        self.save_button.grid(row = len(refs[ngraph])+3, column = 0, padx=10, pady=5)
        
    def createSPboxes(self, ngraph, sel_SP, enlist, key, index):
            print('enteredSPboxes')
            print(sel_SP.get())
            if sel_SP.get() == 'on':
                print(f'Creating SPsubbox for {key}')
                if key not in self.SPsubboxes:
                    spbox = SPsubbox(ngraph, self.SPframe, enlist, key, index, self.RC)
                    spbox.grid(row=len(enlist[ngraph]), column=index, padx=10, pady=5)  
                    print('Creating SPsubbox')
                    self.SPsubboxes[key] = spbox
                else:
                    print('SPsubbox already exists')
            else:
                print(f'Checkbox for {key} is OFF')
                if key in self.SPsubboxes:
                    print(f'Destroying Subbox {key}')
                    self.SPsubboxes[key].destroy()
                    self.SPsubboxes[key].title.destroy()
                    del self.SPsubboxes[key]

    def savesetts(self, ngraph, sel_color, sel_line, sel_bar, mechs, SPdict):
        print(f"Saving settings for ngraph: {ngraph}")
        print('SPsubbox', SPdict )
        for key in SPdict:
            newdic = {}
            #for c in SPdict[key]:
            if len(self.SPsubboxes[key].RCspdict.keys()) != len(self.RC[self.ngraph]):
                print(f'Length of SPdict ({len(self.SPsubboxes[key].RCspdict.keys())}) doesnt match length of RCoord ({len(self.RC[self.ngraph])})')
            for c in self.RC[self.ngraph]:
                newdic[c] = SPdict[key].RCspdict[c].get()
                print('newdic', newdic)
            mechs[ngraph][key].SPdict = newdic
            print(f'Updated {key}: ref = {newdic}')

        for key in sel_color.keys():
            color_name = sel_color[key].get()
            if key in mechs[ngraph]:
                mechs[ngraph][key].color_name = color_name
                mechs[ngraph][key].color_code = colors_html.get(color_name, None)
                print(f'Updated {key}: color_name = {color_name}, color_code = {colors_html.get(color_name, None)}')
            else:
                print(f'Key {key} not found in mechs[{ngraph}]')

        for key in sel_line.keys():
            if key in mechs[ngraph]:
                mechs[ngraph][key].linestyle = sel_line[key].get()
                print(f'Updated {key}: linestyle = {sel_line[key].get()}')
            else:
                print(f'Key {key} not found in mechs[{ngraph}]')

        for key in sel_bar.keys():
            if key in mechs[ngraph]:
                mechs[ngraph][key].barstyle = sel_bar[key].get()
                print(f'Updated {key}: barstyle = {sel_bar[key].get()}')
            else:
                print(f'Key {key} not found in mechs[{ngraph}]')

class SPsubbox(customtkinter.CTkFrame):
    def __init__(self, ngraph, master, enlist, ref, index, RC):
        super().__init__(master)
        self.RC = RC
        self.ref = ref
        self.RCspdict = {}
        self.title = customtkinter.CTkLabel(master, text=f'SP Settings for {self.ref}')
        self.title.grid(row=0, column=index, padx=10, pady=(10, 0), sticky='we')
        # Create checkboxes for SP
        i = 0
        for c in RC[ngraph]:
            self.RCspdict[c] = tk.StringVar(value='False')
            ckSP = customtkinter.CTkCheckBox(self, text=c,variable=self.RCspdict[c], onvalue='True', offvalue='False')
            ckSP.grid(row = i, column=index, padx=10, pady=(10, 0), sticky='we')
            i += 1

class Mech:
    def __init__(self, color, linestyle, barstyle, SPdict=None, S=None):
        self.SPdict = SPdict
        #if SPdict == None:
        #    self.SPdict = 

        if color in colors_html:
            self.color_name = color
            self.color_code = colors_html[color]
            self.linestyle = linestyle
            self.barstyle = barstyle
            self.ref = S
        else:
            raise ValueError(f"Color '{color}' is not in the list of available colors.")
# Global parameters:

eVtokcal = 23.060541945329334
eVtokJ = 96.48530749925794
kcaltokJ = 4.184

# UNITS
Units_list = ['eV', 'kcal/mol', 'kJ/mol']

Conv_mat = [
    [1, eVtokcal, eVtokJ],
    [1/eVtokcal, 1, kcaltokJ],
    [1/eVtokJ, 1/kcaltokJ, 1]
]

# COLORS
colors_html = {
    'Red': '#FF0000',
    'Green': '#008000',
    'Blue': '#0000FF',
    'Yellow': '#FFFF00',
    'Orange': '#FFA500',
    'Purple': '#800080',
    'Pink': '#FFC0CB',
    'Brown': '#A52A2A',
    'Black': '#000000',
    'White': '#FFFFFF',
    'Gray': '#808080',
    'Cyan': '#00FFFF',
    'Magenta': '#FF00FF',
    'Lime': '#00FF00',
    'Olive': '#808000',
    'Navy': '#000080',
    'Teal': '#008080',
    'Maroon': '#800000',
    'Silver': '#C0C0C0',
    'Gold': '#FFD700'
}

# LINEstyleS
linestyles = ['solid', 'dotted', 'dashed', 'dashdot']

root = App()
customtkinter.set_appearance_mode("dark")

# Run the main event loop
root.mainloop()