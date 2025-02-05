M-PES:
GUI to genereate potential energy diagrams. The input files should be .csv. An example is in the directory.
Code is run by complete.py
2.1 Loading a File
Press the ’Load File’ button and load a file with a suitable format (.csv or .mpes). The
information will be loaded into a table at the bottom of the screen. Activation energies
will also be calculated in situ and displayed inside the table.
2.2 Plot the PED
Pressing the ’Plot PES’ button will generate a PED with the imputed energies. This graph
can then be saved as a png in the desired directory.
2.3 Energy Settings
Energy settings affect the PED data and visualization. The settings are listed as follows:
Normalize Use a reference for the energy of the PED. no ref uses the raw energy values
of the input file. all zero sets the first energy value of every mechanism to zero.
Select the name of a mechanism to make its first value the reference.
Unit Conversion Converts the units. Available units: kcal/mol, eV and kJ/mol.
Energy Type Changes the text on the y axis.
1
Change Title Changes the title name.
Show Activation Energy Shows the activation energy in the PED.
Graphical Settings Open the Graphical Settings menu.
Save as mpes Saves the graph as a mpes file. The mpes file saves the graphical settings
of the graph so that it can be loaded with all the graphical modifications previously
done.
2.4 Graphical Settings
Graphical Settings modify the style of the graph. The settings are listed below:
Color Changes the color of the mechanism lines and points.
Linestyle Changes the style of the lines that connect the points in the mechanism.
Barstyle Changes the style of the bars that represent the energy values of the mechanism.
Single Points Turns on the option of selecting some energy points in the mechanism
and making them single points. These points are not used for the calculation of the
activation energies and are represented as a point in the diagram.
Always click save settings to apply the graphical settings.
2.5 Input Files
Recognized input files are .csv and .mpes. The .csv file must follow the following format:
Title of the PED
Reaction Coordinate: #must be
{eV, kcal/mol, kj/mol} #select unit, A, B, TS1, C #points inside mechanism
Energy: #must be
2 #number of mechanisms
mech1, 3.2914905999999, 3.5742076000000, 1.5069175999999, 1.1574676000000
mech2, 0.16611760000, -0.09370239999998, 1.1744876000000, 1.1574676000000
#name of mechanism and energy values
