#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
Created on 8th June 2017

Extracts some key info from a Conquest_out file
and prints it in a neat summary.

If no command line arguments are given, it looks
for a file called Conquest_out in the current
directory.


@author: kshenton
"""

# ## Summary of Conquest output file

# Class to parse conquest_out files:
class Conquest_out:
    """This is a class for extracting useful information from the Conquest_out file.
        In particular we extract the lattice vectors and a bunch of other things.
        TODO add more things... """
    def __init__(self, filename='Conquest_out'):

        re_maxforce = re.compile("Maximum force \:")
        re_stress = re.compile("Total stress\:")
        re_volume = re.compile("volume of cell")
        re_lattice = re.compile("simulation box", re.M)
        re_positions = re.compile("pseudopotential will be used", re.M)
        re_species = re.compile("Label     Mass \(a.u.\)   Charge \(e\)  NLPF")
        re_toten  = re.compile("\|\* DFT total energy")
        re_gridspacing = re.compile("integration grid spacing along x")
        re_mp_mesh = re.compile("Monkhorst-Pack mesh\:")
        re_xc = re.compile("The functional used will be")

        species = []
        positions = []
        with open(filename) as f:
            for line in f:
                if re_positions.search(line):
                    while True:
                        temp = f.next()
                        # starts with a digit:
                        startswithdigit = re.search(r"^\s*\d", temp)
                        if startswithdigit:

                            positions.append([float(temp.split()[1]),
                                             float(temp.split()[2]),
                                             float(temp.split()[3]),
                                             int(temp.split()[4]) # species
                                            ])
                        else:
                            break
                if re_species.search(line):
                    f.next() # underline -----
                    while True:
                        # I'm sure there's better regex for this...
                        temp = f.next()
                        tempinfo = re.search(r"^\s*\d*(\w)\s*(\w*)\s*(\d*[.]\d*)\s*(\d*[.]\d*)\s*(\d*[.]\d*)\s*(\d*)", temp)
                        endline = re.search(r"^\s*---", temp)
                        if tempinfo:
                            spec =  tempinfo.groups()
                            specdic = {"label": spec[1],
                                       "mass": float(spec[2]),
                                       "charge": float(spec[3]),
                                       "radius": float(spec[4]),
                                       "nsf": int(spec[5])}
                            species.append(specdic)
                        elif not endline:
                            # Something spilled onto the next line
                            continue
                        else:
                            break # out of while loop
                if re_toten.search(line):
                    toten = float(line.split()[-2])
                if re_mp_mesh.search(line):
                    mp_mesh = np.array([int(line.split()[2]),
                                        int(line.split()[3]),
                                        int(line.split()[4])])
                if re_lattice.search(line):
                    lat1=f.next()
                if re_xc.search(line):
                    xc = line.split()[-2:]
                if re_gridspacing.search(line):
                    gridspacingx = float(line.split()[-2])
                    gridspacingy = float(f.next().split()[-2])
                    gridspacingz = float(f.next().split()[-2])
                if re_maxforce.search(line):
                    maxforce = float(line.split()[3].split('(')[0])
                if re_stress.search(line):
                    stress = np.array([float(line.split()[2]) ,float(line.split()[3]), float(line.split()[4])])
        lattice = np.array([float(lat1.split()[2]), float(lat1.split()[5]), float(lat1.split()[8])])
        gridspacing = np.array([gridspacingx, gridspacingy, gridspacingz])

        self.lattice = lattice # box lengths (x, y, z) in a0
        self.volume = lattice.prod() # product of box lengths (because cell is orthorhombic) in a0^3
        self.gridspacing = gridspacing # Grid spacing in [x, y, z] in a0
        self.species = species # Dictionary of species
        self.mp_mesh = mp_mesh # Monkhorst pack mesh size
        self.positions = np.array(positions) # positions of atoms (cartesian, a0)
        self.num_atoms = len(positions) # Total number of atoms
        self.xc = xc # [xc type, authors]
        self.maxforce = maxforce # Max force in Ha/a0
        self.stress = stress # Total stress in Ha
        self.toten = toten # Total DFT energy in eV

# Pretty print out
def prettyprint(out):
    print 80*'='
    print 31*'-',"CONQUEST SUMMARY", 31*'-'
    print 80*'='


    print "\n"+21*"-"+" Atomic species "+21*"-"
    print "+"+56*"-"+"+"
    print "|Label     Mass (a.u.)   Charge (e)   NLPF Rad (a0)   NSF|"
    print "+"+56*"-"+"+"
    for species in out.species:
        print "| {0:2}  {1:14.5f}  {2:12.5f}  {3:12.5f}  {4:6d} |".format(species['label'],
                                                                          species['mass'],
                                                                          species['charge'],
                                                                          species['radius'],
                                                                          species['nsf'])
    print "+"+56*"-"+"+"







    print "\n --- Lattice vectors (bohr): ---"
    print out.lattice*np.eye(3)





    print "\n --- Unit cell volume (bohr^3) ---\n {0:3.3f}".format(out.volume)


    print "\n --- Number of atoms ---\n {0}".format(out.num_atoms)



    print "\n"+20*"-"+" Atomic positions "+20*"-"
    print "+"+56*"-"+"+"
    print "|Species         X              Y             Z          |"
    print "+"+56*"-"+"+"
    for pos in out.positions:
        print "| {0:2}  {1:14.5f}  {2:12.5f}  {3:12.5f}         |".format(int(pos[3]),
                                                                              pos[0],
                                                                              pos[1],
                                                                              pos[2])
    print "+"+56*"-"+"+"





    print "\n --- Grid spacing in x, y, z (bohr) ---\n{0:2.4f}, {1:2.4f}, {2:2.4f}".format(out.gridspacing[0],
                                                                                 out.gridspacing[1],
                                                                                 out.gridspacing[2])





    print "\n --- Monkohorst Pack mesh ---\n {0} x {1} x {2}".format(out.mp_mesh[0],
                                                                     out.mp_mesh[1],
                                                                     out.mp_mesh[2])


    print "\n --- Maximum Force (Ha /bohr) ---\n {0:3.6f}".format(out.maxforce)

    print "\n --- Total Stress (Ha) ---\n {0:3.6f} {1:3.6f} {2:3.6f}".format(out.stress[0],
                                                                             out.stress[1],
                                                                             out.stress[2])

    print "\n --- XC functional & authors ---\n {0}, {1}".format(out.xc[0], out.xc[1])

    print "\n --- DFT Total energy (Ha) ---\n {0:3.6f}".format(out.toten)

if __name__ == "__main__":

    import numpy as np
    import sys
    import re

    # Read the path to the Conquest_out file from the command line
    # If nothing is given, it assumes the file is called
    # Conquest_out

    if len(sys.argv) < 2:
        outfile = "Conquest_out"
    elif len(sys.argv) == 2 :
        outfile = str(sys.argv[1])
    else:
        print "Que?!?!"
    print "\n\nOpening {} as conquest_out file.\n\n".format(outfile)


    out = Conquest_out(outfile)
    prettyprint(out)
