#!/usr/bin/env python
# coding: utf-8

# In[1]:





from __future__ import print_function
import sys
import os
import numpy as np
import Bio.PDB

#from sys import argv
#script, pdbfile, chain1, chain2, chain3 = argv

import argparse
parser = argparse.ArgumentParser(description='Process input file.')
parser.add_argument('--pdbfile', type=str)
parser.add_argument('--mhc_1_chain', type=str)
parser.add_argument('--mhc_2_chain', type=str, default=None)
parser.add_argument('--tcr_a_chain', type=str)
parser.add_argument('--tcr_b_chain', type=str)

reffile1="11111.pdb"
reffile2="22222.pdb"




#import mdtraj as md


# In[2]:


# Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Module with assorted geometrical functions on
macromolecules.
"""

from Bio.PDB import Entity

def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """
    
    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else: # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                            "Structure, Model, Chain, Residue, list of Atoms.")
    
    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    
    for atom in atom_list:
        masses.append(atom.mass)
        
        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead.")
    
    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:       
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]


# In[223]:


# Geometrical parameters of TCR-MHC class I
def tcr_mhci_pos(pdfile, mhc_a='A', tcr_a='D', tcr_b='E',
                 mhc_a_init=1, mhc_a_final=180,
                 tcr_a_init=1, tcr_a_final=109, 
                 tcr_b_init=1, tcr_b_final=116):
# Define residues range to align and center of mass calculations
    mhci_atoms_to_align_com = range(mhc_a_init, mhc_a_final + 1)
# Get the structures
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure("reference", reffile1)
    sample_structure = pdb_parser.get_structure("sample", "./%s"%pdbfile)
# Use the first model in the pdb-files for alignment
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]
    ref_atoms    = []
    sample_atoms = []
# Iterate of all residues in each model in order to define proper atoms
    #Reference structure
    for ref_res in ref_model['%s'%mhc_a]:
        if ref_res.get_id()[1] in mhci_atoms_to_align_com:
            ref_atoms.append(ref_res['CA'])
    # Sample structure
    for sample_res in sample_model['%s'%mhc_a]:
        if sample_res.get_id()[1] in mhci_atoms_to_align_com:
            sample_atoms.append(sample_res['CA'])
# Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
# Save the aligned version of pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure) 
    #io.save("%s_aligned.pdb"%pdbid)
# Calculate CoM of MHCI binding groove    
    mhci_com = center_of_mass(sample_atoms, geometric=True)

# Calculate CoM of vTCR 
    tcr_a_atoms_for_com = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_atoms_for_com = range(tcr_b_init, tcr_b_final + 1)
    tcr_atoms_for_com = []
    for tcr_res in sample_model['%s'%tcr_a]:
        if tcr_res.get_id()[1] in tcr_a_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    for tcr_res in sample_model['%s'%tcr_b]:
        if tcr_res.get_id()[1] in tcr_b_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    #vtcr_com = (center_of_mass(tcr_atoms_for_com))
    vtcr_com = (center_of_mass(tcr_atoms_for_com, geometric=True))

# Geomitrical parameters
    dx, dy, dz = (np.subtract(vtcr_com, mhci_com))
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhci_com))))
    r = round(r, 2)
    theta = np.degrees(np.arctan2(dy, dx))
    theta = round(theta, 2)
    phi = np.degrees(np.arccos(dz / r))
    phi = round(phi, 2)
    print('The Geomitrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
    return r, theta, phi


# In[224]:


#tcr_mhci_pos('1ao7')


# In[227]:


# Geometrical parameters of TCR-MHC class II
def tcr_mhcii_pos(pdbid, mhc_a='A', mhc_b='B', tcr_a='D', tcr_b='E',
                 mhc_a_init=1, mhc_a_final=80, mhc_b_init=1, mhc_b_final=90,
                 tcr_a_init=1, tcr_a_final=109, tcr_b_init=1, tcr_b_final=116):
# Define residues range to align and center of mass calculations
    mhc_a_to_be_aligned = range(mhc_a_init, mhc_a_final + 1)
    mhc_b_to_be_aligned = range(mhc_b_init, mhc_b_final + 1)
# Get the structures
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure("reference", reffile2)
    sample_structure = pdb_parser.get_structure("sample", "./%s"%pdbid)
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]
    ref_atoms = []
    sample_atoms = []
# Iterate of all residues in each model in order to define proper atoms
    #Reference structure
    for mhc_res in ref_model['%s'%mhc_a]:
        if mhc_res.get_id()[1] in mhc_a_to_be_aligned:
            ref_atoms.append(mhc_res['CA'])
    for mhc_res in ref_model['%s'%mhc_b]:
        if mhc_res.get_id()[1] in mhc_b_to_be_aligned:
            ref_atoms.append(mhc_res['CA'])
    # Sample structure
    for sample_res in sample_model['%s'%mhc_a]:
        if sample_res.get_id()[1] in mhc_a_to_be_aligned:
            sample_atoms.append(sample_res['CA'])
    for sample_res in sample_model['%s'%mhc_b]:
        if sample_res.get_id()[1] in mhc_b_to_be_aligned:
            sample_atoms.append(sample_res['CA'])
# Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
# Save the aligned version of pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_structure) 
    io.save("aligned_%s"%pdbfile)
# Calculate CoM of MHCII binding groove     
    mhcii_com = center_of_mass(sample_atoms, geometric=True)
# Calculate CoM of vTCR     
    tcr_a_atoms_for_com = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_atoms_for_com = range(tcr_b_init, tcr_b_final + 1)
    tcr_atoms_for_com = []
    for tcr_res in sample_model['%s'%tcr_a]:
        if tcr_res.get_id()[1] in tcr_a_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    for tcr_res in sample_model['%s'%tcr_b]:
        if tcr_res.get_id()[1] in tcr_b_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    #vtcr_com =center_of_mass(tcr_atoms_for_com)
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)
# Geomitrical parameters
    dx, dy, dz = (np.subtract(vtcr_com, mhcii_com))
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhcii_com))))
    r = round(r, 2)
    theta = np.degrees(np.arctan2(dy, dx))
    theta = round(theta, 2)
    phi = np.degrees(np.arccos(dz / r))
    phi = round(phi, 2)
    print('The Geomitrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
    return r, theta, phi    


# In[228]:


#tcr_mhcii_pos('1FYT')

args = parser.parse_args()
pdbfile = args.pdbfile
mhc_1_chain = args.mhc_1_chain
mhc_2_chain = args.mhc_2_chain
tcr_a_chain = args.tcr_a_chain
tcr_b_chain = args.tcr_b_chain

if mhc_2_chain is None:
    tcr_mhci_pos(pdbfile, mhc_1_chain, tcr_a_chain, tcr_b_chain)
else:
    tcr_mhcii_pos(pdbfile, mhc_1_chain, mhc_2_chain, tcr_a_chain, tcr_b_chain) 

# In[ ]:




