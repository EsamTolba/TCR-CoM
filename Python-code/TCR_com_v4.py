#!/usr/bin/env python
# coding: utf-8

"""Import modules"""
from __future__ import print_function
# from sys import argv
import sys
import os
import numpy as np
import Bio.PDB
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import argparse
# import mdtraj as md

""" Define arguments"""
parser = argparse.ArgumentParser(description = 'TCR-CoM calculates the geometrical parameters (r, theta, phi) of T cell receptors (TCR) on the top of MHC proteins')
parser.add_argument('-pdbid', type = str, help = 'PDB-ID', )
parser.add_argument('-mhc_a', type = str, help = 'ID of MHC alpha chain')
parser.add_argument('-mhc_b', type = str, default = None, help = 'ID of MHC beta chain')
parser.add_argument('-tcr_a', type = str, help = 'ID of TCR alpha chain')
parser.add_argument('-tcr_b', type = str, help = 'ID of TCR beta chain')

"""Define reference structures"""
reffile_path = os.path.dirname(os.path.abspath(__file__))
reffile1=os.path.join(reffile_path, "ref1.pdb")
reffile2=os.path.join(reffile_path, "ref2.pdb")


def add_com_to_pdb(mhci_com, vtcr_com, sample_structure):
    #mhc_com
    mhc_com_chain = 'X'
    sample_structure[0].add(Chain(mhc_com_chain))
    res_id = (' ', 1, ' ')
    new_residue = Residue(res_id, "MCM", ' ') 
    new_atom = Atom("C", mhci_com, 0, 0.0, ' ', "C", 1, "C")
    new_residue.add(new_atom)
    sample_structure[0].child_dict[mhc_com_chain].add(new_residue)
    #tcr com
    tcr_com_chain = 'Y'
    sample_structure[0].add(Chain(tcr_com_chain))
    res_id = (' ', 1, ' ')
    new_residue = Residue(res_id, "TCM", ' ') 
    new_atom = Atom("C", vtcr_com, 0, 0.0, ' ', "C", 1, "C")
    new_residue.add(new_atom)
    sample_structure[0].child_dict[tcr_com_chain].add(new_residue)
    return sample_structure

"""Center of mass function"""
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
    else:  # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                         "Structure, Model, Chain, Residue, list of Atoms.")

    masses = []
    positions = [[], [], []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead")

    if geometric:
        return [sum(coord_list) / len(masses) for coord_list in positions]
    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index] * atom_mass)
            w_pos[1].append(positions[1][atom_index] * atom_mass)
            w_pos[2].append(positions[2][atom_index] * atom_mass)

        return [sum(coord_list) / sum(masses) for coord_list in w_pos]


"""Geometrical parameters of TCR-MHC class I"""
def tcr_mhci_pos(pdbid, mhc_a='A', tcr_a='D', tcr_b='E',
                 mhc_a_init=1, mhc_a_final=179,
                 tcr_a_init=1, tcr_a_final=109,
                 tcr_b_init=1, tcr_b_final=116):
    # Define residues range to align and center of mass calculations
    mhci_atoms_to_align_com = range(mhc_a_init, mhc_a_final + 1)
    # Get the structures
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reffile1)
    sample_structure = pdb_parser.get_structure("sample", "./%s.pdb" % pdbid)
    # Use the first model in the pdb-files for alignment
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    ref_atoms = []
    sample_atoms = []
    # Iterate of all residues in each model in order to define proper atoms
    # Reference structure
    for ref_res in ref_model['%s' % mhc_a]:
        if ref_res.get_id()[1] in mhci_atoms_to_align_com:
            ref_atoms.append(ref_res['CA'])
    # Sample structure
    for sample_res in sample_model['%s' % mhc_a]:
        if sample_res.get_id()[1] in mhci_atoms_to_align_com:
            sample_atoms.append(sample_res['CA'])
        # Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCI binding groove
    mhci_com = center_of_mass(sample_atoms, geometric=True)

    # Calculate CoM of vTCR
    tcr_a_atoms_for_com = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_atoms_for_com = range(tcr_b_init, tcr_b_final + 1)
    tcr_atoms_for_com = []
    for tcr_res in sample_model['%s' % tcr_a]:
        if tcr_res.get_id()[1] in tcr_a_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    for tcr_res in sample_model['%s' % tcr_b]:
        if tcr_res.get_id()[1] in tcr_b_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    vtcr_com = (center_of_mass(tcr_atoms_for_com, geometric=True))

    # Save the aligned version of pdb
    ali_structure = add_com_to_pdb(mhci_com, vtcr_com, sample_structure)
    io = Bio.PDB.PDBIO()
    io.set_structure(ali_structure)
    io.save("%s_aligned.pdb"%pdbid)

    # Geometrical parameters
    dx, dy, dz = (np.subtract(vtcr_com, mhci_com))
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhci_com))))
    r = round(r, 2)
    theta = np.degrees(np.arctan2(dy, dx))
    theta = round(theta, 2)
    phi = np.degrees(np.arccos(dz / r))
    phi = round(phi, 2)
    print('The Geomitrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
    return r, theta, phi


"""Geometrical parameters of TCR-MHC class II"""
def tcr_mhcii_pos(pdbid, mhc_a='A', mhc_b='B', tcr_a='D', tcr_b='E',
                  mhc_a_init=1, mhc_a_final=80, mhc_b_init=1, mhc_b_final=90,
                  tcr_a_init=1, tcr_a_final=109, tcr_b_init=1, tcr_b_final=116):
    # Define residues range to align and center of mass calculations
    mhc_a_to_be_aligned = range(mhc_a_init, mhc_a_final + 1)
    mhc_b_to_be_aligned = range(mhc_b_init, mhc_b_final + 1)
    # Get the structures
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reffile2)
    sample_structure = pdb_parser.get_structure("sample", "./%s.pdb" % pdbid)
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    ref_atoms = []
    sample_atoms = []
    # Iterate of all residues in each model in order to define proper atoms
    # Reference structure
    for mhc_res in ref_model['%s' % mhc_a]:
        if mhc_res.get_id()[1] in mhc_a_to_be_aligned:
            ref_atoms.append(mhc_res['CA'])
    for mhc_res in ref_model['%s' % mhc_b]:
        if mhc_res.get_id()[1] in mhc_b_to_be_aligned:
            ref_atoms.append(mhc_res['CA'])
    # Sample structure
    for sample_res in sample_model['%s' % mhc_a]:
        if sample_res.get_id()[1] in mhc_a_to_be_aligned:
            sample_atoms.append(sample_res['CA'])
    for sample_res in sample_model['%s' % mhc_b]:
        if sample_res.get_id()[1] in mhc_b_to_be_aligned:
            sample_atoms.append(sample_res['CA'])
        # Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCII binding groove
    mhcii_com = center_of_mass(sample_atoms, geometric=True)
    # Calculate CoM of vTCR
    tcr_a_atoms_for_com = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_atoms_for_com = range(tcr_b_init, tcr_b_final + 1)
    tcr_atoms_for_com = []
    for tcr_res in sample_model['%s' % tcr_a]:
        if tcr_res.get_id()[1] in tcr_a_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    for tcr_res in sample_model['%s' % tcr_b]:
        if tcr_res.get_id()[1] in tcr_b_atoms_for_com:
            tcr_atoms_for_com.append(tcr_res['CA'])
    # vtcr_com =center_of_mass(tcr_atoms_for_com)
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)

    # Save the aligned version of pdb
    ali_structure = add_com_to_pdb(mhcii_com, vtcr_com, sample_structure)
    io = Bio.PDB.PDBIO()
    io.set_structure(ali_structure)
    io.save("%s_aligned.pdb"%pdbid)

    # Geometrical parameters
    dx, dy, dz = (np.subtract(vtcr_com, mhcii_com))
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhcii_com))))
    r = round(r, 2)
    theta = np.degrees(np.arctan2(dy, dx))
    theta = round(theta, 2)
    phi = np.degrees(np.arccos(dz / r))
    phi = round(phi, 2)
    print('The Geometrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
    return r, theta, phi

args = parser.parse_args()
pdbid = args.pdbid
mhc_a = args.mhc_a
mhc_b = args.mhc_b
tcr_a = args.tcr_a
tcr_b = args.tcr_b


"""Function to check inputs"""
def is_chain_in_pdb(pdbid, input_chain_id):
    from Bio.PDB import PDBParser
    import sys
    pdb_chain_ids_list = []
    structure = PDBParser().get_structure('%s'%pdbid, '%s.pdb'%pdbid)
    model = structure[0]
    [pdb_chain_ids_list.append(str(chain.get_id())) for chain in model]
    if input_chain_id not in pdb_chain_ids_list: sys.exit('Chain "%s" is not found in "%s.pdb"!!!!'%(input_chain_id, pdbid))
def number_of_chains(pdbid):
    from Bio.PDB import PDBParser
    import sys
    pdb_chain_ids_list = []
    structure = PDBParser().get_structure('%s' % pdbid, '%s.pdb' % pdbid)
    model = structure[0]
    counter = 0
    for chain in model: counter+=1
    return (counter)


input_chain_IDs = list(filter(None, [mhc_a, mhc_b, tcr_a, tcr_b]))
[is_chain_in_pdb(pdbid, input_chain_id) for input_chain_id in input_chain_IDs]
if number_of_chains(pdbid) not in (4,5): sys.exit('"%s.pdb" contains unexpected number of chains!!!!'%(pdbid))
tcr_mhci_pos(pdbid, mhc_a, tcr_a, tcr_b) if mhc_b is None else tcr_mhcii_pos(pdbid, mhc_a, mhc_b, tcr_a, tcr_b)
