#!/usr/bin/env python
## - Reference-dependent approach
"""Import modules"""
from __future__ import print_function
import sys
import os
import argparse
import numpy as np
import Bio.PDB
from Bio.PDB import Entity, Chain, Residue, Atom, PDBParser
from Bio.Cluster import pca
from itertools import chain
import time

# from sklearn.decomposition import PCA
# import mdtraj as md

"""
Define arguments
"""
parser = argparse.ArgumentParser(
    description="TCR-CoM calculates the geometrical parameters (r, theta, phi) of T cell receptors (TCR) on the top of MHC proteins"
)
parser.add_argument("-pdbid", type=str, help="PDB-ID")
parser.add_argument("-mhc_a", type=str, help="ID of MHC alpha chain")
parser.add_argument("-mhc_b", type=str, default=None, help="ID of MHC beta chain")
parser.add_argument("-tcr_a", type=str, help="ID of TCR alpha chain")
parser.add_argument("-tcr_b", type=str, help="ID of TCR beta chain")
parser.add_argument("-output_pdb", type=str, help="Output_PDB_structure")

args = parser.parse_args()
pdbid = args.pdbid
if pdbid.endswith(".pdb"):
    pdbid = pdbid.split(".")[0]
else:
    pdbid = pdbid
mhc_a = args.mhc_a
mhc_b = args.mhc_b
tcr_a = args.tcr_a
tcr_b = args.tcr_b
output_pdb = args.output_pdb

"""
Logging info
"""
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
output_name = "%s_tcr_com_%s.log" % (current_time, pdbid)
output_file = open(output_name, "a")
sys.stdout = output_file

"""
Define reference structures
"""
reffile1 = "./dependancies/ref_files/ref1.pdb"
reffile2 = "./dependancies/ref_files/ref2.pdb"

"""
Function to check inputs
"""


def add_com_to_pdb(mhc_com, vtcr_com, sample_structure):
    """
    Function to add pseudoatoms at MHC-CoM, TCR-CoM, and XYZ axise to the output PDB file
    """
    # mhc_com
    mhc_com_chain = "X"
    sample_structure.add(Chain.Chain(mhc_com_chain))
    res_id = (" ", 1, " ")
    new_residue = Residue.Residue(res_id, "MCM", " ")
    new_atom = Atom.Atom("C", mhc_com, 0, 0.0, " ", "C", 1, "C")
    new_residue.add(new_atom)
    sample_structure.child_dict[mhc_com_chain].add(new_residue)
    # tcr com
    tcr_com_chain = "Y"
    sample_structure.add(Chain.Chain(tcr_com_chain))
    res_id = (" ", 1, " ")
    new_residue = Residue.Residue(res_id, "TCM", " ")
    new_atom = Atom.Atom("C", vtcr_com, 0, 0.0, " ", "C", 1, "C")
    new_residue.add(new_atom)
    sample_structure.child_dict[tcr_com_chain].add(new_residue)
    # X,Y,Z atoms
    pos = [[50, 0, 0], [0, 50, 0], [0, 0, 50]]
    resn = ["X", "Y", "Z"]
    xyz_chain = "Z"
    sample_structure.add(Chain.Chain(xyz_chain))
    for i in [0, 1, 2]:
        res_id = (" ", i + 1, " ")
        new_residue = Residue.Residue(res_id, resn[i], " ")
        new_atom = Atom.Atom("O", pos[i], 0, 0.0, " ", "O", 1, "O")
        new_residue.add(new_atom)
        sample_structure.child_dict[xyz_chain].add(new_residue)
    return sample_structure


def center_of_mass(entity, geometric=False):
    # Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
    # This code is part of the Biopython distribution and governed by its
    # license.  Please see the LICENSE file that should have been included
    # as part of this package.
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """

    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, "__iter__") and [x for x in entity if x.level == "A"]:
        atom_list = entity
    else:  # Some other weirdo object
        raise ValueError(
            "Center of Mass can only be calculated from the following objects:\n"
            "Structure, Model, Chain, Residue, list of Atoms."
        )

    masses = []
    positions = [[], [], []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if "ukn" in set(masses) and not geometric:
        raise ValueError(
            "Some Atoms don't have an element assigned.\n"
            "Try adding them manually or calculate the geometrical center of mass instead"
        )

    if geometric:
        return [sum(coord_list) / len(masses) for coord_list in positions]
    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index] * atom_mass)
            w_pos[1].append(positions[1][atom_index] * atom_mass)
            w_pos[2].append(positions[2][atom_index] * atom_mass)

        return [sum(coord_list) / sum(masses) for coord_list in w_pos]


def fetch_atoms(model, selection="A", atom_bounds=[1, 180]):
    """
    Function to fetch atoms from the defined "atom_bound" in "selection" of the "model"
    """
    selection = selection.upper()
    if not isinstance(atom_bounds, (list, tuple)) and len(atom_bounds) > 0:
        raise ValueError("expected non-empty list or tuple, got {}".format(atom_bounds))
    # making sure its a list of lists
    if not isinstance(atom_bounds[0], (list, tuple)):
        atom_bounds = [atom_bounds]
    if not all([len(b) == 2 for b in atom_bounds]):
        raise ValueError(
            "All bounds must be providing one upper "
            "and one lower bound, got {}".format(atom_bounds)
        )
    if not isinstance(selection, (tuple, list)):
        selection = [selection]
    result = []
    for sel in selection:
        for ref_res in model["%s" % sel]:
            resid = ref_res.get_id()[1]
            in_bounds = False
            for bounds in atom_bounds:
                in_bounds |= bounds[0] <= resid and resid <= bounds[1]
            if in_bounds:
                result.append(ref_res["CA"])
    return result


def fetch_entity(model, fetch_atoms=True, selection="A", res_ids=range(1, 180)):
    """
    Function to fetch atoms/resids from the defined "resid_bounds" in "selection" of the "model"
    """
    selection = selection.upper()
    # fetch atoms
    if fetch_atoms is True:
        result = []
        for sel in selection:
            for sample_res in model["%s" % sel]:
                resid = sample_res.get_id()[1]
                if resid in res_ids:
                    result.append(sample_res["CA"])
    # fetch_residues_indeces
    elif fetch_atoms is False:
        result = []
        for sel in selection:
            for sample_res in model["%s" % sel]:
                resid = sample_res.get_id()[1]
                if resid in res_ids:
                    result.append(resid)
    return result


def apply_transformation_to_atoms(model, rotmat, transvec):
    """
    Function to translate/rotate the model by the defined translation vector and rotation matrix
    """
    for chain in model:
        for res in chain:
            for atom in res:
                atom.transform(rotmat, transvec)


def is_chain_in_pdb(pdbid, input_chain_id):
    """
    Function to check if the input_chain_id is in pdb_chain_ids_list
    """

    pdb_chain_ids_list = []
    structure = PDBParser().get_structure("%s" % pdbid, "%s.pdb" % pdbid)
    model = structure[0]
    for chain in model:
        pdb_chain_ids_list.append(str(chain.get_id()))
    return input_chain_id in pdb_chain_ids_list


def number_of_chains(pdbid):

    pdb_chain_ids_list = []
    structure = PDBParser().get_structure("%s" % pdbid, "%s.pdb" % pdbid)
    model = structure[0]
    counter = 0
    for chain in model:
        counter += 1
    return counter


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v is None:
        return True
    if v.lower() in ("TRUE", "True", "yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("FALSE", "False", "no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


"""
Geometrical parameters of TCR-MHC class I
"""


def tcr_mhci_geometrical_parameters(
    pdbid,
    mhc_a="A",
    tcr_a="D",
    tcr_b="E",
    persist_structure=True,
    mhc_a_init=1,
    mhc_a_final=179,
    tcr_a_init=1,
    tcr_a_final=109,
    tcr_b_init=1,
    tcr_b_final=116,
):
    """
    ARGUMENTS:
    pdbid = PDB-ID or name of the input PB file (str)
    mhc_a = ID of MHC alpha chain (str)
    tcr_a = ID of TCR alpha chain (str)
    tcr_b = ID of TCR beta chain (str)
    persist_structure = If true, save the processed structure as PDB-file (Boolean)
    """
    #################################################################
    # Define residues range to align and center of mass calculations#
    #################################################################
    mhc_a_resids = range(mhc_a_init, mhc_a_final + 1)
    tcr_a_resids = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_resids = range(tcr_b_init, tcr_b_final + 1)

    ########################################################################################################
    # Import structure, align to reference, and calculate center of mass of CA atoms in MHCI binding groove#
    ########################################################################################################
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reffile1)
    sample_structure = pdb_parser.get_structure("sample", "./%s.pdb" % pdbid)
    # Use the first model in the pdb-files for alignment
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    # Iterate of all residues in each model in order to define proper atoms
    # Sample structure
    sample_resids = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_a, res_ids=mhc_a_resids
    )
    sample_atoms = fetch_entity(
        sample_model, fetch_atoms=True, selection=mhc_a, res_ids=sample_resids
    )
    # Reference structure
    ref_atoms = fetch_entity(
        ref_model, fetch_atoms=True, selection=mhc_a, res_ids=sample_resids
    )

    # Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCI binding groove
    mhci_com = center_of_mass(sample_atoms, geometric=True)
    # Calculate CoM of vTCR
    tcr_atoms_for_com = fetch_entity(
        sample_model, fetch_atoms=True, selection=tcr_a, res_ids=tcr_a_resids
    )
    tcr_atoms_for_com += fetch_entity(
        sample_model, fetch_atoms=True, selection=tcr_b, res_ids=tcr_b_resids
    )
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)
    print("MHC-CoM: ", [round(x, 2) for x in mhci_com])
    print("vTCR-CoM: ", [round(x, 2) for x in vtcr_com])
    #######################
    # Save final structure#
    #######################
    if persist_structure:
        ali_structure = add_com_to_pdb(mhci_com, vtcr_com, sample_model)
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure)
        io.save("%s_aligned.pdb" % pdbid)

    ###################################
    # Calculate geometrical parameters#
    ###################################
    dx, dy, dz = np.subtract(vtcr_com, mhci_com)
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhci_com))))
    theta = np.degrees(np.arctan2(dy, dx))
    phi = np.degrees(np.arccos(dz / r))
    print(
        "The Geomitrical parameters: r = {:.2f}, "
        "theta = {:.2f}, phi = {:.2f}".format(r, theta, phi)
    )
    return r, theta, phi


# tcr_mhci_pos('1ao7', mhc_a='A', tcr_a='D', tcr_b='E')

# -----------
"""
Geometrical parameters of TCR-MHC class II
"""


def tcr_mhcii_geometrical_parameters(
    pdbid,
    mhc_a="A",
    mhc_b="B",
    tcr_a="D",
    tcr_b="E",
    persist_structure=True,
    mhc_a_init=1,
    mhc_a_final=80,
    mhc_b_init=1,
    mhc_b_final=90,
    tcr_a_init=1,
    tcr_a_final=109,
    tcr_b_init=1,
    tcr_b_final=116,
):
    """
    ARGUMENTS:
    pdbid = PDB-ID or name of the input PB file (str)
    mhc_a = ID of MHC alpha chain (str)
    mhc_b = ID of MHC beta chain (str)
    tcr_a = ID of TCR alpha chain (str)
    tcr_b = ID of TCR beta chain (str)
    persist_structure = If true, save the processed structure as PDB-file (Boolean)
    """
    #################################################################
    # Define residues range to align and center of mass calculations#
    #################################################################
    mhc_a_resids = range(mhc_a_init, mhc_a_final + 1)
    mhc_b_resids = range(mhc_b_init, mhc_b_final + 1)
    tcr_a_resids = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_resids = range(tcr_b_init, tcr_b_final + 1)

    #########################################################################################################
    # Import structure, align to reference, and calculate center of mass of CA atoms in MHCII binding groove#
    #########################################################################################################
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reffile2)
    sample_structure = pdb_parser.get_structure("sample", "./%s.pdb" % pdbid)
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    # Iterate of all residues in each model in order to define proper atoms
    # Sample structure
    sample_resids_a = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_a, res_ids=mhc_a_resids
    )
    sample_resids_b = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_b, res_ids=mhc_b_resids
    )

    sample_atoms = fetch_entity(sample_model, selection=mhc_a, res_ids=sample_resids_a)
    sample_atoms += fetch_entity(sample_model, selection=mhc_b, res_ids=sample_resids_b)
    # Reference structure
    ref_atoms = fetch_entity(ref_model, selection=mhc_a, res_ids=sample_resids_a)
    ref_atoms += fetch_entity(ref_model, selection=mhc_b, res_ids=sample_resids_b)

    # Initiate the superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCII binding groove
    mhcii_com = center_of_mass(sample_atoms, geometric=True)
    # Calculate CoM of vTCR
    tcr_atoms_for_com = fetch_entity(
        sample_model, selection=tcr_a, res_ids=tcr_a_resids
    )
    tcr_atoms_for_com += fetch_entity(
        sample_model, selection=tcr_b, res_ids=tcr_b_resids
    )
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)
    print("MHC-CoM: ", [round(x, 2) for x in mhcii_com])
    print("vTCR-CoM: ", [round(x, 2) for x in vtcr_com])
    #######################
    # Save final structure#
    #######################
    if persist_structure:
        ali_structure = add_com_to_pdb(mhcii_com, vtcr_com, sample_model)
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure)
        io.save("%s_aligned.pdb" % pdbid)

    ###################################
    # Calculate geometrical parameters#
    ###################################
    dx, dy, dz = np.subtract(vtcr_com, mhcii_com)
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhcii_com))))
    theta = np.degrees(np.arctan2(dy, dx))
    phi = np.degrees(np.arccos(dz / r))
    print(
        "The Geomitrical parameters: r = {:.2f}, "
        "theta = {:.2f}, phi = {:.2f}".format(r, theta, phi)
    )
    return r, theta, phi


# tcr_mhcii_pos('1FYT', mhc_a='A', mhc_b='B', tcr_a='D', tcr_b='E')

# ------------

# For Python script
print("PDB-file: %s" % pdbid)
print("MHC alpha chain-ID: %s" % mhc_a)
print("MHC beta chain-ID: %s" % mhc_b)
print("TCR alpha chain-ID: %s" % tcr_a)
print("TCR beta chain-ID: %s" % tcr_b)

input_chain_IDs = list(filter(None, [mhc_a, mhc_b, tcr_a, tcr_b]))
input_chain_IDs_upper = [x.upper() for x in input_chain_IDs]

for input_chain_id in input_chain_IDs_upper:
    if not is_chain_in_pdb(pdbid, input_chain_id):
        raise ValueError(
            'Chain "%s" is not found in "%s.pdb"!' % (input_chain_id, pdbid)
        )
n_chains = number_of_chains(pdbid)
if n_chains < 3 or n_chains > 5:
    raise ValueError(
        "The submitted PDB file contains unexpected number of chains! (expected 5 chains)"
    )
if mhc_b is None:
    tcr_mhci_geometrical_parameters(pdbid, mhc_a, tcr_a, tcr_b, str2bool(output_pdb))
else:
    tcr_mhcii_geometrical_parameters(
        pdbid, mhc_a, mhc_b, tcr_a, tcr_b, str2bool(output_pdb)
    )
