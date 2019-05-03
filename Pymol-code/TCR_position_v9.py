#!/usr/python
'''
For questions about running the script or for reporting bugs, please contact either:
Brian Baker (brian-baker[at]nd.edu)
or
Esam Abualrous (e.abualrous[at]fu-berlin.de)

DESCRIPTION

   Calculates the geometrical parameters (r, theta, phi) of TCR on top of MHC proteins

USAGE
For peptide-MHC class I-TCR complex:
tcr_mhci_pos selection [,mhc_a=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=180 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]

For peptide-MHC class II-TCR complex:
tcr_mhcii_pos selection [,mhc_a=None [,mhc_b=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=80 [,mhc_b_init=1 [,mhc_b_final=90 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]]]]

NOTES
- The code requires defining at least 3 chain IDs for MHCI-TCR or 4 chain IDs for MHCII-TCR
- If no range of residues is defined, the code uses the default values (see USAGE)
- Please make sure that the defined chains do not contain any solvent molecules
'''
from __future__ import print_function
from pymol import cmd, stored
import sys
import os
import numpy as np


# Geometrical parameters of TCR-MHC class I
def tcr_mhci_pos(selection, mhc_a=None, tcr_a=None, tcr_b=None,
                 mhc_a_init=1, mhc_a_final=180,
                 tcr_a_init=1, tcr_a_final=109, tcr_b_init=1, tcr_b_final=116, state=None):
    if selection not in cmd.get_object_list('all'):
        print('unknown object value given (%s)'%selection)
    else:
    # MHC class I
        cmd.load(_get_reference_temporary_file(reference=1), 'refi')
        cmd.super('/%s//%s/%s-%s' % (selection, mhc_a, mhc_a_init, mhc_a_final), '/refi//A')
        mhci_com = get_com('/%s//%s/%s-%s' % (selection, mhc_a, mhc_a_init, mhc_a_final))
        cmd.delete('refi')
        # T cell receptor
        vtcr_com = get_com('/%s//%s/%s-%s or /%s//%s/%s-%s' % (
        selection, tcr_a, tcr_a_init, tcr_a_final, selection, tcr_b, tcr_b_init, tcr_b_final))
        # Geomitrical parameters
        dx, dy, dz = (np.subtract(vtcr_com, mhci_com))
        r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhci_com))))
        theta = np.degrees(np.arctan2(dy, dx))
        phi = np.degrees(np.arccos(dz / r))
        print('The Geomitrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
        return r, theta, phi
cmd.extend("tcr_mhci_pos", tcr_mhci_pos)


# Geometrical parameters of TCR-MHC class II
def tcr_mhcii_pos(selection, mhc_a=None, mhc_b=None, tcr_a=None, tcr_b=None,
                  mhc_a_init=1, mhc_a_final=80, mhc_b_init=1, mhc_b_final=90,
                  tcr_a_init=1, tcr_a_final=109, tcr_b_init=1, tcr_b_final=116, state=None):
    if selection not in cmd.get_object_list('all'):
        print('unknown object value given (%s)' % selection)
    else:
        # MHC class II
        cmd.load(_get_reference_temporary_file(reference=2), 'refii')
        cmd.super('/%s//%s/%s-%s or /%s//%s/%s-%s' % (
            selection, mhc_a, mhc_a_init, mhc_a_final, selection, mhc_b, mhc_b_init, mhc_b_final),
                  '/refii//A or /refii//B')
        mhcii_com = get_com('/%s//%s/%s-%s or /%s//%s/%s-%s' % (
            selection, mhc_a, mhc_a_init, mhc_a_final, selection, mhc_b, mhc_b_init, mhc_b_final))
        cmd.delete('refii')
        # T cell receptor
        vtcr_com = get_com('/%s//%s/%s-%s or /%s//%s/%s-%s' % (
            selection, tcr_a, tcr_a_init, tcr_a_final, selection, tcr_b, tcr_b_init, tcr_b_final))
        # Geometrical parameters
        dx, dy, dz = np.subtract(vtcr_com, mhcii_com)
        r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhcii_com))))
        theta = np.degrees(np.arctan2(dy, dx))
        phi = np.degrees(np.arccos(dz / r))
        print('The Geometrical parameters: r = %s, theta = %s, phi = %s' % (r, theta, phi))
        return r, theta, phi
cmd.extend("tcr_mhcii_pos", tcr_mhcii_pos)


def get_com(selection, state=1, mass=None, quiet=1):
    """
    DESCRIPTION

    Calculates the center of mass
    Author: Sean Law
    Michigan State University
    slaw (at) msu . edu
    """
    quiet = int(quiet)

    totmass = 0.0
    if mass is not None and not quiet:
        print("Calculating mass-weighted COM")

    state = int(state)
    model = cmd.get_model(selection, state)
    x, y, z = 0, 0, 0
    for a in model.atom:
        if mass is not None:
            m = a.get_mass()
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
            totmass += m
        else:
            x += a.coord[0]
            y += a.coord[1]
            z += a.coord[2]

    if mass is not None:
        return x / totmass, y / totmass, z / totmass
    else:
        return x / len(model.atom), y / len(model.atom), z / len(model.atom)
# cmd.extend("get_com", get_com)

def _get_reference_temporary_file(reference=None):
    content_1 = \
    """ATOM      1  N   GLY A   1      13.597 -18.573 -11.059  1.00 33.27           N
ATOM      2  CA  GLY A   1      12.983 -17.304 -11.564  1.00 31.48           C
ATOM      3  C   GLY A   1      12.708 -16.314 -10.451  1.00 30.52           C
ATOM      4  O   GLY A   1      12.511 -16.706  -9.295  1.00 29.73           O
ATOM      5  N   SER A   2      12.699 -15.029 -10.807  1.00 29.85           N
ATOM      6  CA  SER A   2      12.416 -13.952  -9.856  1.00 28.39           C
ATOM      7  C   SER A   2      10.960 -13.950  -9.422  1.00 27.50           C
ATOM      8  O   SER A   2      10.092 -14.454 -10.140  1.00 26.91           O
ATOM      9  CB  SER A   2      12.763 -12.593 -10.451  1.00 27.26           C
ATOM     10  OG  SER A   2      14.156 -12.378 -10.426  1.00 28.53           O
ATOM     11  N   HIS A   3      10.714 -13.388  -8.238  1.00 26.22           N
ATOM     12  CA  HIS A   3       9.369 -13.177  -7.721  1.00 25.56           C
ATOM     13  C   HIS A   3       9.273 -11.853  -6.962  1.00 26.29           C
ATOM     14  O   HIS A   3      10.246 -11.404  -6.358  1.00 27.32           O
ATOM     15  CB  HIS A   3       8.946 -14.344  -6.833  1.00 24.38           C
ATOM     16  CG  HIS A   3       8.862 -15.652  -7.555  1.00 24.09           C
ATOM     17  CD2 HIS A   3       7.964 -16.118  -8.455  1.00 20.94           C
ATOM     18  ND1 HIS A   3       9.797 -16.653  -7.393  1.00 21.39           N
ATOM     19  CE1 HIS A   3       9.471 -17.684  -8.156  1.00 22.49           C
ATOM     20  NE2 HIS A   3       8.363 -17.384  -8.810  1.00 21.48           N
ATOM     21  N   SER A   4       8.098 -11.232  -6.994  1.00 26.00           N
ATOM     22  CA  SER A   4       7.892  -9.953  -6.321  1.00 24.61           C
ATOM     23  C   SER A   4       6.498  -9.835  -5.728  1.00 24.50           C
ATOM     24  O   SER A   4       5.537 -10.351  -6.295  1.00 25.36           O
ATOM     25  CB  SER A   4       8.112  -8.809  -7.306  1.00 24.62           C
ATOM     26  OG  SER A   4       7.136  -8.838  -8.330  1.00 21.59           O
ATOM     27  N   MET A   5       6.398  -9.149  -4.591  1.00 23.67           N
ATOM     28  CA  MET A   5       5.109  -8.809  -3.996  1.00 23.97           C
ATOM     29  C   MET A   5       4.778  -7.340  -4.224  1.00 23.84           C
ATOM     30  O   MET A   5       5.643  -6.474  -4.113  1.00 23.95           O
ATOM     31  CB  MET A   5       5.098  -9.108  -2.501  1.00 24.00           C
ATOM     32  CG  MET A   5       3.717  -9.029  -1.867  1.00 23.13           C
ATOM     33  SD  MET A   5       3.607  -9.981  -0.347  1.00 26.74           S
ATOM     34  CE  MET A   5       4.455  -8.899   0.813  1.00 25.36           C
ATOM     35  N   ARG A   6       3.523  -7.068  -4.554  1.00 24.07           N
ATOM     36  CA  ARG A   6       3.071  -5.705  -4.755  1.00 24.22           C
ATOM     37  C   ARG A   6       1.696  -5.498  -4.154  1.00 24.88           C
ATOM     38  O   ARG A   6       0.816  -6.353  -4.274  1.00 24.13           O
ATOM     39  CB  ARG A   6       3.020  -5.351  -6.242  1.00 23.98           C
ATOM     40  CG  ARG A   6       4.368  -5.190  -6.930  1.00 26.86           C
ATOM     41  CD  ARG A   6       5.163  -4.002  -6.398  1.00 31.28           C
ATOM     42  NE  ARG A   6       6.413  -3.848  -7.138  1.00 33.31           N
ATOM     43  CZ  ARG A   6       7.529  -4.534  -6.897  1.00 32.91           C
ATOM     44  NH1 ARG A   6       7.577  -5.438  -5.920  1.00 31.31           N1+
ATOM     45  NH2 ARG A   6       8.602  -4.315  -7.644  1.00 29.19           N
ATOM     46  N   TYR A   7       1.529  -4.352  -3.501  1.00 25.29           N
ATOM     47  CA  TYR A   7       0.209  -3.864  -3.136  1.00 25.26           C
ATOM     48  C   TYR A   7      -0.098  -2.601  -3.919  1.00 24.64           C
ATOM     49  O   TYR A   7       0.745  -1.718  -4.052  1.00 24.56           O
ATOM     50  CB  TYR A   7       0.101  -3.615  -1.632  1.00 24.56           C
ATOM     51  CG  TYR A   7      -0.088  -4.878  -0.836  1.00 23.39           C
ATOM     52  CD1 TYR A   7       1.000  -5.531  -0.257  1.00 23.39           C
ATOM     53  CD2 TYR A   7      -1.357  -5.435  -0.670  1.00 23.64           C
ATOM     54  CE1 TYR A   7       0.829  -6.705   0.476  1.00 20.08           C
ATOM     55  CE2 TYR A   7      -1.538  -6.609   0.053  1.00 20.71           C
ATOM     56  CZ  TYR A   7      -0.442  -7.234   0.627  1.00 23.00           C
ATOM     57  OH  TYR A   7      -0.619  -8.393   1.349  1.00 26.68           O
ATOM     58  N   PHE A   8      -1.309  -2.534  -4.449  1.00 24.74           N
ATOM     59  CA  PHE A   8      -1.748  -1.372  -5.194  1.00 25.32           C
ATOM     60  C   PHE A   8      -2.899  -0.712  -4.461  1.00 25.13           C
ATOM     61  O   PHE A   8      -3.881  -1.374  -4.115  1.00 25.57           O
ATOM     62  CB  PHE A   8      -2.183  -1.764  -6.611  1.00 25.56           C
ATOM     63  CG  PHE A   8      -1.068  -2.279  -7.473  1.00 25.12           C
ATOM     64  CD1 PHE A   8      -0.813  -3.638  -7.561  1.00 25.99           C
ATOM     65  CD2 PHE A   8      -0.286  -1.405  -8.214  1.00 27.52           C
ATOM     66  CE1 PHE A   8       0.218  -4.116  -8.364  1.00 27.66           C
ATOM     67  CE2 PHE A   8       0.741  -1.874  -9.025  1.00 26.10           C
ATOM     68  CZ  PHE A   8       0.995  -3.230  -9.096  1.00 25.35           C
ATOM     69  N   TYR A   9      -2.767   0.590  -4.226  1.00 25.12           N
ATOM     70  CA  TYR A   9      -3.808   1.385  -3.578  1.00 24.83           C
ATOM     71  C   TYR A   9      -4.402   2.400  -4.544  1.00 25.76           C
ATOM     72  O   TYR A   9      -3.685   3.051  -5.311  1.00 25.89           O
ATOM     73  CB  TYR A   9      -3.263   2.138  -2.356  1.00 23.34           C
ATOM     74  CG  TYR A   9      -2.797   1.284  -1.192  1.00 22.51           C
ATOM     75  CD1 TYR A   9      -3.075  -0.085  -1.129  1.00 19.24           C
ATOM     76  CD2 TYR A   9      -2.108   1.861  -0.126  1.00 22.30           C
ATOM     77  CE1 TYR A   9      -2.659  -0.854  -0.045  1.00 17.81           C
ATOM     78  CE2 TYR A   9      -1.689   1.096   0.961  1.00 21.41           C
ATOM     79  CZ  TYR A   9      -1.969  -0.257   0.993  1.00 19.85           C
ATOM     80  OH  TYR A   9      -1.552  -1.012   2.064  1.00 20.25           O
ATOM     81  N   THR A  10      -5.721   2.537  -4.494  1.00 26.49           N
ATOM     82  CA  THR A  10      -6.396   3.624  -5.183  1.00 27.79           C
ATOM     83  C   THR A  10      -7.354   4.306  -4.210  1.00 28.43           C
ATOM     84  O   THR A  10      -8.252   3.668  -3.661  1.00 29.60           O
ATOM     85  CB  THR A  10      -7.130   3.126  -6.447  1.00 27.28           C
ATOM     86  CG2 THR A  10      -7.706   4.289  -7.229  1.00 28.87           C
ATOM     87  OG1 THR A  10      -6.204   2.432  -7.289  1.00 27.42           O
ATOM     88  N   ALA A  11      -7.133   5.596  -3.976  1.00 29.12           N
ATOM     89  CA  ALA A  11      -8.017   6.390  -3.132  1.00 31.36           C
ATOM     90  C   ALA A  11      -8.665   7.480  -3.972  1.00 33.04           C
ATOM     91  O   ALA A  11      -7.989   8.408  -4.417  1.00 32.40           O
ATOM     92  CB  ALA A  11      -7.251   6.995  -1.966  1.00 31.69           C
ATOM     93  N   MET A  12      -9.970   7.350  -4.204  1.00 35.20           N
ATOM     94  CA  MET A  12     -10.698   8.290  -5.051  1.00 37.36           C
ATOM     95  C   MET A  12     -11.666   9.140  -4.232  1.00 38.53           C
ATOM     96  O   MET A  12     -12.537   8.606  -3.540  1.00 39.38           O
ATOM     97  CB  MET A  12     -11.473   7.565  -6.163  1.00 37.55           C
ATOM     98  CG  MET A  12     -10.728   6.436  -6.868  1.00 41.48           C
ATOM     99  SD  MET A  12     -11.269   4.782  -6.352  1.00 47.78           S
ATOM    100  CE  MET A  12     -12.837   4.640  -7.208  1.00 45.62           C
ATOM    101  N   SER A  13     -11.512  10.460  -4.308  1.00 38.84           N
ATOM    102  CA  SER A  13     -12.489  11.375  -3.727  1.00 39.07           C
ATOM    103  C   SER A  13     -13.644  11.535  -4.706  1.00 40.34           C
ATOM    104  O   SER A  13     -13.452  11.455  -5.922  1.00 40.72           O
ATOM    105  CB  SER A  13     -11.862  12.733  -3.411  1.00 38.66           C
ATOM    106  OG  SER A  13     -11.569  13.453  -4.594  1.00 37.50           O
ATOM    107  N   ARG A  14     -14.843  11.751  -4.176  1.00 41.27           N
ATOM    108  CA  ARG A  14     -16.044  11.832  -4.999  1.00 41.30           C
ATOM    109  C   ARG A  14     -17.096  12.744  -4.360  1.00 42.66           C
ATOM    110  O   ARG A  14     -18.180  12.279  -3.989  1.00 42.64           O
ATOM    111  CB  ARG A  14     -16.605  10.426  -5.246  1.00 41.10           C
ATOM    112  CG  ARG A  14     -16.782   9.593  -3.979  1.00 37.94           C
ATOM    113  CD  ARG A  14     -17.106   8.152  -4.295  1.00 34.59           C
ATOM    114  NE  ARG A  14     -17.487   7.429  -3.086  1.00 32.05           N
ATOM    115  CZ  ARG A  14     -17.625   6.109  -3.006  1.00 28.52           C
ATOM    116  NH1 ARG A  14     -17.413   5.350  -4.069  1.00 28.43           N1+
ATOM    117  NH2 ARG A  14     -17.975   5.552  -1.859  1.00 26.88           N
ATOM    118  N   PRO A  15     -16.790  14.052  -4.241  1.00 43.66           N
ATOM    119  CA  PRO A  15     -17.729  14.944  -3.561  1.00 44.38           C
ATOM    120  C   PRO A  15     -19.016  15.081  -4.362  1.00 45.09           C
ATOM    121  O   PRO A  15     -18.994  14.990  -5.593  1.00 45.29           O
ATOM    122  CB  PRO A  15     -16.982  16.281  -3.510  1.00 44.57           C
ATOM    123  CG  PRO A  15     -15.556  15.965  -3.832  1.00 44.44           C
ATOM    124  CD  PRO A  15     -15.610  14.782  -4.735  1.00 43.79           C
ATOM    125  N   GLY A  16     -20.128  15.289  -3.660  1.00 45.58           N
ATOM    126  CA  GLY A  16     -21.447  15.300  -4.286  1.00 45.98           C
ATOM    127  C   GLY A  16     -22.079  13.921  -4.262  1.00 46.36           C
ATOM    128  O   GLY A  16     -23.304  13.787  -4.349  1.00 46.92           O
ATOM    129  N   ARG A  17     -21.233  12.899  -4.144  1.00 46.08           N
ATOM    130  CA  ARG A  17     -21.684  11.514  -3.992  1.00 45.11           C
ATOM    131  C   ARG A  17     -21.497  11.051  -2.545  1.00 44.20           C
ATOM    132  O   ARG A  17     -22.390  10.420  -1.972  1.00 44.30           O
ATOM    133  CB  ARG A  17     -20.966  10.587  -4.984  1.00 45.25           C
ATOM    134  CG  ARG A  17     -21.003  11.096  -6.438  1.00 45.77           C
ATOM    135  CD  ARG A  17     -20.844   9.985  -7.476  1.00 44.77           C
ATOM    136  NE  ARG A  17     -19.511   9.378  -7.471  1.00 43.77           N
ATOM    137  CZ  ARG A  17     -19.003   8.664  -8.473  1.00 43.77           C
ATOM    138  NH1 ARG A  17     -19.704   8.463  -9.586  1.00 41.58           N1+
ATOM    139  NH2 ARG A  17     -17.784   8.152  -8.367  1.00 45.24           N
ATOM    140  N   GLY A  18     -20.349  11.389  -1.956  1.00 43.05           N
ATOM    141  CA  GLY A  18     -20.080  11.108  -0.546  1.00 41.27           C
ATOM    142  C   GLY A  18     -18.605  10.959  -0.217  1.00 39.98           C
ATOM    143  O   GLY A  18     -17.756  11.599  -0.841  1.00 40.37           O
ATOM    144  N   GLU A  19     -18.309  10.114   0.770  1.00 37.98           N
ATOM    145  CA  GLU A  19     -16.940   9.888   1.237  1.00 36.30           C
ATOM    146  C   GLU A  19     -16.110   9.156   0.186  1.00 34.63           C
ATOM    147  O   GLU A  19     -16.651   8.366  -0.587  1.00 33.65           O
ATOM    148  CB  GLU A  19     -16.923   9.093   2.554  1.00 36.91           C
ATOM    149  CG  GLU A  19     -17.978   9.480   3.585  1.00 39.47           C
ATOM    150  CD  GLU A  19     -18.043  10.973   3.844  1.00 41.93           C
ATOM    151  OE1 GLU A  19     -17.294  11.461   4.721  1.00 40.79           O
ATOM    152  OE2 GLU A  19     -18.860  11.649   3.175  1.00 41.63           O1-
ATOM    153  N   PRO A  20     -14.789   9.416   0.159  1.00 34.07           N
ATOM    154  CA  PRO A  20     -13.875   8.760  -0.780  1.00 33.61           C
ATOM    155  C   PRO A  20     -13.889   7.230  -0.720  1.00 32.89           C
ATOM    156  O   PRO A  20     -14.266   6.640   0.295  1.00 32.40           O
ATOM    157  CB  PRO A  20     -12.505   9.278  -0.350  1.00 33.79           C
ATOM    158  CG  PRO A  20     -12.798  10.598   0.286  1.00 35.21           C
ATOM    159  CD  PRO A  20     -14.075  10.375   1.022  1.00 34.30           C
ATOM    160  N   ARG A  21     -13.483   6.612  -1.824  1.00 32.37           N
ATOM    161  CA  ARG A  21     -13.401   5.164  -1.931  1.00 31.77           C
ATOM    162  C   ARG A  21     -11.943   4.726  -1.871  1.00 30.38           C
ATOM    163  O   ARG A  21     -11.065   5.339  -2.488  1.00 30.38           O
ATOM    164  CB  ARG A  21     -14.060   4.697  -3.231  1.00 31.05           C
ATOM    165  CG  ARG A  21     -14.268   3.189  -3.367  1.00 31.90           C
ATOM    166  CD  ARG A  21     -15.033   2.901  -4.663  1.00 33.79           C
ATOM    167  NE  ARG A  21     -14.899   1.527  -5.147  1.00 37.92           N
ATOM    168  CZ  ARG A  21     -15.311   1.108  -6.345  1.00 42.08           C
ATOM    169  NH1 ARG A  21     -15.882   1.952  -7.202  1.00 43.04           N1+
ATOM    170  NH2 ARG A  21     -15.148  -0.162  -6.696  1.00 41.78           N
ATOM    171  N   PHE A  22     -11.691   3.674  -1.103  1.00 28.86           N
ATOM    172  CA  PHE A  22     -10.363   3.097  -1.000  1.00 27.58           C
ATOM    173  C   PHE A  22     -10.371   1.673  -1.559  1.00 27.22           C
ATOM    174  O   PHE A  22     -11.150   0.824  -1.112  1.00 26.25           O
ATOM    175  CB  PHE A  22      -9.895   3.115   0.458  1.00 27.26           C
ATOM    176  CG  PHE A  22      -8.563   2.466   0.675  1.00 27.56           C
ATOM    177  CD1 PHE A  22      -7.385   3.146   0.372  1.00 27.48           C
ATOM    178  CD2 PHE A  22      -8.481   1.171   1.186  1.00 25.52           C
ATOM    179  CE1 PHE A  22      -6.145   2.542   0.566  1.00 29.56           C
ATOM    180  CE2 PHE A  22      -7.246   0.559   1.383  1.00 24.71           C
ATOM    181  CZ  PHE A  22      -6.078   1.243   1.074  1.00 27.57           C
ATOM    182  N   ILE A  23      -9.523   1.421  -2.550  1.00 27.19           N
ATOM    183  CA  ILE A  23      -9.369   0.067  -3.092  1.00 27.41           C
ATOM    184  C   ILE A  23      -7.925  -0.388  -2.967  1.00 25.77           C
ATOM    185  O   ILE A  23      -7.001   0.338  -3.332  1.00 25.97           O
ATOM    186  CB  ILE A  23      -9.862  -0.061  -4.561  1.00 27.47           C
ATOM    187  CG1 ILE A  23     -11.309   0.429  -4.681  1.00 27.69           C
ATOM    188  CG2 ILE A  23      -9.800  -1.524  -5.022  1.00 28.08           C
ATOM    189  CD1 ILE A  23     -11.814   0.541  -6.106  1.00 29.19           C
ATOM    190  N   ALA A  24      -7.746  -1.587  -2.425  1.00 24.83           N
ATOM    191  CA  ALA A  24      -6.423  -2.158  -2.213  1.00 23.63           C
ATOM    192  C   ALA A  24      -6.341  -3.539  -2.832  1.00 22.96           C
ATOM    193  O   ALA A  24      -7.180  -4.399  -2.550  1.00 22.35           O
ATOM    194  CB  ALA A  24      -6.119  -2.233  -0.736  1.00 23.35           C
ATOM    195  N   VAL A  25      -5.339  -3.744  -3.681  1.00 22.18           N
ATOM    196  CA  VAL A  25      -5.100  -5.061  -4.281  1.00 21.66           C
ATOM    197  C   VAL A  25      -3.659  -5.534  -4.085  1.00 21.01           C
ATOM    198  O   VAL A  25      -2.710  -4.763  -4.229  1.00 21.22           O
ATOM    199  CB  VAL A  25      -5.511  -5.130  -5.783  1.00 21.59           C
ATOM    200  CG1 VAL A  25      -7.025  -5.065  -5.925  1.00 21.62           C
ATOM    201  CG2 VAL A  25      -4.861  -4.025  -6.593  1.00 21.50           C
ATOM    202  N   GLY A  26      -3.512  -6.808  -3.734  1.00 20.09           N
ATOM    203  CA  GLY A  26      -2.200  -7.410  -3.524  1.00 18.58           C
ATOM    204  C   GLY A  26      -1.874  -8.405  -4.614  1.00 18.24           C
ATOM    205  O   GLY A  26      -2.732  -9.190  -5.029  1.00 17.90           O
ATOM    206  N   TYR A  27      -0.634  -8.361  -5.085  1.00 18.03           N
ATOM    207  CA  TYR A  27      -0.166  -9.264  -6.121  1.00 19.17           C
ATOM    208  C   TYR A  27       1.132  -9.939  -5.704  1.00 19.49           C
ATOM    209  O   TYR A  27       2.009  -9.307  -5.108  1.00 18.38           O
ATOM    210  CB  TYR A  27       0.075  -8.511  -7.433  1.00 19.94           C
ATOM    211  CG  TYR A  27      -1.160  -8.183  -8.239  1.00 23.36           C
ATOM    212  CD1 TYR A  27      -2.079  -7.236  -7.794  1.00 23.71           C
ATOM    213  CD2 TYR A  27      -1.395  -8.794  -9.474  1.00 26.17           C
ATOM    214  CE1 TYR A  27      -3.213  -6.927  -8.535  1.00 23.14           C
ATOM    215  CE2 TYR A  27      -2.529  -8.484 -10.226  1.00 24.27           C
ATOM    216  CZ  TYR A  27      -3.431  -7.550  -9.746  1.00 22.65           C
ATOM    217  OH  TYR A  27      -4.555  -7.220 -10.476  1.00 25.65           O
ATOM    218  N   VAL A  28       1.230 -11.232  -6.001  1.00 19.87           N
ATOM    219  CA  VAL A  28       2.519 -11.904  -6.105  1.00 20.15           C
ATOM    220  C   VAL A  28       2.701 -12.142  -7.595  1.00 20.97           C
ATOM    221  O   VAL A  28       1.853 -12.779  -8.228  1.00 21.52           O
ATOM    222  CB  VAL A  28       2.561 -13.235  -5.314  1.00 20.89           C
ATOM    223  CG1 VAL A  28       3.687 -14.145  -5.820  1.00 21.49           C
ATOM    224  CG2 VAL A  28       2.712 -12.969  -3.823  1.00 18.47           C
ATOM    225  N   ASP A  29       3.785 -11.600  -8.154  1.00 21.39           N
ATOM    226  CA  ASP A  29       4.046 -11.649  -9.595  1.00 21.91           C
ATOM    227  C   ASP A  29       2.834 -11.132 -10.381  1.00 21.50           C
ATOM    228  O   ASP A  29       2.411  -9.995 -10.190  1.00 21.44           O
ATOM    229  CB  ASP A  29       4.426 -13.071 -10.047  1.00 21.90           C
ATOM    230  CG  ASP A  29       5.579 -13.651  -9.254  1.00 26.15           C
ATOM    231  OD1 ASP A  29       6.629 -12.983  -9.127  1.00 28.08           O
ATOM    232  OD2 ASP A  29       5.430 -14.789  -8.759  1.00 32.66           O1-
ATOM    233  N   ASP A  30       2.282 -11.974 -11.252  1.00 20.76           N
ATOM    234  CA  ASP A  30       1.090 -11.630 -12.014  1.00 20.79           C
ATOM    235  C   ASP A  30      -0.179 -12.285 -11.437  1.00 21.02           C
ATOM    236  O   ASP A  30      -1.199 -12.381 -12.124  1.00 20.54           O
ATOM    237  CB  ASP A  30       1.279 -11.991 -13.491  1.00 20.03           C
ATOM    238  CG  ASP A  30       2.207 -11.029 -14.217  1.00 23.26           C
ATOM    239  OD1 ASP A  30       2.180  -9.818 -13.910  1.00 26.93           O
ATOM    240  OD2 ASP A  30       2.955 -11.481 -15.114  1.00 22.92           O1-
ATOM    241  N   THR A  31      -0.110 -12.720 -10.179  1.00 21.26           N
ATOM    242  CA  THR A  31      -1.263 -13.303  -9.486  1.00 22.90           C
ATOM    243  C   THR A  31      -1.768 -12.411  -8.356  1.00 23.57           C
ATOM    244  O   THR A  31      -1.078 -12.189  -7.353  1.00 24.32           O
ATOM    245  CB  THR A  31      -0.963 -14.709  -8.905  1.00 23.02           C
ATOM    246  CG2 THR A  31      -2.136 -15.206  -8.074  1.00 22.41           C
ATOM    247  OG1 THR A  31      -0.735 -15.638  -9.972  1.00 23.92           O
ATOM    248  N   GLN A  32      -2.982 -11.908  -8.530  1.00 23.36           N
ATOM    249  CA  GLN A  32      -3.658 -11.164  -7.495  1.00 23.78           C
ATOM    250  C   GLN A  32      -4.124 -12.154  -6.440  1.00 24.12           C
ATOM    251  O   GLN A  32      -4.709 -13.192  -6.774  1.00 24.15           O
ATOM    252  CB  GLN A  32      -4.859 -10.427  -8.083  1.00 24.73           C
ATOM    253  CG  GLN A  32      -5.498  -9.430  -7.129  1.00 26.90           C
ATOM    254  CD  GLN A  32      -6.922  -9.086  -7.504  1.00 25.16           C
ATOM    255  NE2 GLN A  32      -7.574  -8.294  -6.660  1.00 24.01           N
ATOM    256  OE1 GLN A  32      -7.439  -9.533  -8.533  1.00 21.43           O
ATOM    257  N   PHE A  33      -3.857 -11.840  -5.170  1.00 24.11           N
ATOM    258  CA  PHE A  33      -4.224 -12.741  -4.073  1.00 22.34           C
ATOM    259  C   PHE A  33      -5.207 -12.174  -3.042  1.00 23.06           C
ATOM    260  O   PHE A  33      -5.908 -12.930  -2.373  1.00 23.45           O
ATOM    261  CB  PHE A  33      -2.988 -13.366  -3.413  1.00 21.25           C
ATOM    262  CG  PHE A  33      -2.087 -12.387  -2.701  1.00 20.18           C
ATOM    263  CD1 PHE A  33      -1.174 -11.611  -3.405  1.00 15.54           C
ATOM    264  CD2 PHE A  33      -2.113 -12.290  -1.310  1.00 19.39           C
ATOM    265  CE1 PHE A  33      -0.333 -10.729  -2.739  1.00 14.07           C
ATOM    266  CE2 PHE A  33      -1.273 -11.414  -0.639  1.00 14.90           C
ATOM    267  CZ  PHE A  33      -0.384 -10.628  -1.358  1.00 16.68           C
ATOM    268  N   VAL A  34      -5.269 -10.849  -2.932  1.00 23.29           N
ATOM    269  CA  VAL A  34      -6.200 -10.188  -2.020  1.00 22.72           C
ATOM    270  C   VAL A  34      -6.761  -8.898  -2.607  1.00 22.74           C
ATOM    271  O   VAL A  34      -6.093  -8.203  -3.377  1.00 21.57           O
ATOM    272  CB  VAL A  34      -5.557  -9.850  -0.631  1.00 23.19           C
ATOM    273  CG1 VAL A  34      -5.356 -11.102   0.207  1.00 25.58           C
ATOM    274  CG2 VAL A  34      -4.251  -9.075  -0.785  1.00 20.78           C
ATOM    275  N   ARG A  35      -7.990  -8.582  -2.219  1.00 23.12           N
ATOM    276  CA  ARG A  35      -8.583  -7.299  -2.543  1.00 24.85           C
ATOM    277  C   ARG A  35      -9.362  -6.747  -1.358  1.00 25.23           C
ATOM    278  O   ARG A  35      -9.874  -7.502  -0.527  1.00 25.23           O
ATOM    279  CB  ARG A  35      -9.462  -7.396  -3.798  1.00 25.60           C
ATOM    280  CG  ARG A  35     -10.816  -8.040  -3.592  1.00 29.39           C
ATOM    281  CD  ARG A  35     -11.943  -7.042  -3.799  1.00 36.40           C
ATOM    282  NE  ARG A  35     -13.257  -7.638  -3.554  1.00 41.23           N
ATOM    283  CZ  ARG A  35     -13.867  -8.498  -4.369  1.00 44.53           C
ATOM    284  NH1 ARG A  35     -13.291  -8.890  -5.500  1.00 45.54           N1+
ATOM    285  NH2 ARG A  35     -15.061  -8.977  -4.044  1.00 46.99           N
ATOM    286  N   PHE A  36      -9.436  -5.422  -1.293  1.00 25.29           N
ATOM    287  CA  PHE A  36     -10.196  -4.733  -0.271  1.00 25.01           C
ATOM    288  C   PHE A  36     -10.808  -3.471  -0.858  1.00 26.06           C
ATOM    289  O   PHE A  36     -10.100  -2.593  -1.354  1.00 26.21           O
ATOM    290  CB  PHE A  36      -9.302  -4.398   0.925  1.00 24.25           C
ATOM    291  CG  PHE A  36     -10.006  -3.664   2.026  1.00 21.78           C
ATOM    292  CD1 PHE A  36     -10.148  -2.279   1.982  1.00 20.11           C
ATOM    293  CD2 PHE A  36     -10.517  -4.355   3.119  1.00 20.93           C
ATOM    294  CE1 PHE A  36     -10.794  -1.594   3.007  1.00 20.44           C
ATOM    295  CE2 PHE A  36     -11.161  -3.677   4.159  1.00 19.81           C
ATOM    296  CZ  PHE A  36     -11.300  -2.295   4.102  1.00 20.71           C
ATOM    297  N   ASP A  37     -12.131  -3.400  -0.794  1.00 27.48           N
ATOM    298  CA  ASP A  37     -12.878  -2.252  -1.277  1.00 29.20           C
ATOM    299  C   ASP A  37     -13.671  -1.686  -0.122  1.00 30.11           C
ATOM    300  O   ASP A  37     -14.271  -2.434   0.650  1.00 30.07           O
ATOM    301  CB  ASP A  37     -13.819  -2.672  -2.413  1.00 29.77           C
ATOM    302  CG  ASP A  37     -14.425  -1.483  -3.159  1.00 30.46           C
ATOM    303  OD1 ASP A  37     -14.576  -0.392  -2.569  1.00 28.93           O
ATOM    304  OD2 ASP A  37     -14.766  -1.654  -4.347  1.00 32.45           O1-
ATOM    305  N   SER A  38     -13.677  -0.361  -0.011  1.00 32.15           N
ATOM    306  CA  SER A  38     -14.410   0.325   1.050  1.00 34.35           C
ATOM    307  C   SER A  38     -15.917   0.309   0.802  1.00 36.53           C
ATOM    308  O   SER A  38     -16.700   0.470   1.736  1.00 37.67           O
ATOM    309  CB  SER A  38     -13.916   1.763   1.203  1.00 34.01           C
ATOM    310  OG  SER A  38     -13.918   2.439  -0.042  1.00 33.12           O
ATOM    311  N   ASP A  39     -16.310   0.098  -0.454  1.00 38.61           N
ATOM    312  CA  ASP A  39     -17.720   0.122  -0.853  1.00 40.94           C
ATOM    313  C   ASP A  39     -18.422  -1.237  -0.838  1.00 42.38           C
ATOM    314  O   ASP A  39     -19.648  -1.302  -0.966  1.00 42.48           O
ATOM    315  CB  ASP A  39     -17.876   0.768  -2.233  1.00 40.92           C
ATOM    316  CG  ASP A  39     -17.763   2.279  -2.188  1.00 41.64           C
ATOM    317  OD1 ASP A  39     -17.426   2.829  -1.115  1.00 40.97           O
ATOM    318  OD2 ASP A  39     -18.015   2.917  -3.232  1.00 41.90           O1-
ATOM    319  N   ALA A  40     -17.651  -2.312  -0.679  1.00 43.94           N
ATOM    320  CA  ALA A  40     -18.210  -3.663  -0.594  1.00 45.30           C
ATOM    321  C   ALA A  40     -19.208  -3.788   0.564  1.00 46.19           C
ATOM    322  O   ALA A  40     -19.226  -2.950   1.470  1.00 46.36           O
ATOM    323  CB  ALA A  40     -17.091  -4.696  -0.464  1.00 45.09           C
ATOM    324  N   ALA A  41     -20.037  -4.832   0.525  1.00 47.41           N
ATOM    325  CA  ALA A  41     -21.098  -5.045   1.519  1.00 48.12           C
ATOM    326  C   ALA A  41     -20.618  -4.885   2.969  1.00 48.57           C
ATOM    327  O   ALA A  41     -21.145  -4.055   3.713  1.00 49.29           O
ATOM    328  CB  ALA A  41     -21.760  -6.402   1.310  1.00 47.88           C
ATOM    329  N   SER A  42     -19.622  -5.678   3.359  1.00 48.68           N
ATOM    330  CA  SER A  42     -19.008  -5.565   4.683  1.00 48.33           C
ATOM    331  C   SER A  42     -17.488  -5.502   4.542  1.00 47.65           C
ATOM    332  O   SER A  42     -16.839  -6.541   4.392  1.00 48.31           O
ATOM    333  CB  SER A  42     -19.426  -6.732   5.585  1.00 48.96           C
ATOM    334  OG  SER A  42     -19.012  -7.976   5.044  1.00 49.22           O
ATOM    335  N   PRO A  43     -16.917  -4.277   4.580  1.00 46.60           N
ATOM    336  CA  PRO A  43     -15.489  -4.040   4.330  1.00 45.10           C
ATOM    337  C   PRO A  43     -14.562  -4.992   5.085  1.00 43.23           C
ATOM    338  O   PRO A  43     -14.320  -4.819   6.280  1.00 43.44           O
ATOM    339  CB  PRO A  43     -15.288  -2.594   4.801  1.00 45.11           C
ATOM    340  CG  PRO A  43     -16.606  -1.956   4.566  1.00 46.04           C
ATOM    341  CD  PRO A  43     -17.634  -3.018   4.869  1.00 46.35           C
ATOM    342  N   ARG A  44     -14.073  -6.003   4.374  1.00 41.38           N
ATOM    343  CA  ARG A  44     -13.104  -6.956   4.907  1.00 39.89           C
ATOM    344  C   ARG A  44     -12.233  -7.438   3.756  1.00 37.24           C
ATOM    345  O   ARG A  44     -12.667  -7.438   2.605  1.00 36.75           O
ATOM    346  CB  ARG A  44     -13.810  -8.137   5.595  1.00 40.73           C
ATOM    347  CG  ARG A  44     -14.209  -9.265   4.650  1.00 44.88           C
ATOM    348  CD  ARG A  44     -15.582  -9.819   4.958  1.00 48.54           C
ATOM    349  NE  ARG A  44     -16.229 -10.312   3.742  1.00 50.09           N
ATOM    350  CZ  ARG A  44     -17.511 -10.659   3.650  1.00 52.34           C
ATOM    351  NH1 ARG A  44     -18.313 -10.577   4.708  1.00 51.68           N1+
ATOM    352  NH2 ARG A  44     -17.995 -11.093   2.493  1.00 52.62           N
ATOM    353  N   THR A  45     -11.004  -7.834   4.066  1.00 35.50           N
ATOM    354  CA  THR A  45     -10.082  -8.361   3.060  1.00 34.01           C
ATOM    355  C   THR A  45     -10.589  -9.702   2.522  1.00 33.32           C
ATOM    356  O   THR A  45     -10.948 -10.585   3.298  1.00 32.88           O
ATOM    357  CB  THR A  45      -8.672  -8.525   3.643  1.00 33.42           C
ATOM    358  CG2 THR A  45      -7.657  -8.854   2.547  1.00 34.08           C
ATOM    359  OG1 THR A  45      -8.292  -7.312   4.299  1.00 28.89           O
ATOM    360  N   GLU A  46     -10.623  -9.833   1.196  1.00 33.10           N
ATOM    361  CA  GLU A  46     -11.125 -11.044   0.530  1.00 33.78           C
ATOM    362  C   GLU A  46     -10.007 -11.772  -0.228  1.00 33.20           C
ATOM    363  O   GLU A  46      -9.192 -11.126  -0.889  1.00 33.08           O
ATOM    364  CB  GLU A  46     -12.267 -10.701  -0.438  1.00 33.18           C
ATOM    365  CG  GLU A  46     -13.563 -10.214   0.221  1.00 35.33           C
ATOM    366  CD  GLU A  46     -14.585  -9.674  -0.788  1.00 34.96           C
ATOM    367  OE1 GLU A  46     -14.968  -8.490  -0.667  1.00 36.66           O
ATOM    368  OE2 GLU A  46     -14.999 -10.425  -1.702  1.00 34.72           O1-
ATOM    369  N   PRO A  47      -9.970 -13.119  -0.133  1.00 32.40           N
ATOM    370  CA  PRO A  47      -9.016 -13.948  -0.879  1.00 31.62           C
ATOM    371  C   PRO A  47      -9.287 -13.937  -2.384  1.00 31.73           C
ATOM    372  O   PRO A  47     -10.447 -13.945  -2.802  1.00 32.43           O
ATOM    373  CB  PRO A  47      -9.258 -15.355  -0.314  1.00 31.70           C
ATOM    374  CG  PRO A  47     -10.652 -15.326   0.213  1.00 30.52           C
ATOM    375  CD  PRO A  47     -10.854 -13.935   0.724  1.00 32.69           C
ATOM    376  N   ARG A  48      -8.223 -13.913  -3.186  1.00 31.14           N
ATOM    377  CA  ARG A  48      -8.339 -13.919  -4.651  1.00 30.45           C
ATOM    378  C   ARG A  48      -7.426 -14.978  -5.283  1.00 30.13           C
ATOM    379  O   ARG A  48      -7.363 -15.111  -6.504  1.00 29.53           O
ATOM    380  CB  ARG A  48      -8.023 -12.534  -5.235  1.00 30.40           C
ATOM    381  CG  ARG A  48      -8.802 -11.364  -4.633  1.00 32.75           C
ATOM    382  CD  ARG A  48     -10.278 -11.417  -4.979  1.00 38.59           C
ATOM    383  NE  ARG A  48     -10.522 -11.116  -6.387  1.00 42.93           N
ATOM    384  CZ  ARG A  48     -11.582 -11.531  -7.075  1.00 43.48           C
ATOM    385  NH1 ARG A  48     -12.509 -12.285  -6.494  1.00 45.15           N1+
ATOM    386  NH2 ARG A  48     -11.710 -11.201  -8.355  1.00 43.77           N
ATOM    387  N   ALA A  49      -6.712 -15.718  -4.437  1.00 30.36           N
ATOM    388  CA  ALA A  49      -5.899 -16.847  -4.876  1.00 30.50           C
ATOM    389  C   ALA A  49      -6.215 -18.046  -3.998  1.00 30.55           C
ATOM    390  O   ALA A  49      -6.554 -17.875  -2.828  1.00 31.25           O
ATOM    391  CB  ALA A  49      -4.422 -16.506  -4.805  1.00 31.01           C
ATOM    392  N   PRO A  50      -6.107 -19.268  -4.556  1.00 30.56           N
ATOM    393  CA  PRO A  50      -6.477 -20.456  -3.777  1.00 29.51           C
ATOM    394  C   PRO A  50      -5.576 -20.701  -2.566  1.00 28.52           C
ATOM    395  O   PRO A  50      -6.076 -21.059  -1.502  1.00 29.77           O
ATOM    396  CB  PRO A  50      -6.346 -21.597  -4.787  1.00 28.56           C
ATOM    397  CG  PRO A  50      -5.392 -21.090  -5.812  1.00 30.81           C
ATOM    398  CD  PRO A  50      -5.646 -19.620  -5.911  1.00 29.70           C
ATOM    399  N   TRP A  51      -4.270 -20.489  -2.730  1.00 27.63           N
ATOM    400  CA  TRP A  51      -3.274 -20.818  -1.700  1.00 27.07           C
ATOM    401  C   TRP A  51      -3.353 -19.973  -0.427  1.00 27.34           C
ATOM    402  O   TRP A  51      -2.906 -20.409   0.630  1.00 27.39           O
ATOM    403  CB  TRP A  51      -1.852 -20.801  -2.277  1.00 26.15           C
ATOM    404  CG  TRP A  51      -1.614 -19.700  -3.260  1.00 25.42           C
ATOM    405  CD1 TRP A  51      -1.760 -19.770  -4.617  1.00 23.63           C
ATOM    406  CD2 TRP A  51      -1.203 -18.364  -2.966  1.00 23.36           C
ATOM    407  CE2 TRP A  51      -1.113 -17.676  -4.196  1.00 26.05           C
ATOM    408  CE3 TRP A  51      -0.897 -17.680  -1.783  1.00 23.15           C
ATOM    409  NE1 TRP A  51      -1.454 -18.559  -5.188  1.00 24.91           N
ATOM    410  CZ2 TRP A  51      -0.734 -16.329  -4.278  1.00 25.14           C
ATOM    411  CZ3 TRP A  51      -0.514 -16.343  -1.862  1.00 23.97           C
ATOM    412  CH2 TRP A  51      -0.434 -15.685  -3.104  1.00 25.60           C
ATOM    413  N   ILE A  52      -3.920 -18.772  -0.530  1.00 28.16           N
ATOM    414  CA  ILE A  52      -4.046 -17.882   0.623  1.00 28.68           C
ATOM    415  C   ILE A  52      -5.206 -18.273   1.532  1.00 30.22           C
ATOM    416  O   ILE A  52      -5.299 -17.793   2.660  1.00 31.86           O
ATOM    417  CB  ILE A  52      -4.185 -16.389   0.200  1.00 28.54           C
ATOM    418  CG1 ILE A  52      -3.853 -15.455   1.365  1.00 24.78           C
ATOM    419  CG2 ILE A  52      -5.578 -16.091  -0.346  1.00 26.92           C
ATOM    420  CD1 ILE A  52      -2.386 -15.177   1.525  1.00 20.24           C
ATOM    421  N   GLU A  53      -6.091 -19.140   1.049  1.00 31.46           N
ATOM    422  CA  GLU A  53      -7.253 -19.559   1.841  1.00 33.93           C
ATOM    423  C   GLU A  53      -6.880 -20.421   3.054  1.00 35.24           C
ATOM    424  O   GLU A  53      -7.725 -20.699   3.911  1.00 36.65           O
ATOM    425  CB  GLU A  53      -8.299 -20.251   0.956  1.00 34.09           C
ATOM    426  CG  GLU A  53      -9.349 -19.288   0.412  1.00 37.24           C
ATOM    427  CD  GLU A  53      -9.973 -19.742  -0.895  1.00 43.66           C
ATOM    428  OE1 GLU A  53     -10.377 -20.922  -1.006  1.00 45.61           O
ATOM    429  OE2 GLU A  53     -10.069 -18.900  -1.814  1.00 47.77           O1-
ATOM    430  N   GLN A  54      -5.612 -20.824   3.127  1.00 35.31           N
ATOM    431  CA  GLN A  54      -5.103 -21.586   4.261  1.00 35.46           C
ATOM    432  C   GLN A  54      -4.742 -20.689   5.454  1.00 35.02           C
ATOM    433  O   GLN A  54      -4.514 -21.177   6.563  1.00 33.78           O
ATOM    434  CB  GLN A  54      -3.905 -22.443   3.835  1.00 35.58           C
ATOM    435  CG  GLN A  54      -2.601 -21.682   3.637  1.00 35.64           C
ATOM    436  CD  GLN A  54      -1.442 -22.604   3.297  1.00 36.66           C
ATOM    437  NE2 GLN A  54      -1.411 -23.096   2.059  1.00 36.85           N
ATOM    438  OE1 GLN A  54      -0.591 -22.883   4.144  1.00 37.10           O
ATOM    439  N   GLU A  55      -4.690 -19.381   5.216  1.00 35.21           N
ATOM    440  CA  GLU A  55      -4.366 -18.419   6.264  1.00 35.76           C
ATOM    441  C   GLU A  55      -5.537 -18.284   7.226  1.00 36.56           C
ATOM    442  O   GLU A  55      -6.688 -18.187   6.799  1.00 37.13           O
ATOM    443  CB  GLU A  55      -4.009 -17.060   5.661  1.00 34.93           C
ATOM    444  CG  GLU A  55      -2.768 -17.066   4.771  1.00 34.75           C
ATOM    445  CD  GLU A  55      -1.460 -17.232   5.539  1.00 34.11           C
ATOM    446  OE1 GLU A  55      -1.428 -16.971   6.764  1.00 34.74           O
ATOM    447  OE2 GLU A  55      -0.458 -17.615   4.902  1.00 25.99           O1-
ATOM    448  N   GLY A  56      -5.233 -18.274   8.522  1.00 36.36           N
ATOM    449  CA  GLY A  56      -6.257 -18.322   9.560  1.00 36.79           C
ATOM    450  C   GLY A  56      -7.013 -17.026   9.781  1.00 37.71           C
ATOM    451  O   GLY A  56      -6.666 -15.992   9.202  1.00 38.23           O
ATOM    452  N   PRO A  57      -8.058 -17.070  10.631  1.00 37.46           N
ATOM    453  CA  PRO A  57      -8.879 -15.904  10.994  1.00 37.24           C
ATOM    454  C   PRO A  57      -8.052 -14.723  11.497  1.00 36.77           C
ATOM    455  O   PRO A  57      -8.428 -13.568  11.284  1.00 36.60           O
ATOM    456  CB  PRO A  57      -9.763 -16.441  12.123  1.00 37.45           C
ATOM    457  CG  PRO A  57      -9.853 -17.911  11.860  1.00 38.61           C
ATOM    458  CD  PRO A  57      -8.518 -18.299  11.305  1.00 37.41           C
ATOM    459  N   GLU A  58      -6.937 -15.026  12.157  1.00 36.84           N
ATOM    460  CA  GLU A  58      -6.005 -14.025  12.654  1.00 36.06           C
ATOM    461  C   GLU A  58      -5.432 -13.211  11.497  1.00 35.01           C
ATOM    462  O   GLU A  58      -5.188 -12.013  11.635  1.00 35.35           O
ATOM    463  CB  GLU A  58      -4.878 -14.707  13.433  1.00 36.77           C
ATOM    464  CG  GLU A  58      -3.953 -13.764  14.202  1.00 40.63           C
ATOM    465  CD  GLU A  58      -4.508 -13.320  15.559  1.00 45.32           C
ATOM    466  OE1 GLU A  58      -4.046 -12.276  16.064  1.00 47.53           O
ATOM    467  OE2 GLU A  58      -5.391 -14.007  16.125  1.00 43.04           O1-
ATOM    468  N   TYR A  59      -5.232 -13.868  10.356  1.00 33.26           N
ATOM    469  CA  TYR A  59      -4.700 -13.204   9.169  1.00 32.05           C
ATOM    470  C   TYR A  59      -5.699 -12.225   8.559  1.00 31.19           C
ATOM    471  O   TYR A  59      -5.325 -11.115   8.171  1.00 31.36           O
ATOM    472  CB  TYR A  59      -4.244 -14.227   8.118  1.00 31.73           C
ATOM    473  CG  TYR A  59      -3.733 -13.608   6.832  1.00 30.99           C
ATOM    474  CD1 TYR A  59      -2.415 -13.171   6.725  1.00 31.67           C
ATOM    475  CD2 TYR A  59      -4.568 -13.456   5.723  1.00 29.97           C
ATOM    476  CE1 TYR A  59      -1.938 -12.596   5.554  1.00 33.08           C
ATOM    477  CE2 TYR A  59      -4.098 -12.884   4.541  1.00 31.06           C
ATOM    478  CZ  TYR A  59      -2.785 -12.455   4.468  1.00 30.48           C
ATOM    479  OH  TYR A  59      -2.308 -11.890   3.313  1.00 28.73           O
ATOM    480  N   TRP A  60      -6.960 -12.639   8.468  1.00 30.23           N
ATOM    481  CA  TRP A  60      -7.977 -11.820   7.822  1.00 29.76           C
ATOM    482  C   TRP A  60      -8.420 -10.661   8.709  1.00 30.35           C
ATOM    483  O   TRP A  60      -8.681  -9.572   8.205  1.00 31.42           O
ATOM    484  CB  TRP A  60      -9.171 -12.659   7.364  1.00 29.00           C
ATOM    485  CG  TRP A  60      -8.800 -13.808   6.463  1.00 29.04           C
ATOM    486  CD1 TRP A  60      -8.926 -15.138   6.744  1.00 30.77           C
ATOM    487  CD2 TRP A  60      -8.231 -13.730   5.145  1.00 28.74           C
ATOM    488  CE2 TRP A  60      -8.047 -15.056   4.692  1.00 27.20           C
ATOM    489  CE3 TRP A  60      -7.865 -12.670   4.303  1.00 29.73           C
ATOM    490  NE1 TRP A  60      -8.478 -15.892   5.686  1.00 29.46           N
ATOM    491  CZ2 TRP A  60      -7.514 -15.352   3.436  1.00 26.80           C
ATOM    492  CZ3 TRP A  60      -7.336 -12.965   3.051  1.00 30.18           C
ATOM    493  CH2 TRP A  60      -7.163 -14.297   2.632  1.00 30.43           C
ATOM    494  N   ASP A  61      -8.485 -10.898  10.020  1.00 29.81           N
ATOM    495  CA  ASP A  61      -8.789  -9.845  10.995  1.00 29.42           C
ATOM    496  C   ASP A  61      -7.764  -8.711  10.945  1.00 28.93           C
ATOM    497  O   ASP A  61      -8.131  -7.529  10.976  1.00 29.67           O
ATOM    498  CB  ASP A  61      -8.854 -10.420  12.418  1.00 29.77           C
ATOM    499  CG  ASP A  61     -10.074 -11.317  12.647  1.00 32.03           C
ATOM    500  OD1 ASP A  61     -10.787 -11.661  11.675  1.00 34.38           O
ATOM    501  OD2 ASP A  61     -10.312 -11.688  13.816  1.00 30.58           O1-
ATOM    502  N   ARG A  62      -6.485  -9.080  10.862  1.00 27.18           N
ATOM    503  CA  ARG A  62      -5.385  -8.122  10.873  1.00 25.80           C
ATOM    504  C   ARG A  62      -5.339  -7.316   9.583  1.00 26.26           C
ATOM    505  O   ARG A  62      -5.203  -6.091   9.618  1.00 25.85           O
ATOM    506  CB  ARG A  62      -4.051  -8.836  11.109  1.00 24.72           C
ATOM    507  CG  ARG A  62      -2.833  -7.929  11.030  1.00 23.03           C
ATOM    508  CD  ARG A  62      -1.522  -8.691  11.194  1.00 25.43           C
ATOM    509  NE  ARG A  62      -1.244  -9.568  10.056  1.00 25.90           N
ATOM    510  CZ  ARG A  62      -1.481 -10.877  10.037  1.00 23.38           C
ATOM    511  NH1 ARG A  62      -2.006 -11.487  11.094  1.00 24.21           N1+
ATOM    512  NH2 ARG A  62      -1.196 -11.576   8.956  1.00 22.27           N
ATOM    513  N   ASN A  63      -5.449  -8.010   8.451  1.00 26.68           N
ATOM    514  CA  ASN A  63      -5.538  -7.351   7.152  1.00 26.21           C
ATOM    515  C   ASN A  63      -6.711  -6.381   7.119  1.00 25.23           C
ATOM    516  O   ASN A  63      -6.537  -5.216   6.768  1.00 24.83           O
ATOM    517  CB  ASN A  63      -5.651  -8.377   6.017  1.00 27.16           C
ATOM    518  CG  ASN A  63      -4.305  -8.970   5.617  1.00 27.94           C
ATOM    519  ND2 ASN A  63      -3.572  -9.507   6.591  1.00 30.71           N
ATOM    520  OD1 ASN A  63      -3.938  -8.956   4.441  1.00 27.00           O
ATOM    521  N   THR A  64      -7.890  -6.865   7.515  1.00 25.18           N
ATOM    522  CA  THR A  64      -9.100  -6.043   7.635  1.00 25.84           C
ATOM    523  C   THR A  64      -8.831  -4.750   8.409  1.00 26.80           C
ATOM    524  O   THR A  64      -9.123  -3.659   7.918  1.00 26.92           O
ATOM    525  CB  THR A  64     -10.249  -6.826   8.320  1.00 25.68           C
ATOM    526  CG2 THR A  64     -11.477  -5.945   8.508  1.00 26.74           C
ATOM    527  OG1 THR A  64     -10.608  -7.961   7.518  1.00 26.60           O
ATOM    528  N   GLN A  65      -8.267  -4.890   9.609  1.00 27.53           N
ATOM    529  CA  GLN A  65      -7.930  -3.756  10.466  1.00 28.32           C
ATOM    530  C   GLN A  65      -7.006  -2.774   9.754  1.00 28.41           C
ATOM    531  O   GLN A  65      -7.296  -1.570   9.701  1.00 27.77           O
ATOM    532  CB  GLN A  65      -7.275  -4.241  11.764  1.00 28.13           C
ATOM    533  CG  GLN A  65      -6.677  -3.127  12.625  1.00 29.52           C
ATOM    534  CD  GLN A  65      -5.747  -3.633  13.722  1.00 30.94           C
ATOM    535  NE2 GLN A  65      -6.048  -4.808  14.272  1.00 33.33           N
ATOM    536  OE1 GLN A  65      -4.770  -2.964  14.072  1.00 35.99           O
ATOM    537  N   ILE A  66      -5.900  -3.299   9.217  1.00 27.76           N
ATOM    538  CA  ILE A  66      -4.893  -2.486   8.543  1.00 26.88           C
ATOM    539  C   ILE A  66      -5.508  -1.682   7.406  1.00 26.67           C
ATOM    540  O   ILE A  66      -5.193  -0.504   7.252  1.00 28.33           O
ATOM    541  CB  ILE A  66      -3.686  -3.328   8.045  1.00 27.51           C
ATOM    542  CG1 ILE A  66      -2.546  -3.300   9.064  1.00 28.16           C
ATOM    543  CG2 ILE A  66      -3.127  -2.779   6.745  1.00 26.37           C
ATOM    544  CD1 ILE A  66      -2.704  -4.245  10.221  1.00 32.88           C
ATOM    545  N   PHE A  67      -6.401  -2.303   6.635  1.00 25.51           N
ATOM    546  CA  PHE A  67      -7.041  -1.609   5.513  1.00 25.45           C
ATOM    547  C   PHE A  67      -8.121  -0.634   5.957  1.00 26.73           C
ATOM    548  O   PHE A  67      -8.317   0.405   5.319  1.00 28.03           O
ATOM    549  CB  PHE A  67      -7.586  -2.580   4.464  1.00 23.66           C
ATOM    550  CG  PHE A  67      -6.517  -3.330   3.729  1.00 22.36           C
ATOM    551  CD1 PHE A  67      -5.558  -2.651   2.980  1.00 16.88           C
ATOM    552  CD2 PHE A  67      -6.459  -4.720   3.793  1.00 21.61           C
ATOM    553  CE1 PHE A  67      -4.559  -3.340   2.312  1.00 16.73           C
ATOM    554  CE2 PHE A  67      -5.463  -5.423   3.121  1.00 22.36           C
ATOM    555  CZ  PHE A  67      -4.510  -4.730   2.377  1.00 21.04           C
ATOM    556  N   LYS A  68      -8.810  -0.963   7.049  1.00 26.59           N
ATOM    557  CA  LYS A  68      -9.823  -0.074   7.600  1.00 27.00           C
ATOM    558  C   LYS A  68      -9.189   1.204   8.142  1.00 27.19           C
ATOM    559  O   LYS A  68      -9.775   2.283   8.068  1.00 27.52           O
ATOM    560  CB  LYS A  68     -10.620  -0.778   8.699  1.00 27.82           C
ATOM    561  CG  LYS A  68     -11.727  -1.679   8.186  1.00 28.28           C
ATOM    562  CD  LYS A  68     -12.568  -2.224   9.325  1.00 32.89           C
ATOM    563  CE  LYS A  68     -13.506  -1.162   9.902  1.00 36.38           C
ATOM    564  NZ  LYS A  68     -14.299  -1.692  11.055  1.00 39.35           N1+
ATOM    565  N   THR A  69      -7.980   1.069   8.671  1.00 27.88           N
ATOM    566  CA  THR A  69      -7.252   2.184   9.255  1.00 27.74           C
ATOM    567  C   THR A  69      -6.722   3.074   8.130  1.00 28.32           C
ATOM    568  O   THR A  69      -6.640   4.300   8.269  1.00 28.57           O
ATOM    569  CB  THR A  69      -6.126   1.665  10.198  1.00 27.28           C
ATOM    570  CG2 THR A  69      -4.728   1.922   9.647  1.00 28.09           C
ATOM    571  OG1 THR A  69      -6.264   2.259  11.491  1.00 28.11           O
ATOM    572  N   ASN A  70      -6.395   2.447   7.003  1.00 28.02           N
ATOM    573  CA  ASN A  70      -5.958   3.175   5.817  1.00 28.13           C
ATOM    574  C   ASN A  70      -7.100   3.897   5.107  1.00 28.06           C
ATOM    575  O   ASN A  70      -6.887   4.929   4.468  1.00 27.23           O
ATOM    576  CB  ASN A  70      -5.216   2.240   4.862  1.00 28.15           C
ATOM    577  CG  ASN A  70      -3.921   1.715   5.456  1.00 25.76           C
ATOM    578  ND2 ASN A  70      -3.334   0.716   4.813  1.00 24.23           N
ATOM    579  OD1 ASN A  70      -3.457   2.203   6.487  1.00 22.81           O
ATOM    580  N   THR A  71      -8.309   3.351   5.238  1.00 29.00           N
ATOM    581  CA  THR A  71      -9.519   3.981   4.711  1.00 29.80           C
ATOM    582  C   THR A  71      -9.631   5.415   5.217  1.00 30.46           C
ATOM    583  O   THR A  71      -9.777   6.342   4.424  1.00 31.78           O
ATOM    584  CB  THR A  71     -10.785   3.199   5.101  1.00 29.73           C
ATOM    585  CG2 THR A  71     -11.995   3.696   4.320  1.00 31.54           C
ATOM    586  OG1 THR A  71     -10.591   1.809   4.826  1.00 31.65           O
ATOM    587  N   GLN A  72      -9.534   5.590   6.532  1.00 30.77           N
ATOM    588  CA  GLN A  72      -9.626   6.909   7.150  1.00 30.90           C
ATOM    589  C   GLN A  72      -8.427   7.798   6.832  1.00 30.81           C
ATOM    590  O   GLN A  72      -8.572   9.011   6.667  1.00 31.02           O
ATOM    591  CB  GLN A  72      -9.803   6.774   8.660  1.00 31.54           C
ATOM    592  CG  GLN A  72     -11.260   6.803   9.122  1.00 35.00           C
ATOM    593  CD  GLN A  72     -12.061   5.589   8.674  1.00 36.71           C
ATOM    594  NE2 GLN A  72     -11.758   4.429   9.250  1.00 36.82           N
ATOM    595  OE1 GLN A  72     -12.939   5.697   7.815  1.00 37.78           O
ATOM    596  N   THR A  73      -7.251   7.187   6.730  1.00 30.95           N
ATOM    597  CA  THR A  73      -5.997   7.921   6.551  1.00 30.35           C
ATOM    598  C   THR A  73      -5.827   8.476   5.140  1.00 30.58           C
ATOM    599  O   THR A  73      -5.361   9.603   4.974  1.00 29.55           O
ATOM    600  CB  THR A  73      -4.799   7.044   6.935  1.00 29.93           C
ATOM    601  CG2 THR A  73      -3.502   7.858   6.997  1.00 30.64           C
ATOM    602  OG1 THR A  73      -5.056   6.467   8.215  1.00 27.80           O
ATOM    603  N   TYR A  74      -6.200   7.686   4.134  1.00 31.91           N
ATOM    604  CA  TYR A  74      -6.169   8.154   2.747  1.00 32.63           C
ATOM    605  C   TYR A  74      -7.291   9.150   2.461  1.00 33.32           C
ATOM    606  O   TYR A  74      -7.136  10.029   1.619  1.00 33.33           O
ATOM    607  CB  TYR A  74      -6.167   6.983   1.757  1.00 32.69           C
ATOM    608  CG  TYR A  74      -4.813   6.313   1.668  1.00 32.93           C
ATOM    609  CD1 TYR A  74      -4.479   5.259   2.519  1.00 33.46           C
ATOM    610  CD2 TYR A  74      -3.855   6.750   0.754  1.00 33.04           C
ATOM    611  CE1 TYR A  74      -3.235   4.651   2.459  1.00 33.33           C
ATOM    612  CE2 TYR A  74      -2.603   6.148   0.684  1.00 33.48           C
ATOM    613  CZ  TYR A  74      -2.301   5.100   1.542  1.00 34.33           C
ATOM    614  OH  TYR A  74      -1.067   4.498   1.486  1.00 35.05           O
ATOM    615  N   ARG A  75      -8.408   9.015   3.176  1.00 34.45           N
ATOM    616  CA  ARG A  75      -9.462  10.029   3.170  1.00 36.56           C
ATOM    617  C   ARG A  75      -8.911  11.350   3.704  1.00 36.25           C
ATOM    618  O   ARG A  75      -9.117  12.405   3.101  1.00 36.64           O
ATOM    619  CB  ARG A  75     -10.662   9.578   4.012  1.00 36.21           C
ATOM    620  CG  ARG A  75     -11.604   8.624   3.295  1.00 40.24           C
ATOM    621  CD  ARG A  75     -12.692   8.100   4.227  1.00 39.77           C
ATOM    622  NE  ARG A  75     -13.602   7.171   3.549  1.00 46.46           N
ATOM    623  CZ  ARG A  75     -14.641   6.569   4.130  1.00 50.07           C
ATOM    624  NH1 ARG A  75     -14.919   6.788   5.412  1.00 48.77           N1+
ATOM    625  NH2 ARG A  75     -15.407   5.743   3.425  1.00 50.66           N
ATOM    626  N   GLU A  76      -8.204  11.270   4.832  1.00 35.90           N
ATOM    627  CA  GLU A  76      -7.510  12.403   5.433  1.00 35.46           C
ATOM    628  C   GLU A  76      -6.491  12.985   4.468  1.00 34.11           C
ATOM    629  O   GLU A  76      -6.401  14.203   4.312  1.00 34.25           O
ATOM    630  CB  GLU A  76      -6.813  11.963   6.723  1.00 35.53           C
ATOM    631  CG  GLU A  76      -5.892  13.004   7.351  1.00 36.73           C
ATOM    632  CD  GLU A  76      -5.486  12.643   8.772  1.00 37.73           C
ATOM    633  OE1 GLU A  76      -6.368  12.219   9.555  1.00 35.41           O
ATOM    634  OE2 GLU A  76      -4.285  12.792   9.106  1.00 41.90           O1-
ATOM    635  N   SER A  77      -5.735  12.102   3.821  1.00 33.47           N
ATOM    636  CA  SER A  77      -4.726  12.498   2.849  1.00 32.59           C
ATOM    637  C   SER A  77      -5.348  13.314   1.724  1.00 33.52           C
ATOM    638  O   SER A  77      -4.796  14.336   1.315  1.00 33.36           O
ATOM    639  CB  SER A  77      -4.028  11.265   2.283  1.00 32.04           C
ATOM    640  OG  SER A  77      -3.360  10.545   3.297  1.00 28.12           O
ATOM    641  N   LEU A  78      -6.512  12.862   1.252  1.00 34.61           N
ATOM    642  CA  LEU A  78      -7.257  13.525   0.177  1.00 35.29           C
ATOM    643  C   LEU A  78      -7.729  14.931   0.529  1.00 35.89           C
ATOM    644  O   LEU A  78      -7.886  15.771  -0.358  1.00 36.30           O
ATOM    645  CB  LEU A  78      -8.441  12.667  -0.270  1.00 35.14           C
ATOM    646  CG  LEU A  78      -8.110  11.462  -1.158  1.00 36.39           C
ATOM    647  CD1 LEU A  78      -9.257  10.470  -1.142  1.00 36.59           C
ATOM    648  CD2 LEU A  78      -7.773  11.872  -2.585  1.00 34.22           C
ATOM    649  N   ARG A  79      -7.960  15.188   1.813  1.00 36.63           N
ATOM    650  CA  ARG A  79      -8.230  16.546   2.271  1.00 37.25           C
ATOM    651  C   ARG A  79      -6.937  17.354   2.314  1.00 36.85           C
ATOM    652  O   ARG A  79      -6.906  18.506   1.877  1.00 36.98           O
ATOM    653  CB  ARG A  79      -8.910  16.560   3.641  1.00 37.82           C
ATOM    654  CG  ARG A  79     -10.420  16.420   3.592  1.00 41.57           C
ATOM    655  CD  ARG A  79     -11.101  17.154   4.761  1.00 46.94           C
ATOM    656  NE  ARG A  79     -10.616  16.714   6.071  1.00 50.72           N
ATOM    657  CZ  ARG A  79     -11.034  15.620   6.713  1.00 52.05           C
ATOM    658  NH1 ARG A  79     -11.956  14.827   6.176  1.00 50.16           N1+
ATOM    659  NH2 ARG A  79     -10.520  15.313   7.898  1.00 53.62           N
ATOM    660  N   ASN A  80      -5.873  16.746   2.834  1.00 36.13           N
ATOM    661  CA  ASN A  80      -4.584  17.426   2.941  1.00 36.65           C
ATOM    662  C   ASN A  80      -4.081  17.888   1.581  1.00 36.39           C
ATOM    663  O   ASN A  80      -3.612  19.014   1.441  1.00 36.30           O
ATOM    664  CB  ASN A  80      -3.533  16.538   3.617  1.00 36.85           C
ATOM    665  CG  ASN A  80      -3.968  16.037   4.988  1.00 37.40           C
ATOM    666  ND2 ASN A  80      -4.889  16.751   5.628  1.00 37.43           N
ATOM    667  OD1 ASN A  80      -3.478  15.012   5.462  1.00 36.36           O
ATOM    668  N   LEU A  81      -4.196  17.015   0.583  1.00 36.84           N
ATOM    669  CA  LEU A  81      -3.796  17.355  -0.777  1.00 37.37           C
ATOM    670  C   LEU A  81      -4.672  18.450  -1.366  1.00 36.80           C
ATOM    671  O   LEU A  81      -4.155  19.433  -1.908  1.00 35.94           O
ATOM    672  CB  LEU A  81      -3.771  16.116  -1.679  1.00 37.75           C
ATOM    673  CG  LEU A  81      -2.457  15.328  -1.611  1.00 40.54           C
ATOM    674  CD1 LEU A  81      -2.558  14.020  -2.369  1.00 44.76           C
ATOM    675  CD2 LEU A  81      -1.321  16.150  -2.166  1.00 40.05           C
ATOM    676  N   ARG A  82      -5.989  18.292  -1.234  1.00 36.79           N
ATOM    677  CA  ARG A  82      -6.940  19.332  -1.637  1.00 37.04           C
ATOM    678  C   ARG A  82      -6.667  20.658  -0.926  1.00 37.03           C
ATOM    679  O   ARG A  82      -7.143  21.706  -1.358  1.00 37.64           O
ATOM    680  CB  ARG A  82      -8.390  18.881  -1.419  1.00 36.23           C
ATOM    681  CG  ARG A  82      -9.062  18.367  -2.688  1.00 35.81           C
ATOM    682  CD  ARG A  82     -10.531  18.006  -2.481  1.00 37.73           C
ATOM    683  NE  ARG A  82     -10.712  16.696  -1.849  1.00 40.49           N
ATOM    684  CZ  ARG A  82     -11.252  16.505  -0.644  1.00 39.74           C
ATOM    685  NH1 ARG A  82     -11.678  17.543   0.072  1.00 39.03           N1+
ATOM    686  NH2 ARG A  82     -11.374  15.274  -0.157  1.00 30.72           N
ATOM    687  N   GLY A  83      -5.887  20.597   0.151  1.00 37.27           N
ATOM    688  CA  GLY A  83      -5.461  21.785   0.886  1.00 37.93           C
ATOM    689  C   GLY A  83      -4.185  22.396   0.334  1.00 37.78           C
ATOM    690  O   GLY A  83      -4.072  23.615   0.241  1.00 38.32           O
ATOM    691  N   TYR A  84      -3.226  21.549  -0.032  1.00 37.95           N
ATOM    692  CA  TYR A  84      -1.958  22.015  -0.596  1.00 38.47           C
ATOM    693  C   TYR A  84      -2.133  22.637  -1.981  1.00 39.47           C
ATOM    694  O   TYR A  84      -1.556  23.688  -2.273  1.00 40.36           O
ATOM    695  CB  TYR A  84      -0.930  20.880  -0.696  1.00 38.01           C
ATOM    696  CG  TYR A  84      -0.760  19.999   0.527  1.00 36.67           C
ATOM    697  CD1 TYR A  84      -0.739  20.536   1.819  1.00 32.41           C
ATOM    698  CD2 TYR A  84      -0.571  18.625   0.381  1.00 34.90           C
ATOM    699  CE1 TYR A  84      -0.568  19.714   2.935  1.00 33.35           C
ATOM    700  CE2 TYR A  84      -0.395  17.799   1.487  1.00 35.56           C
ATOM    701  CZ  TYR A  84      -0.395  18.346   2.759  1.00 35.03           C
ATOM    702  OH  TYR A  84      -0.216  17.519   3.842  1.00 35.46           O
ATOM    703  N   TYR A  85      -2.914  21.976  -2.833  1.00 39.79           N
ATOM    704  CA  TYR A  85      -3.151  22.454  -4.196  1.00 40.26           C
ATOM    705  C   TYR A  85      -4.271  23.489  -4.261  1.00 40.61           C
ATOM    706  O   TYR A  85      -4.552  24.040  -5.332  1.00 41.23           O
ATOM    707  CB  TYR A  85      -3.448  21.285  -5.138  1.00 40.17           C
ATOM    708  CG  TYR A  85      -2.245  20.411  -5.421  1.00 41.17           C
ATOM    709  CD1 TYR A  85      -1.258  20.814  -6.326  1.00 40.31           C
ATOM    710  CD2 TYR A  85      -2.090  19.183  -4.784  1.00 42.54           C
ATOM    711  CE1 TYR A  85      -0.150  20.011  -6.583  1.00 41.45           C
ATOM    712  CE2 TYR A  85      -0.990  18.373  -5.035  1.00 41.83           C
ATOM    713  CZ  TYR A  85      -0.027  18.789  -5.934  1.00 41.79           C
ATOM    714  OH  TYR A  85       1.060  17.981  -6.168  1.00 39.89           O
ATOM    715  N   ASN A  86      -4.891  23.750  -3.110  1.00 40.00           N
ATOM    716  CA  ASN A  86      -5.997  24.704  -2.984  1.00 40.14           C
ATOM    717  C   ASN A  86      -7.190  24.330  -3.865  1.00 39.88           C
ATOM    718  O   ASN A  86      -7.888  25.194  -4.396  1.00 40.66           O
ATOM    719  CB  ASN A  86      -5.519  26.144  -3.255  1.00 40.63           C
ATOM    720  CG  ASN A  86      -6.491  27.199  -2.736  1.00 39.62           C
ATOM    721  ND2 ASN A  86      -6.339  28.426  -3.217  1.00 40.57           N
ATOM    722  OD1 ASN A  86      -7.364  26.914  -1.917  1.00 35.98           O
ATOM    723  N   GLN A  87      -7.419  23.029  -4.009  1.00 39.58           N
ATOM    724  CA  GLN A  87      -8.476  22.527  -4.874  1.00 39.69           C
ATOM    725  C   GLN A  87      -9.833  22.601  -4.192  1.00 40.47           C
ATOM    726  O   GLN A  87      -9.958  22.334  -2.996  1.00 40.10           O
ATOM    727  CB  GLN A  87      -8.179  21.094  -5.317  1.00 38.77           C
ATOM    728  CG  GLN A  87      -7.017  20.976  -6.291  1.00 37.27           C
ATOM    729  CD  GLN A  87      -6.466  19.562  -6.402  1.00 36.28           C
ATOM    730  NE2 GLN A  87      -6.041  19.190  -7.600  1.00 33.12           N
ATOM    731  OE1 GLN A  87      -6.408  18.821  -5.419  1.00 37.82           O
ATOM    732  N   SER A  88     -10.842  22.973  -4.971  1.00 41.98           N
ATOM    733  CA  SER A  88     -12.218  23.026  -4.495  1.00 43.70           C
ATOM    734  C   SER A  88     -12.807  21.624  -4.442  1.00 44.70           C
ATOM    735  O   SER A  88     -12.270  20.687  -5.038  1.00 45.84           O
ATOM    736  CB  SER A  88     -13.063  23.908  -5.411  1.00 43.59           C
ATOM    737  OG  SER A  88     -13.020  23.436  -6.746  1.00 47.29           O
ATOM    738  N   GLU A  89     -13.920  21.487  -3.732  1.00 45.41           N
ATOM    739  CA  GLU A  89     -14.586  20.201  -3.594  1.00 45.41           C
ATOM    740  C   GLU A  89     -15.505  19.973  -4.794  1.00 45.71           C
ATOM    741  O   GLU A  89     -16.462  19.201  -4.722  1.00 46.79           O
ATOM    742  CB  GLU A  89     -15.366  20.122  -2.270  1.00 45.56           C
ATOM    743  CG  GLU A  89     -14.943  21.116  -1.168  1.00 48.37           C
ATOM    744  CD  GLU A  89     -13.436  21.169  -0.903  1.00 50.97           C
ATOM    745  OE1 GLU A  89     -12.779  20.104  -0.877  1.00 54.54           O
ATOM    746  OE2 GLU A  89     -12.911  22.288  -0.715  1.00 50.84           O1-
ATOM    747  N   ALA A  90     -15.206  20.654  -5.897  1.00 45.23           N
ATOM    748  CA  ALA A  90     -15.972  20.518  -7.134  1.00 44.69           C
ATOM    749  C   ALA A  90     -15.432  19.388  -8.008  1.00 44.43           C
ATOM    750  O   ALA A  90     -16.077  18.984  -8.980  1.00 45.26           O
ATOM    751  CB  ALA A  90     -15.975  21.838  -7.906  1.00 44.38           C
ATOM    752  N   GLY A  91     -14.253  18.882  -7.659  1.00 43.20           N
ATOM    753  CA  GLY A  91     -13.607  17.840  -8.450  1.00 42.14           C
ATOM    754  C   GLY A  91     -13.277  16.554  -7.710  1.00 40.73           C
ATOM    755  O   GLY A  91     -12.956  16.561  -6.517  1.00 40.26           O
ATOM    756  N   SER A  92     -13.362  15.446  -8.441  1.00 39.27           N
ATOM    757  CA  SER A  92     -12.926  14.147  -7.950  1.00 37.86           C
ATOM    758  C   SER A  92     -11.418  14.026  -8.133  1.00 36.54           C
ATOM    759  O   SER A  92     -10.864  14.478  -9.135  1.00 35.66           O
ATOM    760  CB  SER A  92     -13.651  13.026  -8.691  1.00 38.04           C
ATOM    761  OG  SER A  92     -15.037  13.052  -8.402  1.00 37.15           O
ATOM    762  N   HIS A  93     -10.758  13.423  -7.151  1.00 35.65           N
ATOM    763  CA  HIS A  93      -9.298  13.336  -7.137  1.00 34.29           C
ATOM    764  C   HIS A  93      -8.815  11.965  -6.652  1.00 32.23           C
ATOM    765  O   HIS A  93      -9.490  11.303  -5.863  1.00 31.55           O
ATOM    766  CB  HIS A  93      -8.716  14.458  -6.277  1.00 34.80           C
ATOM    767  CG  HIS A  93      -8.916  15.826  -6.854  1.00 36.88           C
ATOM    768  CD2 HIS A  93      -8.090  16.619  -7.577  1.00 38.28           C
ATOM    769  ND1 HIS A  93     -10.091  16.532  -6.708  1.00 38.72           N
ATOM    770  CE1 HIS A  93      -9.979  17.700  -7.316  1.00 39.20           C
ATOM    771  NE2 HIS A  93      -8.775  17.777  -7.853  1.00 38.93           N
ATOM    772  N   ILE A  94      -7.640  11.555  -7.125  1.00 30.28           N
ATOM    773  CA  ILE A  94      -7.151  10.192  -6.929  1.00 28.19           C
ATOM    774  C   ILE A  94      -5.724  10.141  -6.387  1.00 27.54           C
ATOM    775  O   ILE A  94      -4.796  10.695  -6.990  1.00 27.67           O
ATOM    776  CB  ILE A  94      -7.236   9.372  -8.253  1.00 27.93           C
ATOM    777  CG1 ILE A  94      -8.691   9.260  -8.724  1.00 24.79           C
ATOM    778  CG2 ILE A  94      -6.619   7.976  -8.080  1.00 30.04           C
ATOM    779  CD1 ILE A  94      -8.846   8.999 -10.204  1.00 19.52           C
ATOM    780  N   ILE A  95      -5.561   9.472  -5.249  1.00 26.04           N
ATOM    781  CA  ILE A  95      -4.238   9.103  -4.744  1.00 25.77           C
ATOM    782  C   ILE A  95      -3.938   7.646  -5.104  1.00 24.48           C
ATOM    783  O   ILE A  95      -4.709   6.746  -4.774  1.00 24.72           O
ATOM    784  CB  ILE A  95      -4.120   9.281  -3.210  1.00 25.42           C
ATOM    785  CG1 ILE A  95      -4.225  10.755  -2.827  1.00 27.63           C
ATOM    786  CG2 ILE A  95      -2.803   8.718  -2.697  1.00 25.83           C
ATOM    787  CD1 ILE A  95      -4.221  11.000  -1.326  1.00 26.97           C
ATOM    788  N   GLN A  96      -2.818   7.422  -5.782  1.00 23.96           N
ATOM    789  CA  GLN A  96      -2.364   6.069  -6.090  1.00 23.35           C
ATOM    790  C   GLN A  96      -1.089   5.735  -5.326  1.00 23.46           C
ATOM    791  O   GLN A  96      -0.238   6.596  -5.127  1.00 22.99           O
ATOM    792  CB  GLN A  96      -2.165   5.888  -7.596  1.00 22.44           C
ATOM    793  CG  GLN A  96      -3.445   5.477  -8.320  1.00 22.88           C
ATOM    794  CD  GLN A  96      -3.487   5.931  -9.761  1.00 19.41           C
ATOM    795  NE2 GLN A  96      -3.299   4.995 -10.680  1.00 19.68           N
ATOM    796  OE1 GLN A  96      -3.696   7.109 -10.046  1.00 21.13           O
ATOM    797  N   ARG A  97      -0.980   4.485  -4.884  1.00 23.45           N
ATOM    798  CA  ARG A  97       0.219   4.017  -4.204  1.00 23.74           C
ATOM    799  C   ARG A  97       0.560   2.592  -4.618  1.00 24.81           C
ATOM    800  O   ARG A  97      -0.323   1.731  -4.712  1.00 25.10           O
ATOM    801  CB  ARG A  97       0.049   4.088  -2.688  1.00 23.64           C
ATOM    802  CG  ARG A  97       1.363   4.106  -1.902  1.00 23.51           C
ATOM    803  CD  ARG A  97       1.186   3.442  -0.554  1.00 24.00           C
ATOM    804  NE  ARG A  97       2.344   3.602   0.316  1.00 25.73           N
ATOM    805  CZ  ARG A  97       2.444   4.517   1.279  1.00 30.23           C
ATOM    806  NH1 ARG A  97       1.454   5.372   1.504  1.00 28.51           N1+
ATOM    807  NH2 ARG A  97       3.546   4.582   2.016  1.00 31.62           N
ATOM    808  N   MET A  98       1.841   2.352  -4.872  1.00 24.77           N
ATOM    809  CA  MET A  98       2.330   0.989  -5.078  1.00 25.78           C
ATOM    810  C   MET A  98       3.654   0.742  -4.350  1.00 24.35           C
ATOM    811  O   MET A  98       4.626   1.487  -4.509  1.00 25.84           O
ATOM    812  CB  MET A  98       2.410   0.624  -6.573  1.00 25.53           C
ATOM    813  CG  MET A  98       3.284   1.532  -7.433  1.00 26.77           C
ATOM    814  SD  MET A  98       3.048   1.263  -9.204  1.00 29.48           S
ATOM    815  CE  MET A  98       4.198   2.471  -9.871  1.00 25.94           C
ATOM    816  N   TYR A  99       3.672  -0.305  -3.538  1.00 22.36           N
ATOM    817  CA  TYR A  99       4.868  -0.685  -2.798  1.00 21.12           C
ATOM    818  C   TYR A  99       5.086  -2.199  -2.818  1.00 20.83           C
ATOM    819  O   TYR A  99       4.169  -2.964  -3.127  1.00 20.62           O
ATOM    820  CB  TYR A  99       4.793  -0.167  -1.352  1.00 20.57           C
ATOM    821  CG  TYR A  99       3.692  -0.781  -0.517  1.00 16.08           C
ATOM    822  CD1 TYR A  99       2.427  -0.192  -0.458  1.00 17.05           C
ATOM    823  CD2 TYR A  99       3.915  -1.939   0.218  1.00 14.70           C
ATOM    824  CE1 TYR A  99       1.411  -0.747   0.309  1.00 18.90           C
ATOM    825  CE2 TYR A  99       2.901  -2.514   0.985  1.00 16.62           C
ATOM    826  CZ  TYR A  99       1.658  -1.909   1.032  1.00 18.89           C
ATOM    827  OH  TYR A  99       0.659  -2.463   1.792  1.00 15.79           O
ATOM    828  N   GLY A 100       6.299  -2.623  -2.480  1.00 20.32           N
ATOM    829  CA  GLY A 100       6.605  -4.042  -2.391  1.00 20.49           C
ATOM    830  C   GLY A 100       8.072  -4.357  -2.592  1.00 20.44           C
ATOM    831  O   GLY A 100       8.885  -3.466  -2.814  1.00 19.31           O
ATOM    832  N   CYS A 101       8.398  -5.642  -2.509  1.00 21.69           N
ATOM    833  CA  CYS A 101       9.777  -6.107  -2.601  1.00 22.98           C
ATOM    834  C   CYS A 101       9.952  -7.106  -3.740  1.00 22.89           C
ATOM    835  O   CYS A 101       9.045  -7.875  -4.047  1.00 23.10           O
ATOM    836  CB  CYS A 101      10.216  -6.727  -1.273  1.00 21.75           C
ATOM    837  SG  CYS A 101       9.039  -7.926  -0.576  1.00 32.48           S
ATOM    838  N   ASP A 102      11.118  -7.071  -4.376  1.00 23.87           N
ATOM    839  CA  ASP A 102      11.464  -8.031  -5.420  1.00 24.34           C
ATOM    840  C   ASP A 102      12.475  -9.009  -4.860  1.00 24.26           C
ATOM    841  O   ASP A 102      13.326  -8.629  -4.055  1.00 24.33           O
ATOM    842  CB  ASP A 102      12.081  -7.337  -6.641  1.00 25.28           C
ATOM    843  CG  ASP A 102      11.315  -6.100  -7.084  1.00 25.67           C
ATOM    844  OD1 ASP A 102      11.548  -5.648  -8.220  1.00 29.39           O
ATOM    845  OD2 ASP A 102      10.497  -5.566  -6.315  1.00 30.83           O1-
ATOM    846  N   LEU A 103      12.380 -10.265  -5.282  1.00 24.63           N
ATOM    847  CA  LEU A 103      13.395 -11.260  -4.957  1.00 25.47           C
ATOM    848  C   LEU A 103      14.039 -11.760  -6.237  1.00 26.88           C
ATOM    849  O   LEU A 103      13.392 -11.807  -7.287  1.00 27.67           O
ATOM    850  CB  LEU A 103      12.810 -12.435  -4.166  1.00 24.83           C
ATOM    851  CG  LEU A 103      12.186 -12.220  -2.778  1.00 23.77           C
ATOM    852  CD1 LEU A 103      11.833 -13.563  -2.151  1.00 23.06           C
ATOM    853  CD2 LEU A 103      13.065 -11.429  -1.843  1.00 18.87           C
ATOM    854  N   GLY A 104      15.316 -12.128  -6.146  1.00 27.25           N
ATOM    855  CA  GLY A 104      16.047 -12.661  -7.292  1.00 27.71           C
ATOM    856  C   GLY A 104      15.728 -14.118  -7.552  1.00 29.03           C
ATOM    857  O   GLY A 104      15.064 -14.762  -6.736  1.00 28.79           O
ATOM    858  N   PRO A 105      16.200 -14.655  -8.696  1.00 30.48           N
ATOM    859  CA  PRO A 105      16.041 -16.087  -8.987  1.00 30.89           C
ATOM    860  C   PRO A 105      16.800 -16.954  -7.987  1.00 31.41           C
ATOM    861  O   PRO A 105      17.065 -18.123  -8.248  1.00 32.49           O
ATOM    862  CB  PRO A 105      16.649 -16.233 -10.387  1.00 30.50           C
ATOM    863  CG  PRO A 105      17.536 -15.051 -10.556  1.00 30.31           C
ATOM    864  CD  PRO A 105      16.909 -13.946  -9.780  1.00 30.08           C
ATOM    865  N   ASP A 106      17.114 -16.365  -6.839  1.00 32.47           N
ATOM    866  CA  ASP A 106      17.978 -16.961  -5.827  1.00 32.66           C
ATOM    867  C   ASP A 106      17.287 -16.989  -4.469  1.00 31.89           C
ATOM    868  O   ASP A 106      17.724 -17.698  -3.563  1.00 32.13           O
ATOM    869  CB  ASP A 106      19.267 -16.139  -5.718  1.00 33.11           C
ATOM    870  CG  ASP A 106      19.000 -14.637  -5.570  1.00 34.33           C
ATOM    871  OD1 ASP A 106      19.825 -13.841  -6.037  1.00 38.30           O
ATOM    872  OD2 ASP A 106      17.968 -14.240  -4.992  1.00 37.92           O1-
ATOM    873  N   GLY A 107      16.219 -16.202  -4.338  1.00 31.17           N
ATOM    874  CA  GLY A 107      15.498 -16.058  -3.073  1.00 31.23           C
ATOM    875  C   GLY A 107      15.867 -14.788  -2.323  1.00 31.52           C
ATOM    876  O   GLY A 107      15.088 -14.299  -1.503  1.00 30.97           O
ATOM    877  N   ARG A 108      17.056 -14.258  -2.604  1.00 31.50           N
ATOM    878  CA  ARG A 108      17.555 -13.053  -1.953  1.00 32.39           C
ATOM    879  C   ARG A 108      16.839 -11.802  -2.436  1.00 33.01           C
ATOM    880  O   ARG A 108      16.415 -11.732  -3.590  1.00 34.18           O
ATOM    881  CB  ARG A 108      19.060 -12.904  -2.185  1.00 33.15           C
ATOM    882  CG  ARG A 108      19.930 -13.542  -1.115  1.00 32.82           C
ATOM    883  CD  ARG A 108      20.160 -15.012  -1.372  1.00 32.19           C
ATOM    884  NE  ARG A 108      20.500 -15.721  -0.141  1.00 33.35           N
ATOM    885  CZ  ARG A 108      20.876 -16.995  -0.090  1.00 33.96           C
ATOM    886  NH1 ARG A 108      20.975 -17.710  -1.202  1.00 36.02           N1+
ATOM    887  NH2 ARG A 108      21.152 -17.556   1.076  1.00 35.85           N
ATOM    888  N   LEU A 109      16.719 -10.820  -1.546  1.00 32.96           N
ATOM    889  CA  LEU A 109      16.087  -9.542  -1.856  1.00 32.99           C
ATOM    890  C   LEU A 109      16.917  -8.725  -2.846  1.00 33.35           C
ATOM    891  O   LEU A 109      18.108  -8.495  -2.629  1.00 34.29           O
ATOM    892  CB  LEU A 109      15.851  -8.748  -0.566  1.00 33.01           C
ATOM    893  CG  LEU A 109      15.441  -7.271  -0.633  1.00 34.72           C
ATOM    894  CD1 LEU A 109      14.025  -7.084  -1.180  1.00 32.65           C
ATOM    895  CD2 LEU A 109      15.570  -6.631   0.744  1.00 33.59           C
ATOM    896  N   LEU A 110      16.269  -8.304  -3.931  1.00 33.12           N
ATOM    897  CA  LEU A 110      16.879  -7.456  -4.959  1.00 32.42           C
ATOM    898  C   LEU A 110      16.740  -5.970  -4.625  1.00 32.69           C
ATOM    899  O   LEU A 110      17.739  -5.245  -4.598  1.00 33.12           O
ATOM    900  CB  LEU A 110      16.266  -7.756  -6.331  1.00 32.41           C
ATOM    901  CG  LEU A 110      17.052  -8.532  -7.399  1.00 31.09           C
ATOM    902  CD1 LEU A 110      17.943  -9.611  -6.824  1.00 32.44           C
ATOM    903  CD2 LEU A 110      16.092  -9.130  -8.413  1.00 31.55           C
ATOM    904  N   ARG A 111      15.504  -5.525  -4.393  1.00 32.20           N
ATOM    905  CA  ARG A 111      15.221  -4.163  -3.925  1.00 31.94           C
ATOM    906  C   ARG A 111      13.784  -3.971  -3.427  1.00 31.24           C
ATOM    907  O   ARG A 111      12.932  -4.850  -3.589  1.00 30.81           O
ATOM    908  CB  ARG A 111      15.580  -3.111  -4.982  1.00 31.76           C
ATOM    909  CG  ARG A 111      14.835  -3.208  -6.306  1.00 34.76           C
ATOM    910  CD  ARG A 111      15.242  -2.066  -7.242  1.00 34.07           C
ATOM    911  NE  ARG A 111      15.081  -0.756  -6.604  1.00 36.26           N
ATOM    912  CZ  ARG A 111      15.157   0.413  -7.234  1.00 35.91           C
ATOM    913  NH1 ARG A 111      15.396   0.468  -8.540  1.00 34.15           N1+
ATOM    914  NH2 ARG A 111      14.987   1.539  -6.547  1.00 34.35           N
ATOM    915  N   GLY A 112      13.541  -2.822  -2.797  1.00 30.89           N
ATOM    916  CA  GLY A 112      12.208  -2.430  -2.344  1.00 29.82           C
ATOM    917  C   GLY A 112      11.663  -1.258  -3.144  1.00 29.01           C
ATOM    918  O   GLY A 112      12.421  -0.544  -3.805  1.00 29.81           O
ATOM    919  N   HIS A 113      10.347  -1.065  -3.091  1.00 27.67           N
ATOM    920  CA  HIS A 113       9.693   0.027  -3.803  1.00 27.53           C
ATOM    921  C   HIS A 113       8.572   0.651  -2.977  1.00 26.77           C
ATOM    922  O   HIS A 113       7.845  -0.054  -2.282  1.00 24.99           O
ATOM    923  CB  HIS A 113       9.130  -0.464  -5.143  1.00 28.35           C
ATOM    924  CG  HIS A 113      10.150  -1.102  -6.033  1.00 31.39           C
ATOM    925  CD2 HIS A 113      10.509  -2.399  -6.191  1.00 33.24           C
ATOM    926  ND1 HIS A 113      10.946  -0.378  -6.895  1.00 31.59           N
ATOM    927  CE1 HIS A 113      11.747  -1.201  -7.547  1.00 33.28           C
ATOM    928  NE2 HIS A 113      11.501  -2.433  -7.138  1.00 34.63           N
ATOM    929  N   ASP A 114       8.451   1.977  -3.060  1.00 27.22           N
ATOM    930  CA  ASP A 114       7.312   2.710  -2.488  1.00 28.37           C
ATOM    931  C   ASP A 114       7.043   3.988  -3.273  1.00 27.79           C
ATOM    932  O   ASP A 114       7.745   4.989  -3.101  1.00 28.30           O
ATOM    933  CB  ASP A 114       7.538   3.041  -1.004  1.00 29.58           C
ATOM    934  CG  ASP A 114       6.406   3.884  -0.407  1.00 33.95           C
ATOM    935  OD1 ASP A 114       5.233   3.699  -0.798  1.00 40.68           O
ATOM    936  OD2 ASP A 114       6.686   4.740   0.454  1.00 37.68           O1-
ATOM    937  N   GLN A 115       6.024   3.955  -4.127  1.00 27.34           N
ATOM    938  CA  GLN A 115       5.686   5.116  -4.954  1.00 27.46           C
ATOM    939  C   GLN A 115       4.235   5.535  -4.805  1.00 27.39           C
ATOM    940  O   GLN A 115       3.352   4.692  -4.638  1.00 28.29           O
ATOM    941  CB  GLN A 115       5.984   4.840  -6.424  1.00 27.51           C
ATOM    942  CG  GLN A 115       7.456   4.736  -6.762  1.00 25.61           C
ATOM    943  CD  GLN A 115       7.697   3.816  -7.936  1.00 26.60           C
ATOM    944  NE2 GLN A 115       7.986   4.397  -9.091  1.00 29.87           N
ATOM    945  OE1 GLN A 115       7.616   2.595  -7.808  1.00 25.82           O
ATOM    946  N   SER A 116       4.005   6.843  -4.872  1.00 27.29           N
ATOM    947  CA  SER A 116       2.659   7.417  -4.799  1.00 27.12           C
ATOM    948  C   SER A 116       2.413   8.461  -5.896  1.00 26.72           C
ATOM    949  O   SER A 116       3.346   9.101  -6.384  1.00 26.82           O
ATOM    950  CB  SER A 116       2.416   8.039  -3.421  1.00 27.53           C
ATOM    951  OG  SER A 116       2.774   7.149  -2.375  1.00 26.87           O
ATOM    952  N   ALA A 117       1.149   8.624  -6.273  1.00 26.85           N
ATOM    953  CA  ALA A 117       0.763   9.556  -7.330  1.00 26.63           C
ATOM    954  C   ALA A 117      -0.529  10.279  -6.991  1.00 26.70           C
ATOM    955  O   ALA A 117      -1.392   9.733  -6.302  1.00 25.99           O
ATOM    956  CB  ALA A 117       0.612   8.820  -8.657  1.00 26.50           C
ATOM    957  N   TYR A 118      -0.652  11.508  -7.487  1.00 27.29           N
ATOM    958  CA  TYR A 118      -1.880  12.285  -7.357  1.00 27.50           C
ATOM    959  C   TYR A 118      -2.354  12.707  -8.734  1.00 28.51           C
ATOM    960  O   TYR A 118      -1.580  13.242  -9.530  1.00 28.64           O
ATOM    961  CB  TYR A 118      -1.663  13.506  -6.463  1.00 27.43           C
ATOM    962  CG  TYR A 118      -2.931  14.236  -6.059  1.00 28.51           C
ATOM    963  CD1 TYR A 118      -4.049  13.535  -5.592  1.00 27.19           C
ATOM    964  CD2 TYR A 118      -3.005  15.630  -6.112  1.00 26.21           C
ATOM    965  CE1 TYR A 118      -5.209  14.197  -5.210  1.00 24.53           C
ATOM    966  CE2 TYR A 118      -4.163  16.304  -5.726  1.00 25.90           C
ATOM    967  CZ  TYR A 118      -5.261  15.579  -5.278  1.00 26.76           C
ATOM    968  OH  TYR A 118      -6.412  16.233  -4.888  1.00 27.04           O
ATOM    969  N   ASP A 119      -3.631  12.448  -9.004  1.00 29.69           N
ATOM    970  CA  ASP A 119      -4.256  12.701 -10.302  1.00 29.79           C
ATOM    971  C   ASP A 119      -3.442  12.130 -11.462  1.00 30.74           C
ATOM    972  O   ASP A 119      -3.359  12.737 -12.533  1.00 32.97           O
ATOM    973  CB  ASP A 119      -4.554  14.200 -10.498  1.00 29.76           C
ATOM    974  CG  ASP A 119      -5.613  14.730  -9.526  1.00 30.62           C
ATOM    975  OD1 ASP A 119      -6.466  13.942  -9.057  1.00 33.07           O
ATOM    976  OD2 ASP A 119      -5.595  15.944  -9.231  1.00 30.36           O1-
ATOM    977  N   GLY A 120      -2.840  10.964 -11.241  1.00 30.92           N
ATOM    978  CA  GLY A 120      -2.116  10.239 -12.288  1.00 30.72           C
ATOM    979  C   GLY A 120      -0.671  10.642 -12.503  1.00 31.55           C
ATOM    980  O   GLY A 120      -0.005  10.115 -13.396  1.00 32.31           O
ATOM    981  N   LYS A 121      -0.181  11.568 -11.683  1.00 31.88           N
ATOM    982  CA  LYS A 121       1.167  12.109 -11.824  1.00 32.52           C
ATOM    983  C   LYS A 121       2.031  11.788 -10.611  1.00 31.94           C
ATOM    984  O   LYS A 121       1.557  11.868  -9.475  1.00 32.41           O
ATOM    985  CB  LYS A 121       1.098  13.625 -11.983  1.00 32.53           C
ATOM    986  CG  LYS A 121       0.748  14.114 -13.379  1.00 35.99           C
ATOM    987  CD  LYS A 121       0.768  15.641 -13.434  1.00 34.68           C
ATOM    988  CE  LYS A 121       2.127  16.204 -13.022  1.00 36.94           C
ATOM    989  NZ  LYS A 121       1.993  17.517 -12.321  1.00 39.11           N1+
ATOM    990  N   ASP A 122       3.301  11.456 -10.850  1.00 31.14           N
ATOM    991  CA  ASP A 122       4.262  11.190  -9.774  1.00 29.31           C
ATOM    992  C   ASP A 122       4.159  12.264  -8.702  1.00 28.83           C
ATOM    993  O   ASP A 122       4.034  13.452  -9.013  1.00 28.97           O
ATOM    994  CB  ASP A 122       5.695  11.140 -10.311  1.00 28.84           C
ATOM    995  CG  ASP A 122       6.047   9.807 -10.958  1.00 30.61           C
ATOM    996  OD1 ASP A 122       5.199   8.889 -11.007  1.00 30.36           O
ATOM    997  OD2 ASP A 122       7.199   9.677 -11.423  1.00 33.81           O1-
ATOM    998  N   TYR A 123       4.196  11.841  -7.443  1.00 27.99           N
ATOM    999  CA  TYR A 123       4.018  12.754  -6.323  1.00 27.92           C
ATOM   1000  C   TYR A 123       5.169  12.641  -5.333  1.00 29.11           C
ATOM   1001  O   TYR A 123       5.931  13.586  -5.143  1.00 29.96           O
ATOM   1002  CB  TYR A 123       2.664  12.506  -5.647  1.00 26.77           C
ATOM   1003  CG  TYR A 123       2.404  13.391  -4.457  1.00 24.05           C
ATOM   1004  CD1 TYR A 123       2.010  14.718  -4.620  1.00 22.57           C
ATOM   1005  CD2 TYR A 123       2.546  12.901  -3.166  1.00 25.92           C
ATOM   1006  CE1 TYR A 123       1.772  15.530  -3.525  1.00 20.89           C
ATOM   1007  CE2 TYR A 123       2.314  13.706  -2.062  1.00 24.57           C
ATOM   1008  CZ  TYR A 123       1.927  15.013  -2.249  1.00 23.08           C
ATOM   1009  OH  TYR A 123       1.696  15.795  -1.148  1.00 25.28           O
ATOM   1010  N   ILE A 124       5.300  11.478  -4.708  1.00 30.18           N
ATOM   1011  CA  ILE A 124       6.429  11.214  -3.827  1.00 30.97           C
ATOM   1012  C   ILE A 124       6.848   9.751  -3.973  1.00 31.43           C
ATOM   1013  O   ILE A 124       6.014   8.874  -4.224  1.00 31.23           O
ATOM   1014  CB  ILE A 124       6.103  11.594  -2.351  1.00 31.62           C
ATOM   1015  CG1 ILE A 124       7.386  11.753  -1.522  1.00 32.29           C
ATOM   1016  CG2 ILE A 124       5.117  10.600  -1.725  1.00 33.19           C
ATOM   1017  CD1 ILE A 124       7.186  12.429  -0.170  1.00 30.17           C
ATOM   1018  N   ALA A 125       8.145   9.502  -3.843  1.00 31.94           N
ATOM   1019  CA  ALA A 125       8.687   8.156  -3.977  1.00 33.03           C
ATOM   1020  C   ALA A 125       9.846   7.945  -3.022  1.00 33.72           C
ATOM   1021  O   ALA A 125      10.673   8.839  -2.836  1.00 34.31           O
ATOM   1022  CB  ALA A 125       9.134   7.900  -5.413  1.00 32.68           C
ATOM   1023  N   LEU A 126       9.898   6.762  -2.417  1.00 34.17           N
ATOM   1024  CA  LEU A 126      11.040   6.369  -1.601  1.00 35.25           C
ATOM   1025  C   LEU A 126      12.250   6.137  -2.500  1.00 36.07           C
ATOM   1026  O   LEU A 126      12.110   5.612  -3.609  1.00 37.21           O
ATOM   1027  CB  LEU A 126      10.722   5.108  -0.794  1.00 34.47           C
ATOM   1028  CG  LEU A 126      11.638   4.743   0.381  1.00 35.79           C
ATOM   1029  CD1 LEU A 126      11.690   5.856   1.423  1.00 33.07           C
ATOM   1030  CD2 LEU A 126      11.184   3.434   1.029  1.00 35.53           C
ATOM   1031  N   ASN A 127      13.424   6.542  -2.021  1.00 36.41           N
ATOM   1032  CA  ASN A 127      14.670   6.374  -2.769  1.00 37.74           C
ATOM   1033  C   ASN A 127      15.244   4.968  -2.626  1.00 39.42           C
ATOM   1034  O   ASN A 127      14.735   4.161  -1.840  1.00 40.37           O
ATOM   1035  CB  ASN A 127      15.701   7.419  -2.334  1.00 37.42           C
ATOM   1036  CG  ASN A 127      15.269   8.840  -2.654  1.00 37.39           C
ATOM   1037  ND2 ASN A 127      15.956   9.810  -2.070  1.00 37.09           N
ATOM   1038  OD1 ASN A 127      14.333   9.063  -3.421  1.00 38.26           O
ATOM   1039  N   GLU A 128      16.306   4.687  -3.385  1.00 40.63           N
ATOM   1040  CA  GLU A 128      16.975   3.381  -3.367  1.00 41.78           C
ATOM   1041  C   GLU A 128      17.625   3.072  -2.012  1.00 40.73           C
ATOM   1042  O   GLU A 128      17.788   1.906  -1.652  1.00 41.64           O
ATOM   1043  CB  GLU A 128      17.994   3.271  -4.514  1.00 41.82           C
ATOM   1044  CG  GLU A 128      18.756   1.938  -4.573  1.00 45.28           C
ATOM   1045  CD  GLU A 128      19.320   1.599  -5.956  1.00 44.65           C
ATOM   1046  OE1 GLU A 128      18.845   2.162  -6.971  1.00 46.00           O
ATOM   1047  OE2 GLU A 128      20.232   0.741  -6.022  1.00 47.35           O1-
ATOM   1048  N   ASP A 129      17.971   4.112  -1.256  1.00 39.25           N
ATOM   1049  CA  ASP A 129      18.573   3.922   0.067  1.00 38.67           C
ATOM   1050  C   ASP A 129      17.567   3.523   1.154  1.00 37.95           C
ATOM   1051  O   ASP A 129      17.959   3.014   2.205  1.00 38.50           O
ATOM   1052  CB  ASP A 129      19.388   5.155   0.497  1.00 38.49           C
ATOM   1053  CG  ASP A 129      18.545   6.415   0.644  1.00 38.51           C
ATOM   1054  OD1 ASP A 129      17.305   6.321   0.772  1.00 39.58           O
ATOM   1055  OD2 ASP A 129      19.132   7.516   0.636  1.00 38.14           O1-
ATOM   1056  N   LEU A 130      16.281   3.767   0.897  1.00 36.30           N
ATOM   1057  CA  LEU A 130      15.197   3.501   1.855  1.00 34.98           C
ATOM   1058  C   LEU A 130      15.288   4.396   3.104  1.00 34.40           C
ATOM   1059  O   LEU A 130      14.730   4.086   4.159  1.00 33.68           O
ATOM   1060  CB  LEU A 130      15.119   2.005   2.226  1.00 34.87           C
ATOM   1061  CG  LEU A 130      15.145   0.926   1.127  1.00 34.54           C
ATOM   1062  CD1 LEU A 130      14.909  -0.452   1.727  1.00 33.77           C
ATOM   1063  CD2 LEU A 130      14.139   1.191   0.001  1.00 33.45           C
ATOM   1064  N   SER A 131      15.985   5.520   2.958  1.00 34.32           N
ATOM   1065  CA  SER A 131      16.152   6.494   4.031  1.00 33.54           C
ATOM   1066  C   SER A 131      15.587   7.864   3.653  1.00 33.56           C
ATOM   1067  O   SER A 131      15.402   8.721   4.522  1.00 33.54           O
ATOM   1068  CB  SER A 131      17.631   6.632   4.390  1.00 32.81           C
ATOM   1069  OG  SER A 131      18.214   5.363   4.618  1.00 33.26           O
ATOM   1070  N   SER A 132      15.308   8.065   2.364  1.00 32.92           N
ATOM   1071  CA  SER A 132      14.921   9.388   1.860  1.00 32.24           C
ATOM   1072  C   SER A 132      13.846   9.369   0.769  1.00 32.50           C
ATOM   1073  O   SER A 132      13.554   8.327   0.172  1.00 32.13           O
ATOM   1074  CB  SER A 132      16.159  10.133   1.360  1.00 31.61           C
ATOM   1075  OG  SER A 132      16.978   9.274   0.593  1.00 30.95           O
ATOM   1076  N   TRP A 133      13.273  10.545   0.516  1.00 32.22           N
ATOM   1077  CA  TRP A 133      12.188  10.710  -0.441  1.00 31.56           C
ATOM   1078  C   TRP A 133      12.549  11.671  -1.560  1.00 31.42           C
ATOM   1079  O   TRP A 133      13.253  12.660  -1.346  1.00 30.33           O
ATOM   1080  CB  TRP A 133      10.947  11.251   0.266  1.00 31.81           C
ATOM   1081  CG  TRP A 133      10.584  10.512   1.510  1.00 31.83           C
ATOM   1082  CD1 TRP A 133      10.961  10.818   2.787  1.00 30.15           C
ATOM   1083  CD2 TRP A 133       9.773   9.340   1.599  1.00 29.89           C
ATOM   1084  CE2 TRP A 133       9.696   8.989   2.964  1.00 29.50           C
ATOM   1085  CE3 TRP A 133       9.098   8.553   0.657  1.00 30.84           C
ATOM   1086  NE1 TRP A 133      10.433   9.905   3.666  1.00 30.12           N
ATOM   1087  CZ2 TRP A 133       8.975   7.882   3.411  1.00 30.80           C
ATOM   1088  CZ3 TRP A 133       8.379   7.451   1.100  1.00 30.26           C
ATOM   1089  CH2 TRP A 133       8.323   7.130   2.466  1.00 32.00           C
ATOM   1090  N   THR A 134      12.042  11.372  -2.752  1.00 32.38           N
ATOM   1091  CA  THR A 134      12.043  12.310  -3.865  1.00 33.52           C
ATOM   1092  C   THR A 134      10.628  12.868  -4.028  1.00 33.83           C
ATOM   1093  O   THR A 134       9.680  12.111  -4.229  1.00 34.49           O
ATOM   1094  CB  THR A 134      12.494  11.628  -5.164  1.00 33.61           C
ATOM   1095  CG2 THR A 134      12.369  12.569  -6.348  1.00 35.09           C
ATOM   1096  OG1 THR A 134      13.860  11.218  -5.037  1.00 35.36           O
ATOM   1097  N   ALA A 135      10.491  14.187  -3.910  1.00 34.29           N
ATOM   1098  CA  ALA A 135       9.202  14.854  -4.101  1.00 34.58           C
ATOM   1099  C   ALA A 135       9.134  15.458  -5.493  1.00 35.85           C
ATOM   1100  O   ALA A 135      10.126  16.004  -5.989  1.00 36.56           O
ATOM   1101  CB  ALA A 135       8.991  15.920  -3.050  1.00 34.10           C
ATOM   1102  N   ALA A 136       7.963  15.358  -6.119  1.00 36.18           N
ATOM   1103  CA  ALA A 136       7.777  15.811  -7.497  1.00 36.44           C
ATOM   1104  C   ALA A 136       7.522  17.316  -7.591  1.00 37.18           C
ATOM   1105  O   ALA A 136       8.148  18.001  -8.402  1.00 38.05           O
ATOM   1106  CB  ALA A 136       6.663  15.028  -8.169  1.00 36.15           C
ATOM   1107  N   ASP A 137       6.604  17.824  -6.772  1.00 36.94           N
ATOM   1108  CA  ASP A 137       6.377  19.266  -6.664  1.00 36.91           C
ATOM   1109  C   ASP A 137       6.364  19.739  -5.206  1.00 36.43           C
ATOM   1110  O   ASP A 137       6.534  18.936  -4.284  1.00 36.02           O
ATOM   1111  CB  ASP A 137       5.108  19.696  -7.418  1.00 37.49           C
ATOM   1112  CG  ASP A 137       3.843  18.992  -6.928  1.00 38.86           C
ATOM   1113  OD1 ASP A 137       3.862  18.301  -5.885  1.00 40.80           O
ATOM   1114  OD2 ASP A 137       2.805  19.142  -7.602  1.00 40.44           O1-
ATOM   1115  N   THR A 138       6.150  21.038  -5.007  1.00 36.25           N
ATOM   1116  CA  THR A 138       6.219  21.654  -3.673  1.00 36.64           C
ATOM   1117  C   THR A 138       5.119  21.207  -2.698  1.00 35.13           C
ATOM   1118  O   THR A 138       5.281  21.307  -1.478  1.00 33.78           O
ATOM   1119  CB  THR A 138       6.252  23.193  -3.758  1.00 36.72           C
ATOM   1120  CG2 THR A 138       7.563  23.655  -4.371  1.00 37.20           C
ATOM   1121  OG1 THR A 138       5.165  23.652  -4.571  1.00 40.24           O
ATOM   1122  N   ALA A 139       4.012  20.711  -3.241  1.00 34.63           N
ATOM   1123  CA  ALA A 139       2.963  20.111  -2.423  1.00 34.91           C
ATOM   1124  C   ALA A 139       3.466  18.832  -1.748  1.00 34.50           C
ATOM   1125  O   ALA A 139       3.200  18.599  -0.567  1.00 34.37           O
ATOM   1126  CB  ALA A 139       1.729  19.835  -3.256  1.00 34.84           C
ATOM   1127  N   ALA A 140       4.222  18.031  -2.498  1.00 34.25           N
ATOM   1128  CA  ALA A 140       4.805  16.793  -1.981  1.00 34.05           C
ATOM   1129  C   ALA A 140       5.950  17.039  -1.002  1.00 34.17           C
ATOM   1130  O   ALA A 140       6.312  16.150  -0.227  1.00 35.28           O
ATOM   1131  CB  ALA A 140       5.262  15.913  -3.123  1.00 33.80           C
ATOM   1132  N   GLN A 141       6.510  18.247  -1.043  1.00 33.29           N
ATOM   1133  CA  GLN A 141       7.589  18.649  -0.146  1.00 32.73           C
ATOM   1134  C   GLN A 141       7.066  18.920   1.269  1.00 32.21           C
ATOM   1135  O   GLN A 141       7.816  18.843   2.239  1.00 32.71           O
ATOM   1136  CB  GLN A 141       8.284  19.888  -0.706  1.00 33.25           C
ATOM   1137  CG  GLN A 141       9.788  19.739  -0.890  1.00 36.97           C
ATOM   1138  CD  GLN A 141      10.355  20.710  -1.921  1.00 40.42           C
ATOM   1139  NE2 GLN A 141       9.574  20.999  -2.958  1.00 42.96           N
ATOM   1140  OE1 GLN A 141      11.479  21.193  -1.783  1.00 41.96           O
ATOM   1141  N   ILE A 142       5.777  19.248   1.376  1.00 31.96           N
ATOM   1142  CA  ILE A 142       5.093  19.367   2.670  1.00 30.79           C
ATOM   1143  C   ILE A 142       4.864  17.965   3.247  1.00 29.67           C
ATOM   1144  O   ILE A 142       5.012  17.740   4.451  1.00 29.80           O
ATOM   1145  CB  ILE A 142       3.740  20.131   2.541  1.00 30.84           C
ATOM   1146  CG1 ILE A 142       3.954  21.488   1.853  1.00 32.84           C
ATOM   1147  CG2 ILE A 142       3.086  20.319   3.909  1.00 28.44           C
ATOM   1148  CD1 ILE A 142       2.693  22.109   1.249  1.00 30.93           C
ATOM   1149  N   THR A 143       4.513  17.026   2.370  1.00 28.00           N
ATOM   1150  CA  THR A 143       4.410  15.618   2.727  1.00 26.81           C
ATOM   1151  C   THR A 143       5.758  15.091   3.218  1.00 26.20           C
ATOM   1152  O   THR A 143       5.832  14.403   4.245  1.00 23.84           O
ATOM   1153  CB  THR A 143       3.926  14.786   1.533  1.00 26.45           C
ATOM   1154  CG2 THR A 143       3.725  13.333   1.933  1.00 26.61           C
ATOM   1155  OG1 THR A 143       2.684  15.320   1.057  1.00 26.34           O
ATOM   1156  N   GLN A 144       6.816  15.439   2.486  1.00 26.44           N
ATOM   1157  CA  GLN A 144       8.172  15.012   2.823  1.00 27.66           C
ATOM   1158  C   GLN A 144       8.588  15.457   4.227  1.00 28.21           C
ATOM   1159  O   GLN A 144       9.184  14.675   4.968  1.00 28.03           O
ATOM   1160  CB  GLN A 144       9.177  15.516   1.785  1.00 27.66           C
ATOM   1161  CG  GLN A 144      10.487  14.742   1.798  1.00 28.56           C
ATOM   1162  CD  GLN A 144      11.486  15.216   0.761  1.00 29.05           C
ATOM   1163  NE2 GLN A 144      12.759  15.237   1.146  1.00 26.87           N
ATOM   1164  OE1 GLN A 144      11.127  15.547  -0.379  1.00 32.30           O
ATOM   1165  N   ARG A 145       8.262  16.703   4.585  1.00 28.44           N
ATOM   1166  CA  ARG A 145       8.556  17.228   5.922  1.00 28.59           C
ATOM   1167  C   ARG A 145       7.780  16.482   7.005  1.00 29.09           C
ATOM   1168  O   ARG A 145       8.320  16.204   8.075  1.00 29.28           O
ATOM   1169  CB  ARG A 145       8.273  18.729   6.000  1.00 28.50           C
ATOM   1170  CG  ARG A 145       9.278  19.588   5.240  1.00 29.35           C
ATOM   1171  CD  ARG A 145       8.817  21.041   5.109  1.00 28.87           C
ATOM   1172  NE  ARG A 145       9.323  21.636   3.874  1.00 30.78           N
ATOM   1173  CZ  ARG A 145       8.568  22.046   2.859  1.00 30.08           C
ATOM   1174  NH1 ARG A 145       7.244  21.965   2.924  1.00 33.13           N1+
ATOM   1175  NH2 ARG A 145       9.142  22.558   1.777  1.00 29.29           N
ATOM   1176  N   LYS A 146       6.518  16.159   6.717  1.00 29.64           N
ATOM   1177  CA  LYS A 146       5.685  15.377   7.631  1.00 29.55           C
ATOM   1178  C   LYS A 146       6.282  13.990   7.861  1.00 30.18           C
ATOM   1179  O   LYS A 146       6.360  13.520   9.000  1.00 30.49           O
ATOM   1180  CB  LYS A 146       4.266  15.231   7.082  1.00 29.38           C
ATOM   1181  CG  LYS A 146       3.366  16.441   7.251  1.00 29.95           C
ATOM   1182  CD  LYS A 146       2.001  16.152   6.639  1.00 31.71           C
ATOM   1183  CE  LYS A 146       0.970  17.209   7.008  1.00 33.49           C
ATOM   1184  NZ  LYS A 146      -0.387  16.845   6.503  1.00 30.61           N1+
ATOM   1185  N   TRP A 147       6.714  13.349   6.775  1.00 30.27           N
ATOM   1186  CA  TRP A 147       7.224  11.980   6.835  1.00 29.38           C
ATOM   1187  C   TRP A 147       8.621  11.871   7.449  1.00 29.53           C
ATOM   1188  O   TRP A 147       8.963  10.841   8.026  1.00 30.00           O
ATOM   1189  CB  TRP A 147       7.165  11.300   5.459  1.00 28.16           C
ATOM   1190  CG  TRP A 147       5.763  10.981   4.985  1.00 26.93           C
ATOM   1191  CD1 TRP A 147       4.593  11.214   5.656  1.00 26.36           C
ATOM   1192  CD2 TRP A 147       5.391  10.347   3.753  1.00 26.76           C
ATOM   1193  CE2 TRP A 147       3.983  10.242   3.744  1.00 25.55           C
ATOM   1194  CE3 TRP A 147       6.111   9.857   2.656  1.00 26.96           C
ATOM   1195  NE1 TRP A 147       3.522  10.785   4.913  1.00 24.13           N
ATOM   1196  CZ2 TRP A 147       3.284   9.671   2.678  1.00 27.64           C
ATOM   1197  CZ3 TRP A 147       5.417   9.285   1.602  1.00 23.54           C
ATOM   1198  CH2 TRP A 147       4.019   9.201   1.618  1.00 26.95           C
ATOM   1199  N   GLU A 148       9.414  12.933   7.333  1.00 29.65           N
ATOM   1200  CA  GLU A 148      10.746  12.974   7.940  1.00 30.37           C
ATOM   1201  C   GLU A 148      10.653  13.163   9.454  1.00 31.30           C
ATOM   1202  O   GLU A 148      11.311  12.452  10.221  1.00 31.74           O
ATOM   1203  CB  GLU A 148      11.594  14.084   7.313  1.00 30.20           C
ATOM   1204  CG  GLU A 148      12.035  13.802   5.881  1.00 28.53           C
ATOM   1205  CD  GLU A 148      12.699  14.994   5.209  1.00 31.12           C
ATOM   1206  OE1 GLU A 148      12.842  16.063   5.850  1.00 31.42           O
ATOM   1207  OE2 GLU A 148      13.086  14.859   4.027  1.00 34.30           O1-
ATOM   1208  N   ALA A 149       9.825  14.121   9.866  1.00 32.05           N
ATOM   1209  CA  ALA A 149       9.558  14.407  11.278  1.00 32.42           C
ATOM   1210  C   ALA A 149       9.052  13.184  12.035  1.00 33.15           C
ATOM   1211  O   ALA A 149       9.387  13.000  13.212  1.00 33.19           O
ATOM   1212  CB  ALA A 149       8.554  15.548  11.400  1.00 32.03           C
ATOM   1213  N   ALA A 150       8.251  12.364  11.347  1.00 33.17           N
ATOM   1214  CA  ALA A 150       7.651  11.149  11.912  1.00 33.32           C
ATOM   1215  C   ALA A 150       8.457   9.886  11.578  1.00 33.91           C
ATOM   1216  O   ALA A 150       8.061   8.771  11.953  1.00 32.73           O
ATOM   1217  CB  ALA A 150       6.201  11.001  11.442  1.00 32.92           C
ATOM   1218  N   ARG A 151       9.578  10.075  10.875  1.00 34.36           N
ATOM   1219  CA  ARG A 151      10.499   8.991  10.516  1.00 34.44           C
ATOM   1220  C   ARG A 151       9.773   7.798   9.887  1.00 33.78           C
ATOM   1221  O   ARG A 151       9.860   6.678  10.385  1.00 34.31           O
ATOM   1222  CB  ARG A 151      11.303   8.551  11.749  1.00 35.40           C
ATOM   1223  CG  ARG A 151      12.510   9.425  12.080  1.00 36.14           C
ATOM   1224  CD  ARG A 151      12.862   9.420  13.576  1.00 39.98           C
ATOM   1225  NE  ARG A 151      12.504   8.170  14.254  1.00 43.15           N
ATOM   1226  CZ  ARG A 151      11.758   8.093  15.356  1.00 42.82           C
ATOM   1227  NH1 ARG A 151      11.299   9.193  15.941  1.00 40.00           N1+
ATOM   1228  NH2 ARG A 151      11.487   6.906  15.885  1.00 42.45           N
ATOM   1229  N   VAL A 152       9.053   8.049   8.795  1.00 33.35           N
ATOM   1230  CA  VAL A 152       8.237   7.018   8.150  1.00 32.18           C
ATOM   1231  C   VAL A 152       9.073   6.077   7.288  1.00 30.66           C
ATOM   1232  O   VAL A 152       8.724   4.913   7.129  1.00 30.16           O
ATOM   1233  CB  VAL A 152       7.077   7.629   7.314  1.00 32.86           C
ATOM   1234  CG1 VAL A 152       6.097   6.544   6.869  1.00 32.51           C
ATOM   1235  CG2 VAL A 152       6.336   8.675   8.117  1.00 32.99           C
ATOM   1236  N   ALA A 153      10.175   6.588   6.741  1.00 30.39           N
ATOM   1237  CA  ALA A 153      11.122   5.777   5.962  1.00 29.67           C
ATOM   1238  C   ALA A 153      11.641   4.598   6.781  1.00 30.24           C
ATOM   1239  O   ALA A 153      12.063   3.581   6.223  1.00 31.55           O
ATOM   1240  CB  ALA A 153      12.282   6.631   5.471  1.00 28.25           C
ATOM   1241  N   GLU A 154      11.604   4.749   8.103  1.00 29.44           N
ATOM   1242  CA  GLU A 154      11.978   3.692   9.028  1.00 30.04           C
ATOM   1243  C   GLU A 154      10.989   2.537   8.954  1.00 30.04           C
ATOM   1244  O   GLU A 154      11.391   1.380   8.809  1.00 30.05           O
ATOM   1245  CB  GLU A 154      12.022   4.223  10.458  1.00 29.94           C
ATOM   1246  CG  GLU A 154      13.071   5.286  10.727  1.00 29.82           C
ATOM   1247  CD  GLU A 154      13.246   5.561  12.213  1.00 32.10           C
ATOM   1248  OE1 GLU A 154      12.284   5.330  12.986  1.00 33.23           O
ATOM   1249  OE2 GLU A 154      14.347   6.009  12.609  1.00 34.60           O1-
ATOM   1250  N   GLN A 155       9.700   2.865   9.062  1.00 30.56           N
ATOM   1251  CA  GLN A 155       8.613   1.888   8.981  1.00 30.85           C
ATOM   1252  C   GLN A 155       8.699   1.067   7.701  1.00 31.76           C
ATOM   1253  O   GLN A 155       8.567  -0.163   7.731  1.00 32.77           O
ATOM   1254  CB  GLN A 155       7.246   2.585   9.044  1.00 31.01           C
ATOM   1255  CG  GLN A 155       6.796   2.991  10.435  1.00 30.83           C
ATOM   1256  CD  GLN A 155       5.605   3.939  10.402  1.00 30.06           C
ATOM   1257  NE2 GLN A 155       5.881   5.230  10.549  1.00 30.08           N
ATOM   1258  OE1 GLN A 155       4.457   3.518  10.259  1.00 28.93           O
ATOM   1259  N   LEU A 156       8.920   1.756   6.583  1.00 31.64           N
ATOM   1260  CA  LEU A 156       8.987   1.122   5.266  1.00 31.41           C
ATOM   1261  C   LEU A 156      10.175   0.180   5.118  1.00 31.64           C
ATOM   1262  O   LEU A 156      10.013  -0.941   4.642  1.00 32.48           O
ATOM   1263  CB  LEU A 156       9.013   2.176   4.158  1.00 30.97           C
ATOM   1264  CG  LEU A 156       7.686   2.857   3.834  1.00 30.66           C
ATOM   1265  CD1 LEU A 156       7.974   4.149   3.134  1.00 32.06           C
ATOM   1266  CD2 LEU A 156       6.814   1.965   2.970  1.00 32.05           C
ATOM   1267  N   ARG A 157      11.360   0.638   5.519  1.00 30.69           N
ATOM   1268  CA  ARG A 157      12.553  -0.203   5.485  1.00 30.66           C
ATOM   1269  C   ARG A 157      12.386  -1.468   6.334  1.00 30.59           C
ATOM   1270  O   ARG A 157      12.899  -2.531   5.984  1.00 31.06           O
ATOM   1271  CB  ARG A 157      13.791   0.572   5.935  1.00 30.72           C
ATOM   1272  CG  ARG A 157      15.065  -0.191   5.678  1.00 28.09           C
ATOM   1273  CD  ARG A 157      16.265   0.413   6.347  1.00 29.67           C
ATOM   1274  NE  ARG A 157      17.355  -0.560   6.348  1.00 33.46           N
ATOM   1275  CZ  ARG A 157      18.643  -0.264   6.477  1.00 34.24           C
ATOM   1276  NH1 ARG A 157      19.036   0.999   6.619  1.00 32.02           N1+
ATOM   1277  NH2 ARG A 157      19.542  -1.240   6.458  1.00 33.11           N
ATOM   1278  N   ALA A 158      11.666  -1.341   7.446  1.00 29.85           N
ATOM   1279  CA  ALA A 158      11.355  -2.475   8.310  1.00 28.97           C
ATOM   1280  C   ALA A 158      10.390  -3.426   7.614  1.00 27.92           C
ATOM   1281  O   ALA A 158      10.373  -4.625   7.892  1.00 28.52           O
ATOM   1282  CB  ALA A 158      10.770  -1.993   9.623  1.00 28.42           C
ATOM   1283  N   TYR A 159       9.590  -2.882   6.708  1.00 26.19           N
ATOM   1284  CA  TYR A 159       8.663  -3.688   5.935  1.00 24.57           C
ATOM   1285  C   TYR A 159       9.340  -4.306   4.702  1.00 24.11           C
ATOM   1286  O   TYR A 159       9.330  -5.525   4.531  1.00 23.48           O
ATOM   1287  CB  TYR A 159       7.438  -2.859   5.542  1.00 23.26           C
ATOM   1288  CG  TYR A 159       6.523  -3.571   4.586  1.00 22.23           C
ATOM   1289  CD1 TYR A 159       5.577  -4.488   5.043  1.00 21.31           C
ATOM   1290  CD2 TYR A 159       6.610  -3.341   3.217  1.00 20.62           C
ATOM   1291  CE1 TYR A 159       4.742  -5.147   4.158  1.00 21.62           C
ATOM   1292  CE2 TYR A 159       5.787  -3.997   2.327  1.00 18.56           C
ATOM   1293  CZ  TYR A 159       4.857  -4.897   2.800  1.00 21.36           C
ATOM   1294  OH  TYR A 159       4.042  -5.543   1.902  1.00 24.38           O
ATOM   1295  N   LEU A 160       9.932  -3.463   3.859  1.00 23.55           N
ATOM   1296  CA  LEU A 160      10.526  -3.908   2.593  1.00 24.89           C
ATOM   1297  C   LEU A 160      11.678  -4.912   2.749  1.00 26.03           C
ATOM   1298  O   LEU A 160      11.772  -5.875   1.984  1.00 26.19           O
ATOM   1299  CB  LEU A 160      10.957  -2.708   1.741  1.00 24.50           C
ATOM   1300  CG  LEU A 160       9.860  -1.699   1.370  1.00 23.45           C
ATOM   1301  CD1 LEU A 160      10.468  -0.423   0.818  1.00 21.04           C
ATOM   1302  CD2 LEU A 160       8.851  -2.293   0.386  1.00 20.85           C
ATOM   1303  N   GLU A 161      12.535  -4.694   3.743  1.00 26.43           N
ATOM   1304  CA  GLU A 161      13.644  -5.609   4.022  1.00 26.99           C
ATOM   1305  C   GLU A 161      13.280  -6.702   5.037  1.00 27.05           C
ATOM   1306  O   GLU A 161      14.059  -7.641   5.252  1.00 26.61           O
ATOM   1307  CB  GLU A 161      14.873  -4.832   4.503  1.00 27.17           C
ATOM   1308  CG  GLU A 161      15.501  -3.942   3.441  1.00 32.01           C
ATOM   1309  CD  GLU A 161      16.844  -3.384   3.859  1.00 36.68           C
ATOM   1310  OE1 GLU A 161      16.957  -2.886   4.998  1.00 39.26           O
ATOM   1311  OE2 GLU A 161      17.792  -3.437   3.044  1.00 43.10           O1-
ATOM   1312  N   GLY A 162      12.102  -6.579   5.649  1.00 26.92           N
ATOM   1313  CA  GLY A 162      11.695  -7.460   6.741  1.00 25.60           C
ATOM   1314  C   GLY A 162      10.436  -8.250   6.450  1.00 25.61           C
ATOM   1315  O   GLY A 162      10.513  -9.360   5.924  1.00 25.21           O
ATOM   1316  N   LEU A 163       9.280  -7.677   6.795  1.00 25.18           N
ATOM   1317  CA  LEU A 163       7.993  -8.363   6.658  1.00 26.11           C
ATOM   1318  C   LEU A 163       7.673  -8.726   5.215  1.00 26.81           C
ATOM   1319  O   LEU A 163       7.178  -9.819   4.938  1.00 27.75           O
ATOM   1320  CB  LEU A 163       6.853  -7.507   7.213  1.00 26.97           C
ATOM   1321  CG  LEU A 163       6.693  -7.320   8.720  1.00 27.89           C
ATOM   1322  CD1 LEU A 163       5.685  -6.218   8.981  1.00 24.82           C
ATOM   1323  CD2 LEU A 163       6.277  -8.619   9.405  1.00 26.50           C
ATOM   1324  N   CYS A 164       7.943  -7.797   4.305  1.00 26.22           N
ATOM   1325  CA  CYS A 164       7.709  -8.019   2.885  1.00 26.56           C
ATOM   1326  C   CYS A 164       8.412  -9.286   2.385  1.00 26.59           C
ATOM   1327  O   CYS A 164       7.810 -10.084   1.660  1.00 27.43           O
ATOM   1328  CB  CYS A 164       8.156  -6.798   2.083  1.00 26.40           C
ATOM   1329  SG  CYS A 164       7.629  -6.778   0.372  1.00 25.62           S
ATOM   1330  N   VAL A 165       9.673  -9.467   2.787  1.00 24.89           N
ATOM   1331  CA  VAL A 165      10.481 -10.606   2.354  1.00 24.22           C
ATOM   1332  C   VAL A 165       9.950 -11.907   2.952  1.00 24.86           C
ATOM   1333  O   VAL A 165       9.669 -12.859   2.227  1.00 25.15           O
ATOM   1334  CB  VAL A 165      11.989 -10.425   2.708  1.00 24.13           C
ATOM   1335  CG1 VAL A 165      12.784 -11.676   2.379  1.00 22.47           C
ATOM   1336  CG2 VAL A 165      12.586  -9.217   1.984  1.00 23.96           C
ATOM   1337  N   GLU A 166       9.808 -11.925   4.275  1.00 25.66           N
ATOM   1338  CA  GLU A 166       9.379 -13.107   5.014  1.00 25.96           C
ATOM   1339  C   GLU A 166       8.031 -13.640   4.550  1.00 26.45           C
ATOM   1340  O   GLU A 166       7.880 -14.847   4.345  1.00 27.89           O
ATOM   1341  CB  GLU A 166       9.339 -12.814   6.518  1.00 26.58           C
ATOM   1342  CG  GLU A 166       8.841 -13.981   7.360  1.00 28.42           C
ATOM   1343  CD  GLU A 166       9.380 -13.954   8.771  1.00 32.99           C
ATOM   1344  OE1 GLU A 166      10.286 -14.763   9.061  1.00 35.69           O
ATOM   1345  OE2 GLU A 166       8.907 -13.125   9.584  1.00 34.64           O1-
ATOM   1346  N   TRP A 167       7.059 -12.741   4.398  1.00 26.19           N
ATOM   1347  CA  TRP A 167       5.709 -13.115   3.978  1.00 24.95           C
ATOM   1348  C   TRP A 167       5.652 -13.595   2.533  1.00 25.64           C
ATOM   1349  O   TRP A 167       4.953 -14.564   2.224  1.00 26.04           O
ATOM   1350  CB  TRP A 167       4.721 -11.969   4.209  1.00 23.39           C
ATOM   1351  CG  TRP A 167       4.162 -11.992   5.585  1.00 21.24           C
ATOM   1352  CD1 TRP A 167       4.499 -11.172   6.621  1.00 20.33           C
ATOM   1353  CD2 TRP A 167       3.187 -12.909   6.098  1.00 17.54           C
ATOM   1354  CE2 TRP A 167       2.976 -12.576   7.453  1.00 16.25           C
ATOM   1355  CE3 TRP A 167       2.477 -13.979   5.545  1.00 16.45           C
ATOM   1356  NE1 TRP A 167       3.786 -11.513   7.745  1.00 20.19           N
ATOM   1357  CZ2 TRP A 167       2.078 -13.271   8.265  1.00 19.11           C
ATOM   1358  CZ3 TRP A 167       1.579 -14.668   6.347  1.00 23.14           C
ATOM   1359  CH2 TRP A 167       1.387 -14.308   7.700  1.00 22.06           C
ATOM   1360  N   LEU A 168       6.388 -12.914   1.657  1.00 26.10           N
ATOM   1361  CA  LEU A 168       6.495 -13.321   0.261  1.00 26.82           C
ATOM   1362  C   LEU A 168       6.996 -14.760   0.177  1.00 26.95           C
ATOM   1363  O   LEU A 168       6.418 -15.584  -0.542  1.00 27.96           O
ATOM   1364  CB  LEU A 168       7.416 -12.376  -0.519  1.00 26.95           C
ATOM   1365  CG  LEU A 168       7.641 -12.594  -2.021  1.00 26.18           C
ATOM   1366  CD1 LEU A 168       6.350 -12.508  -2.818  1.00 24.72           C
ATOM   1367  CD2 LEU A 168       8.645 -11.584  -2.546  1.00 26.26           C
ATOM   1368  N   ARG A 169       8.046 -15.064   0.935  1.00 26.09           N
ATOM   1369  CA  ARG A 169       8.593 -16.421   0.986  1.00 25.86           C
ATOM   1370  C   ARG A 169       7.538 -17.430   1.431  1.00 25.99           C
ATOM   1371  O   ARG A 169       7.392 -18.485   0.817  1.00 27.34           O
ATOM   1372  CB  ARG A 169       9.843 -16.483   1.867  1.00 24.66           C
ATOM   1373  CG  ARG A 169      11.094 -16.048   1.135  1.00 23.79           C
ATOM   1374  CD  ARG A 169      12.220 -15.651   2.076  1.00 26.03           C
ATOM   1375  NE  ARG A 169      13.387 -15.189   1.322  1.00 23.14           N
ATOM   1376  CZ  ARG A 169      14.492 -14.675   1.855  1.00 20.58           C
ATOM   1377  NH1 ARG A 169      14.620 -14.538   3.174  1.00 16.78           N1+
ATOM   1378  NH2 ARG A 169      15.476 -14.291   1.056  1.00 19.83           N
ATOM   1379  N   ARG A 170       6.787 -17.087   2.474  1.00 25.67           N
ATOM   1380  CA  ARG A 170       5.677 -17.920   2.927  1.00 25.38           C
ATOM   1381  C   ARG A 170       4.658 -18.180   1.802  1.00 24.53           C
ATOM   1382  O   ARG A 170       4.282 -19.331   1.560  1.00 24.60           O
ATOM   1383  CB  ARG A 170       5.023 -17.328   4.188  1.00 24.99           C
ATOM   1384  CG  ARG A 170       3.610 -17.825   4.448  1.00 27.80           C
ATOM   1385  CD  ARG A 170       3.251 -17.941   5.927  1.00 30.22           C
ATOM   1386  NE  ARG A 170       2.173 -18.922   6.081  1.00 34.10           N
ATOM   1387  CZ  ARG A 170       1.315 -18.985   7.096  1.00 33.14           C
ATOM   1388  NH1 ARG A 170       1.377 -18.122   8.104  1.00 34.30           N1+
ATOM   1389  NH2 ARG A 170       0.380 -19.925   7.093  1.00 33.56           N
ATOM   1390  N   TYR A 171       4.237 -17.122   1.108  1.00 23.50           N
ATOM   1391  CA  TYR A 171       3.274 -17.248   0.005  1.00 22.45           C
ATOM   1392  C   TYR A 171       3.823 -18.111  -1.122  1.00 23.53           C
ATOM   1393  O   TYR A 171       3.102 -18.937  -1.687  1.00 22.90           O
ATOM   1394  CB  TYR A 171       2.886 -15.879  -0.565  1.00 21.08           C
ATOM   1395  CG  TYR A 171       2.224 -14.919   0.405  1.00 20.96           C
ATOM   1396  CD1 TYR A 171       1.443 -15.377   1.472  1.00 15.88           C
ATOM   1397  CD2 TYR A 171       2.354 -13.542   0.231  1.00 19.60           C
ATOM   1398  CE1 TYR A 171       0.843 -14.493   2.346  1.00 14.28           C
ATOM   1399  CE2 TYR A 171       1.745 -12.651   1.099  1.00 18.24           C
ATOM   1400  CZ  TYR A 171       0.994 -13.128   2.153  1.00 17.49           C
ATOM   1401  OH  TYR A 171       0.394 -12.229   3.009  1.00 16.70           O
ATOM   1402  N   LEU A 172       5.102 -17.904  -1.446  1.00 24.44           N
ATOM   1403  CA  LEU A 172       5.774 -18.667  -2.493  1.00 24.99           C
ATOM   1404  C   LEU A 172       5.836 -20.165  -2.201  1.00 26.97           C
ATOM   1405  O   LEU A 172       5.735 -20.967  -3.120  1.00 28.37           O
ATOM   1406  CB  LEU A 172       7.183 -18.136  -2.734  1.00 24.32           C
ATOM   1407  CG  LEU A 172       7.349 -16.759  -3.368  1.00 22.68           C
ATOM   1408  CD1 LEU A 172       8.806 -16.359  -3.306  1.00 18.48           C
ATOM   1409  CD2 LEU A 172       6.839 -16.731  -4.798  1.00 20.30           C
ATOM   1410  N   GLU A 173       6.008 -20.545  -0.937  1.00 28.96           N
ATOM   1411  CA  GLU A 173       5.963 -21.963  -0.574  1.00 31.68           C
ATOM   1412  C   GLU A 173       4.531 -22.497  -0.549  1.00 33.10           C
ATOM   1413  O   GLU A 173       4.238 -23.517  -1.179  1.00 33.89           O
ATOM   1414  CB  GLU A 173       6.661 -22.229   0.762  1.00 31.42           C
ATOM   1415  CG  GLU A 173       6.622 -23.689   1.213  1.00 34.88           C
ATOM   1416  CD  GLU A 173       7.338 -24.643   0.261  1.00 39.69           C
ATOM   1417  OE1 GLU A 173       8.379 -24.260  -0.319  1.00 42.76           O
ATOM   1418  OE2 GLU A 173       6.864 -25.787   0.106  1.00 40.89           O1-
ATOM   1419  N   ASN A 174       3.648 -21.802   0.168  1.00 33.74           N
ATOM   1420  CA  ASN A 174       2.247 -22.209   0.289  1.00 34.57           C
ATOM   1421  C   ASN A 174       1.549 -22.532  -1.036  1.00 34.95           C
ATOM   1422  O   ASN A 174       0.687 -23.415  -1.083  1.00 34.56           O
ATOM   1423  CB  ASN A 174       1.446 -21.170   1.077  1.00 34.72           C
ATOM   1424  CG  ASN A 174       1.666 -21.273   2.577  1.00 37.05           C
ATOM   1425  ND2 ASN A 174       1.429 -20.170   3.287  1.00 33.61           N
ATOM   1426  OD1 ASN A 174       2.039 -22.333   3.096  1.00 39.92           O
ATOM   1427  N   GLY A 175       1.930 -21.823  -2.099  1.00 35.19           N
ATOM   1428  CA  GLY A 175       1.365 -22.051  -3.427  1.00 35.63           C
ATOM   1429  C   GLY A 175       2.406 -22.218  -4.516  1.00 36.52           C
ATOM   1430  O   GLY A 175       2.204 -21.765  -5.645  1.00 36.38           O
ATOM   1431  N   LYS A 176       3.508 -22.888  -4.173  1.00 37.37           N
ATOM   1432  CA  LYS A 176       4.651 -23.094  -5.079  1.00 37.74           C
ATOM   1433  C   LYS A 176       4.275 -23.695  -6.431  1.00 37.16           C
ATOM   1434  O   LYS A 176       4.894 -23.395  -7.451  1.00 36.71           O
ATOM   1435  CB  LYS A 176       5.727 -23.963  -4.408  1.00 38.18           C
ATOM   1436  CG  LYS A 176       5.240 -25.319  -3.880  1.00 39.46           C
ATOM   1437  CD  LYS A 176       6.376 -26.108  -3.226  1.00 38.74           C
ATOM   1438  CE  LYS A 176       5.881 -27.424  -2.630  1.00 41.99           C
ATOM   1439  NZ  LYS A 176       5.013 -27.234  -1.427  1.00 41.38           N1+
ATOM   1440  N   GLU A 177       3.250 -24.534  -6.431  1.00 36.75           N
ATOM   1441  CA  GLU A 177       2.887 -25.270  -7.627  1.00 36.10           C
ATOM   1442  C   GLU A 177       2.041 -24.452  -8.617  1.00 35.21           C
ATOM   1443  O   GLU A 177       1.744 -24.915  -9.720  1.00 35.65           O
ATOM   1444  CB  GLU A 177       2.258 -26.615  -7.235  1.00 36.55           C
ATOM   1445  CG  GLU A 177       3.324 -27.622  -6.775  1.00 38.48           C
ATOM   1446  CD  GLU A 177       2.835 -28.632  -5.747  1.00 41.90           C
ATOM   1447  OE1 GLU A 177       1.941 -28.306  -4.935  1.00 42.13           O
ATOM   1448  OE2 GLU A 177       3.372 -29.762  -5.739  1.00 45.06           O1-
ATOM   1449  N   THR A 178       1.678 -23.228  -8.228  1.00 33.36           N
ATOM   1450  CA  THR A 178       1.062 -22.275  -9.160  1.00 31.41           C
ATOM   1451  C   THR A 178       1.860 -20.973  -9.266  1.00 30.83           C
ATOM   1452  O   THR A 178       1.730 -20.239 -10.243  1.00 30.89           O
ATOM   1453  CB  THR A 178      -0.418 -21.963  -8.816  1.00 31.03           C
ATOM   1454  CG2 THR A 178      -1.257 -23.239  -8.803  1.00 30.77           C
ATOM   1455  OG1 THR A 178      -0.499 -21.330  -7.537  1.00 29.53           O
ATOM   1456  N   LEU A 179       2.704 -20.703  -8.274  1.00 30.40           N
ATOM   1457  CA  LEU A 179       3.491 -19.474  -8.260  1.00 30.33           C
ATOM   1458  C   LEU A 179       4.930 -19.647  -8.743  1.00 30.78           C
ATOM   1459  O   LEU A 179       5.555 -18.671  -9.170  1.00 30.58           O
ATOM   1460  CB  LEU A 179       3.486 -18.841  -6.865  1.00 29.80           C
ATOM   1461  CG  LEU A 179       2.144 -18.462  -6.237  1.00 30.18           C
ATOM   1462  CD1 LEU A 179       2.372 -17.976  -4.813  1.00 32.94           C
ATOM   1463  CD2 LEU A 179       1.400 -17.412  -7.059  1.00 29.50           C
ATOM   1464  N   GLN A 180       5.457 -20.870  -8.666  1.00 31.33           N
ATOM   1465  CA  GLN A 180       6.865 -21.122  -9.011  1.00 33.69           C
ATOM   1466  C   GLN A 180       7.053 -21.997 -10.253  1.00 34.94           C
ATOM   1467  O   GLN A 180       8.132 -22.014 -10.852  1.00 35.01           O
ATOM   1468  CB  GLN A 180       7.638 -21.697  -7.817  1.00 32.59           C
ATOM   1469  CG  GLN A 180       7.763 -20.735  -6.638  1.00 33.62           C
ATOM   1470  CD  GLN A 180       8.668 -21.261  -5.536  1.00 35.06           C
ATOM   1471  NE2 GLN A 180       8.109 -21.427  -4.341  1.00 36.83           N
ATOM   1472  OE1 GLN A 180       9.858 -21.500  -5.749  1.00 36.97           O
TER    1473      GLN A 180
ATOM   1474  N   GLU C   1       1.037  -9.700   2.999  1.00 23.48           N
ATOM   1475  CA  GLU C   1       0.949  -8.914   4.259  1.00 22.61           C
ATOM   1476  C   GLU C   1       1.178  -7.444   3.965  1.00 22.03           C
ATOM   1477  O   GLU C   1       2.271  -7.054   3.556  1.00 20.56           O
ATOM   1478  CB  GLU C   1       1.967  -9.418   5.289  1.00 23.16           C
ATOM   1479  CG  GLU C   1       1.916  -8.726   6.660  1.00 23.38           C
ATOM   1480  CD  GLU C   1       0.603  -8.947   7.415  1.00 25.84           C
ATOM   1481  OE1 GLU C   1      -0.175  -9.856   7.039  1.00 23.03           O
ATOM   1482  OE2 GLU C   1       0.353  -8.203   8.391  1.00 22.09           O1-
ATOM   1483  N   PRO C   2       0.148  -6.618   4.194  1.00 22.85           N
ATOM   1484  CA  PRO C   2       0.195  -5.194   3.877  1.00 23.15           C
ATOM   1485  C   PRO C   2       1.147  -4.413   4.777  1.00 24.18           C
ATOM   1486  O   PRO C   2       1.551  -4.903   5.833  1.00 24.30           O
ATOM   1487  CB  PRO C   2      -1.245  -4.737   4.115  1.00 22.64           C
ATOM   1488  CG  PRO C   2      -1.790  -5.705   5.102  1.00 22.94           C
ATOM   1489  CD  PRO C   2      -1.140  -7.013   4.792  1.00 22.97           C
ATOM   1490  N   LEU C   3       1.500  -3.211   4.332  1.00 25.28           N
ATOM   1491  CA  LEU C   3       2.291  -2.271   5.108  1.00 26.09           C
ATOM   1492  C   LEU C   3       1.451  -1.650   6.215  1.00 26.55           C
ATOM   1493  O   LEU C   3       0.480  -0.943   5.932  1.00 27.90           O
ATOM   1494  CB  LEU C   3       2.829  -1.169   4.192  1.00 26.99           C
ATOM   1495  CG  LEU C   3       3.316   0.155   4.797  1.00 28.10           C
ATOM   1496  CD1 LEU C   3       4.736   0.016   5.348  1.00 25.55           C
ATOM   1497  CD2 LEU C   3       3.236   1.254   3.750  1.00 24.40           C
ATOM   1498  N   PRO C   4       1.816  -1.909   7.482  1.00 26.66           N
ATOM   1499  CA  PRO C   4       1.173  -1.239   8.612  1.00 26.67           C
ATOM   1500  C   PRO C   4       1.568   0.230   8.661  1.00 26.15           C
ATOM   1501  O   PRO C   4       2.741   0.554   8.454  1.00 26.32           O
ATOM   1502  CB  PRO C   4       1.747  -1.972   9.835  1.00 26.48           C
ATOM   1503  CG  PRO C   4       2.393  -3.213   9.297  1.00 26.89           C
ATOM   1504  CD  PRO C   4       2.849  -2.857   7.928  1.00 27.01           C
ATOM   1505  N   GLN C   5       0.607   1.111   8.931  1.00 26.00           N
ATOM   1506  CA  GLN C   5       0.901   2.550   9.015  1.00 25.92           C
ATOM   1507  C   GLN C   5       0.209   3.354  10.128  1.00 26.22           C
ATOM   1508  O   GLN C   5       0.834   3.659  11.140  1.00 27.36           O
ATOM   1509  CB  GLN C   5       0.757   3.235   7.649  1.00 26.36           C
ATOM   1510  CG  GLN C   5      -0.283   2.654   6.719  1.00 24.57           C
ATOM   1511  CD  GLN C   5      -0.133   3.174   5.304  1.00 25.72           C
ATOM   1512  NE2 GLN C   5      -0.097   2.264   4.329  1.00 25.51           N
ATOM   1513  OE1 GLN C   5      -0.058   4.379   5.084  1.00 27.99           O
ATOM   1514  N   GLY C   6      -1.055   3.718   9.949  1.00 27.05           N
ATOM   1515  CA  GLY C   6      -1.736   4.603  10.917  1.00 28.24           C
ATOM   1516  C   GLY C   6      -1.389   6.081  10.750  1.00 27.34           C
ATOM   1517  O   GLY C   6      -0.667   6.440   9.828  1.00 25.79           O
ATOM   1518  N   GLN C   7      -1.894   6.936  11.645  1.00 28.21           N
ATOM   1519  CA  GLN C   7      -1.700   8.405  11.551  1.00 28.19           C
ATOM   1520  C   GLN C   7      -0.224   8.758  11.415  1.00 29.90           C
ATOM   1521  O   GLN C   7       0.635   7.934  11.709  1.00 31.32           O
ATOM   1522  CB  GLN C   7      -2.251   9.121  12.788  1.00 28.05           C
ATOM   1523  CG  GLN C   7      -3.610   8.648  13.275  1.00 27.39           C
ATOM   1524  CD  GLN C   7      -3.929   9.153  14.667  1.00 26.68           C
ATOM   1525  NE2 GLN C   7      -4.897  10.055  14.762  1.00 27.33           N
ATOM   1526  OE1 GLN C   7      -3.321   8.729  15.647  1.00 23.09           O
ATOM   1527  N   LEU C   8       0.077   9.980  10.986  1.00 29.70           N
ATOM   1528  CA  LEU C   8       1.468  10.403  10.779  1.00 29.20           C
ATOM   1529  C   LEU C   8       2.048   9.840   9.483  1.00 28.77           C
ATOM   1530  O   LEU C   8       3.137  10.238   9.062  1.00 28.77           O
ATOM   1531  CB  LEU C   8       2.367  10.018  11.965  1.00 29.05           C
ATOM   1532  CG  LEU C   8       2.229  10.735  13.311  1.00 31.68           C
ATOM   1533  CD1 LEU C   8       3.236  10.163  14.303  1.00 30.03           C
ATOM   1534  CD2 LEU C   8       2.417  12.247  13.178  1.00 32.24           C
ATOM   1535  N   THR C   9       1.325   8.913   8.861  1.00 28.36           N
ATOM   1536  CA  THR C   9       1.671   8.426   7.529  1.00 28.88           C
ATOM   1537  C   THR C   9       0.805   9.128   6.482  1.00 28.63           C
ATOM   1538  O   THR C   9       0.830   8.780   5.298  1.00 29.05           O
ATOM   1539  CB  THR C   9       1.500   6.893   7.416  1.00 28.26           C
ATOM   1540  CG2 THR C   9       0.037   6.520   7.275  1.00 29.71           C
ATOM   1541  OG1 THR C   9       2.199   6.422   6.261  1.00 34.11           O
ATOM   1542  N   ALA C  10       0.033  10.112   6.933  1.00 28.97           N
ATOM   1543  CA  ALA C  10      -0.842  10.876   6.055  1.00 28.83           C
ATOM   1544  C   ALA C  10      -0.007  11.743   5.137  1.00 28.88           C
ATOM   1545  O   ALA C  10       1.101  12.147   5.491  1.00 28.88           O
ATOM   1546  CB  ALA C  10      -1.799  11.726   6.860  1.00 27.85           C
ATOM   1547  N   TYR C  11      -0.547  12.004   3.950  1.00 28.98           N
ATOM   1548  CA  TYR C  11       0.103  12.846   2.959  1.00 29.04           C
ATOM   1549  C   TYR C  11       0.097  14.297   3.413  1.00 29.47           C
ATOM   1550  O   TYR C  11      -0.856  14.761   4.040  1.00 28.85           O
ATOM   1551  CB  TYR C  11      -0.589  12.695   1.604  1.00 28.92           C
ATOM   1552  CG  TYR C  11      -0.295  11.372   0.920  1.00 29.82           C
ATOM   1553  CD1 TYR C  11      -0.783  10.170   1.438  1.00 29.13           C
ATOM   1554  CD2 TYR C  11       0.476  11.318  -0.242  1.00 26.06           C
ATOM   1555  CE1 TYR C  11      -0.511   8.956   0.819  1.00 26.78           C
ATOM   1556  CE2 TYR C  11       0.750  10.105  -0.863  1.00 25.72           C
ATOM   1557  CZ  TYR C  11       0.250   8.930  -0.327  1.00 25.83           C
ATOM   1558  OH  TYR C  11       0.510   7.723  -0.931  1.00 28.80           O
ATOM   1559  OXT TYR C  11       1.060  15.029   3.186  1.00 30.62           O1-
END
"""
    content_2 = \
    """ATOM      1  N   ILE A   1      18.718  -7.104 -11.642  1.00134.76           N
ATOM      2  CA  ILE A   1      17.932  -6.554 -12.769  1.00134.76           C
ATOM      3  C   ILE A   1      16.748  -7.425 -13.039  1.00134.76           C
ATOM      4  O   ILE A   1      15.633  -6.936 -13.208  1.00134.76           O
ATOM      5  CB  ILE A   1      18.813  -6.437 -13.993  1.00134.76           C
ATOM      6  CG1 ILE A   1      18.144  -5.603 -15.104  1.00134.76           C
ATOM      7  CG2 ILE A   1      19.249  -7.843 -14.432  1.00134.76           C
ATOM      8  CD1 ILE A   1      16.894  -6.230 -15.722  1.00134.76           C
ATOM      9  N   LYS A   2      16.960  -8.753 -13.061  1.00203.02           N
ATOM     10  CA  LYS A   2      15.874  -9.643 -13.327  1.00203.02           C
ATOM     11  C   LYS A   2      14.934  -9.551 -12.178  1.00203.02           C
ATOM     12  O   LYS A   2      15.341  -9.681 -11.024  1.00203.02           O
ATOM     13  CB  LYS A   2      16.309 -11.114 -13.444  1.00203.02           C
ATOM     14  CG  LYS A   2      15.157 -12.078 -13.734  1.00203.02           C
ATOM     15  CD  LYS A   2      14.578 -11.940 -15.142  1.00203.02           C
ATOM     16  CE  LYS A   2      13.427 -12.908 -15.426  1.00203.02           C
ATOM     17  NZ  LYS A   2      12.935 -12.715 -16.807  1.00203.02           N1+
ATOM     18  N   GLU A   3      13.639  -9.326 -12.469  1.00 98.31           N
ATOM     19  CA  GLU A   3      12.707  -9.211 -11.391  1.00 98.31           C
ATOM     20  C   GLU A   3      11.852 -10.439 -11.368  1.00 98.31           C
ATOM     21  O   GLU A   3      11.127 -10.730 -12.318  1.00 98.31           O
ATOM     22  CB  GLU A   3      11.780  -7.989 -11.500  1.00 98.31           C
ATOM     23  CG  GLU A   3      12.540  -6.668 -11.372  1.00 98.31           C
ATOM     24  CD  GLU A   3      11.539  -5.527 -11.474  1.00 98.31           C
ATOM     25  OE1 GLU A   3      10.314  -5.813 -11.428  1.00 98.31           O
ATOM     26  OE2 GLU A   3      11.988  -4.357 -11.595  1.00 98.31           O1-
ATOM     27  N   GLU A   4      11.985 -11.235 -10.288  1.00 69.01           N
ATOM     28  CA  GLU A   4      11.242 -12.452 -10.136  1.00 69.01           C
ATOM     29  C   GLU A   4       9.800 -12.257  -9.741  1.00 69.01           C
ATOM     30  O   GLU A   4       8.911 -12.806 -10.392  1.00 69.01           O
ATOM     31  CB  GLU A   4      11.881 -13.402  -9.108  1.00 69.01           C
ATOM     32  CG  GLU A   4      11.215 -14.777  -9.054  1.00 69.01           C
ATOM     33  CD  GLU A   4      11.953 -15.621  -8.025  1.00 69.01           C
ATOM     34  OE1 GLU A   4      13.140 -15.962  -8.281  1.00 69.01           O
ATOM     35  OE2 GLU A   4      11.342 -15.935  -6.967  1.00 69.01           O1-
ATOM     36  N   HIS A   5       9.507 -11.461  -8.683  1.00 76.81           N
ATOM     37  CA  HIS A   5       8.128 -11.411  -8.253  1.00 76.81           C
ATOM     38  C   HIS A   5       7.822 -10.105  -7.579  1.00 76.81           C
ATOM     39  O   HIS A   5       8.724  -9.375  -7.171  1.00 76.81           O
ATOM     40  CB  HIS A   5       7.756 -12.499  -7.225  1.00 76.81           C
ATOM     41  CG  HIS A   5       7.915 -13.910  -7.713  1.00 76.81           C
ATOM     42  CD2 HIS A   5       8.784 -14.879  -7.315  1.00 76.81           C
ATOM     43  ND1 HIS A   5       7.138 -14.483  -8.696  1.00 76.81           N
ATOM     44  CE1 HIS A   5       7.572 -15.761  -8.842  1.00 76.81           C
ATOM     45  NE2 HIS A   5       8.570 -16.047  -8.027  1.00 76.81           N
ATOM     46  N   VAL A   6       6.512  -9.771  -7.457  1.00113.97           N
ATOM     47  CA  VAL A   6       6.156  -8.560  -6.771  1.00113.97           C
ATOM     48  C   VAL A   6       4.849  -8.737  -6.050  1.00113.97           C
ATOM     49  O   VAL A   6       3.902  -9.325  -6.570  1.00113.97           O
ATOM     50  CB  VAL A   6       6.012  -7.367  -7.669  1.00113.97           C
ATOM     51  CG1 VAL A   6       4.680  -7.459  -8.427  1.00113.97           C
ATOM     52  CG2 VAL A   6       6.159  -6.103  -6.812  1.00113.97           C
ATOM     53  N   ILE A   7       4.773  -8.208  -4.809  1.00 50.14           N
ATOM     54  CA  ILE A   7       3.565  -8.264  -4.033  1.00 50.14           C
ATOM     55  C   ILE A   7       3.183  -6.844  -3.761  1.00 50.14           C
ATOM     56  O   ILE A   7       4.009  -6.051  -3.311  1.00 50.14           O
ATOM     57  CB  ILE A   7       3.742  -8.940  -2.702  1.00 50.14           C
ATOM     58  CG1 ILE A   7       4.202 -10.393  -2.899  1.00 50.14           C
ATOM     59  CG2 ILE A   7       2.420  -8.819  -1.925  1.00 50.14           C
ATOM     60  CD1 ILE A   7       4.671 -11.065  -1.611  1.00 50.14           C
ATOM     61  N   ILE A   8       1.916  -6.477  -4.043  1.00 44.11           N
ATOM     62  CA  ILE A   8       1.543  -5.105  -3.850  1.00 44.11           C
ATOM     63  C   ILE A   8       0.337  -5.019  -2.969  1.00 44.11           C
ATOM     64  O   ILE A   8      -0.647  -5.735  -3.152  1.00 44.11           O
ATOM     65  CB  ILE A   8       1.202  -4.409  -5.138  1.00 44.11           C
ATOM     66  CG1 ILE A   8       2.422  -4.401  -6.074  1.00 44.11           C
ATOM     67  CG2 ILE A   8       0.674  -3.005  -4.801  1.00 44.11           C
ATOM     68  CD1 ILE A   8       2.092  -3.969  -7.502  1.00 44.11           C
ATOM     69  N   GLN A   9       0.406  -4.122  -1.966  1.00 57.12           N
ATOM     70  CA  GLN A   9      -0.709  -3.843  -1.110  1.00 57.12           C
ATOM     71  C   GLN A   9      -1.278  -2.580  -1.674  1.00 57.12           C
ATOM     72  O   GLN A   9      -0.617  -1.541  -1.656  1.00 57.12           O
ATOM     73  CB  GLN A   9      -0.287  -3.540   0.338  1.00 57.12           C
ATOM     74  CG  GLN A   9      -1.451  -3.237   1.281  1.00 57.12           C
ATOM     75  CD  GLN A   9      -0.861  -2.866   2.635  1.00 57.12           C
ATOM     76  NE2 GLN A   9      -1.340  -1.733   3.216  1.00 57.12           N
ATOM     77  OE1 GLN A   9       0.013  -3.558   3.156  1.00 57.12           O
ATOM     78  N   ALA A  10      -2.516  -2.624  -2.207  1.00 39.66           N
ATOM     79  CA  ALA A  10      -3.023  -1.428  -2.819  1.00 39.66           C
ATOM     80  C   ALA A  10      -4.323  -1.038  -2.192  1.00 39.66           C
ATOM     81  O   ALA A  10      -5.190  -1.876  -1.953  1.00 39.66           O
ATOM     82  CB  ALA A  10      -3.272  -1.575  -4.330  1.00 39.66           C
ATOM     83  N   GLU A  11      -4.484   0.275  -1.918  1.00 60.76           N
ATOM     84  CA  GLU A  11      -5.693   0.766  -1.315  1.00 60.76           C
ATOM     85  C   GLU A  11      -5.979   2.139  -1.849  1.00 60.76           C
ATOM     86  O   GLU A  11      -5.073   2.839  -2.302  1.00 60.76           O
ATOM     87  CB  GLU A  11      -5.599   0.827   0.220  1.00 60.76           C
ATOM     88  CG  GLU A  11      -4.375   1.593   0.728  1.00 60.76           C
ATOM     89  CD  GLU A  11      -4.221   1.302   2.215  1.00 60.76           C
ATOM     90  OE1 GLU A  11      -5.144   1.668   2.991  1.00 60.76           O
ATOM     91  OE2 GLU A  11      -3.181   0.702   2.595  1.00 60.76           O1-
ATOM     92  N   PHE A  12      -7.270   2.550  -1.853  1.00 71.24           N
ATOM     93  CA  PHE A  12      -7.579   3.872  -2.323  1.00 71.24           C
ATOM     94  C   PHE A  12      -8.860   4.355  -1.707  1.00 71.24           C
ATOM     95  O   PHE A  12      -9.642   3.570  -1.168  1.00 71.24           O
ATOM     96  CB  PHE A  12      -7.732   3.984  -3.855  1.00 71.24           C
ATOM     97  CG  PHE A  12      -8.964   3.271  -4.305  1.00 71.24           C
ATOM     98  CD1 PHE A  12     -10.180   3.916  -4.307  1.00 71.24           C
ATOM     99  CD2 PHE A  12      -8.910   1.966  -4.740  1.00 71.24           C
ATOM    100  CE1 PHE A  12     -11.323   3.276  -4.726  1.00 71.24           C
ATOM    101  CE2 PHE A  12     -10.048   1.319  -5.161  1.00 71.24           C
ATOM    102  CZ  PHE A  12     -11.259   1.971  -5.154  1.00 71.24           C
ATOM    103  N   TYR A  13      -9.071   5.691  -1.735  1.00 81.78           N
ATOM    104  CA  TYR A  13     -10.284   6.283  -1.238  1.00 81.78           C
ATOM    105  C   TYR A  13     -10.654   7.340  -2.249  1.00 81.78           C
ATOM    106  O   TYR A  13      -9.812   8.165  -2.606  1.00 81.78           O
ATOM    107  CB  TYR A  13     -10.073   6.967   0.127  1.00 81.78           C
ATOM    108  CG  TYR A  13     -11.388   7.195   0.792  1.00 81.78           C
ATOM    109  CD1 TYR A  13     -12.003   6.154   1.451  1.00 81.78           C
ATOM    110  CD2 TYR A  13     -11.996   8.431   0.784  1.00 81.78           C
ATOM    111  CE1 TYR A  13     -13.212   6.336   2.080  1.00 81.78           C
ATOM    112  CE2 TYR A  13     -13.205   8.617   1.412  1.00 81.78           C
ATOM    113  CZ  TYR A  13     -13.815   7.570   2.059  1.00 81.78           C
ATOM    114  OH  TYR A  13     -15.056   7.759   2.705  1.00 81.78           O
ATOM    115  N   LEU A  14     -11.919   7.339  -2.736  1.00132.50           N
ATOM    116  CA  LEU A  14     -12.382   8.288  -3.725  1.00132.50           C
ATOM    117  C   LEU A  14     -13.379   9.177  -3.029  1.00132.50           C
ATOM    118  O   LEU A  14     -14.138   8.716  -2.180  1.00132.50           O
ATOM    119  CB  LEU A  14     -13.041   7.595  -4.933  1.00132.50           C
ATOM    120  CG  LEU A  14     -13.475   8.513  -6.093  1.00132.50           C
ATOM    121  CD1 LEU A  14     -13.914   7.674  -7.295  1.00132.50           C
ATOM    122  CD2 LEU A  14     -14.581   9.499  -5.695  1.00132.50           C
ATOM    123  N   ASN A  15     -13.369  10.493  -3.342  1.00133.27           N
ATOM    124  CA  ASN A  15     -14.108  11.401  -2.505  1.00133.27           C
ATOM    125  C   ASN A  15     -15.615  11.477  -2.583  1.00133.27           C
ATOM    126  O   ASN A  15     -16.241  11.081  -1.602  1.00133.27           O
ATOM    127  CB  ASN A  15     -13.532  12.828  -2.548  1.00133.27           C
ATOM    128  CG  ASN A  15     -13.861  13.496  -1.220  1.00133.27           C
ATOM    129  ND2 ASN A  15     -12.896  14.301  -0.698  1.00133.27           N
ATOM    130  OD1 ASN A  15     -14.930  13.293  -0.648  1.00133.27           O
ATOM    131  N   PRO A  16     -16.288  11.900  -3.637  1.00105.76           N
ATOM    132  CA  PRO A  16     -17.714  12.067  -3.473  1.00105.76           C
ATOM    133  C   PRO A  16     -18.391  10.761  -3.248  1.00105.76           C
ATOM    134  O   PRO A  16     -19.363  10.689  -2.498  1.00105.76           O
ATOM    135  CB  PRO A  16     -18.187  12.842  -4.696  1.00105.76           C
ATOM    136  CG  PRO A  16     -16.955  13.688  -5.072  1.00105.76           C
ATOM    137  CD  PRO A  16     -15.750  12.856  -4.594  1.00105.76           C
ATOM    138  N   ASP A  17     -17.843   9.710  -3.860  1.00 95.36           N
ATOM    139  CA  ASP A  17     -18.388   8.396  -3.779  1.00 95.36           C
ATOM    140  C   ASP A  17     -18.248   7.888  -2.382  1.00 95.36           C
ATOM    141  O   ASP A  17     -19.069   7.099  -1.918  1.00 95.36           O
ATOM    142  CB  ASP A  17     -17.652   7.421  -4.705  1.00 95.36           C
ATOM    143  CG  ASP A  17     -17.889   7.903  -6.128  1.00 95.36           C
ATOM    144  OD1 ASP A  17     -19.000   7.659  -6.670  1.00 95.36           O
ATOM    145  OD2 ASP A  17     -16.951   8.528  -6.689  1.00 95.36           O1-
ATOM    146  N   GLN A  18     -17.210   8.348  -1.660  1.00 57.89           N
ATOM    147  CA  GLN A  18     -16.968   7.796  -0.359  1.00 57.89           C
ATOM    148  C   GLN A  18     -16.758   6.328  -0.542  1.00 57.89           C
ATOM    149  O   GLN A  18     -17.360   5.507   0.149  1.00 57.89           O
ATOM    150  CB  GLN A  18     -18.126   8.007   0.633  1.00 57.89           C
ATOM    151  CG  GLN A  18     -18.345   9.472   1.018  1.00 57.89           C
ATOM    152  CD  GLN A  18     -19.510   9.536   1.997  1.00 57.89           C
ATOM    153  NE2 GLN A  18     -19.429  10.486   2.967  1.00 57.89           N
ATOM    154  OE1 GLN A  18     -20.460   8.762   1.908  1.00 57.89           O
ATOM    155  N   SER A  19     -15.890   5.971  -1.513  1.00 86.18           N
ATOM    156  CA  SER A  19     -15.592   4.590  -1.766  1.00 86.18           C
ATOM    157  C   SER A  19     -14.139   4.354  -1.489  1.00 86.18           C
ATOM    158  O   SER A  19     -13.293   5.205  -1.753  1.00 86.18           O
ATOM    159  CB  SER A  19     -15.823   4.152  -3.221  1.00 86.18           C
ATOM    160  OG  SER A  19     -17.211   4.120  -3.510  1.00 86.18           O
ATOM    161  N   GLY A  20     -13.814   3.164  -0.946  1.00 28.54           N
ATOM    162  CA  GLY A  20     -12.442   2.851  -0.659  1.00 28.54           C
ATOM    163  C   GLY A  20     -12.261   1.371  -0.801  1.00 28.54           C
ATOM    164  O   GLY A  20     -13.198   0.600  -0.599  1.00 28.54           O
ATOM    165  N   GLU A  21     -11.032   0.933  -1.152  1.00 50.31           N
ATOM    166  CA  GLU A  21     -10.795  -0.475  -1.315  1.00 50.31           C
ATOM    167  C   GLU A  21      -9.416  -0.812  -0.834  1.00 50.31           C
ATOM    168  O   GLU A  21      -8.531   0.042  -0.783  1.00 50.31           O
ATOM    169  CB  GLU A  21     -10.917  -0.944  -2.773  1.00 50.31           C
ATOM    170  CG  GLU A  21     -12.352  -0.872  -3.296  1.00 50.31           C
ATOM    171  CD  GLU A  21     -12.349  -1.296  -4.756  1.00 50.31           C
ATOM    172  OE1 GLU A  21     -11.279  -1.163  -5.406  1.00 50.31           O
ATOM    173  OE2 GLU A  21     -13.413  -1.762  -5.241  1.00 50.31           O1-
ATOM    174  N   PHE A  22      -9.215  -2.090  -0.440  1.00122.00           N
ATOM    175  CA  PHE A  22      -7.950  -2.548   0.067  1.00122.00           C
ATOM    176  C   PHE A  22      -7.734  -3.934  -0.467  1.00122.00           C
ATOM    177  O   PHE A  22      -8.619  -4.784  -0.369  1.00122.00           O
ATOM    178  CB  PHE A  22      -7.987  -2.615   1.609  1.00122.00           C
ATOM    179  CG  PHE A  22      -6.687  -3.074   2.176  1.00122.00           C
ATOM    180  CD1 PHE A  22      -6.362  -4.410   2.205  1.00122.00           C
ATOM    181  CD2 PHE A  22      -5.803  -2.162   2.706  1.00122.00           C
ATOM    182  CE1 PHE A  22      -5.164  -4.825   2.740  1.00122.00           C
ATOM    183  CE2 PHE A  22      -4.605  -2.572   3.242  1.00122.00           C
ATOM    184  CZ  PHE A  22      -4.283  -3.907   3.258  1.00122.00           C
ATOM    185  N   MET A  23      -6.548  -4.213  -1.053  1.00132.60           N
ATOM    186  CA  MET A  23      -6.345  -5.539  -1.569  1.00132.60           C
ATOM    187  C   MET A  23      -4.882  -5.810  -1.760  1.00132.60           C
ATOM    188  O   MET A  23      -4.053  -4.901  -1.770  1.00132.60           O
ATOM    189  CB  MET A  23      -7.083  -5.779  -2.902  1.00132.60           C
ATOM    190  CG  MET A  23      -6.670  -4.828  -4.028  1.00132.60           C
ATOM    191  SD  MET A  23      -5.068  -5.205  -4.795  1.00132.60           S
ATOM    192  CE  MET A  23      -5.716  -6.602  -5.756  1.00132.60           C
ATOM    193  N   PHE A  24      -4.539  -7.111  -1.887  1.00 76.76           N
ATOM    194  CA  PHE A  24      -3.191  -7.558  -2.110  1.00 76.76           C
ATOM    195  C   PHE A  24      -3.099  -8.102  -3.496  1.00 76.76           C
ATOM    196  O   PHE A  24      -4.036  -8.728  -3.992  1.00 76.76           O
ATOM    197  CB  PHE A  24      -2.739  -8.688  -1.166  1.00 76.76           C
ATOM    198  CG  PHE A  24      -2.218  -8.088   0.094  1.00 76.76           C
ATOM    199  CD1 PHE A  24      -3.030  -7.379   0.946  1.00 76.76           C
ATOM    200  CD2 PHE A  24      -0.897  -8.272   0.431  1.00 76.76           C
ATOM    201  CE1 PHE A  24      -2.523  -6.842   2.107  1.00 76.76           C
ATOM    202  CE2 PHE A  24      -0.386  -7.739   1.589  1.00 76.76           C
ATOM    203  CZ  PHE A  24      -1.199  -7.017   2.430  1.00 76.76           C
ATOM    204  N   ASP A  25      -1.952  -7.860  -4.164  1.00 63.70           N
ATOM    205  CA  ASP A  25      -1.789  -8.306  -5.518  1.00 63.70           C
ATOM    206  C   ASP A  25      -0.492  -9.050  -5.632  1.00 63.70           C
ATOM    207  O   ASP A  25       0.502  -8.690  -4.999  1.00 63.70           O
ATOM    208  CB  ASP A  25      -1.762  -7.129  -6.513  1.00 63.70           C
ATOM    209  CG  ASP A  25      -1.993  -7.636  -7.932  1.00 63.70           C
ATOM    210  OD1 ASP A  25      -1.284  -8.586  -8.358  1.00 63.70           O
ATOM    211  OD2 ASP A  25      -2.893  -7.075  -8.612  1.00 63.70           O1-
ATOM    212  N   PHE A  26      -0.490 -10.151  -6.414  1.00159.11           N
ATOM    213  CA  PHE A  26       0.712 -10.896  -6.655  1.00159.11           C
ATOM    214  C   PHE A  26       0.879 -11.024  -8.139  1.00159.11           C
ATOM    215  O   PHE A  26       0.124 -11.734  -8.800  1.00159.11           O
ATOM    216  CB  PHE A  26       0.652 -12.326  -6.105  1.00159.11           C
ATOM    217  CG  PHE A  26       1.948 -12.994  -6.411  1.00159.11           C
ATOM    218  CD1 PHE A  26       2.137 -13.641  -7.612  1.00159.11           C
ATOM    219  CD2 PHE A  26       2.979 -12.967  -5.500  1.00159.11           C
ATOM    220  CE1 PHE A  26       3.334 -14.259  -7.891  1.00159.11           C
ATOM    221  CE2 PHE A  26       4.177 -13.584  -5.772  1.00159.11           C
ATOM    222  CZ  PHE A  26       4.356 -14.234  -6.969  1.00159.11           C
ATOM    223  N   ASP A  27       1.898 -10.345  -8.697  1.00 55.28           N
ATOM    224  CA  ASP A  27       2.194 -10.441 -10.097  1.00 55.28           C
ATOM    225  C   ASP A  27       0.962 -10.199 -10.916  1.00 55.28           C
ATOM    226  O   ASP A  27       0.708 -10.912 -11.887  1.00 55.28           O
ATOM    227  CB  ASP A  27       2.785 -11.805 -10.498  1.00 55.28           C
ATOM    228  CG  ASP A  27       4.187 -11.900  -9.910  1.00 55.28           C
ATOM    229  OD1 ASP A  27       4.693 -10.861  -9.409  1.00 55.28           O
ATOM    230  OD2 ASP A  27       4.771 -13.015  -9.957  1.00 55.28           O1-
ATOM    231  N   GLY A  28       0.158  -9.183 -10.546  1.00 32.08           N
ATOM    232  CA  GLY A  28      -0.964  -8.819 -11.365  1.00 32.08           C
ATOM    233  C   GLY A  28      -2.233  -9.493 -10.931  1.00 32.08           C
ATOM    234  O   GLY A  28      -3.304  -9.152 -11.432  1.00 32.08           O
ATOM    235  N   ASP A  29      -2.177 -10.462  -9.996  1.00 68.40           N
ATOM    236  CA  ASP A  29      -3.404 -11.119  -9.626  1.00 68.40           C
ATOM    237  C   ASP A  29      -3.769 -10.784  -8.212  1.00 68.40           C
ATOM    238  O   ASP A  29      -2.907 -10.606  -7.354  1.00 68.40           O
ATOM    239  CB  ASP A  29      -3.351 -12.650  -9.774  1.00 68.40           C
ATOM    240  CG  ASP A  29      -3.477 -12.972 -11.260  1.00 68.40           C
ATOM    241  OD1 ASP A  29      -4.465 -12.496 -11.878  1.00 68.40           O
ATOM    242  OD2 ASP A  29      -2.589 -13.686 -11.800  1.00 68.40           O1-
ATOM    243  N   GLU A  30      -5.089 -10.670  -7.948  1.00 69.98           N
ATOM    244  CA  GLU A  30      -5.598 -10.322  -6.650  1.00 69.98           C
ATOM    245  C   GLU A  30      -5.501 -11.510  -5.741  1.00 69.98           C
ATOM    246  O   GLU A  30      -6.067 -12.565  -6.031  1.00 69.98           O
ATOM    247  CB  GLU A  30      -7.081  -9.914  -6.687  1.00 69.98           C
ATOM    248  CG  GLU A  30      -7.655  -9.536  -5.321  1.00 69.98           C
ATOM    249  CD  GLU A  30      -9.171  -9.488  -5.443  1.00 69.98           C
ATOM    250  OE1 GLU A  30      -9.670  -9.049  -6.513  1.00 69.98           O
ATOM    251  OE2 GLU A  30      -9.849  -9.906  -4.466  1.00 69.98           O1-
ATOM    252  N   ILE A  31      -4.733 -11.384  -4.636  1.00121.91           N
ATOM    253  CA  ILE A  31      -4.670 -12.430  -3.654  1.00121.91           C
ATOM    254  C   ILE A  31      -5.959 -12.424  -2.888  1.00121.91           C
ATOM    255  O   ILE A  31      -6.630 -13.447  -2.752  1.00121.91           O
ATOM    256  CB  ILE A  31      -3.593 -12.195  -2.633  1.00121.91           C
ATOM    257  CG1 ILE A  31      -2.198 -12.114  -3.279  1.00121.91           C
ATOM    258  CG2 ILE A  31      -3.707 -13.316  -1.586  1.00121.91           C
ATOM    259  CD1 ILE A  31      -1.715 -13.434  -3.869  1.00121.91           C
ATOM    260  N   PHE A  32      -6.330 -11.236  -2.365  1.00 67.51           N
ATOM    261  CA  PHE A  32      -7.545 -11.085  -1.613  1.00 67.51           C
ATOM    262  C   PHE A  32      -7.827  -9.621  -1.499  1.00 67.51           C
ATOM    263  O   PHE A  32      -6.958  -8.788  -1.752  1.00 67.51           O
ATOM    264  CB  PHE A  32      -7.481 -11.666  -0.185  1.00 67.51           C
ATOM    265  CG  PHE A  32      -6.506 -10.897   0.645  1.00 67.51           C
ATOM    266  CD1 PHE A  32      -6.893  -9.743   1.290  1.00 67.51           C
ATOM    267  CD2 PHE A  32      -5.209 -11.335   0.793  1.00 67.51           C
ATOM    268  CE1 PHE A  32      -6.004  -9.033   2.061  1.00 67.51           C
ATOM    269  CE2 PHE A  32      -4.315 -10.628   1.562  1.00 67.51           C
ATOM    270  CZ  PHE A  32      -4.711  -9.475   2.196  1.00 67.51           C
ATOM    271  N   HIS A  33      -9.076  -9.273  -1.134  1.00 60.32           N
ATOM    272  CA  HIS A  33      -9.418  -7.899  -0.912  1.00 60.32           C
ATOM    273  C   HIS A  33     -10.321  -7.897   0.278  1.00 60.32           C
ATOM    274  O   HIS A  33     -10.862  -8.941   0.644  1.00 60.32           O
ATOM    275  CB  HIS A  33     -10.148  -7.214  -2.083  1.00 60.32           C
ATOM    276  CG  HIS A  33     -11.563  -7.668  -2.288  1.00 60.32           C
ATOM    277  CD2 HIS A  33     -12.729  -7.096  -1.879  1.00 60.32           C
ATOM    278  ND1 HIS A  33     -11.924  -8.799  -2.984  1.00 60.32           N
ATOM    279  CE1 HIS A  33     -13.279  -8.854  -2.964  1.00 60.32           C
ATOM    280  NE2 HIS A  33     -13.812  -7.841  -2.305  1.00 60.32           N
ATOM    281  N   VAL A  34     -10.492  -6.730   0.930  1.00 39.78           N
ATOM    282  CA  VAL A  34     -11.311  -6.692   2.110  1.00 39.78           C
ATOM    283  C   VAL A  34     -12.584  -5.975   1.796  1.00 39.78           C
ATOM    284  O   VAL A  34     -12.577  -4.889   1.218  1.00 39.78           O
ATOM    285  CB  VAL A  34     -10.668  -5.971   3.258  1.00 39.78           C
ATOM    286  CG1 VAL A  34     -11.685  -5.874   4.407  1.00 39.78           C
ATOM    287  CG2 VAL A  34      -9.370  -6.707   3.636  1.00 39.78           C
ATOM    288  N   ASP A  35     -13.723  -6.585   2.188  1.00 41.08           N
ATOM    289  CA  ASP A  35     -15.018  -6.006   1.972  1.00 41.08           C
ATOM    290  C   ASP A  35     -15.289  -5.133   3.152  1.00 41.08           C
ATOM    291  O   ASP A  35     -15.381  -5.604   4.282  1.00 41.08           O
ATOM    292  CB  ASP A  35     -16.144  -7.053   1.902  1.00 41.08           C
ATOM    293  CG  ASP A  35     -17.473  -6.345   1.672  1.00 41.08           C
ATOM    294  OD1 ASP A  35     -17.465  -5.100   1.478  1.00 41.08           O
ATOM    295  OD2 ASP A  35     -18.518  -7.047   1.699  1.00 41.08           O1-
ATOM    296  N   MET A  36     -15.415  -3.819   2.906  1.00 59.59           N
ATOM    297  CA  MET A  36     -15.590  -2.895   3.984  1.00 59.59           C
ATOM    298  C   MET A  36     -16.894  -3.122   4.694  1.00 59.59           C
ATOM    299  O   MET A  36     -16.942  -3.057   5.922  1.00 59.59           O
ATOM    300  CB  MET A  36     -15.526  -1.431   3.523  1.00 59.59           C
ATOM    301  CG  MET A  36     -14.168  -1.097   2.903  1.00 59.59           C
ATOM    302  SD  MET A  36     -12.754  -1.378   4.009  1.00 59.59           S
ATOM    303  CE  MET A  36     -11.491  -1.188   2.719  1.00 59.59           C
ATOM    304  N   ALA A  37     -17.994  -3.377   3.954  1.00 43.45           N
ATOM    305  CA  ALA A  37     -19.270  -3.513   4.610  1.00 43.45           C
ATOM    306  C   ALA A  37     -19.310  -4.717   5.506  1.00 43.45           C
ATOM    307  O   ALA A  37     -19.594  -4.612   6.697  1.00 43.45           O
ATOM    308  CB  ALA A  37     -20.432  -3.652   3.613  1.00 43.45           C
ATOM    309  N   LYS A  38     -18.983  -5.895   4.947  1.00115.81           N
ATOM    310  CA  LYS A  38     -18.998  -7.161   5.630  1.00115.81           C
ATOM    311  C   LYS A  38     -17.864  -7.189   6.604  1.00115.81           C
ATOM    312  O   LYS A  38     -17.878  -7.944   7.573  1.00115.81           O
ATOM    313  CB  LYS A  38     -18.776  -8.343   4.671  1.00115.81           C
ATOM    314  CG  LYS A  38     -19.874  -8.518   3.622  1.00115.81           C
ATOM    315  CD  LYS A  38     -21.238  -8.884   4.204  1.00115.81           C
ATOM    316  CE  LYS A  38     -22.319  -9.058   3.135  1.00115.81           C
ATOM    317  NZ  LYS A  38     -23.581  -9.517   3.759  1.00115.81           N1+
ATOM    318  N   LYS A  39     -16.821  -6.383   6.333  1.00 65.49           N
ATOM    319  CA  LYS A  39     -15.638  -6.373   7.144  1.00 65.49           C
ATOM    320  C   LYS A  39     -14.998  -7.723   7.091  1.00 65.49           C
ATOM    321  O   LYS A  39     -14.494  -8.218   8.098  1.00 65.49           O
ATOM    322  CB  LYS A  39     -15.892  -6.050   8.624  1.00 65.49           C
ATOM    323  CG  LYS A  39     -16.345  -4.616   8.891  1.00 65.49           C
ATOM    324  CD  LYS A  39     -16.756  -4.390  10.345  1.00 65.49           C
ATOM    325  CE  LYS A  39     -17.077  -2.933  10.678  1.00 65.49           C
ATOM    326  NZ  LYS A  39     -17.141  -2.766  12.145  1.00 65.49           N1+
ATOM    327  N   GLU A  40     -14.981  -8.358   5.904  1.00 44.44           N
ATOM    328  CA  GLU A  40     -14.397  -9.667   5.826  1.00 44.44           C
ATOM    329  C   GLU A  40     -13.373  -9.705   4.733  1.00 44.44           C
ATOM    330  O   GLU A  40     -13.381  -8.887   3.814  1.00 44.44           O
ATOM    331  CB  GLU A  40     -15.431 -10.769   5.541  1.00 44.44           C
ATOM    332  CG  GLU A  40     -16.162 -10.592   4.208  1.00 44.44           C
ATOM    333  CD  GLU A  40     -17.215 -11.687   4.103  1.00 44.44           C
ATOM    334  OE1 GLU A  40     -17.561 -12.280   5.160  1.00 44.44           O
ATOM    335  OE2 GLU A  40     -17.692 -11.941   2.965  1.00 44.44           O1-
ATOM    336  N   THR A  41     -12.449 -10.682   4.824  1.00 43.70           N
ATOM    337  CA  THR A  41     -11.424 -10.841   3.834  1.00 43.70           C
ATOM    338  C   THR A  41     -11.967 -11.783   2.809  1.00 43.70           C
ATOM    339  O   THR A  41     -12.417 -12.879   3.140  1.00 43.70           O
ATOM    340  CB  THR A  41     -10.162 -11.448   4.378  1.00 43.70           C
ATOM    341  CG2 THR A  41      -9.152 -11.631   3.234  1.00 43.70           C
ATOM    342  OG1 THR A  41      -9.618 -10.620   5.395  1.00 43.70           O
ATOM    343  N   VAL A  42     -11.963 -11.357   1.529  1.00 45.89           N
ATOM    344  CA  VAL A  42     -12.486 -12.187   0.481  1.00 45.89           C
ATOM    345  C   VAL A  42     -11.331 -12.659  -0.348  1.00 45.89           C
ATOM    346  O   VAL A  42     -10.641 -11.866  -0.985  1.00 45.89           O
ATOM    347  CB  VAL A  42     -13.417 -11.442  -0.434  1.00 45.89           C
ATOM    348  CG1 VAL A  42     -13.902 -12.393  -1.542  1.00 45.89           C
ATOM    349  CG2 VAL A  42     -14.544 -10.828   0.413  1.00 45.89           C
ATOM    350  N   TRP A  43     -11.098 -13.985  -0.372  1.00 63.10           N
ATOM    351  CA  TRP A  43      -9.987 -14.502  -1.118  1.00 63.10           C
ATOM    352  C   TRP A  43     -10.400 -14.670  -2.547  1.00 63.10           C
ATOM    353  O   TRP A  43     -11.533 -15.050  -2.842  1.00 63.10           O
ATOM    354  CB  TRP A  43      -9.480 -15.847  -0.571  1.00 63.10           C
ATOM    355  CG  TRP A  43      -8.935 -15.721   0.833  1.00 63.10           C
ATOM    356  CD1 TRP A  43      -9.588 -15.860   2.022  1.00 63.10           C
ATOM    357  CD2 TRP A  43      -7.572 -15.398   1.154  1.00 63.10           C
ATOM    358  CE2 TRP A  43      -7.475 -15.355   2.544  1.00 63.10           C
ATOM    359  CE3 TRP A  43      -6.491 -15.155   0.356  1.00 63.10           C
ATOM    360  NE1 TRP A  43      -8.720 -15.640   3.065  1.00 63.10           N
ATOM    361  CZ2 TRP A  43      -6.290 -15.067   3.161  1.00 63.10           C
ATOM    362  CZ3 TRP A  43      -5.297 -14.866   0.981  1.00 63.10           C
ATOM    363  CH2 TRP A  43      -5.199 -14.823   2.357  1.00 63.10           C
ATOM    364  N   ARG A  44      -9.469 -14.380  -3.478  1.00111.45           N
ATOM    365  CA  ARG A  44      -9.769 -14.468  -4.881  1.00111.45           C
ATOM    366  C   ARG A  44     -10.103 -15.890  -5.197  1.00111.45           C
ATOM    367  O   ARG A  44     -11.069 -16.164  -5.908  1.00111.45           O
ATOM    368  CB  ARG A  44      -8.586 -14.060  -5.775  1.00111.45           C
ATOM    369  CG  ARG A  44      -8.954 -14.032  -7.258  1.00111.45           C
ATOM    370  CD  ARG A  44      -9.942 -12.919  -7.609  1.00111.45           C
ATOM    371  NE  ARG A  44     -10.273 -13.050  -9.054  1.00111.45           N
ATOM    372  CZ  ARG A  44     -11.286 -13.878  -9.445  1.00111.45           C
ATOM    373  NH1 ARG A  44     -11.991 -14.582  -8.512  1.00111.45           N1+
ATOM    374  NH2 ARG A  44     -11.595 -14.002 -10.769  1.00111.45           N
ATOM    375  N   LEU A  45      -9.295 -16.832  -4.672  1.00 57.61           N
ATOM    376  CA  LEU A  45      -9.559 -18.227  -4.866  1.00 57.61           C
ATOM    377  C   LEU A  45      -9.772 -18.779  -3.492  1.00 57.61           C
ATOM    378  O   LEU A  45      -9.018 -18.473  -2.572  1.00 57.61           O
ATOM    379  CB  LEU A  45      -8.383 -18.991  -5.500  1.00 57.61           C
ATOM    380  CG  LEU A  45      -8.035 -18.507  -6.922  1.00 57.61           C
ATOM    381  CD1 LEU A  45      -6.858 -19.298  -7.519  1.00 57.61           C
ATOM    382  CD2 LEU A  45      -9.275 -18.497  -7.829  1.00 57.61           C
ATOM    383  N   GLU A  46     -10.795 -19.634  -3.320  1.00 74.84           N
ATOM    384  CA  GLU A  46     -11.148 -20.122  -2.017  1.00 74.84           C
ATOM    385  C   GLU A  46     -10.002 -20.893  -1.441  1.00 74.84           C
ATOM    386  O   GLU A  46      -9.778 -20.876  -0.231  1.00 74.84           O
ATOM    387  CB  GLU A  46     -12.376 -21.048  -2.018  1.00 74.84           C
ATOM    388  CG  GLU A  46     -12.858 -21.380  -0.604  1.00 74.84           C
ATOM    389  CD  GLU A  46     -14.070 -22.297  -0.701  1.00 74.84           C
ATOM    390  OE1 GLU A  46     -14.222 -22.968  -1.757  1.00 74.84           O
ATOM    391  OE2 GLU A  46     -14.859 -22.340   0.281  1.00 74.84           O1-
ATOM    392  N   GLU A  47      -9.230 -21.575  -2.301  1.00 53.70           N
ATOM    393  CA  GLU A  47      -8.154 -22.408  -1.861  1.00 53.70           C
ATOM    394  C   GLU A  47      -7.192 -21.580  -1.059  1.00 53.70           C
ATOM    395  O   GLU A  47      -6.641 -22.057  -0.070  1.00 53.70           O
ATOM    396  CB  GLU A  47      -7.376 -23.007  -3.044  1.00 53.70           C
ATOM    397  CG  GLU A  47      -6.794 -21.938  -3.972  1.00 53.70           C
ATOM    398  CD  GLU A  47      -6.049 -22.629  -5.105  1.00 53.70           C
ATOM    399  OE1 GLU A  47      -6.000 -23.888  -5.103  1.00 53.70           O
ATOM    400  OE2 GLU A  47      -5.517 -21.906  -5.989  1.00 53.70           O1-
ATOM    401  N   PHE A  48      -6.976 -20.307  -1.437  1.00 83.50           N
ATOM    402  CA  PHE A  48      -6.014 -19.496  -0.741  1.00 83.50           C
ATOM    403  C   PHE A  48      -6.408 -19.396   0.704  1.00 83.50           C
ATOM    404  O   PHE A  48      -5.571 -19.492   1.600  1.00 83.50           O
ATOM    405  CB  PHE A  48      -5.918 -18.044  -1.250  1.00 83.50           C
ATOM    406  CG  PHE A  48      -5.484 -18.021  -2.677  1.00 83.50           C
ATOM    407  CD1 PHE A  48      -4.277 -18.552  -3.062  1.00 83.50           C
ATOM    408  CD2 PHE A  48      -6.271 -17.417  -3.632  1.00 83.50           C
ATOM    409  CE1 PHE A  48      -3.879 -18.515  -4.379  1.00 83.50           C
ATOM    410  CE2 PHE A  48      -5.878 -17.374  -4.950  1.00 83.50           C
ATOM    411  CZ  PHE A  48      -4.681 -17.928  -5.328  1.00 83.50           C
ATOM    412  N   GLY A  49      -7.713 -19.229   0.961  1.00 54.06           N
ATOM    413  CA  GLY A  49      -8.246 -19.019   2.277  1.00 54.06           C
ATOM    414  C   GLY A  49      -7.909 -20.169   3.173  1.00 54.06           C
ATOM    415  O   GLY A  49      -7.937 -20.021   4.390  1.00 54.06           O
ATOM    416  N   ARG A  50      -7.718 -21.384   2.635  1.00175.24           N
ATOM    417  CA  ARG A  50      -7.386 -22.455   3.534  1.00175.24           C
ATOM    418  C   ARG A  50      -5.996 -22.296   4.085  1.00175.24           C
ATOM    419  O   ARG A  50      -5.754 -22.613   5.248  1.00175.24           O
ATOM    420  CB  ARG A  50      -7.578 -23.859   2.930  1.00175.24           C
ATOM    421  CG  ARG A  50      -6.765 -24.157   1.675  1.00175.24           C
ATOM    422  CD  ARG A  50      -7.484 -25.142   0.751  1.00175.24           C
ATOM    423  NE  ARG A  50      -8.759 -24.478   0.353  1.00175.24           N
ATOM    424  CZ  ARG A  50      -9.729 -25.163  -0.321  1.00175.24           C
ATOM    425  NH1 ARG A  50      -9.541 -26.474  -0.651  1.00175.24           N1+
ATOM    426  NH2 ARG A  50     -10.893 -24.535  -0.664  1.00175.24           N
ATOM    427  N   PHE A  51      -5.043 -21.847   3.243  1.00 99.79           N
ATOM    428  CA  PHE A  51      -3.653 -21.683   3.584  1.00 99.79           C
ATOM    429  C   PHE A  51      -3.318 -20.449   4.387  1.00 99.79           C
ATOM    430  O   PHE A  51      -2.345 -20.470   5.142  1.00 99.79           O
ATOM    431  CB  PHE A  51      -2.739 -21.716   2.350  1.00 99.79           C
ATOM    432  CG  PHE A  51      -2.962 -23.062   1.756  1.00 99.79           C
ATOM    433  CD1 PHE A  51      -2.588 -24.192   2.444  1.00 99.79           C
ATOM    434  CD2 PHE A  51      -3.535 -23.195   0.513  1.00 99.79           C
ATOM    435  CE1 PHE A  51      -2.792 -25.439   1.905  1.00 99.79           C
ATOM    436  CE2 PHE A  51      -3.741 -24.440  -0.033  1.00 99.79           C
ATOM    437  CZ  PHE A  51      -3.371 -25.565   0.665  1.00 99.79           C
ATOM    438  N   ALA A  52      -4.038 -19.316   4.218  1.00 53.89           N
ATOM    439  CA  ALA A  52      -3.601 -18.120   4.896  1.00 53.89           C
ATOM    440  C   ALA A  52      -4.758 -17.336   5.442  1.00 53.89           C
ATOM    441  O   ALA A  52      -5.921 -17.675   5.227  1.00 53.89           O
ATOM    442  CB  ALA A  52      -2.805 -17.173   3.984  1.00 53.89           C
ATOM    443  N   SER A  53      -4.445 -16.264   6.214  1.00 85.96           N
ATOM    444  CA  SER A  53      -5.467 -15.450   6.814  1.00 85.96           C
ATOM    445  C   SER A  53      -5.013 -14.015   6.799  1.00 85.96           C
ATOM    446  O   SER A  53      -3.861 -13.722   6.485  1.00 85.96           O
ATOM    447  CB  SER A  53      -5.729 -15.834   8.279  1.00 85.96           C
ATOM    448  OG  SER A  53      -6.742 -15.012   8.836  1.00 85.96           O
ATOM    449  N   PHE A  54      -5.940 -13.071   7.090  1.00 71.95           N
ATOM    450  CA  PHE A  54      -5.614 -11.670   7.137  1.00 71.95           C
ATOM    451  C   PHE A  54      -6.625 -10.989   8.015  1.00 71.95           C
ATOM    452  O   PHE A  54      -7.826 -11.232   7.896  1.00 71.95           O
ATOM    453  CB  PHE A  54      -5.666 -11.018   5.740  1.00 71.95           C
ATOM    454  CG  PHE A  54      -5.325  -9.568   5.828  1.00 71.95           C
ATOM    455  CD1 PHE A  54      -4.022  -9.150   5.973  1.00 71.95           C
ATOM    456  CD2 PHE A  54      -6.320  -8.620   5.739  1.00 71.95           C
ATOM    457  CE1 PHE A  54      -3.721  -7.809   6.047  1.00 71.95           C
ATOM    458  CE2 PHE A  54      -6.024  -7.279   5.810  1.00 71.95           C
ATOM    459  CZ  PHE A  54      -4.722  -6.871   5.966  1.00 71.95           C
ATOM    460  N   GLU A  55      -6.162 -10.089   8.911  1.00 93.00           N
ATOM    461  CA  GLU A  55      -7.069  -9.412   9.798  1.00 93.00           C
ATOM    462  C   GLU A  55      -7.660  -8.249   9.063  1.00 93.00           C
ATOM    463  O   GLU A  55      -7.011  -7.225   8.854  1.00 93.00           O
ATOM    464  CB  GLU A  55      -6.393  -8.877  11.073  1.00 93.00           C
ATOM    465  CG  GLU A  55      -7.356  -8.158  12.018  1.00 93.00           C
ATOM    466  CD  GLU A  55      -8.383  -9.167  12.511  1.00 93.00           C
ATOM    467  OE1 GLU A  55      -7.977 -10.309  12.852  1.00 93.00           O
ATOM    468  OE2 GLU A  55      -9.589  -8.808  12.547  1.00 93.00           O1-
ATOM    469  N   ALA A  56      -8.950  -8.385   8.696  1.00 41.29           N
ATOM    470  CA  ALA A  56      -9.684  -7.436   7.902  1.00 41.29           C
ATOM    471  C   ALA A  56      -9.801  -6.115   8.598  1.00 41.29           C
ATOM    472  O   ALA A  56      -9.725  -5.067   7.960  1.00 41.29           O
ATOM    473  CB  ALA A  56     -11.110  -7.910   7.582  1.00 41.29           C
ATOM    474  N   GLN A  57      -9.963  -6.130   9.931  1.00 58.06           N
ATOM    475  CA  GLN A  57     -10.195  -4.931  10.685  1.00 58.06           C
ATOM    476  C   GLN A  57      -9.071  -3.988  10.393  1.00 58.06           C
ATOM    477  O   GLN A  57      -9.272  -2.778  10.290  1.00 58.06           O
ATOM    478  CB  GLN A  57     -10.191  -5.207  12.197  1.00 58.06           C
ATOM    479  CG  GLN A  57     -11.310  -6.152  12.634  1.00 58.06           C
ATOM    480  CD  GLN A  57     -11.043  -6.563  14.072  1.00 58.06           C
ATOM    481  NE2 GLN A  57     -12.072  -7.159  14.733  1.00 58.06           N
ATOM    482  OE1 GLN A  57      -9.945  -6.371  14.593  1.00 58.06           O
ATOM    483  N   GLY A  58      -7.860  -4.535  10.207  1.00 35.92           N
ATOM    484  CA  GLY A  58      -6.688  -3.745   9.975  1.00 35.92           C
ATOM    485  C   GLY A  58      -6.870  -2.899   8.755  1.00 35.92           C
ATOM    486  O   GLY A  58      -6.348  -1.787   8.687  1.00 35.92           O
ATOM    487  N   ALA A  59      -7.555  -3.430   7.726  1.00 37.75           N
ATOM    488  CA  ALA A  59      -7.773  -2.678   6.523  1.00 37.75           C
ATOM    489  C   ALA A  59      -8.654  -1.507   6.824  1.00 37.75           C
ATOM    490  O   ALA A  59      -8.415  -0.393   6.356  1.00 37.75           O
ATOM    491  CB  ALA A  59      -8.462  -3.502   5.420  1.00 37.75           C
ATOM    492  N   LEU A  60      -9.685  -1.725   7.657  1.00 64.01           N
ATOM    493  CA  LEU A  60     -10.628  -0.684   7.931  1.00 64.01           C
ATOM    494  C   LEU A  60      -9.871   0.451   8.542  1.00 64.01           C
ATOM    495  O   LEU A  60     -10.195   1.618   8.323  1.00 64.01           O
ATOM    496  CB  LEU A  60     -11.723  -1.115   8.917  1.00 64.01           C
ATOM    497  CG  LEU A  60     -12.495  -2.363   8.443  1.00 64.01           C
ATOM    498  CD1 LEU A  60     -13.780  -2.565   9.254  1.00 64.01           C
ATOM    499  CD2 LEU A  60     -12.726  -2.349   6.926  1.00 64.01           C
ATOM    500  N   ALA A  61      -8.845   0.133   9.351  1.00 36.23           N
ATOM    501  CA  ALA A  61      -8.065   1.143  10.003  1.00 36.23           C
ATOM    502  C   ALA A  61      -7.326   1.993   9.008  1.00 36.23           C
ATOM    503  O   ALA A  61      -7.304   3.218   9.130  1.00 36.23           O
ATOM    504  CB  ALA A  61      -7.019   0.552  10.963  1.00 36.23           C
ATOM    505  N   ASN A  62      -6.700   1.372   7.990  1.00 53.69           N
ATOM    506  CA  ASN A  62      -5.926   2.110   7.027  1.00 53.69           C
ATOM    507  C   ASN A  62      -6.827   3.028   6.264  1.00 53.69           C
ATOM    508  O   ASN A  62      -6.481   4.179   6.000  1.00 53.69           O
ATOM    509  CB  ASN A  62      -5.217   1.212   5.995  1.00 53.69           C
ATOM    510  CG  ASN A  62      -4.020   0.537   6.654  1.00 53.69           C
ATOM    511  ND2 ASN A  62      -3.501  -0.540   6.007  1.00 53.69           N
ATOM    512  OD1 ASN A  62      -3.546   0.968   7.704  1.00 53.69           O
ATOM    513  N   ILE A  63      -8.027   2.541   5.909  1.00133.11           N
ATOM    514  CA  ILE A  63      -8.929   3.305   5.102  1.00133.11           C
ATOM    515  C   ILE A  63      -9.314   4.541   5.859  1.00133.11           C
ATOM    516  O   ILE A  63      -9.476   5.610   5.274  1.00133.11           O
ATOM    517  CB  ILE A  63     -10.167   2.543   4.713  1.00133.11           C
ATOM    518  CG1 ILE A  63     -10.760   3.134   3.426  1.00133.11           C
ATOM    519  CG2 ILE A  63     -11.153   2.558   5.890  1.00133.11           C
ATOM    520  CD1 ILE A  63      -9.871   2.903   2.201  1.00133.11           C
ATOM    521  N   ALA A  64      -9.480   4.429   7.189  1.00 32.77           N
ATOM    522  CA  ALA A  64      -9.875   5.567   7.967  1.00 32.77           C
ATOM    523  C   ALA A  64      -8.825   6.625   7.832  1.00 32.77           C
ATOM    524  O   ALA A  64      -9.140   7.809   7.733  1.00 32.77           O
ATOM    525  CB  ALA A  64     -10.020   5.243   9.464  1.00 32.77           C
ATOM    526  N   VAL A  65      -7.539   6.225   7.853  1.00 38.33           N
ATOM    527  CA  VAL A  65      -6.453   7.154   7.715  1.00 38.33           C
ATOM    528  C   VAL A  65      -6.436   7.716   6.327  1.00 38.33           C
ATOM    529  O   VAL A  65      -6.126   8.888   6.125  1.00 38.33           O
ATOM    530  CB  VAL A  65      -5.118   6.523   7.981  1.00 38.33           C
ATOM    531  CG1 VAL A  65      -4.019   7.565   7.711  1.00 38.33           C
ATOM    532  CG2 VAL A  65      -5.122   5.982   9.421  1.00 38.33           C
ATOM    533  N   ASP A  66      -6.753   6.879   5.320  1.00 39.99           N
ATOM    534  CA  ASP A  66      -6.735   7.316   3.952  1.00 39.99           C
ATOM    535  C   ASP A  66      -7.777   8.374   3.782  1.00 39.99           C
ATOM    536  O   ASP A  66      -7.595   9.314   3.009  1.00 39.99           O
ATOM    537  CB  ASP A  66      -7.049   6.187   2.956  1.00 39.99           C
ATOM    538  CG  ASP A  66      -5.874   5.219   2.943  1.00 39.99           C
ATOM    539  OD1 ASP A  66      -4.711   5.679   3.098  1.00 39.99           O
ATOM    540  OD2 ASP A  66      -6.130   3.997   2.782  1.00 39.99           O1-
ATOM    541  N   LYS A  67      -8.911   8.243   4.494  1.00 62.81           N
ATOM    542  CA  LYS A  67      -9.955   9.216   4.352  1.00 62.81           C
ATOM    543  C   LYS A  67      -9.408  10.544   4.766  1.00 62.81           C
ATOM    544  O   LYS A  67      -9.596  11.547   4.080  1.00 62.81           O
ATOM    545  CB  LYS A  67     -11.167   8.954   5.262  1.00 62.81           C
ATOM    546  CG  LYS A  67     -11.939   7.675   4.944  1.00 62.81           C
ATOM    547  CD  LYS A  67     -12.964   7.317   6.021  1.00 62.81           C
ATOM    548  CE  LYS A  67     -13.747   6.041   5.719  1.00 62.81           C
ATOM    549  NZ  LYS A  67     -12.925   4.857   6.049  1.00 62.81           N1+
ATOM    550  N   ALA A  68      -8.695  10.576   5.906  1.00 29.76           N
ATOM    551  CA  ALA A  68      -8.159  11.800   6.428  1.00 29.76           C
ATOM    552  C   ALA A  68      -7.148  12.357   5.470  1.00 29.76           C
ATOM    553  O   ALA A  68      -7.089  13.567   5.258  1.00 29.76           O
ATOM    554  CB  ALA A  68      -7.457  11.608   7.782  1.00 29.76           C
ATOM    555  N   ASN A  69      -6.315  11.482   4.875  1.00 42.92           N
ATOM    556  CA  ASN A  69      -5.266  11.908   3.988  1.00 42.92           C
ATOM    557  C   ASN A  69      -5.832  12.529   2.748  1.00 42.92           C
ATOM    558  O   ASN A  69      -5.288  13.504   2.235  1.00 42.92           O
ATOM    559  CB  ASN A  69      -4.348  10.759   3.543  1.00 42.92           C
ATOM    560  CG  ASN A  69      -3.492  10.372   4.738  1.00 42.92           C
ATOM    561  ND2 ASN A  69      -3.185   9.053   4.865  1.00 42.92           N
ATOM    562  OD1 ASN A  69      -3.105  11.222   5.538  1.00 42.92           O
ATOM    563  N   LEU A  70      -6.942  11.976   2.228  1.00 45.91           N
ATOM    564  CA  LEU A  70      -7.517  12.475   1.011  1.00 45.91           C
ATOM    565  C   LEU A  70      -7.925  13.902   1.210  1.00 45.91           C
ATOM    566  O   LEU A  70      -7.714  14.739   0.334  1.00 45.91           O
ATOM    567  CB  LEU A  70      -8.772  11.693   0.590  1.00 45.91           C
ATOM    568  CG  LEU A  70      -9.456  12.246  -0.673  1.00 45.91           C
ATOM    569  CD1 LEU A  70      -8.517  12.197  -1.886  1.00 45.91           C
ATOM    570  CD2 LEU A  70     -10.797  11.541  -0.930  1.00 45.91           C
ATOM    571  N   GLU A  71      -8.523  14.212   2.374  1.00 33.07           N
ATOM    572  CA  GLU A  71      -8.981  15.546   2.639  1.00 33.07           C
ATOM    573  C   GLU A  71      -7.809  16.478   2.702  1.00 33.07           C
ATOM    574  O   GLU A  71      -7.867  17.594   2.191  1.00 33.07           O
ATOM    575  CB  GLU A  71      -9.734  15.668   3.973  1.00 33.07           C
ATOM    576  CG  GLU A  71     -11.066  14.915   3.984  1.00 33.07           C
ATOM    577  CD  GLU A  71     -11.708  15.110   5.350  1.00 33.07           C
ATOM    578  OE1 GLU A  71     -11.721  16.272   5.838  1.00 33.07           O
ATOM    579  OE2 GLU A  71     -12.196  14.100   5.923  1.00 33.07           O1-
ATOM    580  N   ILE A  72      -6.704  16.033   3.331  1.00 88.56           N
ATOM    581  CA  ILE A  72      -5.551  16.871   3.505  1.00 88.56           C
ATOM    582  C   ILE A  72      -5.006  17.224   2.155  1.00 88.56           C
ATOM    583  O   ILE A  72      -4.708  18.385   1.884  1.00 88.56           O
ATOM    584  CB  ILE A  72      -4.437  16.165   4.228  1.00 88.56           C
ATOM    585  CG1 ILE A  72      -4.869  15.720   5.636  1.00 88.56           C
ATOM    586  CG2 ILE A  72      -3.214  17.098   4.225  1.00 88.56           C
ATOM    587  CD1 ILE A  72      -5.195  16.878   6.576  1.00 88.56           C
ATOM    588  N   MET A  73      -4.857  16.215   1.274  1.00 55.78           N
ATOM    589  CA  MET A  73      -4.312  16.435  -0.036  1.00 55.78           C
ATOM    590  C   MET A  73      -5.221  17.275  -0.869  1.00 55.78           C
ATOM    591  O   MET A  73      -4.759  18.097  -1.657  1.00 55.78           O
ATOM    592  CB  MET A  73      -4.011  15.149  -0.822  1.00 55.78           C
ATOM    593  CG  MET A  73      -2.652  14.539  -0.479  1.00 55.78           C
ATOM    594  SD  MET A  73      -1.256  15.598  -0.974  1.00 55.78           S
ATOM    595  CE  MET A  73      -0.008  14.279  -0.993  1.00 55.78           C
ATOM    596  N   THR A  74      -6.544  17.092  -0.730  1.00 32.93           N
ATOM    597  CA  THR A  74      -7.436  17.840  -1.563  1.00 32.93           C
ATOM    598  C   THR A  74      -7.208  19.296  -1.300  1.00 32.93           C
ATOM    599  O   THR A  74      -7.088  20.092  -2.230  1.00 32.93           O
ATOM    600  CB  THR A  74      -8.874  17.529  -1.276  1.00 32.93           C
ATOM    601  CG2 THR A  74      -9.753  18.349  -2.236  1.00 32.93           C
ATOM    602  OG1 THR A  74      -9.119  16.143  -1.460  1.00 32.93           O
ATOM    603  N   LYS A  75      -7.126  19.673  -0.012  1.00110.73           N
ATOM    604  CA  LYS A  75      -6.938  21.038   0.391  1.00110.73           C
ATOM    605  C   LYS A  75      -5.590  21.522  -0.052  1.00110.73           C
ATOM    606  O   LYS A  75      -5.451  22.652  -0.518  1.00110.73           O
ATOM    607  CB  LYS A  75      -7.018  21.190   1.918  1.00110.73           C
ATOM    608  CG  LYS A  75      -6.922  22.633   2.403  1.00110.73           C
ATOM    609  CD  LYS A  75      -8.134  23.494   2.048  1.00110.73           C
ATOM    610  CE  LYS A  75      -8.047  24.908   2.622  1.00110.73           C
ATOM    611  NZ  LYS A  75      -9.220  25.704   2.197  1.00110.73           N1+
ATOM    612  N   ARG A  76      -4.557  20.666   0.069  1.00 66.61           N
ATOM    613  CA  ARG A  76      -3.220  21.066  -0.262  1.00 66.61           C
ATOM    614  C   ARG A  76      -3.199  21.447  -1.712  1.00 66.61           C
ATOM    615  O   ARG A  76      -2.570  22.432  -2.095  1.00 66.61           O
ATOM    616  CB  ARG A  76      -2.206  19.922  -0.050  1.00 66.61           C
ATOM    617  CG  ARG A  76      -0.737  20.335  -0.188  1.00 66.61           C
ATOM    618  CD  ARG A  76       0.252  19.181   0.025  1.00 66.61           C
ATOM    619  NE  ARG A  76       0.140  18.724   1.440  1.00 66.61           N
ATOM    620  CZ  ARG A  76       0.773  17.584   1.847  1.00 66.61           C
ATOM    621  NH1 ARG A  76       1.501  16.852   0.954  1.00 66.61           N1+
ATOM    622  NH2 ARG A  76       0.679  17.169   3.144  1.00 66.61           N
ATOM    623  N   SER A  77      -3.888  20.649  -2.551  1.00 51.16           N
ATOM    624  CA  SER A  77      -3.979  20.856  -3.971  1.00 51.16           C
ATOM    625  C   SER A  77      -4.968  21.940  -4.280  1.00 51.16           C
ATOM    626  O   SER A  77      -5.166  22.283  -5.447  1.00 51.16           O
ATOM    627  CB  SER A  77      -4.406  19.589  -4.735  1.00 51.16           C
ATOM    628  OG  SER A  77      -5.729  19.216  -4.377  1.00 51.16           O
ATOM    629  N   ASN A  78      -5.613  22.519  -3.251  1.00 55.10           N
ATOM    630  CA  ASN A  78      -6.599  23.536  -3.480  1.00 55.10           C
ATOM    631  C   ASN A  78      -7.708  22.994  -4.331  1.00 55.10           C
ATOM    632  O   ASN A  78      -8.069  23.567  -5.357  1.00 55.10           O
ATOM    633  CB  ASN A  78      -6.033  24.804  -4.142  1.00 55.10           C
ATOM    634  CG  ASN A  78      -5.169  25.512  -3.106  1.00 55.10           C
ATOM    635  ND2 ASN A  78      -4.179  26.313  -3.583  1.00 55.10           N
ATOM    636  OD1 ASN A  78      -5.376  25.362  -1.902  1.00 55.10           O
ATOM    637  N   TYR A  79      -8.272  21.853  -3.895  1.00 68.78           N
ATOM    638  CA  TYR A  79      -9.417  21.239  -4.503  1.00 68.78           C
ATOM    639  C   TYR A  79      -9.247  21.084  -5.978  1.00 68.78           C
ATOM    640  O   TYR A  79     -10.073  21.550  -6.762  1.00 68.78           O
ATOM    641  CB  TYR A  79     -10.724  21.995  -4.208  1.00 68.78           C
ATOM    642  CG  TYR A  79     -10.915  21.907  -2.732  1.00 68.78           C
ATOM    643  CD1 TYR A  79     -10.239  22.762  -1.893  1.00 68.78           C
ATOM    644  CD2 TYR A  79     -11.763  20.971  -2.184  1.00 68.78           C
ATOM    645  CE1 TYR A  79     -10.402  22.687  -0.531  1.00 68.78           C
ATOM    646  CE2 TYR A  79     -11.931  20.891  -0.821  1.00 68.78           C
ATOM    647  CZ  TYR A  79     -11.249  21.751   0.008  1.00 68.78           C
ATOM    648  OH  TYR A  79     -11.415  21.673   1.408  1.00 68.78           O
ATOM    649  N   THR A  80      -8.165  20.401  -6.399  1.00 38.67           N
ATOM    650  CA  THR A  80      -7.972  20.155  -7.797  1.00 38.67           C
ATOM    651  C   THR A  80      -8.526  18.784  -8.062  1.00 38.67           C
ATOM    652  O   THR A  80      -8.137  17.806  -7.427  1.00 38.67           O
ATOM    653  CB  THR A  80      -6.526  20.145  -8.205  1.00 38.67           C
ATOM    654  CG2 THR A  80      -6.435  19.882  -9.718  1.00 38.67           C
ATOM    655  OG1 THR A  80      -5.922  21.391  -7.890  1.00 38.67           O
ATOM    656  N   PRO A  81      -9.457  18.709  -8.977  1.00148.26           N
ATOM    657  CA  PRO A  81     -10.096  17.452  -9.283  1.00148.26           C
ATOM    658  C   PRO A  81      -9.259  16.523 -10.105  1.00148.26           C
ATOM    659  O   PRO A  81      -8.337  16.969 -10.782  1.00148.26           O
ATOM    660  CB  PRO A  81     -11.443  17.801  -9.925  1.00148.26           C
ATOM    661  CG  PRO A  81     -11.359  19.308 -10.227  1.00148.26           C
ATOM    662  CD  PRO A  81     -10.351  19.830  -9.196  1.00148.26           C
TER     663      PRO A  81
ATOM    664  N   GLY B   1     -21.276  18.637 -19.979  1.00 22.15           N
ATOM    665  CA  GLY B   1     -21.600  19.538 -18.850  1.00 22.15           C
ATOM    666  C   GLY B   1     -20.665  19.291 -17.719  1.00 22.15           C
ATOM    667  O   GLY B   1     -20.043  18.231 -17.619  1.00 22.15           O
ATOM    668  N   ASP B   2     -20.550  20.268 -16.800  1.00191.18           N
ATOM    669  CA  ASP B   2     -19.648  20.068 -15.705  1.00191.18           C
ATOM    670  C   ASP B   2     -20.239  19.034 -14.801  1.00191.18           C
ATOM    671  O   ASP B   2     -21.396  19.115 -14.387  1.00191.18           O
ATOM    672  CB  ASP B   2     -19.359  21.336 -14.883  1.00191.18           C
ATOM    673  CG  ASP B   2     -18.060  21.128 -14.109  1.00191.18           C
ATOM    674  OD1 ASP B   2     -17.480  20.012 -14.202  1.00191.18           O
ATOM    675  OD2 ASP B   2     -17.630  22.090 -13.415  1.00191.18           O1-
ATOM    676  N   THR B   3     -19.420  18.006 -14.508  1.00 85.87           N
ATOM    677  CA  THR B   3     -19.814  16.940 -13.639  1.00 85.87           C
ATOM    678  C   THR B   3     -19.288  17.287 -12.279  1.00 85.87           C
ATOM    679  O   THR B   3     -18.647  18.323 -12.117  1.00 85.87           O
ATOM    680  CB  THR B   3     -19.268  15.603 -14.076  1.00 85.87           C
ATOM    681  CG2 THR B   3     -17.729  15.618 -14.002  1.00 85.87           C
ATOM    682  OG1 THR B   3     -19.793  14.563 -13.263  1.00 85.87           O
ATOM    683  N   ARG B   4     -19.583  16.452 -11.255  1.00131.71           N
ATOM    684  CA  ARG B   4     -19.102  16.690  -9.912  1.00131.71           C
ATOM    685  C   ARG B   4     -17.712  16.104  -9.753  1.00131.71           C
ATOM    686  O   ARG B   4     -17.468  14.960 -10.132  1.00131.71           O
ATOM    687  CB  ARG B   4     -20.111  16.216  -8.827  1.00131.71           C
ATOM    688  CG  ARG B   4     -20.611  14.758  -8.857  1.00131.71           C
ATOM    689  CD  ARG B   4     -20.879  14.153 -10.245  1.00131.71           C
ATOM    690  NE  ARG B   4     -21.571  12.840 -10.081  1.00131.71           N
ATOM    691  CZ  ARG B   4     -21.281  11.860 -10.979  1.00131.71           C
ATOM    692  NH1 ARG B   4     -20.386  12.124 -11.979  1.00131.71           N1+
ATOM    693  NH2 ARG B   4     -21.834  10.619 -10.839  1.00131.71           N
ATOM    694  N   PRO B   5     -16.817  16.917  -9.199  1.00 70.70           N
ATOM    695  CA  PRO B   5     -15.412  16.625  -9.119  1.00 70.70           C
ATOM    696  C   PRO B   5     -15.078  15.507  -8.188  1.00 70.70           C
ATOM    697  O   PRO B   5     -15.734  15.354  -7.159  1.00 70.70           O
ATOM    698  CB  PRO B   5     -14.733  17.949  -8.781  1.00 70.70           C
ATOM    699  CG  PRO B   5     -15.677  19.002  -9.394  1.00 70.70           C
ATOM    700  CD  PRO B   5     -17.067  18.342  -9.352  1.00 70.70           C
ATOM    701  N   ARG B   6     -14.054  14.705  -8.543  1.00117.23           N
ATOM    702  CA  ARG B   6     -13.674  13.598  -7.718  1.00117.23           C
ATOM    703  C   ARG B   6     -12.258  13.799  -7.285  1.00117.23           C
ATOM    704  O   ARG B   6     -11.464  14.433  -7.979  1.00117.23           O
ATOM    705  CB  ARG B   6     -13.737  12.238  -8.432  1.00117.23           C
ATOM    706  CG  ARG B   6     -15.161  11.753  -8.707  1.00117.23           C
ATOM    707  CD  ARG B   6     -16.063  11.848  -7.474  1.00117.23           C
ATOM    708  NE  ARG B   6     -17.260  10.987  -7.689  1.00117.23           N
ATOM    709  CZ  ARG B   6     -18.410  11.480  -8.235  1.00117.23           C
ATOM    710  NH1 ARG B   6     -18.455  12.756  -8.713  1.00117.23           N1+
ATOM    711  NH2 ARG B   6     -19.521  10.689  -8.288  1.00117.23           N
ATOM    712  N   PHE B   7     -11.922  13.279  -6.088  1.00 59.16           N
ATOM    713  CA  PHE B   7     -10.586  13.388  -5.586  1.00 59.16           C
ATOM    714  C   PHE B   7     -10.185  11.996  -5.208  1.00 59.16           C
ATOM    715  O   PHE B   7     -10.953  11.282  -4.565  1.00 59.16           O
ATOM    716  CB  PHE B   7     -10.506  14.267  -4.329  1.00 59.16           C
ATOM    717  CG  PHE B   7     -11.117  15.574  -4.712  1.00 59.16           C
ATOM    718  CD1 PHE B   7     -12.478  15.751  -4.606  1.00 59.16           C
ATOM    719  CD2 PHE B   7     -10.347  16.610  -5.185  1.00 59.16           C
ATOM    720  CE1 PHE B   7     -13.064  16.944  -4.955  1.00 59.16           C
ATOM    721  CE2 PHE B   7     -10.928  17.807  -5.538  1.00 59.16           C
ATOM    722  CZ  PHE B   7     -12.287  17.976  -5.424  1.00 59.16           C
ATOM    723  N   LEU B   8      -8.970  11.558  -5.600  1.00 62.84           N
ATOM    724  CA  LEU B   8      -8.619  10.194  -5.319  1.00 62.84           C
ATOM    725  C   LEU B   8      -7.301  10.118  -4.612  1.00 62.84           C
ATOM    726  O   LEU B   8      -6.348  10.820  -4.951  1.00 62.84           O
ATOM    727  CB  LEU B   8      -8.517   9.339  -6.592  1.00 62.84           C
ATOM    728  CG  LEU B   8      -8.163   7.863  -6.348  1.00 62.84           C
ATOM    729  CD1 LEU B   8      -9.267   7.148  -5.552  1.00 62.84           C
ATOM    730  CD2 LEU B   8      -7.832   7.154  -7.671  1.00 62.84           C
ATOM    731  N   TRP B   9      -7.235   9.243  -3.586  1.00 66.84           N
ATOM    732  CA  TRP B   9      -6.029   8.997  -2.843  1.00 66.84           C
ATOM    733  C   TRP B   9      -5.710   7.542  -3.006  1.00 66.84           C
ATOM    734  O   TRP B   9      -6.567   6.687  -2.785  1.00 66.84           O
ATOM    735  CB  TRP B   9      -6.165   9.287  -1.338  1.00 66.84           C
ATOM    736  CG  TRP B   9      -4.953   8.901  -0.523  1.00 66.84           C
ATOM    737  CD1 TRP B   9      -4.753   7.793   0.246  1.00 66.84           C
ATOM    738  CD2 TRP B   9      -3.740   9.668  -0.453  1.00 66.84           C
ATOM    739  CE2 TRP B   9      -2.860   8.970   0.373  1.00 66.84           C
ATOM    740  CE3 TRP B   9      -3.390  10.853  -1.032  1.00 66.84           C
ATOM    741  NE1 TRP B   9      -3.495   7.822   0.797  1.00 66.84           N
ATOM    742  CZ2 TRP B   9      -1.609   9.451   0.634  1.00 66.84           C
ATOM    743  CZ3 TRP B   9      -2.128  11.336  -0.764  1.00 66.84           C
ATOM    744  CH2 TRP B   9      -1.255  10.648   0.052  1.00 66.84           C
ATOM    745  N   GLN B  10      -4.460   7.227  -3.410  1.00 48.95           N
ATOM    746  CA  GLN B  10      -4.073   5.859  -3.624  1.00 48.95           C
ATOM    747  C   GLN B  10      -2.793   5.609  -2.889  1.00 48.95           C
ATOM    748  O   GLN B  10      -1.908   6.462  -2.847  1.00 48.95           O
ATOM    749  CB  GLN B  10      -3.815   5.528  -5.106  1.00 48.95           C
ATOM    750  CG  GLN B  10      -5.070   5.620  -5.978  1.00 48.95           C
ATOM    751  CD  GLN B  10      -4.675   5.369  -7.427  1.00 48.95           C
ATOM    752  NE2 GLN B  10      -4.506   4.073  -7.800  1.00 48.95           N
ATOM    753  OE1 GLN B  10      -4.515   6.304  -8.209  1.00 48.95           O
ATOM    754  N   LEU B  11      -2.673   4.416  -2.272  1.00 61.31           N
ATOM    755  CA  LEU B  11      -1.477   4.108  -1.546  1.00 61.31           C
ATOM    756  C   LEU B  11      -1.052   2.729  -1.939  1.00 61.31           C
ATOM    757  O   LEU B  11      -1.873   1.814  -2.011  1.00 61.31           O
ATOM    758  CB  LEU B  11      -1.694   4.141  -0.022  1.00 61.31           C
ATOM    759  CG  LEU B  11      -0.405   4.029   0.811  1.00 61.31           C
ATOM    760  CD1 LEU B  11       0.516   5.233   0.555  1.00 61.31           C
ATOM    761  CD2 LEU B  11      -0.725   3.842   2.304  1.00 61.31           C
ATOM    762  N   LYS B  12       0.249   2.543  -2.237  1.00139.59           N
ATOM    763  CA  LYS B  12       0.671   1.224  -2.596  1.00139.59           C
ATOM    764  C   LYS B  12       1.987   0.910  -1.965  1.00139.59           C
ATOM    765  O   LYS B  12       2.892   1.742  -1.916  1.00139.59           O
ATOM    766  CB  LYS B  12       0.723   0.987  -4.117  1.00139.59           C
ATOM    767  CG  LYS B  12       1.449   2.068  -4.910  1.00139.59           C
ATOM    768  CD  LYS B  12       1.498   1.758  -6.407  1.00139.59           C
ATOM    769  CE  LYS B  12       0.108   1.628  -7.040  1.00139.59           C
ATOM    770  NZ  LYS B  12       0.231   1.336  -8.486  1.00139.59           N1+
ATOM    771  N   PHE B  13       2.089  -0.317  -1.413  1.00 56.29           N
ATOM    772  CA  PHE B  13       3.310  -0.799  -0.839  1.00 56.29           C
ATOM    773  C   PHE B  13       3.738  -1.902  -1.745  1.00 56.29           C
ATOM    774  O   PHE B  13       3.030  -2.898  -1.893  1.00 56.29           O
ATOM    775  CB  PHE B  13       3.140  -1.419   0.560  1.00 56.29           C
ATOM    776  CG  PHE B  13       2.803  -0.326   1.513  1.00 56.29           C
ATOM    777  CD1 PHE B  13       3.800   0.451   2.057  1.00 56.29           C
ATOM    778  CD2 PHE B  13       1.497  -0.080   1.869  1.00 56.29           C
ATOM    779  CE1 PHE B  13       3.499   1.461   2.939  1.00 56.29           C
ATOM    780  CE2 PHE B  13       1.190   0.930   2.750  1.00 56.29           C
ATOM    781  CZ  PHE B  13       2.192   1.703   3.287  1.00 56.29           C
ATOM    782  N   GLU B  14       4.916  -1.760  -2.378  1.00 43.85           N
ATOM    783  CA  GLU B  14       5.314  -2.770  -3.310  1.00 43.85           C
ATOM    784  C   GLU B  14       6.528  -3.469  -2.786  1.00 43.85           C
ATOM    785  O   GLU B  14       7.521  -2.838  -2.430  1.00 43.85           O
ATOM    786  CB  GLU B  14       5.668  -2.203  -4.695  1.00 43.85           C
ATOM    787  CG  GLU B  14       4.470  -1.571  -5.408  1.00 43.85           C
ATOM    788  CD  GLU B  14       4.946  -1.021  -6.745  1.00 43.85           C
ATOM    789  OE1 GLU B  14       6.159  -1.171  -7.052  1.00 43.85           O
ATOM    790  OE2 GLU B  14       4.103  -0.438  -7.477  1.00 43.85           O1-
ATOM    791  N   CYS B  15       6.468  -4.815  -2.721  1.00 41.54           N
ATOM    792  CA  CYS B  15       7.605  -5.571  -2.284  1.00 41.54           C
ATOM    793  C   CYS B  15       8.129  -6.287  -3.487  1.00 41.54           C
ATOM    794  O   CYS B  15       7.415  -7.075  -4.108  1.00 41.54           O
ATOM    795  CB  CYS B  15       7.283  -6.653  -1.232  1.00 41.54           C
ATOM    796  SG  CYS B  15       6.775  -5.970   0.374  1.00 41.54           S
ATOM    797  N   HIS B  16       9.399  -6.026  -3.851  1.00 57.99           N
ATOM    798  CA  HIS B  16       9.959  -6.679  -5.001  1.00 57.99           C
ATOM    799  C   HIS B  16      10.948  -7.695  -4.522  1.00 57.99           C
ATOM    800  O   HIS B  16      11.849  -7.386  -3.743  1.00 57.99           O
ATOM    801  CB  HIS B  16      10.677  -5.717  -5.961  1.00 57.99           C
ATOM    802  CG  HIS B  16       9.721  -4.773  -6.625  1.00 57.99           C
ATOM    803  CD2 HIS B  16       9.237  -3.574  -6.199  1.00 57.99           C
ATOM    804  ND1 HIS B  16       9.133  -5.003  -7.849  1.00 57.99           N
ATOM    805  CE1 HIS B  16       8.326  -3.941  -8.101  1.00 57.99           C
ATOM    806  NE2 HIS B  16       8.356  -3.048  -7.127  1.00 57.99           N
ATOM    807  N   PHE B  17      10.793  -8.951  -4.989  1.00 60.78           N
ATOM    808  CA  PHE B  17      11.660 -10.014  -4.559  1.00 60.78           C
ATOM    809  C   PHE B  17      12.481 -10.448  -5.729  1.00 60.78           C
ATOM    810  O   PHE B  17      11.965 -10.626  -6.831  1.00 60.78           O
ATOM    811  CB  PHE B  17      10.900 -11.277  -4.110  1.00 60.78           C
ATOM    812  CG  PHE B  17      10.057 -10.955  -2.924  1.00 60.78           C
ATOM    813  CD1 PHE B  17       8.858 -10.298  -3.082  1.00 60.78           C
ATOM    814  CD2 PHE B  17      10.448 -11.327  -1.658  1.00 60.78           C
ATOM    815  CE1 PHE B  17       8.069 -10.003  -1.995  1.00 60.78           C
ATOM    816  CE2 PHE B  17       9.663 -11.035  -0.566  1.00 60.78           C
ATOM    817  CZ  PHE B  17       8.473 -10.370  -0.734  1.00 60.78           C
ATOM    818  N   PHE B  18      13.798 -10.627  -5.509  1.00 70.74           N
ATOM    819  CA  PHE B  18      14.667 -11.090  -6.547  1.00 70.74           C
ATOM    820  C   PHE B  18      15.313 -12.327  -6.012  1.00 70.74           C
ATOM    821  O   PHE B  18      15.880 -12.305  -4.920  1.00 70.74           O
ATOM    822  CB  PHE B  18      15.802 -10.092  -6.846  1.00 70.74           C
ATOM    823  CG  PHE B  18      15.171  -8.805  -7.259  1.00 70.74           C
ATOM    824  CD1 PHE B  18      14.683  -7.933  -6.313  1.00 70.74           C
ATOM    825  CD2 PHE B  18      15.071  -8.462  -8.586  1.00 70.74           C
ATOM    826  CE1 PHE B  18      14.101  -6.744  -6.687  1.00 70.74           C
ATOM    827  CE2 PHE B  18      14.491  -7.275  -8.969  1.00 70.74           C
ATOM    828  CZ  PHE B  18      14.002  -6.413  -8.016  1.00 70.74           C
ATOM    829  N   ASN B  19      15.270 -13.435  -6.781  1.00 65.93           N
ATOM    830  CA  ASN B  19      15.884 -14.655  -6.329  1.00 65.93           C
ATOM    831  C   ASN B  19      15.398 -14.995  -4.950  1.00 65.93           C
ATOM    832  O   ASN B  19      16.168 -15.006  -3.990  1.00 65.93           O
ATOM    833  CB  ASN B  19      17.421 -14.557  -6.296  1.00 65.93           C
ATOM    834  CG  ASN B  19      18.001 -15.885  -5.827  1.00 65.93           C
ATOM    835  ND2 ASN B  19      19.046 -15.810  -4.960  1.00 65.93           N
ATOM    836  OD1 ASN B  19      17.527 -16.955  -6.205  1.00 65.93           O
ATOM    837  N   GLY B  20      14.085 -15.257  -4.807  1.00 44.74           N
ATOM    838  CA  GLY B  20      13.586 -15.607  -3.509  1.00 44.74           C
ATOM    839  C   GLY B  20      13.618 -14.372  -2.659  1.00 44.74           C
ATOM    840  O   GLY B  20      13.191 -13.299  -3.074  1.00 44.74           O
ATOM    841  N   THR B  21      14.018 -14.555  -1.390  1.00102.38           N
ATOM    842  CA  THR B  21      14.181 -13.559  -0.367  1.00102.38           C
ATOM    843  C   THR B  21      15.541 -12.915  -0.370  1.00102.38           C
ATOM    844  O   THR B  21      15.778 -12.005   0.421  1.00102.38           O
ATOM    845  CB  THR B  21      13.896 -14.077   1.011  1.00102.38           C
ATOM    846  CG2 THR B  21      12.435 -14.563   1.039  1.00102.38           C
ATOM    847  OG1 THR B  21      14.781 -15.136   1.344  1.00102.38           O
ATOM    848  N   GLU B  22      16.499 -13.408  -1.180  1.00103.69           N
ATOM    849  CA  GLU B  22      17.851 -12.921  -1.091  1.00103.69           C
ATOM    850  C   GLU B  22      17.896 -11.429  -1.219  1.00103.69           C
ATOM    851  O   GLU B  22      18.579 -10.767  -0.441  1.00103.69           O
ATOM    852  CB  GLU B  22      18.767 -13.464  -2.200  1.00103.69           C
ATOM    853  CG  GLU B  22      20.189 -12.914  -2.084  1.00103.69           C
ATOM    854  CD  GLU B  22      20.990 -13.353  -3.300  1.00103.69           C
ATOM    855  OE1 GLU B  22      21.144 -14.588  -3.499  1.00103.69           O
ATOM    856  OE2 GLU B  22      21.461 -12.454  -4.045  1.00103.69           O1-
ATOM    857  N   ARG B  23      17.186 -10.853  -2.205  1.00 85.52           N
ATOM    858  CA  ARG B  23      17.198  -9.425  -2.330  1.00 85.52           C
ATOM    859  C   ARG B  23      15.775  -8.963  -2.280  1.00 85.52           C
ATOM    860  O   ARG B  23      14.903  -9.575  -2.896  1.00 85.52           O
ATOM    861  CB  ARG B  23      17.796  -8.936  -3.658  1.00 85.52           C
ATOM    862  CG  ARG B  23      17.803  -7.416  -3.783  1.00 85.52           C
ATOM    863  CD  ARG B  23      18.189  -6.911  -5.172  1.00 85.52           C
ATOM    864  NE  ARG B  23      17.856  -5.461  -5.197  1.00 85.52           N
ATOM    865  CZ  ARG B  23      17.514  -4.852  -6.367  1.00 85.52           C
ATOM    866  NH1 ARG B  23      17.536  -5.554  -7.538  1.00 85.52           N1+
ATOM    867  NH2 ARG B  23      17.135  -3.542  -6.362  1.00 85.52           N
ATOM    868  N   VAL B  24      15.492  -7.878  -1.527  1.00 62.02           N
ATOM    869  CA  VAL B  24      14.129  -7.431  -1.445  1.00 62.02           C
ATOM    870  C   VAL B  24      14.128  -5.937  -1.404  1.00 62.02           C
ATOM    871  O   VAL B  24      14.925  -5.334  -0.686  1.00 62.02           O
ATOM    872  CB  VAL B  24      13.445  -7.890  -0.189  1.00 62.02           C
ATOM    873  CG1 VAL B  24      12.006  -7.343  -0.178  1.00 62.02           C
ATOM    874  CG2 VAL B  24      13.538  -9.422  -0.102  1.00 62.02           C
ATOM    875  N   ARG B  25      13.234  -5.287  -2.179  1.00208.86           N
ATOM    876  CA  ARG B  25      13.185  -3.865  -2.045  1.00208.86           C
ATOM    877  C   ARG B  25      11.753  -3.479  -1.840  1.00208.86           C
ATOM    878  O   ARG B  25      10.850  -3.999  -2.496  1.00208.86           O
ATOM    879  CB  ARG B  25      13.818  -3.064  -3.205  1.00208.86           C
ATOM    880  CG  ARG B  25      13.144  -3.158  -4.568  1.00208.86           C
ATOM    881  CD  ARG B  25      12.066  -2.100  -4.806  1.00208.86           C
ATOM    882  NE  ARG B  25      12.770  -0.822  -5.110  1.00208.86           N
ATOM    883  CZ  ARG B  25      12.371  -0.053  -6.166  1.00208.86           C
ATOM    884  NH1 ARG B  25      11.326  -0.452  -6.947  1.00208.86           N1+
ATOM    885  NH2 ARG B  25      13.024   1.111  -6.450  1.00208.86           N
ATOM    886  N   LEU B  26      11.515  -2.566  -0.875  1.00 51.20           N
ATOM    887  CA  LEU B  26      10.186  -2.147  -0.531  1.00 51.20           C
ATOM    888  C   LEU B  26       9.982  -0.752  -1.024  1.00 51.20           C
ATOM    889  O   LEU B  26      10.837   0.114  -0.842  1.00 51.20           O
ATOM    890  CB  LEU B  26       9.919  -2.125   0.984  1.00 51.20           C
ATOM    891  CG  LEU B  26       8.502  -1.645   1.349  1.00 51.20           C
ATOM    892  CD1 LEU B  26       7.427  -2.603   0.816  1.00 51.20           C
ATOM    893  CD2 LEU B  26       8.375  -1.379   2.859  1.00 51.20           C
ATOM    894  N   LEU B  27       8.826  -0.507  -1.676  1.00 63.95           N
ATOM    895  CA  LEU B  27       8.550   0.800  -2.201  1.00 63.95           C
ATOM    896  C   LEU B  27       7.190   1.219  -1.714  1.00 63.95           C
ATOM    897  O   LEU B  27       6.186   0.573  -2.008  1.00 63.95           O
ATOM    898  CB  LEU B  27       8.526   0.789  -3.741  1.00 63.95           C
ATOM    899  CG  LEU B  27       8.419   2.172  -4.403  1.00 63.95           C
ATOM    900  CD1 LEU B  27       9.673   3.014  -4.120  1.00 63.95           C
ATOM    901  CD2 LEU B  27       8.123   2.039  -5.906  1.00 63.95           C
ATOM    902  N   GLU B  28       7.121   2.329  -0.952  1.00 55.58           N
ATOM    903  CA  GLU B  28       5.874   2.823  -0.433  1.00 55.58           C
ATOM    904  C   GLU B  28       5.501   4.013  -1.259  1.00 55.58           C
ATOM    905  O   GLU B  28       6.286   4.954  -1.354  1.00 55.58           O
ATOM    906  CB  GLU B  28       6.003   3.316   1.021  1.00 55.58           C
ATOM    907  CG  GLU B  28       4.709   3.856   1.633  1.00 55.58           C
ATOM    908  CD  GLU B  28       5.036   4.384   3.025  1.00 55.58           C
ATOM    909  OE1 GLU B  28       5.951   3.813   3.678  1.00 55.58           O
ATOM    910  OE2 GLU B  28       4.381   5.371   3.453  1.00 55.58           O1-
ATOM    911  N   ARG B  29       4.305   4.022  -1.895  1.00206.40           N
ATOM    912  CA  ARG B  29       4.010   5.204  -2.660  1.00206.40           C
ATOM    913  C   ARG B  29       2.647   5.774  -2.420  1.00206.40           C
ATOM    914  O   ARG B  29       1.653   5.063  -2.283  1.00206.40           O
ATOM    915  CB  ARG B  29       4.358   5.138  -4.161  1.00206.40           C
ATOM    916  CG  ARG B  29       3.812   3.963  -4.958  1.00206.40           C
ATOM    917  CD  ARG B  29       4.897   2.938  -5.288  1.00206.40           C
ATOM    918  NE  ARG B  29       4.449   2.205  -6.504  1.00206.40           N
ATOM    919  CZ  ARG B  29       4.686   2.755  -7.730  1.00206.40           C
ATOM    920  NH1 ARG B  29       5.323   3.961  -7.823  1.00206.40           N1+
ATOM    921  NH2 ARG B  29       4.280   2.107  -8.861  1.00206.40           N
ATOM    922  N   CYS B  30       2.605   7.126  -2.341  1.00 41.78           N
ATOM    923  CA  CYS B  30       1.402   7.874  -2.110  1.00 41.78           C
ATOM    924  C   CYS B  30       1.021   8.510  -3.414  1.00 41.78           C
ATOM    925  O   CYS B  30       1.849   9.160  -4.055  1.00 41.78           O
ATOM    926  CB  CYS B  30       1.581   9.020  -1.099  1.00 41.78           C
ATOM    927  SG  CYS B  30       2.029   8.432   0.562  1.00 41.78           S
ATOM    928  N   ILE B  31      -0.247   8.335  -3.847  1.00 89.08           N
ATOM    929  CA  ILE B  31      -0.647   8.884  -5.116  1.00 89.08           C
ATOM    930  C   ILE B  31      -1.897   9.700  -4.953  1.00 89.08           C
ATOM    931  O   ILE B  31      -2.899   9.228  -4.418  1.00 89.08           O
ATOM    932  CB  ILE B  31      -0.946   7.830  -6.146  1.00 89.08           C
ATOM    933  CG1 ILE B  31       0.291   6.953  -6.396  1.00 89.08           C
ATOM    934  CG2 ILE B  31      -1.456   8.533  -7.414  1.00 89.08           C
ATOM    935  CD1 ILE B  31       0.636   6.032  -5.227  1.00 89.08           C
ATOM    936  N   TYR B  32      -1.868  10.956  -5.449  1.00 62.13           N
ATOM    937  CA  TYR B  32      -3.006  11.831  -5.390  1.00 62.13           C
ATOM    938  C   TYR B  32      -3.424  12.097  -6.807  1.00 62.13           C
ATOM    939  O   TYR B  32      -2.632  12.571  -7.619  1.00 62.13           O
ATOM    940  CB  TYR B  32      -2.662  13.181  -4.727  1.00 62.13           C
ATOM    941  CG  TYR B  32      -3.846  14.086  -4.762  1.00 62.13           C
ATOM    942  CD1 TYR B  32      -4.896  13.909  -3.892  1.00 62.13           C
ATOM    943  CD2 TYR B  32      -3.892  15.129  -5.658  1.00 62.13           C
ATOM    944  CE1 TYR B  32      -5.982  14.753  -3.925  1.00 62.13           C
ATOM    945  CE2 TYR B  32      -4.972  15.978  -5.697  1.00 62.13           C
ATOM    946  CZ  TYR B  32      -6.021  15.789  -4.829  1.00 62.13           C
ATOM    947  OH  TYR B  32      -7.133  16.656  -4.865  1.00 62.13           O
ATOM    948  N   ASN B  33      -4.688  11.772  -7.148  1.00 67.70           N
ATOM    949  CA  ASN B  33      -5.215  12.010  -8.470  1.00 67.70           C
ATOM    950  C   ASN B  33      -4.274  11.467  -9.511  1.00 67.70           C
ATOM    951  O   ASN B  33      -3.964  12.131 -10.498  1.00 67.70           O
ATOM    952  CB  ASN B  33      -5.447  13.503  -8.769  1.00 67.70           C
ATOM    953  CG  ASN B  33      -6.559  14.014  -7.863  1.00 67.70           C
ATOM    954  ND2 ASN B  33      -6.777  15.354  -7.868  1.00 67.70           N
ATOM    955  OD1 ASN B  33      -7.224  13.239  -7.177  1.00 67.70           O
ATOM    956  N   GLN B  34      -3.814  10.224  -9.299  1.00138.96           N
ATOM    957  CA  GLN B  34      -2.936   9.418 -10.106  1.00138.96           C
ATOM    958  C   GLN B  34      -1.584  10.038 -10.284  1.00138.96           C
ATOM    959  O   GLN B  34      -0.842   9.658 -11.189  1.00138.96           O
ATOM    960  CB  GLN B  34      -3.489   8.970 -11.476  1.00138.96           C
ATOM    961  CG  GLN B  34      -3.679  10.053 -12.532  1.00138.96           C
ATOM    962  CD  GLN B  34      -4.132   9.349 -13.806  1.00138.96           C
ATOM    963  NE2 GLN B  34      -3.255   8.468 -14.358  1.00138.96           N
ATOM    964  OE1 GLN B  34      -5.240   9.569 -14.291  1.00138.96           O
ATOM    965  N   GLU B  35      -1.192  10.969  -9.393  1.00 49.32           N
ATOM    966  CA  GLU B  35       0.123  11.538  -9.484  1.00 49.32           C
ATOM    967  C   GLU B  35       0.836  11.181  -8.212  1.00 49.32           C
ATOM    968  O   GLU B  35       0.407  11.566  -7.125  1.00 49.32           O
ATOM    969  CB  GLU B  35       0.090  13.075  -9.581  1.00 49.32           C
ATOM    970  CG  GLU B  35       1.464  13.738  -9.682  1.00 49.32           C
ATOM    971  CD  GLU B  35       1.234  15.240  -9.770  1.00 49.32           C
ATOM    972  OE1 GLU B  35       0.638  15.687 -10.785  1.00 49.32           O
ATOM    973  OE2 GLU B  35       1.642  15.961  -8.820  1.00 49.32           O1-
ATOM    974  N   GLU B  36       1.954  10.432  -8.313  1.00 73.33           N
ATOM    975  CA  GLU B  36       2.680   9.997  -7.146  1.00 73.33           C
ATOM    976  C   GLU B  36       3.368  11.189  -6.543  1.00 73.33           C
ATOM    977  O   GLU B  36       4.181  11.838  -7.198  1.00 73.33           O
ATOM    978  CB  GLU B  36       3.754   8.947  -7.489  1.00 73.33           C
ATOM    979  CG  GLU B  36       4.449   8.309  -6.283  1.00 73.33           C
ATOM    980  CD  GLU B  36       5.467   7.303  -6.808  1.00 73.33           C
ATOM    981  OE1 GLU B  36       5.786   7.363  -8.024  1.00 73.33           O
ATOM    982  OE2 GLU B  36       5.942   6.461  -5.999  1.00 73.33           O1-
ATOM    983  N   SER B  37       2.994  11.543  -5.293  1.00 73.54           N
ATOM    984  CA  SER B  37       3.559  12.657  -4.576  1.00 73.54           C
ATOM    985  C   SER B  37       4.830  12.343  -3.825  1.00 73.54           C
ATOM    986  O   SER B  37       5.806  13.083  -3.937  1.00 73.54           O
ATOM    987  CB  SER B  37       2.563  13.271  -3.576  1.00 73.54           C
ATOM    988  OG  SER B  37       2.183  12.307  -2.605  1.00 73.54           O
ATOM    989  N   VAL B  38       4.855  11.258  -3.012  1.00137.16           N
ATOM    990  CA  VAL B  38       6.031  10.993  -2.213  1.00137.16           C
ATOM    991  C   VAL B  38       6.220   9.513  -2.086  1.00137.16           C
ATOM    992  O   VAL B  38       5.264   8.745  -2.186  1.00137.16           O
ATOM    993  CB  VAL B  38       5.932  11.411  -0.773  1.00137.16           C
ATOM    994  CG1 VAL B  38       5.567  12.892  -0.675  1.00137.16           C
ATOM    995  CG2 VAL B  38       4.964  10.459  -0.054  1.00137.16           C
ATOM    996  N   ARG B  39       7.473   9.062  -1.848  1.00121.82           N
ATOM    997  CA  ARG B  39       7.620   7.653  -1.627  1.00121.82           C
ATOM    998  C   ARG B  39       8.840   7.335  -0.823  1.00121.82           C
ATOM    999  O   ARG B  39       9.770   8.133  -0.704  1.00121.82           O
ATOM   1000  CB  ARG B  39       7.612   6.787  -2.895  1.00121.82           C
ATOM   1001  CG  ARG B  39       8.711   7.047  -3.918  1.00121.82           C
ATOM   1002  CD  ARG B  39       8.544   6.100  -5.105  1.00121.82           C
ATOM   1003  NE  ARG B  39       9.614   6.376  -6.096  1.00121.82           N
ATOM   1004  CZ  ARG B  39       9.493   5.871  -7.356  1.00121.82           C
ATOM   1005  NH1 ARG B  39       8.360   5.201  -7.717  1.00121.82           N1+
ATOM   1006  NH2 ARG B  39      10.513   6.020  -8.250  1.00121.82           N
ATOM   1007  N   PHE B  40       8.820   6.134  -0.208  1.00 60.62           N
ATOM   1008  CA  PHE B  40       9.910   5.617   0.568  1.00 60.62           C
ATOM   1009  C   PHE B  40      10.435   4.438  -0.185  1.00 60.62           C
ATOM   1010  O   PHE B  40       9.700   3.492  -0.468  1.00 60.62           O
ATOM   1011  CB  PHE B  40       9.465   5.123   1.959  1.00 60.62           C
ATOM   1012  CG  PHE B  40      10.632   4.523   2.672  1.00 60.62           C
ATOM   1013  CD1 PHE B  40      11.537   5.312   3.343  1.00 60.62           C
ATOM   1014  CD2 PHE B  40      10.809   3.157   2.677  1.00 60.62           C
ATOM   1015  CE1 PHE B  40      12.605   4.749   4.001  1.00 60.62           C
ATOM   1016  CE2 PHE B  40      11.875   2.589   3.333  1.00 60.62           C
ATOM   1017  CZ  PHE B  40      12.777   3.385   3.996  1.00 60.62           C
ATOM   1018  N   ASP B  41      11.731   4.480  -0.549  1.00 57.38           N
ATOM   1019  CA  ASP B  41      12.350   3.396  -1.256  1.00 57.38           C
ATOM   1020  C   ASP B  41      13.341   2.788  -0.307  1.00 57.38           C
ATOM   1021  O   ASP B  41      14.208   3.474   0.228  1.00 57.38           O
ATOM   1022  CB  ASP B  41      13.102   3.854  -2.520  1.00 57.38           C
ATOM   1023  CG  ASP B  41      13.462   2.637  -3.362  1.00 57.38           C
ATOM   1024  OD1 ASP B  41      13.424   1.500  -2.820  1.00 57.38           O
ATOM   1025  OD2 ASP B  41      13.776   2.834  -4.566  1.00 57.38           O1-
ATOM   1026  N   SER B  42      13.249   1.462  -0.095  1.00 89.04           N
ATOM   1027  CA  SER B  42      14.075   0.785   0.867  1.00 89.04           C
ATOM   1028  C   SER B  42      15.515   0.983   0.519  1.00 89.04           C
ATOM   1029  O   SER B  42      16.375   1.028   1.396  1.00 89.04           O
ATOM   1030  CB  SER B  42      13.850  -0.737   0.870  1.00 89.04           C
ATOM   1031  OG  SER B  42      12.499  -1.034   1.173  1.00 89.04           O
ATOM   1032  N   ASP B  43      15.827   1.079  -0.781  1.00 52.56           N
ATOM   1033  CA  ASP B  43      17.199   1.232  -1.162  1.00 52.56           C
ATOM   1034  C   ASP B  43      17.701   2.556  -0.659  1.00 52.56           C
ATOM   1035  O   ASP B  43      18.823   2.644  -0.160  1.00 52.56           O
ATOM   1036  CB  ASP B  43      17.378   1.186  -2.690  1.00 52.56           C
ATOM   1037  CG  ASP B  43      17.017  -0.225  -3.144  1.00 52.56           C
ATOM   1038  OD1 ASP B  43      17.013  -1.138  -2.276  1.00 52.56           O
ATOM   1039  OD2 ASP B  43      16.739  -0.409  -4.360  1.00 52.56           O1-
ATOM   1040  N   VAL B  44      16.880   3.619  -0.802  1.00 43.12           N
ATOM   1041  CA  VAL B  44      17.248   4.955  -0.408  1.00 43.12           C
ATOM   1042  C   VAL B  44      17.351   5.074   1.085  1.00 43.12           C
ATOM   1043  O   VAL B  44      18.330   5.615   1.597  1.00 43.12           O
ATOM   1044  CB  VAL B  44      16.266   5.983  -0.880  1.00 43.12           C
ATOM   1045  CG1 VAL B  44      16.729   7.367  -0.393  1.00 43.12           C
ATOM   1046  CG2 VAL B  44      16.141   5.866  -2.409  1.00 43.12           C
ATOM   1047  N   GLY B  45      16.349   4.565   1.830  1.00 33.97           N
ATOM   1048  CA  GLY B  45      16.428   4.628   3.263  1.00 33.97           C
ATOM   1049  C   GLY B  45      15.524   5.687   3.840  1.00 33.97           C
ATOM   1050  O   GLY B  45      15.350   5.729   5.059  1.00 33.97           O
ATOM   1051  N   GLU B  46      14.936   6.583   3.020  1.00 60.20           N
ATOM   1052  CA  GLU B  46      14.057   7.570   3.594  1.00 60.20           C
ATOM   1053  C   GLU B  46      13.117   8.070   2.531  1.00 60.20           C
ATOM   1054  O   GLU B  46      13.282   7.775   1.350  1.00 60.20           O
ATOM   1055  CB  GLU B  46      14.793   8.774   4.212  1.00 60.20           C
ATOM   1056  CG  GLU B  46      15.657   9.562   3.229  1.00 60.20           C
ATOM   1057  CD  GLU B  46      16.301  10.702   4.008  1.00 60.20           C
ATOM   1058  OE1 GLU B  46      16.169  10.713   5.261  1.00 60.20           O
ATOM   1059  OE2 GLU B  46      16.934  11.577   3.360  1.00 60.20           O1-
ATOM   1060  N   TYR B  47      12.088   8.844   2.946  1.00 72.87           N
ATOM   1061  CA  TYR B  47      11.069   9.363   2.070  1.00 72.87           C
ATOM   1062  C   TYR B  47      11.626  10.460   1.218  1.00 72.87           C
ATOM   1063  O   TYR B  47      12.475  11.234   1.656  1.00 72.87           O
ATOM   1064  CB  TYR B  47       9.862   9.970   2.809  1.00 72.87           C
ATOM   1065  CG  TYR B  47       9.076   8.881   3.452  1.00 72.87           C
ATOM   1066  CD1 TYR B  47       8.069   8.251   2.758  1.00 72.87           C
ATOM   1067  CD2 TYR B  47       9.342   8.492   4.745  1.00 72.87           C
ATOM   1068  CE1 TYR B  47       7.336   7.246   3.343  1.00 72.87           C
ATOM   1069  CE2 TYR B  47       8.611   7.488   5.336  1.00 72.87           C
ATOM   1070  CZ  TYR B  47       7.606   6.864   4.634  1.00 72.87           C
ATOM   1071  OH  TYR B  47       6.852   5.833   5.235  1.00 72.87           O
ATOM   1072  N   ARG B  48      11.174  10.514  -0.054  1.00116.50           N
ATOM   1073  CA  ARG B  48      11.585  11.546  -0.967  1.00116.50           C
ATOM   1074  C   ARG B  48      10.352  12.016  -1.686  1.00116.50           C
ATOM   1075  O   ARG B  48       9.439  11.234  -1.949  1.00116.50           O
ATOM   1076  CB  ARG B  48      12.553  11.047  -2.055  1.00116.50           C
ATOM   1077  CG  ARG B  48      13.858  10.447  -1.523  1.00116.50           C
ATOM   1078  CD  ARG B  48      14.929  11.476  -1.159  1.00116.50           C
ATOM   1079  NE  ARG B  48      16.178  10.721  -0.852  1.00116.50           N
ATOM   1080  CZ  ARG B  48      17.057  11.194   0.080  1.00116.50           C
ATOM   1081  NH1 ARG B  48      16.786  12.353   0.746  1.00116.50           N1+
ATOM   1082  NH2 ARG B  48      18.198  10.499   0.358  1.00116.50           N
ATOM   1083  N   ALA B  49      10.274  13.326  -2.000  1.00 51.03           N
ATOM   1084  CA  ALA B  49       9.124  13.820  -2.705  1.00 51.03           C
ATOM   1085  C   ALA B  49       9.306  13.570  -4.174  1.00 51.03           C
ATOM   1086  O   ALA B  49      10.384  13.793  -4.722  1.00 51.03           O
ATOM   1087  CB  ALA B  49       8.899  15.331  -2.518  1.00 51.03           C
ATOM   1088  N   VAL B  50       8.268  13.007  -4.824  1.00 56.68           N
ATOM   1089  CA  VAL B  50       8.231  12.822  -6.249  1.00 56.68           C
ATOM   1090  C   VAL B  50       7.869  14.122  -6.905  1.00 56.68           C
ATOM   1091  O   VAL B  50       8.459  14.514  -7.910  1.00 56.68           O
ATOM   1092  CB  VAL B  50       7.222  11.789  -6.659  1.00 56.68           C
ATOM   1093  CG1 VAL B  50       7.175  11.718  -8.196  1.00 56.68           C
ATOM   1094  CG2 VAL B  50       7.604  10.459  -5.985  1.00 56.68           C
ATOM   1095  N   THR B  51       6.865  14.818  -6.330  1.00 43.36           N
ATOM   1096  CA  THR B  51       6.377  16.068  -6.841  1.00 43.36           C
ATOM   1097  C   THR B  51       6.424  17.019  -5.690  1.00 43.36           C
ATOM   1098  O   THR B  51       6.624  16.611  -4.549  1.00 43.36           O
ATOM   1099  CB  THR B  51       4.954  16.012  -7.312  1.00 43.36           C
ATOM   1100  CG2 THR B  51       4.841  14.951  -8.421  1.00 43.36           C
ATOM   1101  OG1 THR B  51       4.092  15.684  -6.231  1.00 43.36           O
ATOM   1102  N   GLU B  52       6.229  18.322  -5.966  1.00 82.41           N
ATOM   1103  CA  GLU B  52       6.343  19.341  -4.963  1.00 82.41           C
ATOM   1104  C   GLU B  52       5.299  19.155  -3.903  1.00 82.41           C
ATOM   1105  O   GLU B  52       5.528  19.492  -2.742  1.00 82.41           O
ATOM   1106  CB  GLU B  52       6.207  20.765  -5.527  1.00 82.41           C
ATOM   1107  CG  GLU B  52       6.559  21.843  -4.502  1.00 82.41           C
ATOM   1108  CD  GLU B  52       8.046  21.720  -4.192  1.00 82.41           C
ATOM   1109  OE1 GLU B  52       8.814  21.324  -5.110  1.00 82.41           O
ATOM   1110  OE2 GLU B  52       8.434  22.015  -3.030  1.00 82.41           O1-
ATOM   1111  N   LEU B  53       4.121  18.610  -4.264  1.00 67.76           N
ATOM   1112  CA  LEU B  53       3.040  18.458  -3.326  1.00 67.76           C
ATOM   1113  C   LEU B  53       3.458  17.585  -2.184  1.00 67.76           C
ATOM   1114  O   LEU B  53       2.990  17.759  -1.060  1.00 67.76           O
ATOM   1115  CB  LEU B  53       1.784  17.804  -3.934  1.00 67.76           C
ATOM   1116  CG  LEU B  53       0.970  18.721  -4.863  1.00 67.76           C
ATOM   1117  CD1 LEU B  53      -0.262  17.989  -5.425  1.00 67.76           C
ATOM   1118  CD2 LEU B  53       0.592  20.028  -4.148  1.00 67.76           C
ATOM   1119  N   GLY B  54       4.300  16.583  -2.470  1.00 60.69           N
ATOM   1120  CA  GLY B  54       4.791  15.591  -1.554  1.00 60.69           C
ATOM   1121  C   GLY B  54       5.770  16.089  -0.526  1.00 60.69           C
ATOM   1122  O   GLY B  54       5.957  15.436   0.498  1.00 60.69           O
ATOM   1123  N   ARG B  55       6.473  17.208  -0.780  1.00140.74           N
ATOM   1124  CA  ARG B  55       7.575  17.594   0.063  1.00140.74           C
ATOM   1125  C   ARG B  55       7.188  17.738   1.501  1.00140.74           C
ATOM   1126  O   ARG B  55       7.959  17.312   2.361  1.00140.74           O
ATOM   1127  CB  ARG B  55       8.255  18.901  -0.375  1.00140.74           C
ATOM   1128  CG  ARG B  55       9.025  18.727  -1.680  1.00140.74           C
ATOM   1129  CD  ARG B  55      10.462  19.242  -1.637  1.00140.74           C
ATOM   1130  NE  ARG B  55      10.417  20.728  -1.618  1.00140.74           N
ATOM   1131  CZ  ARG B  55      11.496  21.428  -2.072  1.00140.74           C
ATOM   1132  NH1 ARG B  55      12.579  20.756  -2.560  1.00140.74           N1+
ATOM   1133  NH2 ARG B  55      11.491  22.793  -2.042  1.00140.74           N
ATOM   1134  N   PRO B  56       6.065  18.299   1.845  1.00 70.69           N
ATOM   1135  CA  PRO B  56       5.754  18.456   3.237  1.00 70.69           C
ATOM   1136  C   PRO B  56       5.720  17.130   3.930  1.00 70.69           C
ATOM   1137  O   PRO B  56       6.143  17.046   5.082  1.00 70.69           O
ATOM   1138  CB  PRO B  56       4.448  19.245   3.269  1.00 70.69           C
ATOM   1139  CG  PRO B  56       4.511  20.089   1.976  1.00 70.69           C
ATOM   1140  CD  PRO B  56       5.367  19.256   1.004  1.00 70.69           C
ATOM   1141  N   ASP B  57       5.218  16.089   3.248  1.00 51.41           N
ATOM   1142  CA  ASP B  57       5.131  14.775   3.812  1.00 51.41           C
ATOM   1143  C   ASP B  57       6.504  14.194   3.973  1.00 51.41           C
ATOM   1144  O   ASP B  57       6.809  13.568   4.987  1.00 51.41           O
ATOM   1145  CB  ASP B  57       4.321  13.826   2.915  1.00 51.41           C
ATOM   1146  CG  ASP B  57       2.897  14.366   2.845  1.00 51.41           C
ATOM   1147  OD1 ASP B  57       2.373  14.807   3.904  1.00 51.41           O
ATOM   1148  OD2 ASP B  57       2.323  14.361   1.724  1.00 51.41           O1-
ATOM   1149  N   ALA B  58       7.380  14.379   2.967  1.00 31.35           N
ATOM   1150  CA  ALA B  58       8.686  13.788   3.042  1.00 31.35           C
ATOM   1151  C   ALA B  58       9.433  14.382   4.192  1.00 31.35           C
ATOM   1152  O   ALA B  58      10.051  13.672   4.984  1.00 31.35           O
ATOM   1153  CB  ALA B  58       9.520  14.023   1.770  1.00 31.35           C
ATOM   1154  N   GLU B  59       9.364  15.716   4.332  1.00 44.87           N
ATOM   1155  CA  GLU B  59      10.105  16.393   5.354  1.00 44.87           C
ATOM   1156  C   GLU B  59       9.608  15.964   6.701  1.00 44.87           C
ATOM   1157  O   GLU B  59      10.396  15.688   7.603  1.00 44.87           O
ATOM   1158  CB  GLU B  59       9.968  17.920   5.211  1.00 44.87           C
ATOM   1159  CG  GLU B  59      10.588  18.429   3.903  1.00 44.87           C
ATOM   1160  CD  GLU B  59      10.145  19.864   3.660  1.00 44.87           C
ATOM   1161  OE1 GLU B  59       9.315  20.374   4.457  1.00 44.87           O
ATOM   1162  OE2 GLU B  59      10.630  20.468   2.666  1.00 44.87           O1-
ATOM   1163  N   TYR B  60       8.275  15.886   6.863  1.00 58.79           N
ATOM   1164  CA  TYR B  60       7.676  15.560   8.127  1.00 58.79           C
ATOM   1165  C   TYR B  60       7.991  14.154   8.546  1.00 58.79           C
ATOM   1166  O   TYR B  60       8.433  13.918   9.669  1.00 58.79           O
ATOM   1167  CB  TYR B  60       6.148  15.728   8.075  1.00 58.79           C
ATOM   1168  CG  TYR B  60       5.564  15.384   9.401  1.00 58.79           C
ATOM   1169  CD1 TYR B  60       5.732  16.218  10.481  1.00 58.79           C
ATOM   1170  CD2 TYR B  60       4.823  14.235   9.554  1.00 58.79           C
ATOM   1171  CE1 TYR B  60       5.183  15.903  11.702  1.00 58.79           C
ATOM   1172  CE2 TYR B  60       4.271  13.915  10.771  1.00 58.79           C
ATOM   1173  CZ  TYR B  60       4.452  14.749  11.848  1.00 58.79           C
ATOM   1174  OH  TYR B  60       3.884  14.421  13.098  1.00 58.79           O
ATOM   1175  N   TRP B  61       7.796  13.181   7.640  1.00 74.73           N
ATOM   1176  CA  TRP B  61       7.988  11.793   7.968  1.00 74.73           C
ATOM   1177  C   TRP B  61       9.428  11.503   8.256  1.00 74.73           C
ATOM   1178  O   TRP B  61       9.747  10.765   9.188  1.00 74.73           O
ATOM   1179  CB  TRP B  61       7.483  10.867   6.846  1.00 74.73           C
ATOM   1180  CG  TRP B  61       5.975  10.879   6.746  1.00 74.73           C
ATOM   1181  CD1 TRP B  61       5.063  11.195   7.712  1.00 74.73           C
ATOM   1182  CD2 TRP B  61       5.218  10.568   5.564  1.00 74.73           C
ATOM   1183  CE2 TRP B  61       3.869  10.710   5.889  1.00 74.73           C
ATOM   1184  CE3 TRP B  61       5.614  10.197   4.311  1.00 74.73           C
ATOM   1185  NE1 TRP B  61       3.788  11.091   7.210  1.00 74.73           N
ATOM   1186  CZ2 TRP B  61       2.891  10.481   4.961  1.00 74.73           C
ATOM   1187  CZ3 TRP B  61       4.626   9.967   3.379  1.00 74.73           C
ATOM   1188  CH2 TRP B  61       3.290  10.106   3.696  1.00 74.73           C
ATOM   1189  N   ASN B  62      10.349  12.087   7.473  1.00 78.42           N
ATOM   1190  CA  ASN B  62      11.744  11.800   7.653  1.00 78.42           C
ATOM   1191  C   ASN B  62      12.196  12.276   9.000  1.00 78.42           C
ATOM   1192  O   ASN B  62      13.174  11.771   9.549  1.00 78.42           O
ATOM   1193  CB  ASN B  62      12.637  12.408   6.564  1.00 78.42           C
ATOM   1194  CG  ASN B  62      12.399  11.561   5.320  1.00 78.42           C
ATOM   1195  ND2 ASN B  62      12.585  10.221   5.455  1.00 78.42           N
ATOM   1196  OD1 ASN B  62      12.045  12.070   4.261  1.00 78.42           O
ATOM   1197  N   SER B  63      11.508  13.284   9.562  1.00 85.43           N
ATOM   1198  CA  SER B  63      11.874  13.865  10.824  1.00 85.43           C
ATOM   1199  C   SER B  63      11.732  12.878  11.948  1.00 85.43           C
ATOM   1200  O   SER B  63      12.418  13.010  12.960  1.00 85.43           O
ATOM   1201  CB  SER B  63      10.981  15.061  11.189  1.00 85.43           C
ATOM   1202  OG  SER B  63      11.054  16.051  10.174  1.00 85.43           O
ATOM   1203  N   GLN B  64      10.846  11.867  11.822  1.00 74.57           N
ATOM   1204  CA  GLN B  64      10.581  10.992  12.936  1.00 74.57           C
ATOM   1205  C   GLN B  64      11.400   9.738  12.840  1.00 74.57           C
ATOM   1206  O   GLN B  64      11.196   8.892  11.971  1.00 74.57           O
ATOM   1207  CB  GLN B  64       9.096  10.627  13.006  1.00 74.57           C
ATOM   1208  CG  GLN B  64       8.242  11.873  13.242  1.00 74.57           C
ATOM   1209  CD  GLN B  64       6.818  11.564  12.818  1.00 74.57           C
ATOM   1210  NE2 GLN B  64       5.872  12.475  13.166  1.00 74.57           N
ATOM   1211  OE1 GLN B  64       6.556  10.549  12.177  1.00 74.57           O
ATOM   1212  N   LYS B  65      12.319   9.567  13.810  1.00 68.75           N
ATOM   1213  CA  LYS B  65      13.267   8.491  13.821  1.00 68.75           C
ATOM   1214  C   LYS B  65      12.583   7.160  13.917  1.00 68.75           C
ATOM   1215  O   LYS B  65      13.021   6.195  13.295  1.00 68.75           O
ATOM   1216  CB  LYS B  65      14.262   8.586  14.989  1.00 68.75           C
ATOM   1217  CG  LYS B  65      15.454   7.639  14.843  1.00 68.75           C
ATOM   1218  CD  LYS B  65      16.611   7.958  15.793  1.00 68.75           C
ATOM   1219  CE  LYS B  65      17.881   7.160  15.496  1.00 68.75           C
ATOM   1220  NZ  LYS B  65      18.537   7.695  14.281  1.00 68.75           N1+
ATOM   1221  N   ASP B  66      11.509   7.060  14.719  1.00 38.58           N
ATOM   1222  CA  ASP B  66      10.838   5.803  14.912  1.00 38.58           C
ATOM   1223  C   ASP B  66      10.203   5.346  13.635  1.00 38.58           C
ATOM   1224  O   ASP B  66      10.232   4.159  13.310  1.00 38.58           O
ATOM   1225  CB  ASP B  66       9.732   5.887  15.976  1.00 38.58           C
ATOM   1226  CG  ASP B  66      10.407   6.066  17.325  1.00 38.58           C
ATOM   1227  OD1 ASP B  66      11.438   5.383  17.570  1.00 38.58           O
ATOM   1228  OD2 ASP B  66       9.910   6.900  18.128  1.00 38.58           O1-
ATOM   1229  N   LEU B  67       9.607   6.282  12.880  1.00 54.76           N
ATOM   1230  CA  LEU B  67       8.916   5.956  11.662  1.00 54.76           C
ATOM   1231  C   LEU B  67       9.886   5.380  10.679  1.00 54.76           C
ATOM   1232  O   LEU B  67       9.608   4.364  10.042  1.00 54.76           O
ATOM   1233  CB  LEU B  67       8.275   7.202  11.029  1.00 54.76           C
ATOM   1234  CG  LEU B  67       7.598   6.961   9.668  1.00 54.76           C
ATOM   1235  CD1 LEU B  67       6.466   5.929   9.770  1.00 54.76           C
ATOM   1236  CD2 LEU B  67       7.132   8.289   9.052  1.00 54.76           C
ATOM   1237  N   LEU B  68      11.063   6.013  10.544  1.00 47.95           N
ATOM   1238  CA  LEU B  68      12.054   5.571   9.606  1.00 47.95           C
ATOM   1239  C   LEU B  68      12.541   4.219  10.023  1.00 47.95           C
ATOM   1240  O   LEU B  68      12.831   3.375   9.178  1.00 47.95           O
ATOM   1241  CB  LEU B  68      13.267   6.517   9.519  1.00 47.95           C
ATOM   1242  CG  LEU B  68      12.927   7.905   8.938  1.00 47.95           C
ATOM   1243  CD1 LEU B  68      14.172   8.806   8.870  1.00 47.95           C
ATOM   1244  CD2 LEU B  68      12.212   7.785   7.584  1.00 47.95           C
ATOM   1245  N   GLU B  69      12.644   3.985  11.347  1.00 79.45           N
ATOM   1246  CA  GLU B  69      13.164   2.741  11.839  1.00 79.45           C
ATOM   1247  C   GLU B  69      12.286   1.611  11.400  1.00 79.45           C
ATOM   1248  O   GLU B  69      12.776   0.592  10.918  1.00 79.45           O
ATOM   1249  CB  GLU B  69      13.242   2.661  13.375  1.00 79.45           C
ATOM   1250  CG  GLU B  69      14.335   3.530  13.995  1.00 79.45           C
ATOM   1251  CD  GLU B  69      14.452   3.140  15.462  1.00 79.45           C
ATOM   1252  OE1 GLU B  69      14.331   1.921  15.758  1.00 79.45           O
ATOM   1253  OE2 GLU B  69      14.670   4.051  16.305  1.00 79.45           O1-
ATOM   1254  N   GLN B  70      10.958   1.774  11.541  1.00121.42           N
ATOM   1255  CA  GLN B  70      10.005   0.759  11.190  1.00121.42           C
ATOM   1256  C   GLN B  70      10.079   0.444   9.727  1.00121.42           C
ATOM   1257  O   GLN B  70      10.128  -0.723   9.339  1.00121.42           O
ATOM   1258  CB  GLN B  70       8.570   1.231  11.504  1.00121.42           C
ATOM   1259  CG  GLN B  70       7.456   0.669  10.610  1.00121.42           C
ATOM   1260  CD  GLN B  70       7.507  -0.847  10.547  1.00121.42           C
ATOM   1261  NE2 GLN B  70       6.497  -1.455   9.869  1.00121.42           N
ATOM   1262  OE1 GLN B  70       8.422  -1.477  11.071  1.00121.42           O
ATOM   1263  N   ARG B  71      10.124   1.481   8.877  1.00105.45           N
ATOM   1264  CA  ARG B  71      10.098   1.309   7.450  1.00105.45           C
ATOM   1265  C   ARG B  71      11.311   0.575   6.968  1.00105.45           C
ATOM   1266  O   ARG B  71      11.226  -0.246   6.056  1.00105.45           O
ATOM   1267  CB  ARG B  71      10.054   2.638   6.678  1.00105.45           C
ATOM   1268  CG  ARG B  71       8.734   3.394   6.831  1.00105.45           C
ATOM   1269  CD  ARG B  71       7.510   2.543   6.492  1.00105.45           C
ATOM   1270  NE  ARG B  71       6.339   3.458   6.421  1.00105.45           N
ATOM   1271  CZ  ARG B  71       5.098   3.016   6.782  1.00105.45           C
ATOM   1272  NH1 ARG B  71       4.940   1.761   7.293  1.00105.45           N1+
ATOM   1273  NH2 ARG B  71       4.015   3.832   6.628  1.00105.45           N
ATOM   1274  N   ARG B  72      12.484   0.870   7.550  1.00 76.66           N
ATOM   1275  CA  ARG B  72      13.704   0.261   7.102  1.00 76.66           C
ATOM   1276  C   ARG B  72      13.643  -1.217   7.355  1.00 76.66           C
ATOM   1277  O   ARG B  72      14.202  -2.006   6.598  1.00 76.66           O
ATOM   1278  CB  ARG B  72      14.937   0.858   7.803  1.00 76.66           C
ATOM   1279  CG  ARG B  72      15.148   2.324   7.414  1.00 76.66           C
ATOM   1280  CD  ARG B  72      16.348   3.004   8.073  1.00 76.66           C
ATOM   1281  NE  ARG B  72      16.385   4.404   7.562  1.00 76.66           N
ATOM   1282  CZ  ARG B  72      17.043   5.378   8.257  1.00 76.66           C
ATOM   1283  NH1 ARG B  72      17.668   5.077   9.434  1.00 76.66           N1+
ATOM   1284  NH2 ARG B  72      17.077   6.653   7.774  1.00 76.66           N
ATOM   1285  N   ALA B  73      13.016  -1.603   8.478  1.00 72.75           N
ATOM   1286  CA  ALA B  73      12.786  -2.939   8.963  1.00 72.75           C
ATOM   1287  C   ALA B  73      11.754  -3.687   8.169  1.00 72.75           C
ATOM   1288  O   ALA B  73      11.735  -4.917   8.188  1.00 72.75           O
ATOM   1289  CB  ALA B  73      12.324  -2.959  10.431  1.00 72.75           C
ATOM   1290  N   ALA B  74      10.838  -2.960   7.503  1.00 58.23           N
ATOM   1291  CA  ALA B  74       9.676  -3.496   6.845  1.00 58.23           C
ATOM   1292  C   ALA B  74      10.055  -4.534   5.841  1.00 58.23           C
ATOM   1293  O   ALA B  74       9.280  -5.457   5.597  1.00 58.23           O
ATOM   1294  CB  ALA B  74       8.861  -2.418   6.108  1.00 58.23           C
ATOM   1295  N   VAL B  75      11.234  -4.410   5.210  1.00115.95           N
ATOM   1296  CA  VAL B  75      11.598  -5.374   4.219  1.00115.95           C
ATOM   1297  C   VAL B  75      11.566  -6.738   4.851  1.00115.95           C
ATOM   1298  O   VAL B  75      11.164  -7.706   4.206  1.00115.95           O
ATOM   1299  CB  VAL B  75      12.957  -5.127   3.621  1.00115.95           C
ATOM   1300  CG1 VAL B  75      14.057  -5.536   4.614  1.00115.95           C
ATOM   1301  CG2 VAL B  75      13.026  -5.858   2.273  1.00115.95           C
ATOM   1302  N   ASP B  76      12.057  -6.862   6.105  1.00 73.18           N
ATOM   1303  CA  ASP B  76      12.002  -8.097   6.841  1.00 73.18           C
ATOM   1304  C   ASP B  76      10.655  -8.380   7.449  1.00 73.18           C
ATOM   1305  O   ASP B  76      10.107  -9.470   7.295  1.00 73.18           O
ATOM   1306  CB  ASP B  76      13.024  -8.131   7.989  1.00 73.18           C
ATOM   1307  CG  ASP B  76      14.410  -8.198   7.367  1.00 73.18           C
ATOM   1308  OD1 ASP B  76      14.606  -9.034   6.446  1.00 73.18           O
ATOM   1309  OD2 ASP B  76      15.288  -7.401   7.793  1.00 73.18           O1-
ATOM   1310  N   THR B  77      10.096  -7.395   8.183  1.00 76.61           N
ATOM   1311  CA  THR B  77       8.878  -7.590   8.923  1.00 76.61           C
ATOM   1312  C   THR B  77       7.695  -7.738   8.020  1.00 76.61           C
ATOM   1313  O   THR B  77       6.790  -8.508   8.330  1.00 76.61           O
ATOM   1314  CB  THR B  77       8.591  -6.456   9.863  1.00 76.61           C
ATOM   1315  CG2 THR B  77       9.760  -6.330  10.856  1.00 76.61           C
ATOM   1316  OG1 THR B  77       8.428  -5.245   9.140  1.00 76.61           O
ATOM   1317  N   TYR B  78       7.585  -6.866   6.998  1.00122.49           N
ATOM   1318  CA  TYR B  78       6.506  -6.891   6.045  1.00122.49           C
ATOM   1319  C   TYR B  78       6.679  -7.802   4.849  1.00122.49           C
ATOM   1320  O   TYR B  78       5.926  -8.756   4.658  1.00122.49           O
ATOM   1321  CB  TYR B  78       6.263  -5.449   5.557  1.00122.49           C
ATOM   1322  CG  TYR B  78       5.134  -5.334   4.593  1.00122.49           C
ATOM   1323  CD1 TYR B  78       3.834  -5.503   5.010  1.00122.49           C
ATOM   1324  CD2 TYR B  78       5.378  -5.006   3.281  1.00122.49           C
ATOM   1325  CE1 TYR B  78       2.795  -5.373   4.120  1.00122.49           C
ATOM   1326  CE2 TYR B  78       4.342  -4.874   2.386  1.00122.49           C
ATOM   1327  CZ  TYR B  78       3.048  -5.060   2.805  1.00122.49           C
ATOM   1328  OH  TYR B  78       1.980  -4.925   1.893  1.00122.49           O
ATOM   1329  N   CYS B  79       7.720  -7.533   4.026  1.00 68.86           N
ATOM   1330  CA  CYS B  79       7.872  -8.193   2.751  1.00 68.86           C
ATOM   1331  C   CYS B  79       8.210  -9.647   2.864  1.00 68.86           C
ATOM   1332  O   CYS B  79       7.446 -10.509   2.432  1.00 68.86           O
ATOM   1333  CB  CYS B  79       8.937  -7.529   1.856  1.00 68.86           C
ATOM   1334  SG  CYS B  79       8.511  -5.812   1.420  1.00 68.86           S
ATOM   1335  N   ARG B  80       9.337  -9.966   3.524  1.00 99.78           N
ATOM   1336  CA  ARG B  80       9.789 -11.325   3.575  1.00 99.78           C
ATOM   1337  C   ARG B  80       8.745 -12.135   4.260  1.00 99.78           C
ATOM   1338  O   ARG B  80       8.538 -13.301   3.932  1.00 99.78           O
ATOM   1339  CB  ARG B  80      11.127 -11.508   4.311  1.00 99.78           C
ATOM   1340  CG  ARG B  80      12.324 -11.048   3.476  1.00 99.78           C
ATOM   1341  CD  ARG B  80      13.677 -11.269   4.151  1.00 99.78           C
ATOM   1342  NE  ARG B  80      14.730 -10.871   3.177  1.00 99.78           N
ATOM   1343  CZ  ARG B  80      15.093  -9.560   3.059  1.00 99.78           C
ATOM   1344  NH1 ARG B  80      14.467  -8.609   3.812  1.00 99.78           N1+
ATOM   1345  NH2 ARG B  80      16.081  -9.200   2.189  1.00 99.78           N
ATOM   1346  N   HIS B  81       8.063 -11.536   5.246  1.00 47.65           N
ATOM   1347  CA  HIS B  81       7.065 -12.253   5.979  1.00 47.65           C
ATOM   1348  C   HIS B  81       5.941 -12.686   5.077  1.00 47.65           C
ATOM   1349  O   HIS B  81       5.532 -13.845   5.110  1.00 47.65           O
ATOM   1350  CB  HIS B  81       6.438 -11.415   7.104  1.00 47.65           C
ATOM   1351  CG  HIS B  81       5.401 -12.179   7.872  1.00 47.65           C
ATOM   1352  CD2 HIS B  81       4.046 -12.172   7.753  1.00 47.65           C
ATOM   1353  ND1 HIS B  81       5.686 -13.079   8.875  1.00 47.65           N
ATOM   1354  CE1 HIS B  81       4.497 -13.566   9.311  1.00 47.65           C
ATOM   1355  NE2 HIS B  81       3.472 -13.044   8.660  1.00 47.65           N
ATOM   1356  N   ASN B  82       5.414 -11.770   4.243  1.00 51.16           N
ATOM   1357  CA  ASN B  82       4.282 -12.057   3.402  1.00 51.16           C
ATOM   1358  C   ASN B  82       4.645 -13.047   2.338  1.00 51.16           C
ATOM   1359  O   ASN B  82       3.810 -13.846   1.917  1.00 51.16           O
ATOM   1360  CB  ASN B  82       3.721 -10.809   2.698  1.00 51.16           C
ATOM   1361  CG  ASN B  82       3.016  -9.961   3.747  1.00 51.16           C
ATOM   1362  ND2 ASN B  82       2.997  -8.618   3.531  1.00 51.16           N
ATOM   1363  OD1 ASN B  82       2.490 -10.478   4.731  1.00 51.16           O
ATOM   1364  N   TYR B  83       5.899 -13.008   1.855  1.00 59.52           N
ATOM   1365  CA  TYR B  83       6.335 -13.874   0.798  1.00 59.52           C
ATOM   1366  C   TYR B  83       6.189 -15.287   1.262  1.00 59.52           C
ATOM   1367  O   TYR B  83       5.661 -16.140   0.548  1.00 59.52           O
ATOM   1368  CB  TYR B  83       7.822 -13.635   0.462  1.00 59.52           C
ATOM   1369  CG  TYR B  83       8.255 -14.544  -0.638  1.00 59.52           C
ATOM   1370  CD1 TYR B  83       8.609 -15.848  -0.372  1.00 59.52           C
ATOM   1371  CD2 TYR B  83       8.324 -14.086  -1.933  1.00 59.52           C
ATOM   1372  CE1 TYR B  83       9.016 -16.682  -1.387  1.00 59.52           C
ATOM   1373  CE2 TYR B  83       8.729 -14.916  -2.952  1.00 59.52           C
ATOM   1374  CZ  TYR B  83       9.077 -16.218  -2.679  1.00 59.52           C
ATOM   1375  OH  TYR B  83       9.495 -17.073  -3.720  1.00 59.52           O
ATOM   1376  N   GLY B  84       6.634 -15.562   2.500  1.00 31.78           N
ATOM   1377  CA  GLY B  84       6.602 -16.881   3.054  1.00 31.78           C
ATOM   1378  C   GLY B  84       5.188 -17.360   3.186  1.00 31.78           C
ATOM   1379  O   GLY B  84       4.914 -18.548   3.027  1.00 31.78           O
ATOM   1380  N   VAL B  85       4.252 -16.453   3.522  1.00 55.28           N
ATOM   1381  CA  VAL B  85       2.894 -16.846   3.779  1.00 55.28           C
ATOM   1382  C   VAL B  85       2.241 -17.460   2.567  1.00 55.28           C
ATOM   1383  O   VAL B  85       1.621 -18.518   2.664  1.00 55.28           O
ATOM   1384  CB  VAL B  85       2.044 -15.688   4.214  1.00 55.28           C
ATOM   1385  CG1 VAL B  85       0.593 -16.173   4.380  1.00 55.28           C
ATOM   1386  CG2 VAL B  85       2.652 -15.090   5.494  1.00 55.28           C
ATOM   1387  N   GLY B  86       2.324 -16.759   1.419  1.00 73.17           N
ATOM   1388  CA  GLY B  86       1.767 -17.063   0.121  1.00 73.17           C
ATOM   1389  C   GLY B  86       2.503 -18.097  -0.680  1.00 73.17           C
ATOM   1390  O   GLY B  86       1.950 -18.632  -1.640  1.00 73.17           O
ATOM   1391  N   GLU B  87       3.778 -18.373  -0.351  1.00 69.04           N
ATOM   1392  CA  GLU B  87       4.660 -19.110  -1.217  1.00 69.04           C
ATOM   1393  C   GLU B  87       4.047 -20.382  -1.724  1.00 69.04           C
ATOM   1394  O   GLU B  87       4.189 -20.693  -2.906  1.00 69.04           O
ATOM   1395  CB  GLU B  87       5.986 -19.474  -0.527  1.00 69.04           C
ATOM   1396  CG  GLU B  87       7.017 -20.093  -1.474  1.00 69.04           C
ATOM   1397  CD  GLU B  87       8.274 -20.394  -0.674  1.00 69.04           C
ATOM   1398  OE1 GLU B  87       8.716 -19.500   0.097  1.00 69.04           O
ATOM   1399  OE2 GLU B  87       8.809 -21.524  -0.821  1.00 69.04           O1-
ATOM   1400  N   SER B  88       3.333 -21.148  -0.882  1.00 95.08           N
ATOM   1401  CA  SER B  88       2.794 -22.393  -1.351  1.00 95.08           C
ATOM   1402  C   SER B  88       1.896 -22.147  -2.528  1.00 95.08           C
ATOM   1403  O   SER B  88       2.070 -22.740  -3.588  1.00 95.08           O
ATOM   1404  CB  SER B  88       1.915 -23.092  -0.299  1.00 95.08           C
ATOM   1405  OG  SER B  88       2.666 -23.364   0.873  1.00 95.08           O
ATOM   1406  N   PHE B  89       0.870 -21.299  -2.342  1.00123.26           N
ATOM   1407  CA  PHE B  89      -0.098 -21.024  -3.365  1.00123.26           C
ATOM   1408  C   PHE B  89       0.327 -20.065  -4.445  1.00123.26           C
ATOM   1409  O   PHE B  89      -0.230 -20.113  -5.539  1.00123.26           O
ATOM   1410  CB  PHE B  89      -1.484 -20.633  -2.808  1.00123.26           C
ATOM   1411  CG  PHE B  89      -1.376 -19.487  -1.862  1.00123.26           C
ATOM   1412  CD1 PHE B  89      -1.447 -18.188  -2.307  1.00123.26           C
ATOM   1413  CD2 PHE B  89      -1.215 -19.724  -0.516  1.00123.26           C
ATOM   1414  CE1 PHE B  89      -1.352 -17.141  -1.422  1.00123.26           C
ATOM   1415  CE2 PHE B  89      -1.120 -18.683   0.376  1.00123.26           C
ATOM   1416  CZ  PHE B  89      -1.189 -17.388  -0.078  1.00123.26           C
ATOM   1417  N   THR B  90       1.153 -19.045  -4.132  1.00178.77           N
ATOM   1418  CA  THR B  90       1.535 -18.099  -5.153  1.00178.77           C
ATOM   1419  C   THR B  90       2.678 -18.472  -6.054  1.00178.77           C
ATOM   1420  O   THR B  90       2.497 -18.820  -7.219  1.00178.77           O
ATOM   1421  CB  THR B  90       1.921 -16.805  -4.527  1.00178.77           C
ATOM   1422  CG2 THR B  90       2.232 -15.829  -5.652  1.00178.77           C
ATOM   1423  OG1 THR B  90       0.858 -16.300  -3.738  1.00178.77           O
ATOM   1424  N   VAL B  91       3.907 -18.466  -5.483  1.00150.68           N
ATOM   1425  CA  VAL B  91       5.104 -18.621  -6.266  1.00150.68           C
ATOM   1426  C   VAL B  91       5.141 -19.999  -6.829  1.00150.68           C
ATOM   1427  O   VAL B  91       5.635 -20.223  -7.931  1.00150.68           O
ATOM   1428  CB  VAL B  91       6.365 -18.319  -5.512  1.00150.68           C
ATOM   1429  CG1 VAL B  91       6.373 -16.816  -5.188  1.00150.68           C
ATOM   1430  CG2 VAL B  91       6.428 -19.193  -4.260  1.00150.68           C
ATOM   1431  N   GLN B  92       4.662 -20.962  -6.034  1.00 90.49           N
ATOM   1432  CA  GLN B  92       4.594 -22.363  -6.322  1.00 90.49           C
ATOM   1433  C   GLN B  92       3.486 -22.761  -7.257  1.00 90.49           C
ATOM   1434  O   GLN B  92       3.585 -23.817  -7.874  1.00 90.49           O
ATOM   1435  CB  GLN B  92       4.463 -23.204  -5.045  1.00 90.49           C
ATOM   1436  CG  GLN B  92       5.699 -23.079  -4.153  1.00 90.49           C
ATOM   1437  CD  GLN B  92       5.490 -23.935  -2.918  1.00 90.49           C
ATOM   1438  NE2 GLN B  92       5.988 -23.439  -1.753  1.00 90.49           N
ATOM   1439  OE1 GLN B  92       4.904 -25.013  -2.982  1.00 90.49           O
ATOM   1440  N   ARG B  93       2.384 -21.990  -7.365  1.00176.67           N
ATOM   1441  CA  ARG B  93       1.232 -22.494  -8.071  1.00176.67           C
ATOM   1442  C   ARG B  93       1.543 -22.842  -9.491  1.00176.67           C
ATOM   1443  O   ARG B  93       2.294 -22.151 -10.179  1.00176.67           O
ATOM   1444  CB  ARG B  93       0.006 -21.568  -8.048  1.00176.67           C
ATOM   1445  CG  ARG B  93       0.180 -20.249  -8.789  1.00176.67           C
ATOM   1446  CD  ARG B  93      -0.883 -19.224  -8.401  1.00176.67           C
ATOM   1447  NE  ARG B  93      -2.190 -19.938  -8.360  1.00176.67           N
ATOM   1448  CZ  ARG B  93      -2.928 -20.113  -9.492  1.00176.67           C
ATOM   1449  NH1 ARG B  93      -2.470 -19.634 -10.684  1.00176.67           N1+
ATOM   1450  NH2 ARG B  93      -4.122 -20.772  -9.430  1.00176.67           N
ATOM   1451  N   ARG B  94       0.972 -23.977  -9.953  1.00 81.38           N
ATOM   1452  CA  ARG B  94       1.180 -24.413 -11.302  1.00 81.38           C
ATOM   1453  C   ARG B  94      -0.102 -25.013 -11.803  1.00 81.38           C
ATOM   1454  O   ARG B  94      -0.695 -25.873 -11.150  1.00 81.38           O
ATOM   1455  CB  ARG B  94       2.250 -25.510 -11.428  1.00 81.38           C
ATOM   1456  CG  ARG B  94       3.654 -25.063 -11.020  1.00 81.38           C
ATOM   1457  CD  ARG B  94       4.380 -24.263 -12.104  1.00 81.38           C
ATOM   1458  NE  ARG B  94       5.784 -24.057 -11.647  1.00 81.38           N
ATOM   1459  CZ  ARG B  94       6.114 -22.968 -10.893  1.00 81.38           C
ATOM   1460  NH1 ARG B  94       5.155 -22.062 -10.543  1.00 81.38           N1+
ATOM   1461  NH2 ARG B  94       7.407 -22.785 -10.494  1.00 81.38           N
TER    1462      ARG B  94
ATOM   1463  N   PRO C   1      -1.341 -19.042   9.092  1.00 39.38           N
ATOM   1464  CA  PRO C   1      -0.373 -17.925   8.953  1.00 39.38           C
ATOM   1465  C   PRO C   1      -1.139 -16.708   8.551  1.00 39.38           C
ATOM   1466  O   PRO C   1      -2.168 -16.829   7.877  1.00 39.38           O
ATOM   1467  CB  PRO C   1       0.673 -18.370   7.932  1.00 39.38           C
ATOM   1468  CG  PRO C   1       0.004 -19.490   7.133  1.00 39.38           C
ATOM   1469  CD  PRO C   1      -0.959 -20.142   8.135  1.00 39.38           C
ATOM   1470  N   LYS C   2      -0.633 -15.538   8.998  1.00 67.75           N
ATOM   1471  CA  LYS C   2      -1.296 -14.294   8.753  1.00 67.75           C
ATOM   1472  C   LYS C   2      -0.452 -13.389   7.924  1.00 67.75           C
ATOM   1473  O   LYS C   2       0.736 -13.183   8.189  1.00 67.75           O
ATOM   1474  CB  LYS C   2      -1.579 -13.448  10.007  1.00 67.75           C
ATOM   1475  CG  LYS C   2      -2.719 -13.915  10.909  1.00 67.75           C
ATOM   1476  CD  LYS C   2      -2.800 -13.084  12.190  1.00 67.75           C
ATOM   1477  CE  LYS C   2      -2.395 -11.625  11.968  1.00 67.75           C
ATOM   1478  NZ  LYS C   2      -3.598 -10.780  11.801  1.00 67.75           N1+
ATOM   1479  N   TYR C   3      -1.107 -12.794   6.903  1.00 76.27           N
ATOM   1480  CA  TYR C   3      -0.552 -11.704   6.147  1.00 76.27           C
ATOM   1481  C   TYR C   3      -0.504 -10.488   7.016  1.00 76.27           C
ATOM   1482  O   TYR C   3      -1.370 -10.287   7.865  1.00 76.27           O
ATOM   1483  CB  TYR C   3      -1.362 -11.283   4.902  1.00 76.27           C
ATOM   1484  CG  TYR C   3      -1.128 -12.195   3.745  1.00 76.27           C
ATOM   1485  CD1 TYR C   3      -0.030 -12.006   2.936  1.00 76.27           C
ATOM   1486  CD2 TYR C   3      -2.004 -13.215   3.445  1.00 76.27           C
ATOM   1487  CE1 TYR C   3       0.202 -12.824   1.855  1.00 76.27           C
ATOM   1488  CE2 TYR C   3      -1.775 -14.037   2.364  1.00 76.27           C
ATOM   1489  CZ  TYR C   3      -0.673 -13.842   1.568  1.00 76.27           C
ATOM   1490  OH  TYR C   3      -0.437 -14.683   0.458  1.00 76.27           O
ATOM   1491  N   VAL C   4       0.549  -9.659   6.850  1.00116.59           N
ATOM   1492  CA  VAL C   4       0.640  -8.441   7.605  1.00116.59           C
ATOM   1493  C   VAL C   4       0.534  -7.266   6.690  1.00116.59           C
ATOM   1494  O   VAL C   4       1.130  -7.237   5.613  1.00116.59           O
ATOM   1495  CB  VAL C   4       1.904  -8.285   8.398  1.00116.59           C
ATOM   1496  CG1 VAL C   4       1.838  -9.232   9.597  1.00116.59           C
ATOM   1497  CG2 VAL C   4       3.106  -8.580   7.490  1.00116.59           C
ATOM   1498  N   LYS C   5      -0.254  -6.254   7.107  1.00174.60           N
ATOM   1499  CA  LYS C   5      -0.397  -5.096   6.277  1.00174.60           C
ATOM   1500  C   LYS C   5       0.563  -4.060   6.756  1.00174.60           C
ATOM   1501  O   LYS C   5       1.175  -4.197   7.815  1.00174.60           O
ATOM   1502  CB  LYS C   5      -1.806  -4.466   6.221  1.00174.60           C
ATOM   1503  CG  LYS C   5      -2.137  -3.412   7.285  1.00174.60           C
ATOM   1504  CD  LYS C   5      -2.225  -3.913   8.723  1.00174.60           C
ATOM   1505  CE  LYS C   5      -2.801  -2.863   9.678  1.00174.60           C
ATOM   1506  NZ  LYS C   5      -1.823  -1.778   9.908  1.00174.60           N1+
ATOM   1507  N   GLN C   6       0.754  -3.008   5.940  1.00 70.45           N
ATOM   1508  CA  GLN C   6       1.621  -1.930   6.309  1.00 70.45           C
ATOM   1509  C   GLN C   6       0.713  -0.781   6.621  1.00 70.45           C
ATOM   1510  O   GLN C   6      -0.243  -0.521   5.891  1.00 70.45           O
ATOM   1511  CB  GLN C   6       2.576  -1.510   5.178  1.00 70.45           C
ATOM   1512  CG  GLN C   6       3.704  -0.577   5.624  1.00 70.45           C
ATOM   1513  CD  GLN C   6       4.712  -1.397   6.422  1.00 70.45           C
ATOM   1514  NE2 GLN C   6       4.367  -2.684   6.693  1.00 70.45           N
ATOM   1515  OE1 GLN C   6       5.777  -0.907   6.794  1.00 70.45           O
ATOM   1516  N   ASN C   7       0.974  -0.074   7.736  1.00 58.21           N
ATOM   1517  CA  ASN C   7       0.107   0.995   8.152  1.00 58.21           C
ATOM   1518  C   ASN C   7       0.325   2.204   7.304  1.00 58.21           C
ATOM   1519  O   ASN C   7       1.450   2.536   6.932  1.00 58.21           O
ATOM   1520  CB  ASN C   7       0.318   1.413   9.615  1.00 58.21           C
ATOM   1521  CG  ASN C   7       1.742   1.935   9.759  1.00 58.21           C
ATOM   1522  ND2 ASN C   7       1.873   3.242  10.113  1.00 58.21           N
ATOM   1523  OD1 ASN C   7       2.712   1.206   9.560  1.00 58.21           O
ATOM   1524  N   THR C   8      -0.778   2.906   6.981  1.00 38.50           N
ATOM   1525  CA  THR C   8      -0.679   4.105   6.202  1.00 38.50           C
ATOM   1526  C   THR C   8      -0.301   5.195   7.157  1.00 38.50           C
ATOM   1527  O   THR C   8      -0.612   5.125   8.344  1.00 38.50           O
ATOM   1528  CB  THR C   8      -1.974   4.488   5.548  1.00 38.50           C
ATOM   1529  CG2 THR C   8      -2.466   3.291   4.717  1.00 38.50           C
ATOM   1530  OG1 THR C   8      -2.942   4.820   6.532  1.00 38.50           O
ATOM   1531  N   LEU C   9       0.394   6.234   6.659  1.00 49.05           N
ATOM   1532  CA  LEU C   9       0.838   7.292   7.522  1.00 49.05           C
ATOM   1533  C   LEU C   9      -0.011   8.497   7.272  1.00 49.05           C
ATOM   1534  O   LEU C   9      -0.410   8.762   6.139  1.00 49.05           O
ATOM   1535  CB  LEU C   9       2.296   7.704   7.260  1.00 49.05           C
ATOM   1536  CG  LEU C   9       3.309   6.588   7.572  1.00 49.05           C
ATOM   1537  CD1 LEU C   9       4.749   7.046   7.290  1.00 49.05           C
ATOM   1538  CD2 LEU C   9       3.130   6.062   9.004  1.00 49.05           C
ATOM   1539  N   LYS C  10      -0.323   9.258   8.340  1.00 71.09           N
ATOM   1540  CA  LYS C  10      -1.118  10.439   8.172  1.00 71.09           C
ATOM   1541  C   LYS C  10      -0.232  11.488   7.571  1.00 71.09           C
ATOM   1542  O   LYS C  10       0.966  11.544   7.848  1.00 71.09           O
ATOM   1543  CB  LYS C  10      -1.720  10.984   9.481  1.00 71.09           C
ATOM   1544  CG  LYS C  10      -2.791  12.058   9.262  1.00 71.09           C
ATOM   1545  CD  LYS C  10      -3.658  12.314  10.498  1.00 71.09           C
ATOM   1546  CE  LYS C  10      -4.629  11.174  10.814  1.00 71.09           C
ATOM   1547  NZ  LYS C  10      -5.409  11.490  12.032  1.00 71.09           N1+
ATOM   1548  N   LEU C  11      -0.815  12.337   6.707  1.00121.03           N
ATOM   1549  CA  LEU C  11      -0.091  13.371   6.020  1.00121.03           C
ATOM   1550  C   LEU C  11       0.201  14.517   6.931  1.00121.03           C
ATOM   1551  O   LEU C  11      -0.426  14.695   7.973  1.00121.03           O
ATOM   1552  CB  LEU C  11      -0.827  13.905   4.780  1.00121.03           C
ATOM   1553  CG  LEU C  11      -0.627  13.023   3.535  1.00121.03           C
ATOM   1554  CD1 LEU C  11      -0.828  11.537   3.854  1.00121.03           C
ATOM   1555  CD2 LEU C  11      -1.539  13.482   2.389  1.00121.03           C
ATOM   1556  N   ALA C  12       1.222  15.311   6.553  1.00 26.08           N
ATOM   1557  CA  ALA C  12       1.607  16.454   7.324  1.00 26.08           C
ATOM   1558  C   ALA C  12       0.513  17.464   7.211  1.00 26.08           C
ATOM   1559  O   ALA C  12      -0.121  17.601   6.165  1.00 26.08           O
ATOM   1560  CB  ALA C  12       2.913  17.110   6.844  1.00 26.08           C
ATOM   1561  N   THR C  13       0.259  18.191   8.315  1.00 15.40           N
ATOM   1562  CA  THR C  13      -0.767  19.190   8.322  1.00 15.40           C
ATOM   1563  C   THR C  13      -0.080  20.548   8.261  1.00 15.40           C
ATOM   1564  O   THR C  13      -0.628  21.458   7.580  1.00 15.40           O
ATOM   1565  CB  THR C  13      -1.596  19.169   9.572  1.00 15.40           C
ATOM   1566  CG2 THR C  13      -2.653  20.282   9.477  1.00 15.40           C
ATOM   1567  OG1 THR C  13      -2.224  17.904   9.719  1.00 15.40           O
ATOM   1568  OXT THR C  13       0.996  20.697   8.900  1.00 15.40           O1-
END
"""
    import tempfile
    file_name = tempfile.mktemp(prefix='reference', suffix='.pdb')
    with open(file_name, 'w') as fh:
        if reference == 1:
            fh.write(content_1)
        elif reference == 2:
            fh.write(content_2)
        else:
            raise RuntimeError('unknown reference value given (%s)' % reference)
    return file_name
