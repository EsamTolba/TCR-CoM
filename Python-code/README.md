# TCR-CoM
---------
For questions about running the script or for reporting bugs, please contact:<br/>
   *__Brian Baker__(brian-baker[at]nd.edu)*<br/>
or<br/>
   *__Esam Abualrous__(e.abualrous[at]fu-berlin.de)*<br/>

### Description
   Calculates the geometrical parameters (r, theta, phi) of T cell receptors (TCR) on the top of MHC proteins

### Citation
If you find this tool useful, please cite:
```
Nishant K. Singh†, Esam T. Abualrous†, Cory M. Ayres, Ragul Gowthaman, Brian G. Pierce, Frank Noé, Brian M. Baker* "Geometrical characterization of T cell receptor binding modes reveals class-specific binding to maximize access to antigen" Proteins, 2019
```

### Python script
__USAGE__

1- For peptide-MHC class I-TCR complex:
```
python tcr_com_v4.py -pdbid XXX -mhc_chain_a A -tcr_chain_a D -tcr_chain_b E -output_pdb True
```

2- For peptide-MHC class II-TCR complex:

```
python tcr_com_v4.py -pdbid XXX -mhc_chain_a A -mhc_chain_b A -tcr_chain_a D -tcr_chain_b E -output_pdb True
```

**Notes**
---------
- The approach is also available as online tool: https://tcr3d.ibbr.umd.edu/tcr_com
- Refine the submitted PDB file by one of the available PDB-editors to avoid inaccurate results
- The submitted structure of MHCI-TCR or MHCII-TCR should not contain more than five chains
- The code requires defining at least 3 chain IDs for MHCI-TCR or 4 chain IDs for MHCII-TCR
- Please make sure that the defined chains do not contain any non-protein atoms
