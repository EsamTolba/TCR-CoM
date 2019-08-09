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
python tcr_com_ref-dep_v5.py -pdbid XXX -mhc_chain_a A -tcr_chain_a D -tcr_chain_b E -output_pdb True
```

2- For peptide-MHC class II-TCR complex:
```
python tcr_com_ref-dep_v5.py -pdbid XXX -mhc_chain_a A -mhc_chain_b A -tcr_chain_a D -tcr_chain_b E -output_pdb True
```

### Pymol plugin
__USAGE__

1- For peptide-MHC class I-TCR complex:
```
tcr_mhci_pos selection [,mhc_a=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=180 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]
```

2- For peptide-MHC class II-TCR complex:
```
tcr_mhcii_pos selection [,mhc_a=None [,mhc_b=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=80 [,mhc_b_init=1 [,mhc_b_final=90 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]]]]
```

**Notes**
---------
- The approach is available as Python code, Pymol plugin, and online tool: https://tcr3d.ibbr.umd.edu/tcr_com
- The submitted PDB-file should be refined by one of the available PDB-editors to avoid inaccurate results
- Due to differences in implementation, small differences may be observed between the output values of the Python code and the Pymol plugin. 
- The submitted structure of MHCI-TCR or MHCII-TCR should not contain more than five chains
- The code requires defining at least three chain-IDs for MHCI-TCR or four chain-IDs for MHCII-TCR
- Please make sure that the defined chains do not contain any non-protein atoms
