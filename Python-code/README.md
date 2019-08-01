# TCR-CoM
---------
For questions about running the script or for reporting bugs, please contact:<br/>
   *__Brian Baker__(brian-baker[at]nd.edu)*<br/>
or<br/>
   *__Esam Abualrous__(e.abualrous[at]fu-berlin.de)*<br/>
### Description
   Calculates the geometrical parameters (r, theta, phi) of T cell receptors (TCR) on the top of MHC proteins

### Python script
__USAGE__

1- For peptide-MHC class I-TCR complex:
```
python TCR_position_BioPy_v1.py --pdbfile XXX.pdb --mhc_chain_a A --tcr_chain_a D --tcr_chain_b E
``` 

2- For peptide-MHC class II-TCR complex:
```
python TCR_position_BioPy_v1.py --pdbfile XXX.pdb --mhc_chain_a A --mhc_chain_b A --tcr_chain_a D --tcr_chain_b E
```

**Notes**
---------
- The code requires defining at least 3 chain IDs for MHCI-TCR or 4 chain IDs for MHCII-TCR
- If no range of residues is defined, the code uses the default values (see USAGE)
- Please make sure that the defined chains do not contain any solvent molecules
