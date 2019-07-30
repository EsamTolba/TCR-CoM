## TCR-CoM
For questions about running the script or for reporting bugs, please contact either:
Brian Baker (brian-baker[at]nd.edu)
or
Esam Abualrous (e.abualrous[at]fu-berlin.de)

# DESCRIPTION

   Calculates the geometrical parameters (r, theta, phi) of T cell receptors (TCR) on the top of MHC proteins

## Python script
USAGE

## Pymol plugin

For peptide-MHC class I-TCR complex:
tcr_mhci_pos selection [,mhc_a=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=180 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]

For peptide-MHC class II-TCR complex:
tcr_mhcii_pos selection [,mhc_a=None [,mhc_b=None [,tcr_a = None [,tcr_b = None [,mhc_a_init=1 [,mhc_a_final=80 [,mhc_b_init=1 [,mhc_b_final=90 [,tcr_a_init=1 [,tcr_a_final=109 [, tcr_b_init=1 [,tcr_b_final=116 [,state=None ]]]]]]]]]]]]]

## <u>NOTES</u>
- The code requires defining at least 3 chain IDs for MHCI-TCR or 4 chain IDs for MHCII-TCR
- If no range of residues is defined, the code uses the default values (see USAGE)
- Please make sure that the defined chains do not contain any solvent molecules
