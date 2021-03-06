CONFLIST GLC        GLCBK 

NATOM    GLCBK      7

IATOM    GLCBK  N   0
IATOM    GLCBK  H   1
IATOM    GLCBK  CA  2
IATOM    GLCBK 1HA  3
IATOM    GLCBK 2HA  4
IATOM    GLCBK  O   5
IATOM    GLCBK  C   6

ATOMNAME GLCBK    0  N  
ATOMNAME GLCBK    1  H  
ATOMNAME GLCBK    2  CA 
ATOMNAME GLCBK    3 1HA 
ATOMNAME GLCBK    4 2HA 
ATOMNAME GLCBK    5  O  
ATOMNAME GLCBK    6  C  



#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  GLCBK  N   sp2       -1    C   0     CA  0     H
CONNECT  GLCBK  H   s         0     N
CONNECT  GLCBK  CA  sp3       0     N   0     C   0    1HA  0    2HA
CONNECT  GLCBK 1HA  s         0     CA
CONNECT  GLCBK 2HA  s         0     CA 
CONNECT  GLCBK  O   sp2       0     C 
CONNECT  GLCBK  C   sp2       0     CA  0     O   1     N

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   GLC    N   1.55
RADIUS   GLC    H   1.20
RADIUS   GLC    CA  1.70
RADIUS   GLC   1HA  1.20
RADIUS   GLC   2HA  1.20
RADIUS   GLC    O   1.52
RADIUS   GLC    C   1.70
