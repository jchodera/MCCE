CONFLIST ASL        ASLBK ASL01 

NATOM    ASLBK      6
NATOM    ASL01      7

IATOM    ASLBK  N   0
IATOM    ASLBK  H   1
IATOM    ASLBK  CA  2
IATOM    ASLBK  HA  3
IATOM    ASLBK  C   4
IATOM    ASLBK  O   5
IATOM    ASL01  CB  0
IATOM    ASL01 1HB  1
IATOM    ASL01 2HB  2
IATOM    ASL01  CG  3
IATOM    ASL01  OD1 4
IATOM    ASL01  HD1 5
IATOM    ASL01  OD2 6

ATOMNAME ASLBK    0  N  
ATOMNAME ASLBK    1  H  
ATOMNAME ASLBK    2  CA 
ATOMNAME ASLBK    3  HA 
ATOMNAME ASLBK    4  C  
ATOMNAME ASLBK    5  O  
ATOMNAME ASL01    0  CB 
ATOMNAME ASL01    1 1HB 
ATOMNAME ASL01    2 2HB 
ATOMNAME ASL01    3  CG 
ATOMNAME ASL01    4  OD1
ATOMNAME ASL01    5  HD1
ATOMNAME ASL01    6  OD2






#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   ASL01      0
PKA      ASL01      0.0
ELECTRON ASL01      0
EM       ASL01      0.0
RXN      ASL01      -3.51

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  ASLBK  N   sp2       -1    C   0     CA  0     H
CONNECT  ASLBK  H   s         0     N
CONNECT  ASLBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  ASLBK  HA  s         0     CA
CONNECT  ASLBK  C   sp2       0     CA  0     O   1     N
CONNECT  ASLBK  O   sp2       0     C
CONNECT  ASL01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASL01 1HB  s         0     CB
CONNECT  ASL01 2HB  s         0     CB
CONNECT  ASL01  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASL01  OD1 sp3       0     CG  0     HD1
CONNECT  ASL01  HD1 s         0     OD1
CONNECT  ASL01  OD2 sp2       0     CG
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
H_DIHED  ASLBK  1    H    N    CA   C   60
H_DIHED  ASL01  1    HD1  OD1  CG   CB  60

#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
ACCEPTOR ASLBK  1    C  = O
ACCEPTOR ASL01  1    CG - OD2
DONOR    ASL01  1    OD1- HD1
DONOR    ASL02  1    OD2- HD2

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   ASL    N   1.55
RADIUS   ASL    H   1.20
RADIUS   ASL    CA  1.70
RADIUS   ASL    HA  1.20
RADIUS   ASL    C   1.70
RADIUS   ASL    O   1.52
RADIUS   ASL    CB  1.70
RADIUS   ASL   1HB  1.20
RADIUS   ASL   2HB  1.20
RADIUS   ASL    CG  1.70
RADIUS   ASL    OD1 1.52
RADIUS   ASL    HD1 1.20
RADIUS   ASL    OD2 1.52

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
#=========================================================================