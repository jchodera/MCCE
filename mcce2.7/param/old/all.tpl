CONFLIST ALL        ALLBK ALL01 

NATOM    ALLBK      6
NATOM    ALL01      4

IATOM    ALLBK  N   0
IATOM    ALLBK  H   1
IATOM    ALLBK  CA  2
IATOM    ALLBK  HA  3
IATOM    ALLBK  C   4
IATOM    ALLBK  O   5
IATOM    ALL01  CB  0
IATOM    ALL01 1HB  1
IATOM    ALL01 2HB  2
IATOM    ALL01 3HB  3

ATOMNAME ALLBK    0  N  
ATOMNAME ALLBK    1  H  
ATOMNAME ALLBK    2  CA 
ATOMNAME ALLBK    3  HA 
ATOMNAME ALLBK    4  C  
ATOMNAME ALLBK    5  O  
ATOMNAME ALL01    0  CB 
ATOMNAME ALL01    1 1HB 
ATOMNAME ALL01    2 2HB 
ATOMNAME ALL01    3 3HB 









#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   ALL01      0
PKA      ALL01      0.0
ELECTRON ALL01      0
EM       ALL01      0.0
RXN      ALL01      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  ALLBK  N   sp2       -1    C   0     CA  0     H
CONNECT  ALLBK  H   s         0     N
CONNECT  ALLBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  ALLBK  HA  s         0     CA
CONNECT  ALLBK  C   sp2       0     CA  0     O   1     N
CONNECT  ALLBK  O   sp2       0     C
CONNECT  ALL01  CB  sp3       0     CA  0    1HB  0    2HB  0    3HB
CONNECT  ALL01 1HB  s         0     CB
CONNECT  ALL01 2HB  s         0     CB
CONNECT  ALL01 3HB  s         0     CB
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
H_DIHED  ALLBK  1    H    N    CA   C    60
H_DIHED  ALL01  1   1HB   CB   CA   N    60

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   ALL    N   1.55
RADIUS   ALL    H   1.20
RADIUS   ALL    CA  1.70
RADIUS   ALL    HA  1.20
RADIUS   ALL    C   1.70
RADIUS   ALL    O   1.52
RADIUS   ALL    CB  1.70
RADIUS   ALL   1HB  1.20
RADIUS   ALL   2HB  1.20
RADIUS   ALL   3HB  1.20

#4.Rotomer
# None


