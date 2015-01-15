CONFLIST MET        METBK MET01

NATOM    METBK      6
NATOM    MET01      11

IATOM    METBK  N   0
IATOM    METBK  H   1
IATOM    METBK  CA  2
IATOM    METBK  HA  3
IATOM    METBK  C   4
IATOM    METBK  O   5
IATOM    MET01  CB  0
IATOM    MET01 1HB  1
IATOM    MET01 2HB  2
IATOM    MET01  CG  3
IATOM    MET01 1HG  4
IATOM    MET01 2HG  5
IATOM    MET01  SD  6
IATOM    MET01  CE  7
IATOM    MET01 1HE  8
IATOM    MET01 2HE  9
IATOM    MET01 3HE  10

ATOMNAME METBK    0  N
ATOMNAME METBK    1  H
ATOMNAME METBK    2  CA
ATOMNAME METBK    3  HA
ATOMNAME METBK    4  C
ATOMNAME METBK    5  O
ATOMNAME MET01    0  CB
ATOMNAME MET01    1 1HB
ATOMNAME MET01    2 2HB
ATOMNAME MET01    3  CG
ATOMNAME MET01    4 1HG
ATOMNAME MET01    5 2HG
ATOMNAME MET01    6  SD
ATOMNAME MET01    7  CE
ATOMNAME MET01    8 1HE
ATOMNAME MET01    9 2HE
ATOMNAME MET01   10 3HE


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   MET01      0
PKA      MET01      0.0
ELECTRON MET01      0
EM       MET01      0.0
RXN      MET01      -0.1

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  METBK  N   sp2       -1    C   0     CA  0     H
CONNECT  METBK  H   s         0     N
CONNECT  METBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  METBK  HA  s         0     CA
CONNECT  METBK  C   sp2       0     CA  0     O   1     N
CONNECT  METBK  O   sp2       0     C
CONNECT  MET01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  MET01 1HB  s         0     CB
CONNECT  MET01 2HB  s         0     CB
CONNECT  MET01  CG  sp3       0     CB  0     SD  0    1HG  0    2HG
CONNECT  MET01 1HG  s         0     CG
CONNECT  MET01 2HG  s         0     CG
CONNECT  MET01  SD  sp3       0     CG  0     CE
CONNECT  MET01  CE  sp3       0     SD  0    1HE  0    2HE  0    3HE
CONNECT  MET01 1HE  s         0     CE
CONNECT  MET01 2HE  s         0     CE
CONNECT  MET01 3HE  s         0     CE
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I

#3.Atom Parameters: Partial Charges and Radii
CHARGE   METBK  N    -0.350
CHARGE   METBK  H     0.250
CHARGE   METBK  CA    0.100
CHARGE   METBK  C     0.550
CHARGE   METBK  O    -0.550
CHARGE   MET01  CG    0.060
CHARGE   MET01  SD   -0.120
CHARGE   MET01  CE    0.060

# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   MET    N   1.50
RADIUS   MET    H   1.00
RADIUS   MET    CA  2.00
RADIUS   MET    HA  0.00
RADIUS   MET    C   1.70
RADIUS   MET    O   1.40
RADIUS   MET    CB  2.00
RADIUS   MET   1HB  0.00
RADIUS   MET   2HB  0.00
RADIUS   MET    CG  2.00
RADIUS   MET   1HG  0.00
RADIUS   MET   2HG  0.00
RADIUS   MET    SD  1.85
RADIUS   MET    CE  2.00
RADIUS   MET   1HE  0.00
RADIUS   MET   2HE  0.00
RADIUS   MET   3HE  0.00

#=========================================================================
#        GRP   #      BOND     AFFECTED_ATOMS
#123456789012345678901234567890
#-------|---|----|-|---------|----|----|----|----|----|----|----|----|----
ROTAMER  MET   0     CA - CB   CG   SD   CE
ROTAMER  MET   1     CB - CG   SD   CE
ROTAMER  MET   2     CG - SD   CE
#=========================================================================

