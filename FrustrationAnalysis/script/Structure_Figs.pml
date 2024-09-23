reinitialize
load ./pdb/6xc2.pdb
set seq_view, 1
bg_color white
hide all
show cartoon, chain A
color grey80, chain A 

load ./pdb/WT.pdb

align WT, 6xc2 

color lightblue, WT and chain A 
color pink, WT and chain B

#Site 27 (26 in PDB numbering). 5.7 Angstroms between Asn 487 of Spike and Phe 26 of Heavy Chain

show sticks, 6xc2 and Chain A and resi 487 and (not name c+o+n)
show sticks, WT and Chain A and resi 26 and (not name c+o+n)
pseudoatom pi_center, /WT//A/PHE`26/CG+CZ
distance dist1, pi_center, /6xc2//A/ASN`487/ND2

#Site 53
show sticks, WT and chain A and resi 52 and (not name c+o+n)
show sticks, 6xc2 and chain A and resi 473 and (not name c+o+n)
distance dist2, /WT//A/SER`52/OG, /6xc2//A/TYR`473/OH

hide labels

util.cnc all

set dash_color, black

remove (hydro)

set_view (\
    -0.192667186,   -0.160288721,    0.968086123,\
    -0.802904010,   -0.541401327,   -0.249434665,\
     0.564104438,   -0.825334489,   -0.024387654,\
    -0.000449432,    0.000398032,  -48.377693176,\
    19.739419937,  -17.569862366,   39.529762268,\
    33.187007904,   63.567756653,  -20.000000000 )

ray 1200,1000
png ./image/Spike_S53.png



reinitialize

load ./pdb/6xc2.pdb
set seq_view, 1
bg_color white
hide all
show cartoon, chain A
color grey80, chain A 

load ./pdb/F27I.pdb
 
align F27I, 6xc2

color lightblue, F27I and chain A 
color pink, F27I and chain B

show sticks, 6xc2 and Chain A and resi 487 and (not name c+o+n)
show sticks, F27I and Chain A and resi 26 and (not name c+o+n)

hide labels

util.cnc all

set dash_color, black

remove (hydro)

set_view (\
    -0.192667186,   -0.160288721,    0.968086123,\
    -0.802904010,   -0.541401327,   -0.249434665,\
     0.564104438,   -0.825334489,   -0.024387654,\
    -0.000449432,    0.000398032,  -48.377693176,\
    19.739419937,  -17.569862366,   39.529762268,\
    33.187007904,   63.567756653,  -20.000000000 )

ray 1200,1000
png ./image/Spike_I27.png



reinitialize

load ./pdb/6xc2.pdb
set seq_view, 1
bg_color white
hide all
show cartoon, chain A
color grey80, chain A 

load ./pdb/F27L.pdb

align F27L, 6xc2

color lightblue, F27L and chain A 
color pink, F27L and chain B

show sticks, 6xc2 and Chain A and resi 487 and (not name c+o+n)
show sticks, F27L and Chain A and resi 26 and (not name c+o+n)

hide labels

util.cnc all

set dash_color, black

remove (hydro)

set_view (\
    -0.480050117,    0.672069490,    0.563805938,\
    -0.599620461,   -0.720502079,    0.348308682,\
     0.640311360,   -0.170865074,    0.748860538,\
    -0.000484694,   -0.000015028,  -48.399753571,\
    21.457820892,  -14.275853157,   36.177444458,\
    31.078458786,   65.676307678,  -20.000000000 )

ray 1200,1000
png ./image/Spike_L27.png



reinitialize

load ./pdb/6xc2.pdb
set seq_view, 1
bg_color white
hide all
show cartoon, chain A
color grey80, chain A 

load ./pdb/F27V.pdb

align F27V, 6xc2

color lightblue, F27V and chain A 
color pink, F27V and chain B

show sticks, 6xc2 and Chain A and resi 487 and (not name c+o+n)
show sticks, F27V and Chain A and resi 26 and (not name c+o+n)

hide labels

util.cnc all

set dash_color, black

remove (hydro)

set_view (\
    -0.480050117,    0.672069490,    0.563805938,\
    -0.599620461,   -0.720502079,    0.348308682,\
     0.640311360,   -0.170865074,    0.748860538,\
    -0.000484694,   -0.000015028,  -48.399753571,\
    21.457820892,  -14.275853157,   36.177444458,\
    31.078458786,   65.676307678,  -20.000000000 )

ray 1200,1000
png ./image/Spike_V27.png



reinitialize

load ./pdb/6xc2.pdb
set seq_view, 1
bg_color white
hide all
show cartoon, chain A
color grey80, chain A 

load ./pdb/S53P.pdb

align S53P, 6xc2

color lightblue, S53P and chain A 
color pink, S53P and chain B

show sticks, S53P and chain A and resi 52 and (not name c+o+n)
show sticks, S53P and chain A and resi 55 and (not name c+o+n)
show sticks, 6xc2 and chain A and resi 460 and (not name c+o+n)
distance dist3, /S53P//A/SER`55/OG, /6xc2//A/ASN`460/ND2
distance dist4, /S53P//A/SER`55/OG, /6xc2//A/ASN`460/OD1

hide labels

util.cnc all

set dash_color, black

remove (hydro)

set_view (\
    -0.781787813,    0.407293737,    0.472133785,\
    -0.576383471,   -0.183212638,   -0.796366334,\
    -0.237854674,   -0.894721389,    0.377991199,\
     0.000005461,    0.000700896,  -49.111980438,\
     9.148626328,  -24.010543823,   39.925968170,\
    38.333068848,   59.682811737,  -20.000000000 )

ray 1200,1000
png ./image/Spike_P53.png
