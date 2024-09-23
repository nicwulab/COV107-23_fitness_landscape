# For site 53 (site 52 in PDB numbering)

reinitialize
load ./pdb/S53P.pdb
bg_color white
hide all
set dynamic_measures, 0
show cartoon, all
color grey, all
run ./script/Mut_S53P_draw_links.py
run ./script/spectrumany.py
run ./script/data2bfactor.py

draw_links resi 31 and name CA and Chain A, resi 52 and name CA and Chain A, color=red, color2=red, radius=0.05, object_name=31:52_red

zoom all
hide labels
color red, max_frst_wm
color green, min_frst_wm

set_view (\
    -0.483052820,   -0.407135606,    0.775177956,\
    -0.744714379,    0.656650841,   -0.119185410,\
    -0.460497767,   -0.634859622,   -0.620398045,\
     0.000080943,   -0.000085436,  -99.509063721,\
    19.682664871,   31.850175858,   16.156400681,\
  -1048.714477539, 1247.738525391,  -20.000000000 )

hide dashes
hide labels

alter Chain A and resi 1-116, b=0
data2b_res Chain A and resi 1-116, ./data/S53P-FrustIndex.txt
spectrumany b, red gray80 green, Chain A and resi 1-116, minimum=-1.7, maximum=2
show sticks, Chain A and resi 52 and (not name c+o+n)
remove (hydro)

ray 1200,1000
png ./image/Site53-Mut.png
