set orthoscopic, on
bg_color white

# 110, 181, 255
set_color cPLhead, [0.430, 0.710, 1.00]

# 0, 110, 219
set_color cPLtail, [0.000, 0.430, 0.86]

# 36, 255, 36
set_color cTGhead, [0.140, 1.000, 0.140]

# 255, 255, 0
set_color cTGtail, [1.000, 1.000, 0.000]

run ../../vis_sel.py


load step3.gro
hide

cmd.select('PLhead', PLhead)
cmd.select('PLtail', PLtail)
select TIP3, name OH2
deselect

show  spheres, PLhead
color cPLhead, PLhead

show  spheres, PLtail
color cPLtail, PLtail

show  spheres, TIP3
color white,   TIP3

util.performance(100)
rebuild

rotate x, 90
