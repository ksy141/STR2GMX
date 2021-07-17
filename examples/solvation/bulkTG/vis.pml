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


load vis.gro
hide

cmd.select('TGhead', TGhead)
cmd.select('TGtail', TGtail)
select TIP3, name OH2
deselect

show  spheres, TGhead
color cTGhead, TGhead

show  spheres, TGtail
color cTGtail, TGtail

show  spheres, TIP3
color white,   TIP3

util.performance(100)
rebuild

rotate y, 90
