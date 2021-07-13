import MDAnalysis as mda

u  = mda.Universe('step4.gro')
ag = u.select_atoms('prop x > 50')
ag.write('vis_step4.gro')

u  = mda.Universe('step5.gro')
ag = u.select_atoms('prop x > %f' % (u.dimensions[0]/2))
ag.write('vis_step5.gro')

