import MDAnalysis as mda

u  = mda.Universe('step4.gro')
ag = u.select_atoms('prop x > 50')
ag.write('vis_step4.gro')

u  = mda.Universe('step5.gro')
ag = u.select_atoms('prop x > %f' % (u.dimensions[0]/2))
ag.write('vis_step5.gro')

u  = mda.Universe('step6.gro')
ag1 = u.select_atoms('resname POPC TRIO and prop x > %f' % (u.dimensions[0]/2))
ag2 = u.select_atoms('resname TIP3 and prop x > %f' % (u.dimensions[0]/2 + 10))
mda.Merge(ag1, ag2).atoms.write('vis_step6.gro')

