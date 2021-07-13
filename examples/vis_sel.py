HEAD = ['N', 'C12', 'H12A', 'H12B', 
        'C13', 'H13A', 'H13B', 'H13C', 
        'C14', 'H14A', 'H14B', 'H14C', 
        'C15', 'H15A', 'H15B', 'H15C', 
        'C11', 'H11A', 'H11B', 
        'P', 'O13', 'O14', 'O12', 'O11', 
        'C1', 'HA', 'HB', 
        'C2', 'HS', 
        'O21', 'C21', 'O22', 
        'C3',  'HX',  'HY', 
        'O31', 'C31', 'O32', 
        'HN1', 'HN2', 'HN3', 
        'H2', 'O2', 'HO2', 
        'H3', 'O3', 'HO3', 
        'H4', 'O4', 'HO4', 'H5', 
        'O5', 'HO5', 'C16', 
        'H6', 'O6', 'HO6', 'H1']

TAIL = []
TAIL += ['C2%d' %n for n in range(2, 23)]
TAIL += ['H%dR' %n for n in range(2, 23)]
TAIL += ['H%dS' %n for n in range(2, 23)]
TAIL += ['C3%d' %n for n in range(2, 23)]
TAIL += ['H%dX' %n for n in range(2, 23)]
TAIL += ['H%dY' %n for n in range(2, 23)]
TAIL += ['H20T', 'H18Z', 'H91', 'H101', 'H18T', 'H16Z']

for i in range(2, 23):
    TAIL.append('C2' + str(i))
    TAIL.append('H' + str(i) + 'R')

PLhead = 'resname POPC+DOPE+SAPI and name ' + '+'.join(HEAD)
PLtail = 'resname POPC+DOPE+SAPI and name ' + '+'.join(TAIL)

print(PLhead)
print('\n')
print(PLtail)
print('\n')

TGhead = 'resname TRIO+TRIN and name C1+C2+C3+C11+C21+C31+O11+O21+O31+O12+O22+O32+HA+HB+HX+HY+HS'
TTAIL = []
TTAIL += ['C1%d' %n for n in range(2, 19)]
TTAIL += ['H%dA' %n for n in range(2, 19)]
TTAIL += ['H%dB' %n for n in range(2, 19)]
TTAIL += ['C2%d' %n for n in range(2, 19)]
TTAIL += ['H%dR' %n for n in range(2, 19)]
TTAIL += ['H%dS' %n for n in range(2, 19)]
TTAIL += ['C3%d' %n for n in range(2, 19)]
TTAIL += ['H%dX' %n for n in range(2, 19)]
TTAIL += ['H%dY' %n for n in range(2, 19)]
TTAIL += ['H18C', 'H18T', 'H18Z']
TGtail = 'resname TRIO+TRIN and name ' + '+'.join(TTAIL)

print(TGhead)
print('\n')
print(TGtail)
print('\n')

