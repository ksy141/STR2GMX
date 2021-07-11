import numpy as np

elements = ['LP', 'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
            'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc', 
            'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 
            'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 
            'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 
            'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 
            'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 
            'Ra', 'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 
            'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 
            'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

masses = np.array([  0.        ,   1.008     ,   4.002602  ,   6.94      ,
                     9.0121831 ,  10.81      ,  12.011     ,  14.007     ,
                    15.999     ,  18.99840316,  20.1797    ,  22.98976928,
                    24.305     ,  26.9815385 ,  28.085     ,  30.973762  ,
                    32.06      ,  35.45      ,  39.948     ,  39.0983    ,
                    40.078     ,  44.955908  ,  47.867     ,  50.9415    ,
                    51.9961    ,  54.938044  ,  55.845     ,  58.933194  ,
                    58.6934    ,  63.546     ,  65.38      ,  69.723     ,
                    72.63      ,  74.921595  ,  78.971     ,  79.904     ,
                    83.798     ,  85.4678    ,  87.62      ,  88.90584   ,
                    91.224     ,  92.90637   ,  95.95      ,  97.90721   ,
                   101.07      , 102.9055    , 106.42      , 107.8682    ,
                   112.414     , 114.818     , 118.71      , 121.76      ,
                   127.6       , 126.90447   , 131.293     , 132.90545196,
                   137.327     , 138.90547   , 140.116     , 140.90766   ,
                   144.242     , 144.91276   , 150.36      , 151.964     ,
                   157.25      , 158.92535   , 162.5       , 164.93033   ,
                   167.259     , 168.93422   , 173.045     , 174.9668    ,
                   178.49      , 180.94788   , 183.84      , 186.207     ,
                   190.23      , 192.217     , 195.084     , 196.966569  ,
                   200.592     , 204.38      , 207.2       , 208.9804    ,
                   209.        , 210.        , 222.        , 223.        ,
                   226.        , 227.        , 232.0377    , 231.03588   ,
                   238.02891   , 237.        , 244.        , 243.        ,
                   247.        , 247.        , 251.        , 252.        ,
                   257.        , 258.        , 259.        , 262.        ,
                   267.        , 268.        , 271.        , 274.        ,
                   269.        , 276.        , 281.        , 281.        ,
                   285.        , 286.        , 289.        , 288.        ,
                   293.        , 294.        , 294.        ])

