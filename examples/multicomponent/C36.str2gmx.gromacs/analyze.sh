#!/bin/bash

mkdir energy
echo 1  | gmx energy -f step7.edr -o energy/bond.xvg
echo 2  | gmx energy -f step7.edr -o energy/UB.xvg
echo 3  | gmx energy -f step7.edr -o energy/dihe.xvg
echo 4  | gmx energy -f step7.edr -o energy/impr.xvg
echo 5  | gmx energy -f step7.edr -o energy/LJ14.xvg
echo 6  | gmx energy -f step7.edr -o energy/CB14.xvg
echo 7  | gmx energy -f step7.edr -o energy/LJ.xvg
echo 8  | gmx energy -f step7.edr -o energy/CB.xvg
echo 9  | gmx energy -f step7.edr -o energy/CBr.xvg
echo 10 | gmx energy -f step7.edr -o energy/pot.xvg

