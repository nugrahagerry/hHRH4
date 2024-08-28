#target preparation
LoadYOb (MacroDir)\519.yob
SplitObj 1
NiceOriAll
Cell Auto,Extension=7,Shape=Cuboid,Obj 2
FixAll
ForceField NOVA,SetPar=Yes
Boundary Wall
DelObj 2
SaveSce (MacroDir)\hrh4_receptor.sce
Clear
