#target preparation
LoadYOb (MacroDir)\YourModel_minimization.yob
Color Element
SplitObj 1
NiceOriAll
Cell Auto,Extension=7,Shape=Cuboid,Obj 2
FixAll
ForceField NOVA,SetPar=Yes
Boundary Wall
DelObj 2
SaveSce (MacroDir)\hRH4_receptor.sce
Clear
