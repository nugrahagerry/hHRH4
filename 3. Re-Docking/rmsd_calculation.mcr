#rmsd_calculation
for l=001 to 1000
  LoadYOb (MacroDir)\hrh4_ligandref.yob
  LoadYOb (MacroDir)\hrh4_(l)_001.yob
  FreeAll 
  SplitAll Center=No,Keep=ObjNum 
  DelObj 2
  NumberObj 3,First=2
  DelAtom Element H 
  TransferObj 2,1,Local=Match 
  RecordLog (MacroDir)\hrh4_(l).rmsd.log 
  RMSDObj 2,1,Match=AltLoc,Flip=Yes,Unit=Obj
  StopLog 
  Clear