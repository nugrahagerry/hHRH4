#konvertion pdb to sce
for k=00001 to 01000
  LoadPDB (MacroDir)\aktif_344(k).pdb,Center=No
  SplitObj 1
  JoinObj 3-278,1
  NiceOriAll
  Cell Auto,Extension=7,Shape=Cuboid,Obj 2
  FixAll
  ForceField NOVA,SetPar=Yes
  Boundary Wall
  SaveSce (MacroDir)\aktif_344(k)_complex.sce
  Clear

#Calculation
method = 'VINALS'
runs = 1
rmsdmin = 5.0
rigid = 1
for j=00001 to 01000
  LoadSce (MacroDir)\aktif_344(j)_complex.sce
  NameObj 1,receptor
  NameObj 2,ligand
  ForceField AMBER03
  Boundary Wall
  Longrange None
#FixAll

  Experiment Docking
    Method (method)
    ReceptorObj (1)
    LigandObj 2
    Runs (runs)
    ClusterRMSD (rmsdmin)
  # Result file number must be separated with '_', in case MacroTarget also ends with number
    ResultFile temp_001
  # Uncomment below to set the number of energy evaluations (AutoDock ga_num_evals):
  # DockPar ga_num_evals 25000000
  # Uncomment the two lines below to provide your own atom parameters (VdW radii etc.)
  # GridPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # DockPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # If you want to keep all temporary files, uncomment below
  # TmpFileID YourChoice
  Experiment On
  Wait ExpEnd
  # Save a scene with receptor and all ligand conformations
  SaveSce (MacroDir)\aktif_344(j)_result.sce
  # Save a log file with an analysis

  Console Off
  RecordLog (MacroDir)\aktif_344(j)_result.log
  print 'Local docking result analysis'
  print '============================='
  print
  print '(runs) (method) docking runs of the ligand object 2 to the receptor object 1 yielded the following results,'
  print 'sorted by binding energy: [More positive values indicates stronger binding, and negative values mean no binding]'
  print
  print 'Run |Bind.energy[kcal/mol]|Dissoc.constant[pM]| Contacting receptor residues'
  print '----+---------------------+------------------+-----------------------------'
  clusters=0
  for i=001 to runs
    # Ligands have SegmentID C001, C002 etc.
    result = DockingResult i,'Obj 1','Segment C(i)'
    print (result)
  # Count the clusters
    exists = FileSize (MacroTarget)_(i).yob
    if exists
      clusters=clusters+1
    if !(i%10)
      ShowMessage 'Analyzing docking results, (100*i/runs)% completed...'
      Wait 1
  print
  print 'After clustering the (runs) runs, the following (clusters) distinct complex conformations were found:'
  print '[They all differ by at least (rmsdmin) A heavy atom RMSD]' 
  print
  print 'Cluster|Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
  print '-------+---------------------+---------------------+-----------------------------'
  for i=001 to clusters
    DelAll
    LoadYOb (MacroTarget)_(i)
    # Ligands have SegmentID C001, C002 etc.
    result = DockingResult i,'Segment !C???','Segment C???'
    print (result)
  print 
  print 'The results of the (runs) runs have been combined in a single scene at (MacroTarget).sce',Convert=No
  print 'The (clusters) clusters have been saved as:'
  for i=001 to clusters
    print '(MacroTarget)_(i).yob',Convert=No
  StopLog
  HideMessage
  Console On
  LoadSce (MacroDir)\aktif_344(j)_result.sce
  # Exit YASARA if this macro was provided as command line argument in console mode
  if runWithMacro and ConsoleMode
    Exit 

# EXTRACT DOCKING RESULT
# ======================
# Returns a result string. 'num' is a sequential number, 'recsel' selects the receptor, 'ligsel' the ligand.
  def DockingResult num,recsel,ligsel
  # Get the binding energy from the B-factor...
    bindnrg = BFactorAtom (ligsel)
  # ...and the dissociation constant from the atomic property.
    dissconst = PropAtom (ligsel)
    if bindnrg<=0
      dissconst='               None'
    else
      dissconst= 00000000000000.0000+dissconst
  # Get receptor residues that contact the ligand
    reslist() = ListRes (recsel) with distance<4.0 from (ligsel),Format='MOLNAME RESNAME RESNUM'
    if !count reslist
      reslist='No residues in contact'
    else
      reslist=join reslist
    result='(000+num) |         (000000.0000+bindnrg) | (dissconst) | (reslist)'
    return result
  Clear
