# YASARA MACRO
# TOPIC:       5. Structure prediction
# TITLE:       Docking a ligand to a receptor
# REQUIRES:    Structure
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro runs VINA or AutoDock to predict the structure of a ligand-receptor complex. It can also continue a docking run that got interrupted. An analysis log file is written at the end.

# Parameter section - adjust as needed, but NOTE that some changes only take effect
# if you start an entirely new docking job, not if you continue an existing one. 
# =================================================================================
Processors CPUThreads=48,GPU=0
Antialias 0
Console Off

# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget='/home/myname/projects/docking/1sdf'

# Docking method, either AutoDockLGA or VINA
method='VINA'

# Number of docking runs (maximally 999, each run can take up to an hour)
runs=999

# Docking results usually cluster around certain hot spot conformations,
# and the lowest energy complex in each cluster is saved. Two complexes belong to
# different clusters if the ligand RMSD is larger than this minimum [A]:
rmsdmin=5.0

# Set to 1 to keep the ligand completely rigid (alternatively you can provide
# the ligand as a *.yob file and fix certain dihedral angles only).
rigid=0

# A selection of receptor residues whose side-chains should be kept flexible, e.g.
# flexres = 'Lys 91 Leu 100', or (if the receptor is not a monomer)
# flexres = 'Res Lys 91 Mol A or Res Leu 100 Mol B'.
# Alternatively you can provide the receptor as a *.sce or *.yob file with fixed
# atoms, which gives better control (e.g. you can keep only part of a side-chain
# or even a terminal backbone flexible). 
flexres='' 

# Set to 1 to calculate the binding energy spread in each cluster (takes longer)
calcspread=1

# Set to 1 if you also want to save PDB files of all clusters
pdbsaved=0

# YASARA commands follow below, no '=' is used
# --------------------------------------------

# Force field used for charge assignment (AutoDock only, not VINA) and to determine which ligand bonds can rotate.
# Do not change to an in vacuo force field (NOVA), all other force fields should work too (not tested).
ForceField AMBER03

# Uncomment below to set a certain random seed
# RandomSeed 1000

# Normally no change required below this point
# ============================================

# Sanity checks
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"
if runs>999
  RaiseError 'Too many docking runs selected, (runs) would take forever'

structlist='receptor','ligand'

# Load receptor and ligand
Clear
# Do we already have a receptor scene?
scefilename='(MacroTarget)_receptor.sce'
scene = FileSize (scefilename)
if !scene
  # No, load PDB or YOb files of receptor and ligand
  for struct in structlist
    (struct)=0
    for type in 'yob','pdb','sdf'
      filename='(MacroTarget)_(struct).(type)'
      exists = FileSize (filename)
      if exists and (type!='sdf' or struct=='ligand')
        (struct) = Load(type) (filename)
        break
    if not (struct)
      RaiseError 'The (struct) was not found at (filename). Please follow the instructions in the "Recipes" section exactly, especially when setting the macro target'
  # Orient the receptor, create a docking cell that includes the entire receptor 
  # and also has enough empty space around to accomodate the ligand (which has to
  # be removed temporarily, since 'Cell Auto' encloses the entire soup). 
  ligandradius = RadiusObj (ligand)  
  NiceOriObj (receptor)
  RemoveObj (ligand)
  Cell Auto,Extension=(ligandradius*2)
  AddObj (ligand)
else
  # Yes, a scene is present
  LoadSce (scefilename)
  # Verify that the cell is present
  simcell = CountObj SimCell
  if !simcell
    RaiseError 'If you provide a scene, it must contain a simulation cell, but none was found in (scefilename)'
  # The receptor is the object containing the first atom
  receptor = ListObj Atom 1
  if Objects>2
    ShowMessage 'Your scene contains (Objects) objects, while only the receptor and the docking cell are expected.'
    Wait ContinueButton
    HideMessage
  # Load the ligand
  ligand=0
  for type in 'yob','pdb','sdf'
    filename='(MacroTarget)_ligand.(type)'
    exists = FileSize (filename)
    if exists 
      ligand = Load(type) (filename)
      break
  if not ligand
    RaiseError 'The ligand was not found at (filename)'
  ligandradius = RadiusObj (ligand)
# Warn if the cell is too large
x,y,z = Cell
if method!='VINA' and (x>384 or y>384 or z>384)
  # AutoDock now uses dynamic memory allocation for grids up to 1024x1024x1024 (int MAX_GRID_PTS=1025 in autocomm.h)
  # Given the standard grid spacing of 0.375 A, these are 384 A. VINA doesn't have a compile-time limit
  ShowMessage "A cell axis is longer than 384 A, which may reduce docking accuracy. Consider focusing on the active site."
  Wait ContinueButton
  HideMessage
# Rename receptor and ligand objects to make sure they can be identified later
NameObj (receptor),Receptor
NameObj (ligand),Ligand
# Keep selected side-chains flexible
if flexres!=''
  FixObj Receptor
  FreeAtom Res (flexres) and Sidechain Obj Receptor
  FixRes Cys Atom SG with bond to Atom SG or Ala Pro and Obj Receptor
# The segment name C*** is used to tag ligand Conformations, and cannot be used for the receptor
hit = ListAtom Segment C??? Obj Receptor
if hit
  MarkAtom (hit)
  segname = SegAtom (hit)
  RaiseError 'Segment names starting with "C" are unfortunately reserved to identify ligand conformations, please click "Edit > Name > Segment" to rename the receptor segment "(segname)"' 
# The molecule name 'l' is reserved for the ligand
hit = ListAtom Mol l
if hit
  RaiseError 'The receptor must not contain a molecule named "l", please click Edit > Name > Molecule and try again'
# Docking is done without periodic boundaries
Boundary Wall
Longrange None  
# Do not show secondary structure for protein ligands, this slows things down
HideSecStrObj Ligand
ShowObj Ligand
StickObj Ligand
# Move the ligand out of the cell, so that it does not block the view
cellpos() = PosObj SimCell
PosObj (ligand),(cellpos1),(cellpos2+y*0.5+ligandradius),(cellpos3)
# Perform rigid docking if requested
if rigid
  FixObj Ligand
# Perform the docking
Experiment Docking
  Method (method)
  ReceptorObj (receptor)
  LigandObj (ligand)
  Runs (runs)
  ClusterRMSD (rmsdmin)
  # Result file number must be separated with '_', in case MacroTarget also ends with number.
  # If PDB files are preferred, replace .yob with .pdb
  ResultFile (MacroTarget)_001.yob
  # Uncomment below to set the number of energy evaluations in AutoDock (ga_num_evals):
  # DockPar ga_num_evals 25000000
  # Uncomment below to set the search exhaustiveness in Vina:
  # DockPar exhaustiveness 25
  # Uncomment the two lines below to provide your own atom parameters (VdW radii etc.)
  # GridPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # DockPar parameter_file /Path/To/Custom/AD4_parameters.dat
  # Uncomment below to keep temporary files in the current working directory:
  # TmpFileID 1adb_tmp
  # Uncomment below to stop after the setup step (combine with TmpFileID above to use YASARA only for setup)
  # SetupOnly yes
Experiment On
Wait ExpEnd
# Save a scene with receptor and all ligand conformations
SaveSce (MacroTarget)

# Save a log file with an analysis
Console Off
RecordLog (MacroTarget)
print 'Global docking result analysis'
print '=============================='
print
print '(runs) (method) docking runs of the ligand object (ligand) to the receptor object (receptor) yielded the following results,'
print 'sorted by binding energy [more positive energies indicate stronger binding, and negative energies mean no binding]'
print
print 'Run |Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
print '----+---------------------+---------------------+-----------------------------'
clusters=0
for i=001 to runs
  # Ligands have SegmentID C001, C002 etc.
  result = DockingResult i,'Obj (receptor)','Segment C(i)'
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
print '[They all differ by at least (rmsdmin) A heavy atom RMSD after superposing on the receptor]' 
print
print 'Clu |Bind.energy[kcal/mol]|Dissoc. constant [pM]| Contacting receptor residues'
print '----+---------------------+---------------------+-----------------------------'
for i=001 to clusters
  DelAll
  LoadYOb (MacroTarget)_(i)
  # Convert clusters to PDB format if requested
  if pdbsaved
    SavePDB 1,(MacroTarget)_(i)
  # Determine SegmentID of the ligand that seeds this cluster
  segnamelist(i) = SegAtom Segment C???
  # Ligands have SegmentID C001, C002 etc.
  result = DockingResult i,'Segment !C???','Segment C???'
  print (result)
print 
if calcspread
  # Additionally calculate the binding energy spread in each cluster (average and standard deviation)
  print 'While the table above lists the best binding energy in each cluster, it is sometimes'
  print 'helpful to also look at the energy spread [average and standard deviation], the' 
  print 'dissociation constant has been recalculated from the average binding energy:'
  print
  print 'Clu |Members|Bind.energy spread [kcal/mol]|Dissoc. constant [pM]'
  print '----+-------+-----------------------------+---------------------'
  LoadSce (MacroTarget)
  for i=001 to clusters
    # Collect the binding energies of all poses that belong to the current cluster
    ShowMessage 'Calculating energy spread, (100*i/clusters)% completed...'
    Wait 1
    bindnrglist()=0
    lastbindnrg=-1e99
    members=0
    for j=001 to runs
      if !count assigned(j)
        r = RMSDAtom Segment C(j) Element !H Obj Ligand,Segment (segnamelist(i)) Element !H Obj Ligand
        if r<rmsdmin
          # Found a member of this cluster
          bindnrg = BFactorAtom Segment C(j)
          if bindnrg!=lastbindnrg
            # Not a duplicate pose
            members=members+1
            bindnrglist(members)=bindnrg
          lastbindnrg=bindnrg
          assigned(j)=1
    # Calculate the mean dissociation constant: Ki = exp(deltaG/(R*T))
    bindnrg=mean bindnrglist
    if bindnrg<=0
      dissconst='               None'
    else
      # 4184. converts from kcal/mol to J/mol, AvoConst*BoltzConst/JToUnit is R in [J/(K*mol)] 
      dissconst= 00000000000000.0000+exp (-bindnrg*4184./(AvoConst*(BoltzConst/JToUnit)*298.15))*1e12
    # Print the spread
    print '(000+i) |  (000+members)  |    (000000.0000+bindnrg)+-(000000.0000+stddev bindnrglist) | (dissconst)'
print 
rndseed = RandomSeed
print 'The current random seed is (rndseed).'
fof = ForceField
print 'Point charges and dihedral barriers were obtained from the (fof) force field.'
print 'The results of the (runs) runs have been combined in a single scene at (MacroTarget).sce',Convert=No
print 'The (clusters) clusters have been saved as:'
LoadSce (MacroTarget)
for i=001 to clusters
  print '(MacroTarget)_(i).yob',Convert=No
StopLog
HideMessage
Console On
# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit

# EXTRACT DOCKING RESULT
# ======================
# Returns a result string and binding energy, extracted from atomic BFactor and Property.
# 'num' is a sequential number, 'recsel' selects the receptor, 'ligsel' the ligand.
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


