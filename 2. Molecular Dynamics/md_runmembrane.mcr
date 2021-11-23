# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Running a molecular dynamics simulation of a membrane protein with normal or fast speed
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro sets up and runs a simulation of a membrane protein. It scans the protein for secondary structure elements with hydrophobic surface residues, orients it accordingly and embeds it in a membrane of adjustable lipid composition. Finally a 250 ps restrained equilibration simulation is run, which ensures that the membrane can adapt to the newly embedded protein. Then the real simulation starts.

# Parameter section - adjust as needed, but NOTE that some changes only take
# effect if you start an entirely new simulation, not if you continue an existing one. 
# ====================================================================================

# The structure to simulate must be present with a .pdb or .sce extension.
# If a .sce (=YASARA scene) file is present, the membrane and cell must have been added.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget = 'c:\MyProject\1crn'
# Extension of the cell on each side of the protein in the membrane plane (=XZ plane)
# '15' means that the membrane will be 30 A larger than the protein
memextension=15

# Extension of the cell on each side of the protein along the third (water) axis (=Y-axis)
# '10' means that the cell will be 20 A higher than the protein
waterextension=10

# Flag to use a square membrane. This makes sure that also elongated proteins
# embedded in the membrane can rotate freely during very long simulations. If
# only a short simulation is planned, it can be speeded up by setting the flag
# to 0, creating a rectangular membrane that fits the solute more tightly.
square=1

# Membrane composition: The percentages of phosphatidyl-ethanolamine (PEA), phosphatidyl-
# choline (PCH, also known as POPC) and phosphatidyl-serine (PSE), must sum up to 100 in
# each column. All lipids are 1-palmitoyl, 2-oleoyl by default.
# The left column is for the bottom side of the membrane, the right column is for the top side.
# When YASARA shows you the suggested membrane embedding, you need to check that the protein
# orientation matches the membrane composition. If not, flip left and right columns below and
# rerun the macro. Note that PCH has a large headgroup which cannot form hydrogen bonds, and
# thus reduces membrane stability. PEA is the most stable membrane lipid.
# Percentage of phosphatidyl-ethanolamine, bottom and top membrane side
PEApercent=100,100
# Percentage of phosphatidyl-choline, bottom and top membrane side
PCHpercent=  0,  0
# Percentage of phosphatidyl-serine, bottom and top membrane side
PSEpercent=  0,  0
# Or uncomment below to use your own membrane template with 10x10 lipids on each side,
# see membrane simulation recipes for details. (If usermemname='YourChoice', the membrane
# must be saved as yasara/yob/membrane_YourChoice.yob)
#usermemname='YourChoice'
#usermemsize=77.21,73.24

# pH at which the simulation should be run, by default physiological pH 7.4.
ph=7.4

# The ion concentration as a mass fraction, here we use 0.9% NaCl (physiological solution)
ions='Na,Cl,0.9'

# Forcefield to use (this is a YASARA command, so no '=' used)
ForceField AMBER14

# Simulation temperature, which also serves as the random number seed (see Temp command).
# If you increase the temperature significantly by X%, you also need to reduce the timestep by X%
# by changing the 'tslist' that matches your speed below.
temperature='298K'

# Pressure at which the simulation should be run [bar].
pressure=1 

# Cutoff
cutoff=8

# Equilibration period in picoseconds:
# During this initial equilibration phase, the membrane is artificially stabilized
# so that it can repack and cover the solute, while solvent molecules are kept outside.
equiperiod=250

# Delay for animations, 1=maximum speed
delay=100

# The format used to save the trajectories: YASARA 'sim', GROMACS 'xtc' or AMBER 'mdcrd'.
# If you don't pick 'sim', a single *.sim restart file will be saved too, since the other
# two formats don't contain velocities, only positions.
format='sim'

# Duration of the complete simulation, must be longer than equiperiod above.
# Alternatively use e.g. duration=5000 to simulate for 5000 picoseconds
# 'if !count duration' simply checks if variable 'duration' as been defined previously (e.g. by an including macro)
if !count duration
  duration='100000'

# The simulation speed, either 'slow' (2*1 fs timestep), 'normal' (2*1.25 fs timestep) or
# 'fast' (maximize performance with 2*2.5 fs timestep and constraints)
# Do not use 'fast' if you simulate incorrect molecules (that would not be stable in reality) 
# 'if !count speed' simply checks if variable 'speed' as been defined previously (e.g. by an including macro) 
if !count speed
  speed='normal'

# The save interval for snapshots. Normally you don't need more than 500-1000 snapshots
# of your simulation, since that's the resolution limit of a typical figure in a journal.
if speed=='fast'
  # Fast speed, save simulation snapshots every 250000 fs, i.e. 250 ps.
  saveinterval=250000
else  
  # Slow or normal speed, save simulation snapshots every 100000 fs, i.e. 100 ps.
  saveinterval=100000

# Flag if the protein's membrane embedding needs to be confirmed
confirmneeded=1

# Normally no change required below this point
# ============================================

RequireVersion 15.1.1

# Treat all simulation warnings as errors that stop the macro
WarnIsError On

# Membrane simulations are always periodic
Boundary periodic

# The membrane values below should not be changed, since they only describe
# membrane characteristics that originate from the simulation: during the
# simulation, the membrane reaches an equilibrium, which depends only on
# the lipid composition, the embedded protein and the force field parameters,
# but not on the initial setup. It is thus impossible to say 'I want a membrane
# 39 A thick', because if you build such a membrane, it will adapt to its
# inherent thickness during the simulation.  

# Membrane height in A, density should be ~0.861 g/ml 
memheight=43.
# Height of the hydrophobic membrane core
memcoreheight=28.
      
# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

# When run as a macro in text mode, add configuration details to log file
if runWithMacro and ConsoleMode
  Processors

Clear
Console off
SimSteps 1,1
# Do we already have a scene with water?
waterscene = FileSize (MacroTarget)_water.sce
if waterscene
  LoadSce (MacroTarget)_water
else
  # No scene with water present yet
  # Do we have a scene with the protein embedded in the membrane?
  scene = FileSize (MacroTarget).sce
  if scene
    # Yes, protein with membrane is already present
    LoadSce (MacroTarget)
    # Search for the membrane
    c = CountObj Membrane
    if c!=1
      # No membrane object, this must be a user-provided initial scene to prevent modifications like Clean/OptHyd, load later
      scene=0
  if !scene
    # No membrane scene present yet
    # Has the user already oriented the protein inside the membrane interactively?
    oriscene = FileSize (MacroTarget)_ori.sce
    if oriscene
      LoadSce (MacroTarget)_ori
      Unselect
    else
      # No orientation scene present yet
      if ConsoleMode and confirmneeded
        RaiseError "The membrane embedding needs to be confirmed visually. Please start YASARA in graphics mode, then continue in text mode"
      # Load the protein, assuming it's a SCE, PDB or YOB file
      for filetype in 'sce','yob','pdb'
        size = FileSize (MacroTarget).(filetype)
        if size
          break
      if !size
        RaiseError 'Initial structure not found, expected (MacroTarget).pdb or .yob. Make sure to create a project directory and place the structure there in PDB or YOB format'
      # Load structure
      Load(filetype) (MacroTarget)
      DelObj SimCell
      nameclash = CountObj Membrane
      if nameclash
        NameObj all,Solute
      # In case user accidentally provided a YOb file with selected atoms
      Unselect
      # Align object with major axes to minimize cell size
      NiceOriAll
      # As initial guess, assume the major axis is perpendicular to the membrane
      Style Ribbon,Stick
      ShowAtom fixed
      NiceOri Axis1=Y,Axis2=X
      r = Radius
      AutoPos Z=(r*2),Steps=(delay)
      if filetype!='sce'
        # Delete long peptide bonds that bridge gaps in the structure, which tells CleanAll to add ACE/NME
        # capping groups (the structure of the missing residues could also be predicted, see LoadPDB docs).
        DelBond N,C,LenMin=5
        # Delete waters that are not involved in metal binding and not fixed, to help the calculation of binding energies 
        DelRes !fixed Water with 0 arrows to all
        # Make sure that all hydrogens are present etc. and adapt protonation states to chosen pH
        CleanAll
        pH (ph)
        # pH command recreated water hydrogens, free them
        FreeAtom Element H
      # Now search for hydrophobic transmembrane helices and strands
      ShowMessage "Searching for transmembrane elements..."
      Wait 1
      # To find exposed hydrophobic residues, we first need the maximum side-chain surface areas
      SurfPar Probe=1.4,Resolution=3,Radii=Solvation
      namelist()='Ala','Val','Leu','Ile','Met','Phe','Tyr','Trp'
      for name in namelist
        obj = BuildRes (name)
        CleanObj (obj)
        AddEnvObj (obj)
        surfmax(name) = SurfAtom Sidechain Obj (obj),Accessible,Unit=Obj
        DelObj (obj)
      # Get a list of hydrophobic side-chain surface areas
      AddEnv
      reslist() = ListRes (join namelist) with contribution to accessible surface of all
      if count reslist
        # Found hydrophobic residues contributing to the surface
        surflist() = SurfAtom Res (join reslist) Sidechain,Accessible,Unit=Res
        # Get a list of all heavily exposed (>30%) hydrophobic side-chains
        for i=1 to count reslist
          name = NameRes (reslist(i))
          if surflist(i)>surfmax(name)*0.3
            # Create list of exposed hydrophobic residues, explist
            explist(count explist+1)=reslist(i)
      # Prepare to calculate the major protein axis, normal to the membrane plane
      axis=0.,0.,0.
      if count explist>5
        # First search for transmembrane helices, then for beta strands (beta barrels)
        for type in 'helix','strand'
          secstrelements=0
          # Get a list of all helices or strands
          atomlist() = ListAtom SecStr (type),compress=yes
          for i=1 to count atomlist step 2
            # Select CA atom of first and last residue
            first = ListAtom CA Res Atom (atomlist(i))
            last = ListAtom CA Res Atom (atomlist(i)+atomlist(i+1)-1)
            totresidues = CountRes Atom (first)-(last)
            hydresidues = CountRes (join namelist) Atom (first)-(last)
            expresidues = CountRes Atom (first)-(last) and (join explist)
            print '(type) (first)-(last) is (totresidues) residues long, of which (hydresidues) are hydrophobic and (expresidues) are exposed hydrophobic ones' 
            if (type=='helix' and totresidues>16 and hydresidues>7 and expresidues>3) or
               (type=='strand' and totresidues>8 and hydresidues>3 and expresidues>1)
              # Found transmembrane element
              secstrelements=secstrelements+1
              reslist() = ListRes Atom (first) (last),Format='RESName MOLNAME RESNUM'
              ShowMessage 'Found transmembrane (type) (secstrelements) spanning residues (reslist1) to (reslist2)...' 
              MarkAtom (first),(last)
              ColorRes Atom (first)-(last),Yellow
              Wait (delay)
              # Add direction vector to major axis
              px,py,pz,dir() = GroupLine (first)-(last) and CA
              if sum (axis*dir)<0
                # Scalar product is negative, turn direction around
                dir()=-dir
              axis()=axis+dir
              #p1,p2,p3,dir() = GroupLine (first)-(last) and CA,CoordSys=Global
              #ShowArrow Start=Point,(p),End=Point,(p+dir*30)
          if norm axis!=0
            break
        MarkAtom None
        if type=='strand' and secstrelements<5
          ShowMessage "Not enough transmembrane secondary structure elements found that could help to tune the protein orientation."
          axis=0.,0.,0.
          Wait (delay)
      else
        ShowMessage "Not enough exposed hydrophobic residues found. Looking for other features..."
        Wait (delay)
      if not norm axis
        # No transmembrane axis found, look for lipid anchors
        lipid='SmilesMEMBER CH2?CH2CH2CH2CH2?'
        liplist() = ListRes (lipid)
        lipatoms = CountAtom (lipid) 
        for lip in liplist
          id = ListRes (lip),Format="RESNAME RESNUM"
          first,last=SpanAtom Res (lip)
          ShowMessage 'Found transmembrane lipid anchor (id)...' 
          MarkAtom (first),(last)
          ColorAtom (first)-(last),Yellow
          Wait (delay/2)
          px,py,pz,dir() = GroupLine Res (lip) and (lipid)
          if sum (axis*dir)<0
            # Scalar product is negative, turn direction around
            dir()=-dir
          axis()=axis+dir
      if norm axis  
        # Orient the protein such that the axis through the membrane helices is perpenticular to the membrane
        ShowMessage "Orienting transmembrane regions perpendicular to the membrane..."
        oriold() = OriObj 1
        orinew()=oriold
        OriObj 1,0,0,0
        alpha,beta,gamma = OriVec (axis)
        RotateObj 1,Y=(-beta)
        RotateObj 1,Z=(90-gamma)
        # Now rotate around the Y-axis until the user sees most of the protein
        xmax=0
        for i=1 to 30
          x = GroupBox all
          if x>xmax
            xmax=x
            orinew() = OriObj 1
          RotateObj 1,Y=5
        # And rotate to the new orientation
        OriObj 1,(oriold)
        AutoPos X=0,Y=0,Z=(r*3),Steps=(delay),Wait=No
        AutoOriObj 1,(orinew),Steps=(delay)
      residuesmax=0
      if count explist or count liplist
        # Scan vertically for highest density of exposed hydrophobic residues or lipid anchors
        ShowMessage "Scanning protein for transmembrane region..."
        plane = ShowPolygon Points,Yellow,Vertices=4, (-r),0,(-r),  (r),0,(-r), (r),0,(r), (-r),0,(r)
        for i=-r to r step 2
          PosObj (plane),Y=(i),Z=(r*3)
          Wait 1
          if count liplist
            # Count lipid atoms in current transmembrane region
            memlipatoms = CountAtom (lipid) GlobalY>(i-memcoreheight/2) GlobalY<(i+memcoreheight/2)
            residues=(0.+memlipatoms)*count liplist/lipatoms
          else
            # Count exposed hydrophobic residues in current transmembrane region
            residues = CountRes (join explist) GlobalY>(i-memcoreheight/2) GlobalY<(i+memcoreheight/2)
          #print 'PosY=(i), ExpResidues=(residues)'
          if residues>residuesmax
            residuesmax=residues
            memposy=i
        DelObj (plane)
      # First get the cell that is required to enclose the entire protein, ignoring the membrane
      size1,size2,size3,cpos() = GroupBox all,Type=VdW,CoordSys=Global
      if residuesmax
        # Found a membrane embedding
        # Get the cell that is required to enclose the protein
        posmin1()=cpos-size*0.5-waterextension
        posmax1()=cpos+size*0.5+waterextension
        # Now get the cell that is required to enclose the protein membrane region
        size1,size2,size3,cpos() = GroupBox GlobalY>(memposy-memcoreheight/2) GlobalY<(memposy+memcoreheight/2),Type=VdW,CoordSys=Global
        # Along the Y-axis, the size is memheight+waterextension*2
        size2=memheight+waterextension*2-memextension*2
        posmin2()=cpos-size*0.5-memextension
        posmax2()=cpos+size*0.5+memextension
        # Calculate the new overall cell position, which encloses both cells above
        for i=1 to 3
          for f in 'min','max'
            poslist=pos(f)1(i),pos(f)2(i)
            pos(f)3(i)=(f) poslist
        cpos()=(posmin3+posmax3)*0.5
        size()=posmax3-posmin3
        if count explist
          # Color exposed hydrophobic residues orange
          ColorRes (join explist),150
      else
        # No embedding found, put membrane in the middle
        if size2<memheight
          size2=memheight
        size1=size1+memextension*2
        size2=size2+waterextension*2
        size3=size3+memextension*2
        memposy=cpos2
        Color Element
        if ConsoleMode
          RaiseError "This object does not look like a membrane protein. Please run YASARA in graphics mode and embed the protein manually"
        ShowMessage "This object does not look like a membrane protein. Click 'Continue', then place it in the membrane manually..."
        Wait ContinueButton
        confirmneeded=1
      if square
        # Make sure that X- and Z-axis have the same length and we get a square membrane
        if size1>size3
          size3=size1
        else
          size1=size3
      # Create and place cell
      Cell (size)
      PosObj SimCell,(cpos)
      if Structure and filetype!='sce'
        # Optimize the hydrogen-bonding network (more stable trajectories)
        ShowMessage "Optimizing hydrogen bonding network..."
        Wait 1
        OptHydAll
      # Create MembPreview, 'MEMBRANE' text that symbolizes the future membrane
      MakeTextObj Name=MembPreview,Width=(size3*10),Height=(size3)
      Font Arial,Height=45%,Color=Grey,Alpha=80,Depth=(memcoreheight),DepthCol=606060,DepthAlpha=50
      PosText 50%,50%,Justify=Center
      Print 'MEMB'
      PosText 50%,0%,Justify=Center
      Print 'RANE'
      ScaleMesh MembPreview,X=((size1/size3)*0.5)
      Style Ribbon,Stick
      PrintCon
      # Let the membrane preview fly in
      UserInput Off
      ShowMessage "Moving membrane placeholder to suggested position. Exposed hydrophobic residues are colored yellow..."
      PosOriObj MembPreview,Y=-40,Z=-18,Alpha=60
      AutoPosObj MembPreview,Y=(memposy),Z=(cpos3),Steps=(delay*2),Wait=No
      AutoOriObj MembPreview,Alpha=90,Steps=(delay*2)
      UserInput On
      # Let the user tune the embedding
      HUD Obj
      if confirmneeded
        ShowMessage "Move the objects and adjust the cell size to fine-tune the suggested membrane embedding. Click 'Continue' when done."
        Wait ContinueButton
      SaveSce (MacroTarget)_ori
    # Next step is to build the membrane
    # Get the required membrane size
    memsize1,cellsizey,memsize2 = Cell
    if count usermemname
      ShowMessage 'Building a membrane of size (0.0+memsize1) x (0.0+memsize2) A^2 based on user-supplied membrane (usermemname)...'
      membranename='(MacroTarget)_membrane_(0+memsize1)x(0+memsize2)_(usermemname).yob'
    else
      ShowMessage 'Building a membrane of size (0.0+memsize1) x (0.0+memsize2) A^2 with composition (PEApercent)% PEA, (PCHpercent)% PCH, (PSEpercent)%PSE...'
      if PEApercent1==PEApercent2 and PCHpercent1==PCHpercent2 and PSEpercent1==PSEpercent2
        membranename='(MacroTarget)_membrane_(0+memsize1)x(0+memsize2)_pea(PEApercent1)pch(PCHpercent1)pse(PSEpercent1).yob'
      else
        membranename='(MacroTarget)_membrane_(0+memsize1)x(0+memsize2)_pea(PEApercent1)-(PEApercent2)pch(PCHpercent1)-(PCHpercent2)pse(PSEpercent1)-(PSEpercent2).yob'
    membrane = FileSize (membranename)
    if membrane
      # Membrane has been built before
      memobj = LoadYOb (membranename)
    else
      # Membrane must be built now
      AutoPosObj 1,X=-150,Z=100,Steps=(delay),Wait=No
      AutoPosObj MembPreview,X=150,Z=100,Steps=(delay),Wait=No
      AutoPosObj SimCell,X=0,Y=150,Z=100,Steps=(delay)
      DelAll
      # Load the equilibrated membrane snapshot that matches best
      if count usermemname
        memobj = LoadYOb (YASARADir)/yob/membrane_(usermemname)
        memfragsize()=usermemsize
      elif PCHpercent1>=33 and PCHpercent2>=33
        memobj = LoadYOb (YASARADir)/yob/membrane_pea66pch33
        memfragsize=84.82,82.52
      elif PCHpercent1>=20 and PCHpercent2>=20
        memobj = LoadYOb (YASARADir)/yob/membrane_pea80pch20
        memfragsize=85.37,81.02
      else
        memobj = LoadYOb (YASARADir)/yob/membrane_pea
        memfragsize=81.68,79.99
      NameObj 1,Membrane
      PosObj 1,Y=-33,Z=-35
      AutoPosObj 1,Z=77,Steps=(delay)
      AutoOriObj 1,Alpha=90,Steps=(delay*2),Wait=No
      AutoPosObj 1,0,0,(norm memsize*1.5),Steps=(delay*2)
      TransformObj 1
      Cell Auto,Extension=10
      Cell (memfragsize)
      SwitchObj SimCell,off
      if !count usermemname
        # Prepare to mutate membrane residues to get the requested composition.
        # First we search for headgroups that don't have many neighbors and can
        # thus be changed without causing too many bumps. This search is done
        # with periodic boundaries and the original membrane fragment size.
        # Handle both membrane faces separately, each equilibrated membrane fragment consists
        # of 200 residues, residue numbers 1-100 are on side 1, 101-200 on side 2.
        for i=1 to 2
          # Set segment name so that each membrane side can be selected easily afterwards
          siderange='(-99+i*100)-(i*100)'
          SegRes (siderange),MEM(i)
          percentsum=PEApercent(i)+PCHpercent(i)+PSEpercent(i)
          if percentsum!=100
            RaiseError 'The percentages of PEA, PCH and PSE in column (i) sum up to (percentsum), but should sum up to 100'
          for name in 'PSE','PCH'
            lipids = CountRes (name) and (siderange)
            if lipids<(name)percent(i)
              # Determine those PEA residues that can best be mutated, i.e. have
              # the smallest number of neighbor atoms that could cause bumps.
              ShowMessage 'Locating exposed lipids on side (i) that can be mutated easily to (name)...'
              Wait 1
              Sim On
              pealist() = ListRes PEA and (siderange),Format=RESNUM
              for j=1 to count pealist
                if name=='PCH'
                  neighbors = CountAtom all with distance<5 from N Res (pealist(j))
                else
                  neighbors = CountAtom all with distance<5 from PDB 2HA Res (pealist(j)) 
                # Add neighbors*1000 to the residue number so that we can sort easily
                pealist(j)=pealist(j)+neighbors*1000
              Sim Off
              pealist()=sort pealist
              for j=1 to count pealist
                resnum=pealist(j)%1000
                if name=='PCH'
                  # Turn phosphatidyl-ethanolamine into phosphatidyl-choline
                  SwapAtom H with bond to N Res (resnum),C,Update=Yes
                  NameRes (resnum),PCH
                else
                  # Turn phosphatidyl-ethanolamine into phosphatidyl-serine
                  # SwapAtom creates C, with an unknown atom number since hydrogens are rebuilt 
                  SwapAtom PDB 2HA Res (resnum),C,Update=Yes
                  c = ListAtom C Res (resnum) 
                  NameAtom (c),C
                  DelAtom Element H with bond to (c)
                  AddHydAtom (c),2
                  SwapAtom (c+1) (c+2),O
                  Distance (c),all,Bound=Yes,Set=1.231
                  SwapBond (c),(c+1) (c+2),1.5
                  NameAtom (c+1),OT1
                  NameAtom (c+2),OT2
                  NameRes (resnum),PSE
                lipids=lipids+1
                ShowMessage 'Created (name) residue (lipids) of ((name)percent(i)) on side (i)...'
                Wait 1
                if lipids==(name)percent(i)
                  break
          AutoRotate Y=(180./delay)
          Wait (delay)
          AutoRotate Y=0
      OriObj all,0,0,0
      DelObj SimCell
      # Replicate the membrane fragments 
      axislist='X','Y'
      for i=1 to 2
        memrepsize(i)=memfragsize(i)
        while memrepsize(i)<memsize(i)
          ShowMessage 'Replicating membrane fragment to reach a membrane size of (0+memsize1)x(0+memsize2) A^2. (axislist(i)) (memrepsize(i))'
          obj = DuplicateObj (memobj)
          MoveObj (obj),(axislist(i))=(memrepsize(i))
          memrepsize(i)=memrepsize(i)+memfragsize(i)
          Wait 1
        JoinObj Membrane,(memobj)
      CenterObj Membrane
      AutoPosObj Membrane,0,0,Steps=(delay)
      # Keep aspect ratio
      if memrepsize1<memrepsize2
        memrepsize2=memrepsize1*(memsize2/memsize1)
      else
        memrepsize1=memrepsize2*(memsize1/memsize2)
      # Calculate how many lipids we need per side, to 1+ avoids that we round downwards
      lipidsneeded=1+100*(memsize1*memsize2/(memfragsize1*memfragsize2))
      # For each side, iteratively find the region that contains the required number of lipids
      # (We can't do it for both sides together, since we might end up with a different number of lipids per side)
      for i=1 to 2
        print 'Side (i)'
        side='Segment MEM(i)'
        memselmax()=0.5e0*memrepsize
        memselmin()=0.,0.
        lipidsmax = CountRes (side)
        lipidsmin=0
        do
          ShowMessage 'Locating membrane region on side (i) that contains the required number of (lipidsneeded) lipids...'
          Wait 1
          memsel()=(memselmin+memselmax)*0.5
          region='GlobalX>(-memsel1) GlobalX<(memsel1) GlobalY>(-memsel2) GlobalY<(memsel2)'
          lipids = CountRes (side) (region)
          #print 'Min=(memselmin), Max=(memselmax), Sel=(memsel), Lipids=(lipids)'
          if lipids>lipidsneeded
            memselmax()=memsel
            ColorRes (side) and not (region),Grey
          elif lipids<lipidsneeded
            memselmin()=memsel
        while lipids!=lipidsneeded and norm (memselmax-memselmin)>0.1
        # Delete surplus lipids
        DelRes (side) and not (region)
        HideMessage
      NumberRes all,1  
      # Enclose in cell, which is now larger than needed
      Cell Auto,Extension=0
      _,_,z = Cell
      Cell Z=(z+20)
      # Shrink the cell step-by-step while minimizing the lipids
      Longrange None
      Cutoff 5.24
      Temp 298K
      # Compressed membrane size is 5% smaller to close holes faster
      f=0.95
      commemsize()=memsize*f
      TimeStep 2,1
      Sim On
      # Since this is done in vacuo, prevent electrostatic interactions
      ChargeAtom Obj Membrane,0
      shrinkstep=0.6e0
      while 1
        x,y = Cell
        ShowMessage 'Compressing membrane to reach a size of (0+commemsize1)x(0+commemsize2) A^2, current size is (0+x)x(0+y) A^2...'
        if x<=commemsize1 and y<=commemsize2
          break
        Sim Off
        if x>commemsize1
          x=x-shrinkstep
        if y>commemsize2
          y=y-shrinkstep
        Cell X=(x),Y=(y)
        Sim On
        TempCtrl SteepDes
        for i=1 to 500
          Wait 1
          if SpeedMax<4000
            break
        TempCtrl Rescale
        for i=1 to 10
          StabilizeMembrane 'Z',35.,0
          Wait 10,femtoseconds
      # Allow membrane to relax
      ShowMessage 'Running short MD simulation of the compressed membrane to close gaps...'
      Sim On
      TempCtrl Anneal
      Wait 200,femtoseconds
      TempCtrl Rescale
      Temp 298K
      for i=1 to 100
        StabilizeMembrane 'Z',35.,0
        Wait 10,femtoseconds
      TempCtrl Anneal
      Wait 200,femtoseconds
      # Uncompress the membrane
      ScaleObj Membrane,X=(1./f),Y=(1./f)
      # The membrane ends have now been fused, the fusion region is at the periodic boundary.
      # Since the fusion region is not equilibrated and less ideally packed, we shift the
      # membrane by half the cell size to make sure that the worst part ends up in the cell
      # center - where it will be deleted anyway, since it overlaps with the protein.
      MoveAtom all,(memsize*0.5)
      Cell (memsize)
      # Wrap the lipids
      Sim On
      Sim Off
      # Finished
      Color Element
      SaveYOb 1,(membranename)
      AutoPosObj 1,Z=-45,Steps=(delay)
      # Load the orientation scene again and include the membrane
      LoadSce (MacroTarget)_ori
      memobj = LoadYOb (membranename)
    # Replace the 'MEMBRANE' text with the real membrane
    NameObj (memobj),Membrane
    TransferObj Membrane,MembPreview,Local=Keep
    SwapObj MembPreview,Membrane
    DelObj MembPreview
    # Align the cell with the global coordinate system in case it has been rotated by the user
    alpha,beta,gamma = OriObj SimCell
    Rotate Y=(-beta)
    Rotate Z=(-gamma)
    Rotate X=(-alpha)
    # The membrane placeholder and the actual membrane have top and bottom flipped
    RotateObj Membrane,X=180
    # Use the opportunity to remember the local Y-coordinate of the membrane center in the cell
    _,memposy = PosObj Membrane
    _,cellposy = PosObj SimCell
    memposy=memposy-cellposy+cellsizey*0.5
    # Freeze the current orientation
    TransformObj all
    # Make sure that protein and membrane share the same local coordinate system
    # The protein origin 0/0/0 is then at the center of the membrane
    TransferObj 1,Membrane,Local=Fix
    # Get a list of atoms embedded inside the membrane
    atomlist() = ListAtom Obj 1 LocalY>(-memcoreheight/2+2) LocalY<(memcoreheight/2-2)
    if count atomlist
      # Found embedded atoms
      # Calculate the bounding box of the transmembrane region
      transmemsize1,_,transmemsize2 = GroupBox (join atomlist), Type=Nuclear
      #ColorAtom Obj 1 LocalY>(-memcoreheight/2+2) LocalY<(memcoreheight/2-2),magenta
      # The following minimizations will be done quickly
      Cutoff 5.24
      Longrange None
      # Calculate initial XZ-scaling factor required to shrink the membrane region
      # 12 was too much (with old bumpdis 0.75 below), 10 was too little (with new bumpdis 1)
      scalexz=(max transmemsize-11)/max transmemsize
      # scalexz could be <=0 for small peptides now
      if scalexz<0.5
        scalexz=0.5e0
      # Prepare to embed the protein in the membrane. We cannot simply delete bumping
      # lipids: if a single overlap was enough to delete a lipid, we would end up with too few 
      # lipids and large holes between membrane and protein. Instead we first shrink the protein,
      # then delete bumping lipids, and finally grow the protein again during an energy minimization,
      # so that it can slowly push the surrounding lipids away until they smoothly cover the protein.
      #
      # Initialize simulation before the object is scaled, otherwise the incorrect
      # bond lengths prevent YASARA from caching the parameters to be on the safe side
      Sim Init
      # First shrink the protein strongly, in case it has a large pore
      ScaleObj 1,X=(scalexz*0.5),Z=(scalexz*0.5)
      ShowMessage "Deleting lipids that strongly bump into the shrunk protein..."
      Wait 1
      # Delete bumping lipids. 0.75 was too short, yielded lipid almost inside protein
      DelRes Obj Membrane with distance<1 from Obj 1
      ScaleObj 1,X=2,Z=2
      DelRes Obj Membrane with distance<1 from Obj 1
      # During the expansion, we don't want the protein atoms to move
      FixObj 1
      scalexzstart=scalexz
      ShowMessage "Expanding protein to fill the membrane pore..."
      Sim On
      ChargeAtom Obj Membrane,0
      # Iterate while the protein has not reached its original size yet
      while scalexz<1
        # Quick energy minimization
        Sim On
        TempCtrl SteepDes
        for i=1 to 500
          Wait 1
          if SpeedMax<4000
            break
        Temp 298K
        TempCtrl Rescale
        for i=1 to 20
          StabilizeMembrane 'Y',memposy,0
          Wait 10,femtoseconds
        f=1.02
        ScaleObj 1,X=(f),Z=(f)
        ShowMessage 'Expanding protein to fill the membrane pore, (0+(scalexz-scalexzstart)*100/(1.-scalexzstart))% completed...'
        scalexz=scalexz*f
      # Set original size again, since we did not reach 1.0 exactly
      ScaleObj 1,X=(1./scalexz),Z=(1./scalexz)
      ShowMessage "Short energy minimization to optimize membrane geometry..."
      # Minimize before saving, since the pulling caused incorrect bond lengths
      # Final scaling can have caused stress, start with steepest descent
      Sim On
      TempCtrl SteepDes
      for i=1 to 500
        Wait 1
        if SpeedMax<4000
          break
      TempCtrl Anneal
      Wait 400,femtoseconds
      Sim Off
    SaveSce (MacroTarget)
  # Set reliable simulation parameters again
  Cutoff (cutoff)
  Longrange Coulomb
  # Get membrane Y-position in the local coordinate system of the simulation cell
  # (needed in case we have not just created the membrane but loaded it instead).
  Sim On
  _,memposy = PosAtom Obj Membrane,Mean=yes
  Sim Off
  # Fill the cell with water, keep user-provided water fixed
  FreeRes !Water
  Experiment Neutralization
    WaterDensity 0.998
    pH (ph)
    Ions (ions)
    pKaFile (MacroTarget).pka
    #Vacuum LocalY>(memposy-memcoreheight*0.5) LocalY<(memposy+memcoreheight*0.5)
    Speed Fast
  Experiment On
  Wait ExpEnd
  # Delete leftover waters in the membrane region, except those that are far away
  # from the membrane (waters inside the channel of a pore protein) 
  ShowMessage 'Deleting water molecules inside the membrane...'
  Wait 1
  Sim On
  memregion='LocalY>(memposy-memcoreheight*0.5) LocalY<(memposy+memcoreheight*0.5)'
  DelRes HOH Atom O (memregion) !fixed with distance<9 from Obj Membrane
  # Tag those waters that are allowed in the membrane region (usually pore waters)
  # No pulling force will be applied to these waters during the equilibration period
  Sim On
  PropRes all,0
  PropRes HOH Atom O (memregion),100
  Sim Off
  StickRes Water
  # Save scene with water
  SaveSce (MacroTarget)_water
  HideMessage
Wait 1
# Don't keep selected atoms, LoadXTC/LoadMDCRD would load only selected ones
Unselect

# Choose timestep and activate constraints
if speed=='fast'
  # Fast simulation speed
  # Constrain bonds to hydrogens
  FixBond all,Element H
  # Constrain certain bond angles involving hydrogens
  FixHydAngle all
  # Choose a multiple timestep of 2*2.5 = 5 fs
  # For structures with severe errors, 2*2 = 4 fs is safer (tslist=2,2)
  tslist=2,2.5
else
  # Slow or normal simulation speed
  # Remove any constraints
  FreeBond all,all
  FreeAngle all,all,all
  if speed=='slow'
    # Choose a multiple timestep of 2*1.00 = 2.0 fs
    tslist=2,1.0
  else
    # Choose a multiple timestep of 2*1.25 = 2.5 fs
    tslist=2,1.25
    # With this timestep, atoms may get too fast in very rare circumstances (only
    # in a specific protein, only once every few nanoseconds). The command below
    # slows down atoms moving faster than 13000 m/s. Such a 'random collision' every
    # few nanoseconds has no more impact than the random number seed. You can comment
    # it out for most proteins, or use the smaller timestep with speed 'slow' above:
    Brake 13000
# During equilibration update the pairlist every 10 steps
SimSteps Screen=10,Pairlist=10
# Calculate total timestep, we want a float, so tslist2 is on the left side
ts=tslist2*tslist1
# Snapshots are saved every 'savesteps'
savesteps=saveinterval/ts
# Set final simulation parameters
TimeStep (tslist)
Temp (temperature)
Cutoff (cutoff)
Longrange Coulomb
# Make sure all atoms are free to move
FreeAll
# Alread a snapshot/trajectory present?
i=00000
if format=='sim'
  trajectfilename='(MacroTarget)(i).sim'
else  
  trajectfilename='(MacroTarget).(format)'
  restartfilename='(MacroTarget).sim'  
  # Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
  old = FileSize (MacroTarget)(i).xtc
  if old
    RenameFile (MacroTarget)(i).xtc,(trajectfilename)
running = FileSize (trajectfilename)
if not running
  # Perform energy minimization
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  # And now start the real simulation
  Sim On
else
  # Simulation has been running before
  ShowMessage "Simulation has been running before, loading last snapshot..."
  # Switch console off to load the snapshots quickly
  Console Off
  if format=='sim'
    # Find and load the last SIM snapshot
    do
      i=i+1
      found = FileSize (MacroTarget)(i).sim
    while found
    i=i-1
    LoadSim (MacroTarget)(i)
    # Adjust savesteps to save snapshots in the same interval as previously
    if i>0
      t = Time
      savesteps=0+t/(ts*i)
  else
    # Do we have a restart file with atom velocities?
    found = FileSize (restartfilename)
    if found
      # Yes. First determine the savesteps if possible by loading the 2nd XTC/MDCrd snapshot
      last,t = Load(format) (trajectfilename),1
      if !last
        last,t = Load(format) (trajectfilename),2
        savesteps=0+t/ts
      # Then load the restart file
      LoadSim (restartfilename)
    else
      # No restart file found, load the last snapshot in the XTC/MDCrd trajectory
      do
        i=i+1
        last,t = Load(format) (trajectfilename),(i)
        ShowMessage 'Searching (format) trajectory for last snapshot, showing snapshot (i) at (0+t) fs'
        Sim Pause
        Wait 1
      while !last
      savesteps=0+t/(ts*(i-1))
      Sim Continue
HideMessage

# Set temperature and pressure control
TempCtrl Rescale
PressureCtrl Manometer2D,Pressure=(pressure)

# Now the simulation is running, here you can make changes to the force field

# Uncomment to add distance constraints
# AddSpring O Res Lys 80,H Res Glu 84,Len=1.9
AddSpring 1438,6135,Len=2,SFC=100
# Uncomment to modify charges, e.g. let Trp 12 in Mol A lose an electron:
# ChargeRes Trp 12 Mol A,+1

# And finally, make sure that future snapshots are saved
Save(format) (trajectfilename),(savesteps)
if format!='sim'
  # We save an XTC/MDCrd trajectory plus a single Sim restart file
  SaveSim (restartfilename),(savesteps),Number=no

# Get current membrane Y position and disable potentially conflicting drift correction
_,memposy = PosAtom Obj Membrane,Mean=yes
CorrectDrift off
    
# At the beginning of the simulation, the membrane is not ideally packed yet and needs
# about 250ps equilibration time. During this period it is still 'vulnerable' to
# water molecules that are squeezed in. We thereforce keep these waters out.
t = Time
if t<equiperiod*1000
  ShowMessage 'Running (equiperiod)ps of equilibration simulation, preventing water from entering the membrane...'
  while t<equiperiod*1000
    for i=1 to 10
      StabilizeMembrane 'Y',memposy,'Water'
      Wait (10*ts),femtoseconds
    t = Time
  HideMessage

# Membrane has been equilibrated, now keep the membrane protein from diffusing around and crossing periodic boundaries
CorrectDrift on

# After equilibration update the pairlist every 10 (CPU) or 25 (GPU) steps
_,_,gpu = Processors
if gpu
  SimSteps Screen=25,Pairlist=25
else    
  SimSteps Screen=10,Pairlist=10

if duration=='forever'
  Console On
  if ConsoleMode
    # In the console, we need to wait forever to avoid a prompt for user input
    Wait forever
else
  Console Off
  measurements=0
  # Wait for given number of picoseconds
  do
    # Tabulate properties you want to monitor during the simulation,
    # e.g. the speeds and velocity vectors of atoms 4, 5 and 7:
    #Tabulate SpeedAtom 4 5 7
    # Or apply a pulling force, in this example downwards
    #AccelRes GLI 1027,Y=-10000
    # Wait for one screen update
    Wait 1
    measurements=measurements+1
    t = Time
  while t<1000.*duration+1
  # Did we create a table with measurements?
  vallist() = Tab Default
  if count vallist
    # Yes, save the table
    SaveTab default,(MacroTarget)_duringsim,Format=Text,Columns=(count vallist/measurements),Header='Insert your own header here'
  Sim Off
# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit

# STABILIZE A MEMBRANE BY PULLING THE LIPIDS 
# ==========================================
# 'axis' is the membrane axis (normally 'Y' or 'Z'), 'pos' is the
# position of the membrane center along 'axis'. 'water' is either 0
# (no water) or selects the water object to be kept out of the membrane.
def StabilizeMembrane axis,pos,water
  global memheight,memcoreheight
  
  #t=Time
  #ShowMessage 'Stabilizing membrane at (0+t), pos (pos)'
  # Set single atom force f and residue acceleration a
  s = SimSteps
  if s>10
    # Don't permit more than 10 steps per screen update during stabilization, otherwise forces get too large
    s=10
    SimSteps (s)
  f=1200.*s
  a=4000.*s
  # Accelerate complete lipid residues that try to leave the membrane
  AccelRes Obj Membrane Local(axis)>(pos+memheight*0.5),(axis)=(-a)
  AccelRes Obj Membrane Local(axis)<(pos-memheight*0.5),(axis)=(a)
  # Pull up the lipid heads that get sucked into the membrane
  # Replace 'P' with the name of the phosphorus atom in your lipid if needed
  ForceAtom P Segment MEM1 Local(axis)>(pos-14),(axis)=(-f)
  ForceAtom P Segment MEM2 Local(axis)<(pos+14),(axis)=(f)
  # Pull back the lipid tails that get too close to the membrane surface
  # Replace 'C34' and 'C18' with the names of the terminal carbons in the long
  # and short lipid tails
  ForceAtom C34 Segment MEM1 Local(axis)<(pos+2),(axis)=(f)
  ForceAtom C18 Segment MEM1 Local(axis)<(pos-2),(axis)=(f)
  ForceAtom C34 Segment MEM2 Local(axis)>(pos-2),(axis)=-(f)
  ForceAtom C18 Segment MEM2 Local(axis)>(pos+2),(axis)=-(f)
  """
  # Coloring for debugging
  Color Element
  ColorRes Obj Membrane Local(axis)>(pos+memheight*0.5),Magenta
  ColorRes Obj Membrane Local(axis)<(pos-memheight*0.5),Magenta
  ColorAtom P Segment MEM1 Local(axis)>(pos-14),Magenta
  ColorAtom P Segment MEM2 Local(axis)<(pos+14),Magenta
  ColorAtom C34 Segment MEM1 Local(axis)<(pos+2),Yellow
  ColorAtom C18 Segment MEM1 Local(axis)<(pos-2),Yellow
  ColorAtom C34 Segment MEM2 Local(axis)>(pos-2),Yellow
  ColorAtom C18 Segment MEM2 Local(axis)>(pos+2),Yellow
  """
  if water
    # Keep the water out of the membrane (atoms with property 100 are allowed membrane pore waters)
    AccelRes Property!=100 Local(axis)>(pos-memcoreheight*0.5) Local(axis)<(pos) Obj (water),(axis)=(-a)
    AccelRes Property!=100 Local(axis)<(pos+memcoreheight*0.5) Local(axis)>(pos) Obj (water),(axis)=(a)
    # Ions need more force
    AccelRes !Water Local(axis)>(pos-memcoreheight*0.5) Local(axis)<(pos) Obj (water),(axis)=(-a*10)
    AccelRes !Water Local(axis)<(pos+memcoreheight*0.5) Local(axis)>(pos) Obj (water),(axis)=(a*10)
    """
    # Coloring for debugging
    ColorRes Property!=100 Local(axis)>(pos-memcoreheight*0.5) Local(axis)<(pos) Obj (water),Magenta
    ColorRes Property!=100 Local(axis)<(pos+memcoreheight*0.5) Local(axis)>(pos) Obj (water),Magenta
    """
