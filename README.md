# analysis codes for CLAS12 SIDIS analyses
## Tmothy B. Hayward, last updated January 25, 2023


repository for various analysis codes (EventBuilder, fitters, etc.) used for analyzing CLAS12 data at Jefferson Lab. Primarily SIDIS focused. 

A directory with myClara installation (with my coatjava extensions already in the correct directory), the clasQA and the below example scripts can be found on ifarm at  
> /work/clas12/thayward/analysis_tools/ 

Timothy Hayward's  
[Thesis](https://www.jlab.org/Hall-B/general/thesis/THayward_thesis.pdf )  
[Letter](https://arxiv.org/abs/2101.04842 )

--------

included files:  
I. **analysis_fitter.java**  
&nbsp;&nbsp;&nbsp;This is the class for the kinematic fitter I use to build events (takes the detector responses to assign particle ID to tracks and adds them to the event). The idea is to take the CLAS12 EventBuilder as a basis and enhance the PID on top of that. Loops through all particles in REC::Particle bank and sees if they pass the enhanced particle PID cuts (e.g. tightened sampling fraction, fiducial cuts, chi2pid cuts for hadron identification etc.) Start reading around line 700,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"public PhysicsEvent getPhysicsEvent(DataEvent event) {"

IIa. **Inclusive.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'X events (inclusive DIS). 
  
IIb. **Hadron.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.)  for ep -> e'hX events (single hadron SIDIS). 

IIc. **Dihadrons.java**
  Extension (technically written first) of the single hadron case to two hadrons, ep -> e' h1 h2 X, includes additional variables for each hadron, e.g. z1, z2, PT1, PT2, etc.
  
IId. **Trihadron.java**
  Extension of the single hadron case to three hadrons, ep -> e' h1 h2 h3 X, includes additional variables for each hadron, e.g. z1, z2, z3, PT1, PT2, PT3 etc. as well as all permutations of hadron combinations, i.e. Mh12, Mh13, Mh23
  
IIe. **Hadron_RICH_studies.java**  
&nbsp;&nbsp;&nbsp;This is the class used to calculate relevant kinematic variables (Q2, W, Mx, xF, PT, phi_trento, etc.) for ep -> e'hX events (single hadron SIDIS) and gives full access to all the variables in the RICH::Particle bank. 
  
&nbsp;&nbsp;&nbsp;Both classes start with a "channel_test" function that allows for cuts on Q2, W, xF, y, Mx, etc. By default Q2 > 1, W > 2, y < 0.8* are enabled and the others are commented out. If you are sure of the cuts you desire, more can be turned on to significantly increase processing speed. *I'd recommend double checking these.
  
  
III. **hayward_coatjava_extensions.jar**
&nbsp;&nbsp;&nbsp;The distribution version of all my classes if you want to just try running the code as is first. You'll need this one in your myClara/plugins/clas12/lib/services/ directory to use this with CLAS12 Groovy code.

IV. **processing_single_hadrons.groovy**
&nbsp;&nbsp;&nbsp;I prefer to process the clas12 hipo4 files and create text outputs that I can then import into Mathematica, root, etc. It's not the most efficient but it is very nicely universal. Includes a section (commented out by default to increase computation speed) for the printing of revelant RICH variables for hadron PID studies. If no track is present in the RICH (which only exists in one sector; second sector installed for RGC forward) then it prints PID = 0. This script accepts 4 input arguments:  
  1. hipo file directory  
  a directory such as, /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/, that contains hipo4 data files.  
  2. pid for the hadron i.e. 211 for pi+  
  3. output text file name  
  4. number of files to process in the directory

IVb. **processing_dihadrons.groovy**
&nbsp;&nbsp;&nbsp;dihadron version of the above. This script accepts 5 input arguments:  
  1. hipo file directory  
  a directory such as, /cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/, that contains hipo4 data files.  
  2. pid for p1 (hadron 1), i.e. 211 for pi+  
  3. pid for p2 (hadron 2), i.e. -211 for pi-
  4. output text file name  
  5. number of files to process in the directory

V. **processing_mc_single_hadrons.groovy**
&nbsp;&nbsp;&nbsp;same as **processing_single_hadrons.groovy** except for the use with MC. Places generated and reconstructed values next to each other, i.e. rec_Q2 - gen_Q2. For studies of various systematics like bin migration.

Vb. **processing_mc_dihadrons.groovy**
&nbsp;&nbsp;&nbsp;same as **processing_dihadrons.groovy** except for the use with MC.

I have also included processing scripts for the inclusive and trihadron case which should be self explanatory; the descriptions of each column can be found at the end of the groovy scripts. 

--------

The processing scripts (at least for data) rely on the [CLAS12 Quality Assurance](https://github.com/JeffersonLab/clasqaDB) code (developed by Christopher Dilks). You can do a git pull on it and source the following file before running my scripts. At the moment we rely on the QA cut "OkForAsymmetry" which is defined in the QA readme.  
> source clasqaDB/env.csh 
