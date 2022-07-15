# LFA_RFpulseDesign
SAR constrained slice selective arbitrary flip angle parallel transmit RF pulse design for MRI

Start from runSpokes.m under runFromHere folder.

This MATLAB code can design arbitrary flip angle parallel transmission slice selective refocusing and excitation pulses to mitigate transmit field inhomogeneity while constraining SAR and/or RF power. 
A full bloch simulation of the RF pulse is run at every step of optimization to characterize the effect of the RF pulse.
fmincon optimization function is utilized. Jacobian and (approximated Hessian) for the objective function of the optimization
is calculated analytically and provided to fmincon to speed up computation. 
Also, GPU and CPU parallelization is implemented as options to use to speed up the computation time. 
The refocusing and excitation pulse pair is set up to satisfy the CPMG condition.
Windowed sinc subpulses are used for slice selection.


runFromHere: main operating directory

SimData: includes simulated transmit fields (8 channels), off resonance map and
other information such as the region of interest, position in x and y, etc..

mintverse_v10_filiz: includes code to apply the minimum time VERSE algorithm on the RF subpulses

Spokes_LFA_v12_filiz_v2: includes most of the functions used in pulse design


runFromHere

runSpokes.m: main file, change options for the simulation and some of the pulse parameters
determining e.g. pulse duration, pulse shape, etc. from here. You can choose to design 
CP mode pulses or 1-spoke or multispoke pulses using least squares or magnitude least squares

spokes_def_90.txt & spokes_def_180.txt: define some system and pulse parameters here
e.g. the directory for the transmit field and off resonance map data, slice thickness, maximum
gradient strength, slew rate, maximum RF voltage, time discretization step (e.g. 10us), etc.

sarmats: includes text files containing the virtual observation points and the average SAR matrices
obtained from the simulation data the transmit fields were generated from.


Spokes_LFA_v12_filiz_v2

Spokes_LFA.m: main pulse design file which calls designSTA and designLTA

designSTA.m: main file to design small tip angle RF pulses to initialize the large flip angle pulse design

designLTA.m: main file to design arbitrary flip angle RF pulses. Calls 

f0_LFA_mexCudaInterf.m: details of the objective function and the Jacobian matrix calculation. Includes calls to
mex files to utilize CPU or GPU parallelization. The parallelization is implemented across pixels of the ROI.

hessian_LFA_noB0_mexCudaInterf.m: details of the approximated Hessian matrix calculation. The approximation
is that off resonance field is assumed to be 0 simplifying the Bloch simulation rotation matrices from 100s to
a few.
