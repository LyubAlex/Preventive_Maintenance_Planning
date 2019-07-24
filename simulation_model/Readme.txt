The code in PMSMP.cpp allows to generate a MEX function to be used in Matlab environment for reliability indices evaluation such as inherent availability and achieved availability of radio communication devices depending on preventive maintenance periodicity. 
The code is based on a Monte Carlo simulation model of Semi-Markov process corresponded to the operational process of the mentioned facilities. The model takes into account the following exploitation factors: sudden, gradual, latent and fictitious failures, a human factor of service staff and time parameters of maintenance. Consequently, the model considers eight possible states: operable state, misalignment state, non-operable state, preventive maintenance of the operative system, maintenance of the system with misalignment, latent failure, maintenance of system being in latent failure and fictitious failure. This model can be utilized for the analysis of reliability and optimization of maintenance periodicity of radio-electronic and communication facilities whose operational process is assumed to be described by the above-mentioned states. 
The function Preventive Maintenance Semi-Markov Process (PMSMP) is proposed to be used as a multivariant analysis module whose work is stopped automatically when calculations achieve the allowable level of inherent availability set in model input parameters. This automatic stop feature is added in order to use the multivariant analysis module as a part of a CAD system of preventive maintenance schedule radio communication equipment. 
[MxKg, MxKti, epsKg, epsKti, masT, masP, TobL, SFlag] = PMSMP(Tk,Kgdop,P0,[L01 L02 L12],[Tp Tr Ta Ts a1 a2 b1 b2],Pop,eps,dT,Tobmax) 
Input parameters: 
Tk is the simulation time in hours (Tk=10e6 is recommended); 
Kgdop is the allowable level of inherent availability (for example, 0.997); 
P0 is the vector of initial state (for example, for P0=(1,0,0,0,0,0,0,0) is assumed that each run (replication) of the model starts from the operable state S0 with probability equal to 1.0); 
[L01 L02 L12] is the vector of the rates of misalignments, sudden failures, and the failure rate of misaligned system, respectively, hours^(-1); 
[Tp Tr Ta Ts a1 a2 b1 b2] is the vector of maintenance time parameters and parameters characterizing built-in and external diagnostic devices: 
Tp is the testing time, hours; 
Tr is the time for system tuning and configuration, hours; 
Ts is the fault search time, hours; 
Ta is the emergency repair time, hours; 
a1, b1 are the type I (a) and type II (b) errors of built-in diagnostic modules; 
a2, b2 are type I (a) and type II (b) errors of the external diagnostic devices used for maintenance. Pop is the function of service staff errors probability; 
eps is the required estimation accuracy of MxKg and MxKti; 
dT is the step for changing the maintenance periodicity, hours (dT=300 hours is recommended). 
Tobmax is the maximal maintenance periodicity which is considered for reliability indicators calculation. 

Output parameters: 
MxKg is the array of mean values of the inherent availability for each value of periodicity in [dT:dT:TobL];
MxKti is the array of mean values of the achieved availability for each value of periodicity in [dT:dT:TobL];
epsKg is the array  of MxKg estimation errors;
epsKt is the array  of MxKti estimation errors;
masT is the matrix of size MxN of means of the sojourn time in each state of the process. M is the number of values of periodicity in [dT:dT:TobL], N=8 is the number of states of the considered stochastic process
masP is the matrix of size MxN of mean probabilities of staying in each state of the process;
TobL is the maintenance periodicity that corresponds to the allowable level of inherent availability Kgdop
SFlag is the flag that indicates the achievement of  Kgdop within the analyzed interval of periodicity from dT to Tobmax

In this folder Win64_Kgdop (works on 64-bit Windows), you can find the following files:
PMSMP.cpp is the source code of the simulation model
Example.m is the Matlab script for the demonstration
PMSMP_MinGW.mexw64 is the compiled mex-function in Matlab R2016b using MinGW64 v6.3.0. The function also was tested in R2018b with the same installed compiler and it works well.  
PMSMP_MSVC.mexw64 is the compiled mex-function in earlier versions of Matlab using MS Visual C++ compiler, the version is unknown. PMSMP_MSVC works well in R2016b, doesn’t run in R2018b.  
I provide here two already compiled versions because PMSMP_MSVC has better performance, i.e. needs less computational time, in the Example.m file, PMSMP_MSVC performs three times faster the experiment than PMSMP_MinGW  
RND.h is the code of the pseudo-random generator with uniform distribution based on L’Ecuyer algorithm.
If this simulation model seems to be useful for your research and you use it, please, make a reference in your publications as A. Lyubchenko,  “MEX function for multivariant analysis of reliability indices depending on maintenance periodicity of radio communication equipment,” source code, 2012. DOI: 10.5281/zenodo.321435
