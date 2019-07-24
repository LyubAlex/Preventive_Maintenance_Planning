clc
clear all
format long; 
iFig = 0;
%============================================================%
% Input parameters 
%============================================================%
P0 = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]; % P0 is the vector of initial state 
%(for example, for P0=(1,0,0,0,0,0,0,0) is assumed that each run (replication) 
%of the model starts from the operable state S0 with probability equal to 1.0)
Tk = 1e6; % Tk is the simulation time in hours
dT = 100; % dT is the step for changing the maintenance periodicity, hours
eps = 0.5e-4; % eps is the required estimation accuracy of MxKg and MxKti (confidence level = 0.95)

Kgdop = 0.997; % Kgdop is the allowable value of inherent availability
%============================================================%
%Pop is the function of service staff errors probability
TkO = 8760;
Pm = 0.05;
% Pop = [0.0001 0.2*Pm 0.4*Pm 0.6*Pm 0.8*Pm Pm Pm Pm];
Pop = [0.0 Pm Pm Pm Pm Pm Pm];
T   = [0 0.5*TkO 0.6*TkO 0.8*TkO TkO 1.5*TkO 2*TkO];
TA = 0:2*TkO;
Pop = interp1(T,Pop,TA,'linear'); 

iFig = iFig + 1;
figure(iFig);
plot(TA,Pop);
axis ([0 TkO 0.0 Pm*1.05]);
grid on;
%============================================================%
% [L01 L02 L12] is the vector of the rates of misalignments, sudden failures, 
% and failure rate of misaligned system, respectively, hours^(-1)
L01 = 14.15E-6;
L02 = 3.99E-6;
L12 = 1.71E-6;

% [Tp Tr Ta Ts a1 a2 b1 b2] is the vector of maintenance time parameters and 
% parameters characterizing built-in and external diagnostic devices:    
Tp = 0.25; %Tp is the testing time, hours;
Tr = 0.5; %Tr is the time for system tuning and configuration, hours;
Ta = 1.2; %Ts is the failure search time, hours;
Ts = 0.8; %Ta is the emergency repair time, hours; 

a1 = 0.01; % a1, b1 are the type I (a) and type II (b) errors of built-in diagnostic modules; 
a2 = 0.005;
b1 = 0.02; % a2, b2 are type I (a) and type II (b) errors of the external diagnostic devices used for maintenance. 
b2 = 0.05;

Tobmax = 1.5*8760; % Tobmax is the maximal maintenance periodicity which is considered for reliability indicators calculation   


if ['Release R' version('-release')] == 'Release R2016b'
    tic
    [MxKg, MxKti, epsKg, epsKti, masT, masP, TobL, SFlag] = PMSMP_MSVC(Tk, Kgdop, P0,[L01 L02 L12],[Tp Tr Ta Ts a1 a2 b1 b2],Pop,eps,dT,Tobmax);
    toc
elseif ['Release R' version('-release')] == 'Release R2018b'
    tic
    [MxKg, MxKti, epsKg, epsKti, masT, masP, TobL, SFlag] = PMSMP_MinGW(Tk, Kgdop, P0,[L01 L02 L12],[Tp Tr Ta Ts a1 a2 b1 b2],Pop,eps,dT,Tobmax);
    toc
end


iFig = iFig + 1;
figure(iFig);
plot([dT:dT:TobL],MxKti,[dT:dT:TobL], MxKg);
xlabel('Maintenance periodicity, hours')
ylabel('Reliability indicators')
legend('Inherent Availability','Achieved Availability ')
axis ([dT TobL/2 0.998 1]);
grid on;

