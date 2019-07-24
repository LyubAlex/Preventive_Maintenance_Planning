tic
clc
format long; 

P0 = [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]; % the vector of initial state
P = []; % matrix of transition probabilities
L01 = 14.15E-6;
L02 = 3.99E-6;
L12 = 1.71E-6;
Tp = 0.25;
Tr = 0.5;
Ta = 1.2;
Ts = 0.8;
a1 = 0.01;
a2 = 0.005;
b1 = 0.02;
b2 = 0.01;

Tk = 6E5; % simulation time, hours
Ntr =500; % number of runs
TobL = 4000; % max value of maintenance periodicity
dT = 100; % step of the maintenance periodicity
NTk = TobL/dT; % number of periodicity values for calculation

masFSi = zeros(length(P0));
masTSi = zeros(length(P0));

for n = 1:Ntr
n
for intT = dT:dT:TobL 
NS1 = zeros(length(P0));
NSi = zeros(length(P0));
NTSi = zeros(length(P0));

Ti = 0;
T = 0;
I = 0; % index of the current state within a run
%============================================================%
%              Determination of initial state 
%============================================================%
A = 0;
PS = rand; 
    while PS >= A
        
          if (I >= length(P0))
             I = 0;     
          else
             I = I + 1;
             A = A + P0(I); 
          end            
    end
    NS1(I)=NS1(I)+1;   
%============================================================%

while T < Tk
      F = [((-1/(L01+L02))*log(rand)) ((-1/(L12))*log(rand)) (Tp+Ts+Ta) Tp  (Tp+Tr) intT Tp Tp];
    
   if (F(I) <= intT)   
    
        if (T + F(I)>=Tk)
            F(I) = Tk - T;
            T = Tk;
            NTSi(I)=NTSi(I)+F(I);
        else
            NTSi(I)=NTSi(I)+F(I);
            T = T + F(I);  
        end
        
    else 
        F(I) = intT;
            if (T + F(I)>=Tk)
            F(I) = Tk - T;
            T = Tk;
            NTSi(I)=NTSi(I)+F(I);
        else
            NTSi(I)=NTSi(I)+F(I);
            T = T + F(I);  
        end
    end
    NSi(I)=NSi(I)+1;
    P = getmasP(a1,a2,b1,b2,L01,L02,L12,intT); % calculation of the current transition probability matrix
    J = OPRZ(I,P); % get index of the next state       
I = J;
Ti = Ti+1;
end

 masKti(n,intT/dT) = (NTSi(1)+NTSi(2)+NTSi(8))/(NTSi(1)+NTSi(2)+NTSi(3)+NTSi(4)+NTSi(5)+NTSi(6)+NTSi(7)+NTSi(8));
 masKg(n,intT/dT) = (NTSi(1)+NTSi(2)+NTSi(8))/(NTSi(1)+NTSi(2)+NTSi(3)+NTSi(6)+NTSi(8));
end
    for i = 1:intT/dT
        addMasKti = 0;
        addMasKg = 0;
        
           for in = 1:n
               addMasKti=addMasKti+ masKti(in,i);
               addMasKg=addMasKg+ masKg(in,i);         
           end    
        MxKti(i) = addMasKti/n;
        MxKg(i) = addMasKg/n;
 
    end

end 

figure(3)
plot([dT:dT:TobL],MxKti([1:NTk]), [dT:dT:TobL],MxKg([1:NTk]));
axis ([0 TobL 0.999 1]);

