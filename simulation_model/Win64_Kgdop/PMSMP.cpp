#include <math.h>
#include <string.h>
#include <windows.h>
#include <matrix.h>
#include <vector>
#include "mex.h"
#include "RND.h"
#include <vector>

void outputIntT (double pIntT) {
double pr[] = {pIntT};
mxArray *array_ptr; int num_out, num_in; 
mxArray *output_array[1], *input_array[1];

array_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
memcpy(mxGetPr(array_ptr), pr, 1*sizeof(double)); 
        
num_out = 0; 
num_in = 1; 
input_array[0] = array_ptr; 
mexCallMATLAB(num_out, output_array, num_in, input_array, "disp");
}

void mexFunction(int nOut, mxArray*	pOut[], int	nIn, const mxArray*	pIn[])   
{ 

//============================================================//
//Model input parameters
//============================================================//
	int Tk   = mxGetScalar(pIn[0]);
	
	double Kgdop = mxGetScalar(pIn[1]);
	
	double *InP0;
	InP0  = mxGetPr(pIn[2]);

	double *InL;
	InL        = mxGetPr(pIn[3]);
	double L01 = InL[0];
	double L02 = InL[1];
	double L12 = InL[2];

	double *InTimes;
	InTimes    = mxGetPr(pIn[4]);
	double Tp  = InTimes[0];
	double Tr  = InTimes[1];
	double Ta  = InTimes[2];
	double Ts  = InTimes[3];
	double a1  = InTimes[4];
	double a2  = InTimes[5];
	double b1  = InTimes[6];
	double b2  = InTimes[7];

	double *InPop;
	InPop = mxGetPr(pIn[5]);
	
	double eps = mxGetScalar(pIn[6]);
	int dT = mxGetScalar(pIn[7]);

	int Tobmax = mxGetScalar(pIn[8]);
//============================================================//

//============================================================//
//Additional variables of the model
//============================================================//
	std::vector<double> masMxKg; 
	std::vector<double> masMxKti;

	std::vector<double> masEpsKg;
	std::vector<double> masEpsKti;

	std::vector<double> masKgS;
	std::vector<double> masKtiS;
	std::vector<double> masNtr;

	double SKg = 0;
	double SKti = 0;

	double DxKg = 0;
	double DxKti = 0;

	double MxKg = 0;
	double MxKti = 0;

	double NSi[8];
	double NTSi[8];
	double NTSi1[8];
		
	int I,J;
	
	double addMxKg = 0;
	double addMxKti = 0;

	double addDxKg = 0;
	double addDxKti = 0;

	double Kg = 0;
	double Kti = 0;

	double Ntr = 0;
	double kolEps = 0;
	double epsKg, epsKti;

	double P2,PS;

	double (*masT1)[8] = new double[18000/dT][8];
	double (*masP)[8]  = new double[18000/dT][8];
	
	double addMasT[8];
	double addMasP[8];
	double kolTr;

    int intT;
	int TobL = 0;
	int SFlag;

//============================================================//

//============================================================//
//Initilizing of random number generator 
//============================================================//
	long seed = time(NULL);
	Seed(seed);
//============================================================//
	

for (intT = dT; intT < Tk; intT=intT+dT)
{
	kolEps = 0;
	epsKg =0;
	epsKti = 0;
	
	SKg = 0;
	SKti = 0;

	DxKg = 0;
	DxKti = 0;

	MxKg = 0;
	MxKti = 0;

	addMxKg = 0;
	addMxKti = 0;

	addDxKg = 0;
	addDxKti = 0;

	Kg = 0;
	Kti = 0;

	Ntr = 0;
	
	double F01 = 1-exp(-L01*intT);
	double F02 = 1-exp(-L02*intT);
	double F12 = 1-exp(-L12*intT);
	double Fto = InPop[intT];

	double P41 = (1-a2)*(1-Fto);
	double P12 = (1-F02)*F01;
	double P13 = (1-b1)*F02;
	double P23 = (1-b1)*F12;
	double P43 = (1-b2)*Fto;
	double P14 = (1-a1)*(1-F02)*(1-F01);
	double P25 = (1-a1)*(1-F12);
	double P18 = a1*(1-F02)*(1-F01);
	double P28 = a1*(1-F12);
	double P48 = a2*(1-Fto);

	double P[8][8] = { 
			0,  P12,  P13,  P14,  0,   b1*F02,  0,  P18,
			0,   0,   P23,   0,  P25,  b1*F12,  0,  P28,
			1,   0,    0,    0,   0,     0,     0,   0,
		       P41,  0,   P43,   0,   0,   b2*Fto,  0,  P48,
		       P41,  0,   P43,   0,   0,   b2*Fto,  0,  P48,
			0,   0,    0,    0,   0,     0,     1,   0,
		        0,   0,   1-b2,  0,   0,     b2,    0,   0, 
			1,   0,    0,    0,   0,     0,     0,   0 	};

	for (int i = 0; i<8; i++)
		{
			addMasT[i] = 0;
			addMasP[i] = 0;
		}
	do
	{		
			Ntr++;  
			for ( int i = 0; i<8; i++)
				{
					NSi[i] = 0;
					NTSi[i] = 0;
					NTSi1[i] = 0;
					kolTr = 0;
				}
			double T = 0;
			I = 0;
			double A = 0;
			PS = unirand();
			A = A + InP0[0]; 
			while (PS > A)
			{
				I = I + 1;
				A = A + InP0[I];
			}
		
		double tau = 0;
		double tau1 = 0;
		while (T < Tk)
		{
			P2 = unirand();
			J = 0;
			double B = P[I][J];

			while (P2 > B)
			{
						J = J + 1;
						B = B + P[I][J];
			}
			NSi[I] = ++NSi[I];
			double F[] = {-1.0/(L01+L02)*log(P2), (-1.0/(L12)*log(P2)), Tp+Ts+Ta, Tp, Tp+Tr, intT, Tp, Tp};

		if (F[I] < intT)
		{	
					double Fi1 = F[I];					
					if (T + Fi1 >= Tk)
						{
							Fi1 = Tk - T;
							T = Tk;
							NTSi[I] = NTSi[I] + Fi1;
						}
						else
						{
							NTSi[I] = NTSi[I] + Fi1;
							T = T + Fi1;
						}
		}
		else
		{	
					double Fi2 = intT;			
						
						if (T + Fi2 >= Tk)
							{
								Fi2 = Tk - T;
								T = Tk;
								NTSi[I] = NTSi[I] + Fi2;
							}
						else
							{
								NTSi[I] = NTSi[I] + Fi2;
								T = T + Fi2;
							}
		}						
			I = J;
			kolTr++;
		}

		Kti = (NTSi[0]+NTSi[1]+NTSi[7])/(NTSi[0]+NTSi[1]+NTSi[2]+NTSi[3]+NTSi[4]+NTSi[5]+NTSi[6]+NTSi[7]);
		Kg  = (NTSi[0]+NTSi[1]+NTSi[7])/(NTSi[0]+NTSi[1]+NTSi[2]+NTSi[5]+NTSi[7]);

        
		addMxKti = addMxKti + Kti;
		addMxKg  = addMxKg + Kg;
    
        
		addDxKti = addDxKti + pow(Kti,2);
		addDxKg  = addDxKg + pow(Kg,2);

		MxKti= addMxKti/Ntr;
		MxKg = addMxKg/Ntr;

       	for (int i = 0; i < 8; i++)
		{	
					addMasT[i] = addMasT[i] + NTSi[i];
					masT1[intT/dT-1][i] = addMasT[i]/Ntr;

					addMasP[i] = addMasP[i] + NSi[i]/kolTr;
					masP[intT/dT-1][i] = addMasP[i]/Ntr;
		}          
        
                
	if (Ntr > 120)
		{
			DxKti= (addDxKti/Ntr - pow(MxKti,2))*Ntr/(Ntr-1);
			DxKg = (addDxKg/Ntr  - pow(MxKg,2))*Ntr/(Ntr-1);

			SKti = sqrt(DxKti);
			SKg  = sqrt(DxKg);

			epsKti = (SKti/sqrt(Ntr))*1.96;
			epsKg  = (SKg/sqrt(Ntr))*1.96;

			if ( epsKg < eps && epsKti < eps)
				{
					kolEps++;
				}
			else 
				{
					continue;
				}
		}
	else
		{
			continue;
		}

	}
	while ( kolEps != 10 );
	if (intT <= Tobmax)
		{
     
                  
			if (MxKg <= (Kgdop-(Kgdop*0.0001)))
				{	
					SFlag = 0;
					TobL = intT;

					masMxKg.push_back(MxKg);
					masMxKti.push_back(MxKti);
					masEpsKg.push_back(epsKg);
					masEpsKti.push_back(epsKti);			
					break;
				}
			else
				{
					TobL = intT;
					masMxKg.push_back(MxKg);
					masMxKti.push_back(MxKti);
					masEpsKg.push_back(epsKg);
					masEpsKti.push_back(epsKti);
					outputIntT(intT); //output of the current value of the calculated point in Matlab workspace 
					continue;
				}
		}
		else
		{
        outputIntT(Tobmax);
        outputIntT(intT);
			SFlag = 1;
			break;
		}
}
//============================================================//
//Model output parameters
//============================================================//
	double *outMasMxKg;
	pOut[0]=mxCreateDoubleMatrix(1,masMxKg.size(),mxREAL);
	outMasMxKg  = mxGetPr( pOut[0] );

	double *outMasMxKti;
	pOut[1]=mxCreateDoubleMatrix(1,masMxKti.size(),mxREAL);
	outMasMxKti = mxGetPr( pOut[1] );

	double *outMasEpsKg;
	pOut[2]=mxCreateDoubleMatrix(1,masEpsKg.size(),mxREAL);
	outMasEpsKg = mxGetPr( pOut[2] );

	double *outMasEpsKti;
	pOut[3]=mxCreateDoubleMatrix(1,masEpsKti.size(),mxREAL);
	outMasEpsKti = mxGetPr( pOut[3] );
	
	double *outMasT;
	pOut[4]=mxCreateDoubleMatrix(intT/dT,8, mxREAL);
	outMasT = mxGetPr( pOut[4] );

	double *outMasP;
	pOut[5]=mxCreateDoubleMatrix(intT/dT,8, mxREAL);
	outMasP = mxGetPr( pOut[5] );

	double *outMasTobL;
	pOut[6]=mxCreateDoubleMatrix(1,1,mxREAL);
	outMasTobL = mxGetPr( pOut[6] );

	double *outFlag;
	pOut[7]=mxCreateDoubleMatrix(1,1,mxREAL);
	outFlag = mxGetPr( pOut[7] );

	outFlag[0] = SFlag; 
	if (SFlag == 0)
	{
		outMasTobL[0] = TobL;
		for (int j = 0; j < masMxKg.size(); j++)
			{
				outMasMxKg[j]   = masMxKg[j];
				outMasMxKti[j]  = masMxKti[j];

				outMasEpsKg[j]  = masEpsKg[j];
				outMasEpsKti[j] = masEpsKti[j];
			}
	int k=0;
	for (int i = 0; i < 8; i++)
		{
		for (int j =0 ; j<intT/dT; j++)
			{
				outMasT[k] = masT1[j][i];
				outMasP[k] = masP[j][i];
				k++;
			}
		}
	delete []masT1;
	delete []masP; 
    return;
	}
	else
	{
	delete []masT1;
	delete []masP; 
    return;
	}
}
//============================================================//

