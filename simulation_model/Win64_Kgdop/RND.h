#define NTAB 32
#define NWUP 8
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define IM1 2147483563
#define IM2 2147483399
#undef AM
#define AM (1./IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#undef NDIV
#define NDIV (1+IMM1/NTAB)

#include <time.h>
static long dummy;

void Seed(long dum) {dummy=dum;}

float unirand(void) {
 int j;
 long k;
 static long dummy2=123456789;
 static long iy=0;
 static long iv[NTAB];
 float temp;
  if(dummy<=0 || !iy) {
 if(dummy<0) dummy=-dummy;
 else if(dummy==0) dummy=1;
 dummy2=dummy;
 for(j=NTAB+NWUP-1;j>=0;j--) {
 k=dummy/IQ1;
 if((dummy=IA1*(dummy-k*IQ1)-IR1*k)<0) dummy+=IM1;
 if(j<NTAB) iv[j]=dummy;
 }
 iy=iv[0];
 }

 k=dummy/IQ1;
 if((dummy=IA1*(dummy-k*IQ1)-IR1*k)<0) dummy+=IM1;
 k=dummy2/IQ2;
 if((dummy2=IA2*(dummy2-k*IQ2)-IR2*k)<0) dummy2+=IM2;
 iy=iv[j=iy/NDIV]-dummy2;iv[j]=dummy;
 if(iy<1) iy+=IMM1;
 if((temp=AM*iy)>RNMX) return(RNMX);
 else return(temp);
}