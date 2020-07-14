#include <stdlib.h>
#include "header.h"

double ran3(int idum)
{ static int inext, inextp;
  static int iff=0;
  const int MBIG=1000000000,MSEED=161803398,MZ=0;
  const double FAC=(1.0/MBIG);
//  static valarray<int> ma(56);
  static int ma[56];
  int i, ii, k, mj, mk;
  if(idum<0 || iff==0)
  { iff=1;
    mj=labs(MSEED-labs(idum));
    mj%=MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<55;i++) 
    { ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk<int(MZ)) mk+=MBIG;
      mj = ma[ii];
    }
    for(k=1;k<=4;k++)
    { for(i=1;i<=55;i++)
      { ma[i] = ma[i] - ma[1 + ((i+30)%55)];
        if (ma[i]<int(MZ)) ma[i]+=MBIG;
      }
    }
    inext=0;
    inextp=31;
    //idum=1;		//this must be in here, if you don't have ran3(int idum=1) in  header file
    			//and then the reference must be passed 
    			//and function must always be called with ran3(idum)
  }
  if (++inext == 56)   inext  =1;
  if (++inextp ==56)  inextp =1;
  mj = ma[inext] - ma[inextp];
  if (mj<int(MZ)) mj+=MBIG;
  ma[inext]=mj;
  return  mj *FAC;
}


//----------------------------GENERATE RANDOM VELOCITY ACCORDING TO BOLTZMANN DISTRIBUTION------------------------


double RandomVelocity(double beta,double mass)
{
        double v1,v2,rsq,fac;
        double temperature;
        double velocityMB;

        do
        {
                v1=2.0*ran3()-1.0;
                v2=2.0*ran3()-1.0;
                rsq= v1*v1 + v2*v2;
        }while((rsq>=1) || (rsq==0));
        fac = sqrt(-2.0*log(rsq)/rsq);

        //v1*fac and v2*fac are gaussian distributed random numbers (cf. numerical recipies etc.)

        temperature=1.0/beta;

        velocityMB = sqrt(temperature/mass)*v1*fac;

        return velocityMB;
}

//-------------------------------Generate Gaussian random number-------------------------------------------------------

double gauss()
{
        double v1,v2,rsq,fac;
	double grand;
	static int getnew = 0;
	static double savegrand = 0.0;


	if(getnew == 0){
        	do
        	{
                	v1=2.0*ran3()-1.0;
                	v2=2.0*ran3()-1.0;
                	rsq= v1*v1 + v2*v2;
        	}while((rsq>=1) || (rsq==0));
        	fac = sqrt(-2.0*log(rsq)/rsq);

        	//v1*fac and v2*fac are gaussian distributed random numbers (cf. numerical recipies etc.)
		grand = v1*fac;
		savegrand = v2*fac;
		getnew = 1;
	}
	else{
		grand = savegrand;
		getnew = 0;
	}



        return grand;
}


