// checked
#include "math.h"
#include "iostream"
#include "interface.h"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;

class diffusionfluxinterface
{
	public:
	float diffusionfluxvector[5] ;

	diffusionfluxinterface(vector<float>& vectorleftminus, vector<float>& vectorleft, vector<float>& vectorright,
	vector<float>& vectorrightplus, vector<float>& areavectorleft, vector<float>& areavectorright, vector<float>& areavectorrightplus,
	float volumeleftmins, float volumeleft, float volumeright, float volumerightplus, float deltat)
	{
		float thetai[5];
		float thetaiplus[5];
		float gvector[5] ;
		float gvectorplus[5] ;
		float signale ;
		float betainterface[5] ;
		float phiinterface[5] ;
		float Zvectorinterface[5];
		float pshiinterface[5] ;
	
	interface left(vectorleftminus,vectorleft,areavectorleft,volumeleftmins,volumeleft,deltat) ;
	interface right(vectorleft, vectorright, areavectorright,volumeleft,volumeright, deltat) ;
	interface rightplus(vectorright,vectorrightplus,areavectorrightplus,volumeright,volumerightplus,deltat);

// here I am somehow trying to defie the gi and gi+1
for (int i = 0; i < 5; ++i)
{
	if (right.gvectorinterface[i] >= 0)
	{
		signale = 1.0;
	}
	else
	{
		signale = -1.0 ;
	}
	// this is pura jugaad so change karna hai baad me "gvector[i]" ka expression
	float tempB = (left.gvectorinterface[i]*signale)  ;
	float tempA = fabs(right.gvectorinterface[i]) ;
	tempA = min(tempA,tempB) ;
	float temp_zero=0.0;
	tempB = max(temp_zero, tempA) ;
	gvector[i] = signale *  tempB ;

	// gvector[i] = signale * max(0.0, min( fabs(right.gvectorinterface[i]),
	// 		left.gvectorinterface[i]*signale ) ) ;
}

 for (int i = 0; i < 5; ++i)
{
	if (rightplus.gvectorinterface[i] >= 0)
	{
		signale = 1.0;
	}
	else
	{
		signale = -1.0 ;
	}
// this is pura jugaad so change karna hai baad me "gvector[i]" ka expression
	float tempA = fabs(rightplus.gvectorinterface[i]) ;
	float tempB = (right.gvectorinterface[i]*signale)  ;
	tempA = min(tempA,tempB) ;
	float temp_zero=0.0;
	tempB = max(temp_zero, tempA) ;
	gvectorplus[i] = signale * tempB ;
	// gvector[i] = signale * max(0.0, min( fabs(right.gvectorinterface[i]),
	// 		left.gvectorinterface[i]*signale ) ) ;
	
 }
// theta i term defying 
for (int i = 0; i < 5; ++i)
	{
		// here I have doubt about not equal to symbol because it cant be exect 0.00000 so most of the 
		//time we endup choosing theta i = 0.0
		if((fabs(right.alphavectorinterface[i]) + fabs(left.alphavectorinterface[i]))!=0.0)
		{
			thetai[i] = fabs(right.alphavectorinterface[i] - left.alphavectorinterface[i])/
				(fabs(right.alphavectorinterface[i]) + fabs(left.alphavectorinterface[i])) ;
		}
		else
		{
			thetai[i] = 0.0;
		}
	}

// theta i plus 
for (int i = 0; i < 5; ++i)
{
	// here I have doubt about not equal to symbol because it cant be exect 0.00000 so most of the
	//time we endup choosing theta i = 0.0
	if((fabs(rightplus.alphavectorinterface[i])+fabs(right.alphavectorinterface[i]))!=0.0)
	{
		thetaiplus[i]=fabs(rightplus.alphavectorinterface[i]-right.alphavectorinterface[i])/
			(fabs(rightplus.alphavectorinterface[i]) + fabs(right.alphavectorinterface[i])) ;
	}
	else
	{
		thetaiplus[i] = 0.0 ;
	}
}

// beta defying
// before that we nee to define the omega values which are constant 
float omega[5] = {0.25, 1.0, 1.0, 1.0, 0.25} ;
for (int i = 0; i < 5; ++i)
{
	betainterface[i] = 1.0 + omega[i]*(max(thetai[i],thetaiplus[i]));
}
/////////////////////////////////
// phi defying
for (int i = 0; i < 5; ++i)
{
	if (right.alphavectorinterface[i] != 0)
	{
		phiinterface[i] = (gvectorplus[i] - gvector[i]) / right.alphavectorinterface[i] ;
	}
	else
	{
		phiinterface[i] = 0.0 ;
	}
}

// now re-defying the Zinterface
for (int i = 0; i < 5; ++i)
{
	Zvectorinterface[i] = right.Zvectorinterface[i] + betainterface[i]*phiinterface[i] ;
}

// pshiinterface re-defying (for invisid wall flow)
for (int i = 0; i < 5; ++i)
{
	pshiinterface[i] = pow(Zvectorinterface[i],2) + 0.25 ;
}

// re-defying the pshiinterface
// // again we need to specify the deltaf = 0.2
// 	float deltaf = 0.2 ;
// 	for (int i = 0; i < 5; ++i)
// 	{
// 		if (Zvectorinterface[i] >= deltaf)
// 		{
// 			pshiinterface[i] = fabs(Zvectorinterface[i]) ;
// 		}
// 		else
// 		{
// 			pshiinterface[i] = 0.5*(pow(Zvectorinterface[i],2) + pow(deltaf,2))/ deltaf ;
// 		}
// 	}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// finally calculating the interface diffusion term = diffusionfluxvetor
// before that we can calculate tempvector as
float tempvector[5];
for (int i = 0; i < 5; ++i)
{
	tempvector[i] = ( betainterface[i]*(gvector[i]+gvectorplus[i]) -
	 pshiinterface[i]*right.alphavectorinterface[i] )/deltat ;
}
// now the diffusionfluxvector
for (int i = 0; i < 5; ++i)
{
	diffusionfluxvector[i] = 0.0 ;
	for (int l = 0; l < 5; ++l)
	{
		diffusionfluxvector[i] += right.eigenvectormatrix[i][l] * tempvector[l] ;
	}
}
};
	// ~diffusionfluxinterface();
	
};