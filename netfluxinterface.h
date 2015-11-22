#include "math.h"
#include "iostream"
#include "eulerflux.h"
// #include "viscusflux.h"
#include "diffusionfluxinterface.h"
#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;

class netfluxinterface
{
	public:
	float netflux[5] ;
	
	netfluxinterface(vector<float>& vectorleftminus, vector<float>& vectorleft, vector<float>& vectorright,
vector<float>& vectorrightplus, vector<float>& areavectorleft, vector<float>& areavectorright, vector<float>& areavectorrightplus,
float volumeleftmins, float volumeleft, float volumeright, float volumerightplus, float deltat){

eulerflux left(vectorleft);
eulerflux right(vectorright);

float volumeinterfaceleft = (volumeleftmins +volumeleft)/2 ;
float volumeinterfaceright = (volumeleft + volumeright)/2 ;
float volumeinterfacerightplus = (volumeright + volumerightplus)/2 ;

diffusionfluxinterface diffusion(vectorleftminus,vectorleft,vectorright,vectorrightplus,
	areavectorleft,areavectorright,areavectorrightplus,volumeleftmins,volumeleft,
	volumeright,volumerightplus,deltat); 

// Averaged interface Euler flux 
		float hx = areavectorright[0] / volumeinterfaceright ;
		float hy = areavectorright[1] / volumeinterfaceright ;	
		float hz = areavectorright[2] / volumeinterfaceright ;

		float xnetfluxinterface[5] ; 
		float ynetfluxinterface[5] ; 
		float znetfluxinterface[5] ; 
		for (int i = 0; i < 5; ++i)
		{
			xnetfluxinterface[i] = 0.5*(left.xeulerflux[i]+right.xeulerflux[i]) ; 
			ynetfluxinterface[i] = 0.5*(left.yeulerflux[i]+right.yeulerflux[i]) ; 
			znetfluxinterface[i] = 0.5*(left.zeulerflux[i]+right.zeulerflux[i]) ;

			netflux[i] = (xnetfluxinterface[i]*hx + ynetfluxinterface[i]*hy + znetfluxinterface[i]*hz)*volumeinterfaceright + 
			0.5*diffusion.diffusionfluxvector[i] ;
		}

	};

	// ~netfluxinterface();
	
};