// checked 
#include "math.h"
#include "iostream"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;

class eulerflux
{
	public:
	float xeulerflux[5] ;
	float yeulerflux[5] ;
	float zeulerflux[5] ;

	eulerflux(vector<float>& vector)
	{
		float pressure = (gamma -1)*( vector[4] - 0.5*(pow(vector[1],2)+pow(vector[2],2)+
		pow(vector[3],2))/vector[0] ) ;  

		//  Euler flux
		xeulerflux[0] = vector[1] ;
		xeulerflux[1] = (pow(vector[1],2)/vector[0]) + pressure ;     
		xeulerflux[2] = vector[1]*vector[2]/vector[0] ; 
		xeulerflux[3] = vector[1]*vector[3] /vector[0] ;
		xeulerflux[4] = ((vector[4] + pressure)*vector[1]) / vector[0] ;

		yeulerflux[0] = vector[2] ;
		yeulerflux[1] = vector[1]*vector[2]/vector[0] ; 
		yeulerflux[2] = (pow(vector[2],2)/vector[0]) + pressure ;     
		yeulerflux[3] = vector[2]*vector[3] /vector[0] ;
		yeulerflux[4] = ((vector[4] + pressure)*vector[2]) / vector[0] ;

		zeulerflux[0] = vector[3] ;
		zeulerflux[1] = vector[1]*vector[3]/vector[0] ; 
		zeulerflux[2] = vector[2]*vector[3] /vector[0] ;
		zeulerflux[3] = (pow(vector[3],2)/vector[0]) + pressure ;     
		zeulerflux[4] = ((vector[4] + pressure)*vector[3]) / vector[0] ;

	 };
	 // ~eulerflux();	
};