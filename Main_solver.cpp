// assuming that only these three things are given 
// (1)Grid points
// (2)Boundary condition
// (3)Initial condition are given
#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "netfluxinterface.h"
#include "BC.h"
// #include "grid_straight_duct.h"
#include "grid_bump.h"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;
// BC function implements the boundary condition
void BC(vector<vector<vector<vector<float> > > > & variablesvector, vector<vector<vector<vector<float> > > > & y_face_area, int Nx, int Ny, int Nz) ;

// grid function genrates the area vector and volume for the all cells in the domain
void grid(vector<vector<vector<vector<float> > > > & x_face_area , vector<vector<vector<vector<float> > > > & y_face_area,
	vector<vector<vector<vector<float> > > > & z_face_area, vector<vector<vector<float> > >& cell_volume,
	int Nx, int Ny, int Nz,float lenght,float delta) ;


int main ()
{
	time_t start, end ;
	time(&start); // noteing the starting time

	float deltat = 0.00015; // this is for CFL = 0.2
	float TIME = 100000*deltat; 
	int totaltimesteps = floor(TIME/deltat) ;

	float lenght = 25 ; 
	float delta = 1.0 ; // this basically defines the grid size

	int N = floor(lenght/delta) ;
	// extra 4 is added for ghost cell
	int Nx = (3*N + 4);
	int Ny = N+4 ;
	int Nz = 4+4 ; // Because this is 2D-simulation so no need to take large number of grids in z direction 

	// Creating a 4D vector object
	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	// this store previous values of variables (density , three momentum, energy)
	matrix4D variablesvector(Nx,Dim3(Ny,Dim2(Nz,Dim1(5)))); 
	// this store new values of variables(density , three momentum, energy) 
	matrix4D variablevectornew(Nx,Dim3(Ny,Dim2(Nz,Dim1(5)))); 
	
	// this store grid information in the domain
	matrix4D x_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D y_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D z_face_area(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	Dim3 cell_volume(Nx,Dim2(Ny,Dim1(Nz)));

	// Now the grid function has defined the area vector and volume od cells 
	grid(x_face_area,y_face_area,z_face_area,cell_volume,Nx,Ny,Nz,lenght,delta) ;

	// Initial condition(these are just rendom values )
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				variablesvector[i][j][k][0] = 1.16;
				variablesvector[i][j][k][1] = 0 ;
				variablesvector[i][j][k][2] = 0 ;
				variablesvector[i][j][k][3] = 0;
				variablesvector[i][j][k][4] = 271767;

				variablevectornew[i][j][k][0] = 1.16;
				variablevectornew[i][j][k][1] = 0 ;
				variablevectornew[i][j][k][2] = 0 ;
				variablevectornew[i][j][k][3] = 0 ;
				variablevectornew[i][j][k][4] = 271767;
			}
		}
	}

	// time progression
	// this file is opend to store the mass residual at each time step
	ofstream kullu_mass ;
	// kullu_mass.open("mass_residual_SD.csv");
	kullu_mass.open("mass_residual_B.csv");
	kullu_mass <<  "t(secs)" << "," << "mass_residual" <<  endl ;

	for (int t = 0; t < totaltimesteps; ++t)
	{
		// cout << "timestep  = " << t << endl ; 
		// Befor every time step we need to have proper value in the ghost cells 
		// So applying the Boundary condition
		// BC takes  care of Inlet, Exit, y-wall and Z-wall boundary condition
		BC(variablesvector,y_face_area,Nx,Ny,Nz) ; 

		// next timestep variablesvectors calcuation
		for (int i = 1; i < Nx-2; ++i)
		{
			for (int j = 1; j < Ny-2; ++j)
			{
				for (int k = 1; k < Nz-2; ++k)
				{
					//x right interface volume
					float xcellinterfacevolumeright = 0.5*(cell_volume[i][j][k] + cell_volume[i+1][j][k]);

					//y right interface volume
					float ycellinterfacevolumeright = 0.5*(cell_volume[i][j][k] + cell_volume[i][j+1][k]);

					//z right interface volume
					float zcellinterfacevolumeright = 0.5*( cell_volume[i][j][k] + cell_volume[i][j][k+1]);

					// net flux using the class netfluxinterface
					netfluxinterface xrightface(variablesvector[i-1][j][k],variablesvector[i][j][k],variablesvector[i+1][j][k],
							variablesvector[i+2][j][k] ,x_face_area[i-1][j][k] ,x_face_area[i][j][k] ,x_face_area[i+1][j][k],
						cell_volume[i-1][j][k], cell_volume[i][j][k], cell_volume[i+1][j][k], cell_volume[i+2][j][k] ,deltat) ;

					netfluxinterface yrightface(variablesvector[i][j-1][k],variablesvector[i][j][k],variablesvector[i][j+1][k],
							variablesvector[i][j+2][k] ,y_face_area[i][j-1][k] ,y_face_area[i][j][k] ,y_face_area[i][j+1][k],
						cell_volume[i][j-1][k], cell_volume[i][j][k], cell_volume[i][j+1][k], cell_volume[i][j+2][k] ,deltat) ;

					netfluxinterface zrightface(variablesvector[i][j][k-1],variablesvector[i][j][k],variablesvector[i][j][k+1],
							variablesvector[i][j][k+2] ,z_face_area[i][j][k-1] ,z_face_area[i][j][k] ,z_face_area[i][j][k+1],
						cell_volume[i][j][k-1], cell_volume[i][j][k], cell_volume[i][j][k+1], cell_volume[i][j][k+2] ,deltat) ;

					// updating the variablevectornew using flux at the right interfaces
					for (int l = 0; l < 5; ++l)
					{
						variablevectornew[i][j][k][l] -=(deltat/xcellinterfacevolumeright)*(xrightface.netflux[l]);
						variablevectornew[i+1][j][k][l] +=(deltat/xcellinterfacevolumeright)*(xrightface.netflux[l]);

						variablevectornew[i][j][k][l] -=(deltat/ycellinterfacevolumeright)*(yrightface.netflux[l]);
						variablevectornew[i][j+1][k][l] +=(deltat/ycellinterfacevolumeright)*(yrightface.netflux[l]);

						variablevectornew[i][j][k][l] -=(deltat/zcellinterfacevolumeright)*(zrightface.netflux[l]);
						variablevectornew[i][j][k+1][l] +=(deltat/zcellinterfacevolumeright)*(zrightface.netflux[l]);
					}

				}
			}
		}

		// before going to the new timestep update variablesvector by variablevectornew
		for (int i = 2; i < Nx-2; ++i)
		{
			for (int j = 2; j < Ny-2; ++j)
			{
				for (int k = 2; k < Nz-2; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						variablesvector[i][j][k][l] = variablevectornew[i][j][k][l] ;
					}				
				}
			}
		}

		// cout << "Mach at (x,4,8) "  << endl ;
		// // printing the mach at three location so to wheter solution is coming correct or has is diverged
		// for (int i = 1; i < 10; ++i)
		// {
		// 	float pressure = (gamma-1)*( variablesvector[i][8][4][4] - (0.5*( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )/variablesvector[i][8][4][0]) ) ;
		// 	float temperature = pressure/( gasconstant*variablesvector[i][8][4][0] ) ;
		// 	cout <<  ( (sqrt( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )) / ( variablesvector[i][8][4][0] * sqrt(gamma*gasconstant*temperature) ) ) << endl ;
		// }
		// cout<< endl ;

		
		// for (int i = (Nx-10)/2; i < (Nx+10)/2; ++i)
		// {
		// 	float pressure = (gamma-1)*( variablesvector[i][8][4][4] - (0.5*( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )/variablesvector[i][8][4][0]) ) ;
		// 	float temperature = pressure/( gasconstant*variablesvector[i][8][4][0] ) ;
		// 	cout <<  ( (sqrt( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )) / ( variablesvector[i][8][4][0] * sqrt(gamma*gasconstant*temperature) ) ) << endl ;
		// }
		// cout << endl ;


		// for (int i = Nx-12; i < Nx-2; ++i)
		// {
		// 	float pressure = (gamma-1)*( variablesvector[i][8][4][4] - (0.5*( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )/variablesvector[i][8][4][0]) ) ;
		// 	float temperature = pressure/( gasconstant*variablesvector[i][8][4][0] ) ;
		// 	cout <<  ( (sqrt( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
		// 					pow(variablesvector[i][8][4][3],2) )) / ( variablesvector[i][8][4][0] * sqrt(gamma*gasconstant*temperature) ) ) << endl ;
		// }
		
			
		// mass residual calculation after each timestep and writing that into the mass_residual file
		
		float m_differance = 0.0 ;
		float m_in = 0.0 ;
		for (int j = 2; j < Ny-2; ++j)
			{
				for (int k = 2; k < Nz-2; ++k)
				{
					m_differance += (variablesvector[1][j][k][1]-variablesvector[Nx-2][j][k][1]-variablesvector[Nx-2][j][k][2]-
						variablesvector[Nx-2][j][k][3])  ;  
					m_in += variablesvector[1][j][k][1] ; 
				}
			}

			kullu_mass << t*deltat << "," << m_differance/m_in  << endl ;
			cout << m_differance/m_in  << endl ;

			if (t%250 == 0)
			{
				// // storing the premetive variable at few time steps
				// ofstream kullu_premitive ;
				// // kullu_premitive.open("non_dim_parameters_SD.csv");
				// kullu_premitive.open("premitive_variable.csv");
				// kullu_premitive << "timestep" << t << endl ;
				// kullu_premitive << "density" << "," << "density*u" << ","<< "density*v" << "," << "density*w" << "," << "energy"  << endl ;
				// for (int i = 0; i < Nx; ++i)
				// {

				// 	kullu_premitive << variablesvector[i][8][4][0] << "," << variablesvector[i][8][4][1]<<","<<variablesvector[i][8][4][2] << "," <<
				// 	 variablesvector[i][8][4][3] << "," << variablesvector[i][8][4][4] << endl ;
				// }
				
				// storing the velocity in one plane
				ofstream kullu_2D ;
				kullu_2D.open("2D_parameters_B.csv");
				kullu_2D << "density" << "," << "density*u" << ","<< "density*v" << "," << "density*w" << "," << "energy"  << endl ;
				for (int i = 2; i < Nx-2; ++i)
				{
					for (int j = 2; j < Ny-2; ++j)
					{
						kullu_2D << variablesvector[i][j][4][0] << "," << variablesvector[i][j][4][1]<<","<<variablesvector[i][j][4][2] << "," <<
						 variablesvector[i][j][4][3] << "," << variablesvector[i][j][4][4] << endl ;
					}
				}


				// opening the file in write mode(this file contains non-dimentionalized parameters)
				ofstream kullu_non_dim ;
				// kullu_non_dim.open("non_dim_parameters_SD.csv");
				kullu_non_dim.open("non_dim_parameters_B.csv");
				kullu_non_dim << "X(m)" << "," << "density" << ","<< "pressure" << "," << "temperature" << "," << "Mach" << ","  << "net_velocity" << endl ;
				float free_stream_temperature = 288.2 ;
				float free_stream_pressure = 101325 ;
				float free_stream_density = free_stream_pressure /(gasconstant*free_stream_temperature) ;
				for (int i = 0; i < Nx; ++i)
				{
					float density = variablesvector[i][8][4][0]  ;
					float pressure = (gamma-1)*( variablesvector[i][8][4][4] - ( 0.5*( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
									pow(variablesvector[i][8][4][3],2) )/variablesvector[i][8][4][0] ) );
					float temperature = pressure/( gasconstant*density) ;

					float Mach = ( (sqrt( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
									pow(variablesvector[i][8][4][3],2) )) / ( variablesvector[i][8][4][0] * sqrt(gamma*gasconstant*temperature) ) ) ;

					float net_velocity = Mach*sqrt(gamma * gasconstant * temperature) ;
					
					// writeing the data inside the file

					kullu_non_dim << (delta*i)/lenght << "," << density/free_stream_density<<","<<pressure/free_stream_pressure << "," <<
					 temperature/free_stream_temperature << "," << Mach << ","  << net_velocity/(sqrt(gamma*gasconstant*free_stream_temperature)) << endl ;
				}
			}

	} 
	// time progression ends here


	// // opening the file in write mode(this file contains dimentional parameters)
	// ofstream kullu_data ;

	// // this will contain all the dimentional parameters as maintioned below
	// // kullu_data.open("dim_parameters_SD.csv"); 
	// kullu_data.open("dim_parameters_B.csv"); 

	// kullu_data << "X(m)" << "," << "density" << ","<< "pressure" << "," << "temperature" << "," << "Mach" << ","  << "net_velocity" << endl ;
	// for (int i = 0; i < Nx; ++i)
	// {
	// 	float pressure = (gamma-1)*( variablesvector[i][8][4][4] - (0.5*( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
	// 					pow(variablesvector[i][8][4][3],2) )/variablesvector[i][8][4][0]) ) ;

	// 	float temperature = pressure/( gasconstant*variablesvector[i][8][4][0] ) ;
	// 	float Mach = ( (sqrt( pow(variablesvector[i][8][4][1],2) + pow(variablesvector[i][8][4][2],2) + 
	// 					pow(variablesvector[i][8][4][3],2) )) / ( variablesvector[i][8][4][0] * sqrt(gamma*gasconstant*temperature) ) ) ;

	// 	float net_velocity = Mach*sqrt(gamma*gasconstant*temperature) ;

	// 	float density = pressure /(gasconstant*temperature) ;
		
		
	// 	// writeing the data inside the file
	// 	kullu_data << delta*i << "," << density<<","<<pressure << "," << temperature << "," << Mach << ","  << net_velocity << endl ;
	// }

	


	time(&end) ;
	double diff = difftime (end,start);
	cout << "Time taken by the solver in secs = " << diff << endl ;
	return 0;
}

