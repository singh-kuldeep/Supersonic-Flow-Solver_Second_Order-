// Checked 
// this file genrates the grid for the straight duct
#include "iostream"
#include <vector>
#include "math.h"

// Function defines the area vector and cell volumes 
void grid(vector<vector<vector<vector<float> > > > & x_face_area ,vector<vector<vector<vector<float> > > > & y_face_area,
	vector<vector<vector<vector<float> > > > & z_face_area, vector<vector<vector<float> > > & cell_volume,
	int Nx, int Ny, int Nz,float lenght, float delta)  
{
	// Grids for straight duct(just to check wheter the basic code working)
	float deltatx = delta ;
	float deltaty = delta ;
	float deltatz = delta ;

	for (int i = 0; i  < Nx+1; ++i)
	{
		for (int  j = 0;  j < Ny+1; ++j)
		{
			for (int  k = 0;  k < Nz+1; ++k)
			{
				// area vector for x direction cells
				x_face_area[i][j][k][0] = deltaty*deltatz ;   
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				// area vector for y direction cells
				y_face_area[i][j][k][0] = 0 ;   
				y_face_area[i][j][k][1] = deltatz*deltatx ;
				y_face_area[i][j][k][2] = 0 ;

				// area vector for z direction cells
				z_face_area[i][j][k][0] = 0 ;   
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = deltaty*deltatx ;
			}
		}	
	}

	// cell volume 
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j= 0;  j < Ny; ++j)
		{
			for (int  k= 0;  k < Nz; ++k)
			{
				cell_volume[i][j][k] = deltatx*deltaty*deltatz ;
			}
		}	
	}
}