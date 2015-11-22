#include "iostream"
#include <vector>
#include "math.h"
#include <fstream>

// this class calculates the area and volues in the domain

// Function defines the area vector and cell volumes 
void grid( vector<vector<vector<vector<float> > > > & x_face_area ,
		  vector<vector<vector<vector<float> > > > & y_face_area, vector<vector<vector<vector<float> > > > & z_face_area ,
		  vector<vector<vector<float> > > & cell_volume, int Nx, int Ny, int Nz, float lenght, float delta)  
{
	// Grids for bump(this is the first case which will test the scheme)

	// Creating a 4D vector object for grid points
	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	// this store previous values of variables (density , three momentum, energy)
	matrix4D variablesvector(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 
	matrix4D grid_point(Nx+1,Dim3(Ny+1,Dim2(Nz+1,Dim1(3)))); 

	int N = floor(lenght/delta) ; 
	float deltatx = delta ;
	float deltaty = delta ;
	float deltatz = delta ;


	// First defining the grid points
	for (int i =0; i < Nx+1; ++i) 
	{
		for (int  j=0;  j < Ny+1; j++)
		{
			for (int  k=0;  k < Nz+1; k++)
			{
				grid_point[i][j][k][0] = i*deltatx ;   
				grid_point[i][j][k][2] = k*deltatz ;
				if (i < N+2)
				{
					grid_point[i][j][k][1] = j*deltaty ;
				}
				else if( i >= N+2 && i <= floor(3*N/2) + 2 )
				{
					grid_point[i][j][k][1] = j * (deltaty - 0.2*(i-N-2)/Ny)  ; 
				}
				else 
				{
					grid_point[i][j][k][1] = grid_point[3*N+4-i][j][k][1] ;
				}

			}
		}	
	}


	// this file is opend to store the grid poits
	ofstream kullu_grid ;
	kullu_grid.open("grids_Bump_2D.csv");
	kullu_grid <<  "x" << "," << "y" << "," << "z" <<   endl ;
	for (int i = 2; i < Nx-2; ++i)
	{
		for (int j = 2; j < Ny-2; ++j)
		{
			kullu_grid << grid_point[i][j][4][0] << "," << grid_point[i][j][4][1] << endl ; 
		}
	}

	// here comes the area vectors
	for (int i = 0; i  < Nx; ++i)
	{
		for (int  j = 0;  j < Ny; ++j)
		{
			for (int  k = 0;  k < Nz; ++k)
			{
				x_face_area[i][j][k][0] = (grid_point[i][j+1][k][1]-grid_point[i][j][k][1])*deltatz ;
				x_face_area[i][j][k][1] = 0 ;
				x_face_area[i][j][k][2] = 0 ;

				y_face_area[i][j][k][0] = -deltatz*(grid_point[i+1][j][k][1]-grid_point[i][j][k][1]) ;
				y_face_area[i][j][k][1] =  deltatz*(grid_point[i+1][j][k][0]-grid_point[i][j][k][0]) ;
				y_face_area[i][j][k][2] = 0 ;

				z_face_area[i][j][k][0] = 0 ; 
				z_face_area[i][j][k][1] = 0 ;
				z_face_area[i][j][k][2] = 0.5*deltatx*( (grid_point[i][j+1][k][1]-grid_point[i][j][k][1]) + 
					(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1]) ); 
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
				if (i < N+2)
				{
					cell_volume[i][j][k] = deltatx*deltaty*deltatz ; 
				}
				else if( i >= N+2 && i <= floor(3*N/2) + 2 )
				{
					cell_volume[i][j][k] = 0.5 * deltatx * ((grid_point[i][j+1][k][1] - grid_point[i][j][k][1])  + 
					(grid_point[i+1][j+1][k][1]-grid_point[i+1][j][k][1])) * deltatz ; 
				}
				else 
				{
					cell_volume[i][j][k] = cell_volume[3*N+4-i][j][k]  ;
				}
			}
		}	
	}
}

