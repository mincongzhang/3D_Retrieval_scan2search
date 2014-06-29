#include "stdafx.h"
#include "MeshOperation.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <random>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

double candidate_index_array[DATASIZE] = {}; 
using namespace std; // make std:: accessible

/*Add random Gaussian Noise to verteices*/
void AddNoise(double noise_standard_deviation,MyMesh &mesh)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,noise_standard_deviation); //Gaussian distribution: mean value = 0.0

	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it)
	{
		double Pt[3] = {};
		for (int d=0;d<3;d++)
		{
			Pt[d]=*(mesh.point(it).data()+d);
			double randn = distribution(generator);
			if ((randn>=-1.0)&&(randn<=1.0))//Gaussian distribution range [-1.0,1.0]							        
			{
				Pt[d]= Pt[d]*(1.0+randn);
				*(mesh.point(it).data()+d)=float(Pt[d]);
			}
		}
	}
	NOISE_CONTROL = false;
}

/*normalize the model and rasterize to 2R*2R*2R voxel grid*/
void NormalizeMesh(MyMesh &mesh,vector<double> &grid_id_x,vector<double> &grid_id_y,vector<double> &grid_id_z,
				   vector<double> &dist_vector)
{
	double x_max,y_max,z_max,x_min,y_min,z_min;
	FindMaxMin(mesh,x_max,y_max,z_max,x_min,y_min,z_min);

	double distance_x = x_max - x_min;
	double distance_y = y_max - y_min;
	double distance_z = z_max - z_min;
	double max_distance = distance_x;

	if (distance_y > max_distance) max_distance = distance_y;
	if (distance_z > max_distance) max_distance = distance_z;

	/*normalize and rasterize*/
	//initial centroid
	double centroid_x = 0.0,centroid_y = 0.0,centroid_z = 0.0;

	//create grid 
	bool *grid;
	int grid_size = 2*RADIUS*2*RADIUS*2*RADIUS;
	grid = new bool [grid_size];
	memset(grid,false,grid_size*sizeof(bool));

	if(grid!=NULL)
	{
		for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end(); ++v_it)
		{
			//move to positive, normalize to 1 
			double x_normalize = (mesh.point(v_it).data()[0] - x_min)/max_distance;
			double y_normalize = (mesh.point(v_it).data()[1] - y_min)/max_distance;
			double z_normalize = (mesh.point(v_it).data()[2] - z_min)/max_distance;

			*(mesh.point(v_it).data()+0) = x_normalize;
			*(mesh.point(v_it).data()+1) = y_normalize;
			*(mesh.point(v_it).data()+2) = z_normalize;

			//rasterize to 2Rx2Rx2R voxel grid
			int x_rasterize = static_cast<int>(round(x_normalize*(2*RADIUS-1)));
			int y_rasterize = static_cast<int>(round(y_normalize*(2*RADIUS-1)));
			int z_rasterize = static_cast<int>(round(z_normalize*(2*RADIUS-1)));
		
			//If this vertex hasn't been registered
			if(grid[x_rasterize*2*RADIUS*2*RADIUS + y_rasterize*2*RADIUS + z_rasterize]!=true)
			{
				//register
				grid[x_rasterize*2*RADIUS*2*RADIUS + y_rasterize*2*RADIUS + z_rasterize]=true;

				//push back
				grid_id_x.push_back(double(x_rasterize));
				grid_id_y.push_back(double(y_rasterize));
				grid_id_z.push_back(double(z_rasterize));

				//get centroid
				centroid_x += double(x_rasterize);
				centroid_y += double(y_rasterize);
				centroid_z += double(z_rasterize);
			}
		}
	}//end if(grid!=NULL)
	//delete grid
	delete [] grid;
	
	centroid_x/=double(grid_id_x.size());
	centroid_y/=double(grid_id_y.size());
	centroid_z/=double(grid_id_z.size());

	//move centroid to origin and get distance
	float mean_dist = 0.0;
	for (unsigned int grid_iter = 0;grid_iter<grid_id_x.size();grid_iter++)
	{
		//move centroid to origin
		grid_id_x.at(grid_iter) -= centroid_x;
		grid_id_y.at(grid_iter) -= centroid_y;
		grid_id_z.at(grid_iter) -= centroid_z;

		//get distance between vertex and origin
		double current_dist = sqrt(pow(grid_id_x.at(grid_iter),2) + pow(grid_id_y.at(grid_iter),2) + pow(grid_id_z.at(grid_iter),2));
		mean_dist += current_dist;
		dist_vector.push_back(current_dist);
	}
	mean_dist/=double(grid_id_x.size());

	//scale and make the average distance to center of mass is R/2
	float scale_ratio = (double(RADIUS)/2.0)/mean_dist;
	for (unsigned int grid_iter = 0;grid_iter<grid_id_x.size();grid_iter++)
	{
		grid_id_x.at(grid_iter) *= scale_ratio;
		grid_id_y.at(grid_iter) *= scale_ratio;
		grid_id_z.at(grid_iter) *= scale_ratio;
	}

	NORMALIZE_CONTROL = FALSE;
}

/*compute spherical harmonics*/
void ComputeSpharm(vector<double> &grid_id_x,vector<double> &grid_id_y,vector<double> &grid_id_z,vector<double> &dist_vector)
{
	vector<double> phi_vector,theta_vector;
	bool get_polar,get_sorted;
	get_polar = GetPolarCoordinate(grid_id_x,grid_id_y,grid_id_z,dist_vector,phi_vector,theta_vector);
	get_sorted = qsortPolarCoordinate(0,(dist_vector.size()-1),dist_vector,phi_vector,theta_vector);

	//test sorting
	//vector<double> vector1;
	//vector<double> vector2;
	//vector<double> vector3;
	//for(int testid = 7;testid>0;testid--)
	//{
	//	vector1.push_back(double(testid));
	//	vector2.push_back(double(testid));
	//	vector3.push_back(double(testid));
	//}
	//bool testresult = qsortPolarCoordinate(0,(vector1.size()-1),vector1,vector2,vector3);
	//int testend = 0;

	if(get_polar&&get_sorted)
	{
		//TEST	SH
		double SH;
		SH = gsl_sf_legendre_sphPlm(1, 1, cos(M_PI/4));
		//TEST END
	 }

}

void ChooseCandidate(double candidate_index_array[],int candidateIndx)
{
	int index = candidateIndx;

	if(index>=10) index = index-10;

	int CandidateIdx = int(candidate_index_array[DATASIZE-index-1]);
	//number to string
	string CandidateIdx_S = static_cast<ostringstream*>( &(ostringstream() << CandidateIdx) )->str();
	string back_filname = "./MeshData/back/back"+CandidateIdx_S+".obj";
	string seat_filname = "./MeshData/seat/seat"+CandidateIdx_S+".obj";
	string leg_filname = "./MeshData/leg/leg"+CandidateIdx_S+".obj";

	//load candidate mesh 
	MyMesh back_mesh,seat_mesh,leg_mesh;
	OpenMesh::IO::read_mesh(back_mesh, back_filname);
	OpenMesh::IO::read_mesh(seat_mesh, seat_filname);
	OpenMesh::IO::read_mesh(leg_mesh,  leg_filname);
	meshQueue.clear();
	if(candidateIndx<10)
	{
		meshQueue.push_back(seat_mesh);
		meshQueue.push_back(back_mesh);	
		meshQueue.push_back(leg_mesh);
	}
	if(candidateIndx>=10)
	{
		meshQueue.push_back(back_mesh);	
		meshQueue.push_back(seat_mesh);
		meshQueue.push_back(leg_mesh);
	}
}