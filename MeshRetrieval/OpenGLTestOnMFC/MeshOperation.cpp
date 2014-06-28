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
	//TEST	SH
	double SH;
	SH = gsl_sf_legendre_sphPlm(1, 1, cos(M_PI/4));
	//TEST END
	NOISE_CONTROL = false;
}

/*normalize the model and rasterize to 2R*2R*2R voxel grid*/
void NormalizeMesh(MyMesh &mesh,vector<float> &grid_id_x,vector<float> &grid_id_y,vector<float> &grid_id_z,
				   vector<float> &dist_vector)
{
	float x_max,y_max,z_max,x_min,y_min,z_min;
	FindMaxMin(mesh,x_max,y_max,z_max,x_min,y_min,z_min);

	float distance_x = x_max - x_min;
	float distance_y = y_max - y_min;
	float distance_z = z_max - z_min;
	float max_distance = distance_x;

	if (distance_y > max_distance) max_distance = distance_y;
	if (distance_z > max_distance) max_distance = distance_z;

	/*normalize and rasterize*/
	//initial centroid
	float centroid_x = 0.0,centroid_y = 0.0,centroid_z = 0.0;

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
			float x_normalize = (mesh.point(v_it).data()[0] - x_min)/max_distance;
			float y_normalize = (mesh.point(v_it).data()[1] - y_min)/max_distance;
			float z_normalize = (mesh.point(v_it).data()[2] - z_min)/max_distance;

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
				grid_id_x.push_back(float(x_rasterize));
				grid_id_y.push_back(float(y_rasterize));
				grid_id_z.push_back(float(z_rasterize));

				//get centroid
				centroid_x += float(x_rasterize);
				centroid_y += float(y_rasterize);
				centroid_z += float(z_rasterize);
			}
		}
	}//end if(grid!=NULL)
	//delete grid
	delete [] grid;
	
	centroid_x/=float(grid_id_x.size());
	centroid_y/=float(grid_id_y.size());
	centroid_z/=float(grid_id_z.size());

	//move centroid to origin and get distance
	float mean_dist = 0.0;
	for (unsigned int grid_iter = 0;grid_iter<grid_id_x.size();grid_iter++)
	{
		//move centroid to origin
		grid_id_x.at(grid_iter) -= centroid_x;
		grid_id_y.at(grid_iter) -= centroid_y;
		grid_id_z.at(grid_iter) -= centroid_z;

		//get distance between vertex and origin
		float current_dist = sqrt(pow(grid_id_x.at(grid_iter),2) + pow(grid_id_y.at(grid_iter),2) + pow(grid_id_z.at(grid_iter),2));
		mean_dist += current_dist;
		dist_vector.push_back(current_dist);
	}
	mean_dist/=float(grid_id_x.size());

	//scale and make the average distance to center of mass is R/2
	float scale_ratio = (float(RADIUS)/2.0)/mean_dist;
	for (unsigned int grid_iter = 0;grid_iter<grid_id_x.size();grid_iter++)
	{
		grid_id_x.at(grid_iter) *= scale_ratio;
		grid_id_y.at(grid_iter) *= scale_ratio;
		grid_id_z.at(grid_iter) *= scale_ratio;
	}

	//TEST
	//float mean_dist_test = 0.0;
	//for (unsigned int grid_iter = 0;grid_iter<grid_id_x.size();grid_iter++)
	//{
	//	float dist_temp = sqrt(pow(grid_id_x.at(grid_iter),2) + pow(grid_id_y.at(grid_iter),2) + pow(grid_id_z.at(grid_iter),2));
	//	mean_dist_test += dist_temp;
	//}
	//mean_dist_test /= float(grid_id_x.size());
	//int test_done = 1;
	//TEST END

	NORMALIZE_CONTROL = FALSE;
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