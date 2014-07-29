#include "stdafx.h"
#include "MeshOperation.h"
#include "point.h"
#include "vector.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <random>
#include <string>

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

int fillGridLine(Point point1,Point point2,vector<Point> &inter12,bool grid[],vector<Point> &grid_points)
{
	Vector v12 = point2-point1;
	double len12 = v12.norm();
	v12.normalise();

	int n_points = 0;
	for(unsigned int n=0;n<=floor(len12);n++)
	{
		Point temp_p;
		temp_p.x() = round(point1.x()+v12.x()*double(n));
		temp_p.y() = round(point1.y()+v12.y()*double(n));
		temp_p.z() = round(point1.z()+v12.z()*double(n));

		int grid_coordinate	= int(temp_p.x()*2*RADIUS*2*RADIUS + temp_p.y()*2*RADIUS + temp_p.z());
		if(grid[grid_coordinate]!=true)
		{
			grid[grid_coordinate] = true;
			inter12.push_back(temp_p);
			grid_points.push_back(temp_p);
			n_points++;
		}
	}

	int grid_coordinate_p2 = int(point2.x()*2*RADIUS*2*RADIUS + point2.y()*2*RADIUS + point2.z());
	if(grid[grid_coordinate_p2]!=true)
	{
		grid[grid_coordinate_p2] = true;
		inter12.push_back(point2);
		grid_points.push_back(point2);
		n_points++;
	} 

	return n_points;
}

/*write the output SH*/
void WriteSH(string &write_filename,const int &max_r,const int &max_l,double SH_descriptor[])
{
	// open file
	ofstream myfile;
	myfile.open(write_filename);
	//write file header
	if (myfile.is_open())
	{
		//write SH
		//(idx_r-1)*max_l+idx_l
		for(unsigned int idx_r=0;idx_r<max_r;idx_r++)
		{
			for(unsigned int idx_l=0;idx_l<max_l;idx_l++)
			{
				myfile << SH_descriptor[idx_r*max_l+idx_l]<< " ";
			}
			myfile << "\n";
		}
		myfile.close();
	}
}

/*normalize the model */
void NormalizeMesh(MyMesh &mesh)
{
	double x_max,y_max,z_max,x_min,y_min,z_min;
	FindMaxMin(mesh,x_max,y_max,z_max,x_min,y_min,z_min);

	double distance_x = x_max - x_min;
	double distance_y = y_max - y_min;
	double distance_z = z_max - z_min;
	double max_distance = distance_x;

	if (distance_y > max_distance) max_distance = distance_y;
	if (distance_z > max_distance) max_distance = distance_z;

	/*normalize*/
	//initial centroid
	Point centroid(0.0,0.0,0.0);
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end(); ++v_it)
	{
		//move to positive, normalize to 1 
		double x_normalize = (mesh.point(v_it).data()[0] - x_min)/max_distance;
		double y_normalize = (mesh.point(v_it).data()[1] - y_min)/max_distance;
		double z_normalize = (mesh.point(v_it).data()[2] - z_min)/max_distance;

		*(mesh.point(v_it).data()+0) = x_normalize;
		*(mesh.point(v_it).data()+1) = y_normalize;
		*(mesh.point(v_it).data()+2) = z_normalize;

		//get centroid
		centroid.x() += double(x_normalize);
		centroid.y() += double(y_normalize);
		centroid.z() += double(z_normalize);
	}
	centroid.x()/=double(mesh.n_vertices());
	centroid.y()/=double(mesh.n_vertices());
	centroid.z()/=double(mesh.n_vertices());

	//get distance
	float mean_dist = 0.0;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end(); ++v_it)
	{
		//get distance between vertex and origin
		double current_dist = sqrt(pow(*(mesh.point(v_it).data()+0)-centroid.x(),2.0) 
			+ pow(*(mesh.point(v_it).data()+1)-centroid.y(),2.0) + pow(*(mesh.point(v_it).data()+2)-centroid.z(),2.0));
		mean_dist += current_dist;
	}
	mean_dist/=double(mesh.n_vertices());

	//scale and make the average distance to center of mass is R/2
	float scale_ratio = (double(RADIUS)/2.0)/mean_dist;
	for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end(); ++v_it)
	{
		*(mesh.point(v_it).data()+0) *= scale_ratio;
		*(mesh.point(v_it).data()+1) *= scale_ratio;
		*(mesh.point(v_it).data()+2) *= scale_ratio;

		//scale back to range [0,1] for better viewing 
		*(mesh.point(v_it).data()+0) /= double(2*RADIUS);
		*(mesh.point(v_it).data()+1) /= double(2*RADIUS);
		*(mesh.point(v_it).data()+2) /= double(2*RADIUS);
	}

	NORMALIZE_CONTROL = FALSE;
}

/*rasterize the model to 2R*2R*2R voxel grid*/
void RasterizeMesh(MyMesh &mesh,vector<Point> &grid_points,vector<double> &dist_vector)
{
	//create grid 
	bool *grid;
	int grid_size = 2*RADIUS*2*RADIUS*2*RADIUS;
	grid = new bool [grid_size];
	memset(grid,false,grid_size*sizeof(bool));

	for(MyMesh::FaceIter f_it=mesh.faces_begin();f_it!=mesh.faces_end();++f_it)
	{
		vector<Point> face_points;

		for(MyMesh::FaceVertexIter v_it=mesh.fv_iter(f_it);v_it;++v_it)
		{
			Point temp_face_point;
			//scale back from range [0,1] to [0,2R]
			temp_face_point.x() = *(mesh.point(v_it).data()+0)*2*RADIUS;
			temp_face_point.y() = *(mesh.point(v_it).data()+1)*2*RADIUS;
			temp_face_point.z() = *(mesh.point(v_it).data()+2)*2*RADIUS;
			face_points.push_back(temp_face_point);

			//get grid coordinate
			int grid_coordinate	= int(temp_face_point.x()*2*RADIUS*2*RADIUS + temp_face_point.y()*2*RADIUS + temp_face_point.z());

			//discard points out of range
			if(grid_coordinate>(2*RADIUS*2*RADIUS*2*RADIUS-1) || grid_coordinate<0)  continue; 

			//if not registered
			if(grid[grid_coordinate]!=true)
			{
				grid[grid_coordinate] = true;
				grid_points.push_back(temp_face_point);
			}
		}

		//fill the edges for a triangle face
		vector<Point> inter12,inter23,inter13;
		int n_points12 = fillGridLine(face_points.at(0),face_points.at(1),inter12,grid,grid_points);
		int n_points23 = fillGridLine(face_points.at(1),face_points.at(2),inter23,grid,grid_points);
		int n_points13 = fillGridLine(face_points.at(0),face_points.at(2),inter13,grid,grid_points);

		//fill in the face
		if(n_points12>=n_points13 && n_points13>0)
		{
			for(unsigned int n=0;n<n_points13;n++)
			{
				vector<Point> temp_inter;
				fillGridLine(inter12.at(n),inter13.at(n),temp_inter,grid,grid_points);
			}
			for(unsigned int n=n_points13;n<n_points12;n++)
			{
				vector<Point> temp_inter;
				fillGridLine(inter12.at(n),inter13.at(inter13.size()-1),temp_inter,grid,grid_points);
			}
		} 
		else if(n_points13>=n_points12 && n_points12>0)
		{
			for(unsigned int n=0;n<n_points12;n++)
			{
				vector<Point> temp_inter;
				fillGridLine(inter12.at(n),inter13.at(n),temp_inter,grid,grid_points);
			}
			for(unsigned int n=n_points12;n<n_points13;n++)
			{
				vector<Point> temp_inter;
				fillGridLine(inter12.at(inter12.size()-1),inter13.at(n),temp_inter,grid,grid_points);
			}
		}
	}

	//get centroid again
	Point centroid_after(0.0,0.0,0.0);
	for(unsigned int p_it = 0;p_it<grid_points.size();p_it++)
	{
		centroid_after.x()+=grid_points.at(p_it).x();
		centroid_after.y()+=grid_points.at(p_it).y();
		centroid_after.z()+=grid_points.at(p_it).z();

	}
	centroid_after.x()/=double(grid_points.size());
	centroid_after.y()/=double(grid_points.size());
	centroid_after.z()/=double(grid_points.size());

	//move grid to the origin and get distance vector
	for(unsigned int p_it = 0;p_it<grid_points.size();p_it++)
	{
		grid_points.at(p_it).x() -= centroid_after.x();
		grid_points.at(p_it).y() -= centroid_after.y();
		grid_points.at(p_it).z() -= centroid_after.z();

		double current_dist = pow(grid_points.at(p_it).x(),2.0)+pow(grid_points.at(p_it).y(),2.0)+pow(grid_points.at(p_it).z(),2.0);
		current_dist = sqrt(current_dist);
		dist_vector.push_back(current_dist);
	}

	delete [] grid;

	//TEST rotate xyz
	//double rotate_angle = M_PI/6.0;
	//for(unsigned int id=0;id<grid_points.size();id++)
	//{
	//	double x = grid_points.at(id).x();
	//	double y = grid_points.at(id).y();

	//	grid_points.at(id).x() = x*cos(rotate_angle)-y*sin(rotate_angle);
	//	grid_points.at(id).y() = x*sin(rotate_angle)+y*cos(rotate_angle);
	//}

	RASTERIZE_CONTROL = FALSE;
}

/*compute spherical harmonics*/
void ComputeSpharm(vector<Point> &grid_points,vector<double> &dist_vector,string write_filename)
{
	vector<double> phi_vector,theta_vector;	//radian
	bool get_polar,get_sorted;

	get_polar = GetPolarCoordinate(grid_points,dist_vector,phi_vector,theta_vector);
	get_sorted = qsortPolarCoordinate(0,(dist_vector.size()-1),dist_vector,phi_vector,theta_vector);

	if(get_polar&&get_sorted)
	{
		//begin
		int max_l = RADIUS;
		int max_r = RADIUS;
		int idx_n = 0;
		double *SH_descriptor;
		SH_descriptor = new double [max_r*max_l]();

		for(unsigned int idx_r = 1;idx_r<=RADIUS;idx_r++)
		{
			//for each frequency
			for(int idx_l = 0;idx_l<max_l;idx_l++)
			{
				vector<double> 	a_ml_pow;
				//traverse frequency range
				for(int idx_m = -idx_l;idx_m<=idx_l;idx_m++)
				{
					vector<double> Yml_real,Yml_imag;
					double NPml;
					while(dist_vector.at(idx_n)<idx_r)
					{
						//initial Y_ml = N*P(m,l,cos(theta))*exp(i*m*phi)
						//function gsl_sf_legendre_sphPlm returns value of N*P(m,l,cos(theta))
						//cos(x) input should be a radian
						//theta_vector and phi_vector are radians
						if(idx_m>=0)
						{
							NPml = gsl_sf_legendre_sphPlm(idx_l,idx_m,cos(theta_vector.at(idx_n)));
							Yml_real.push_back(NPml*cos(double(idx_m)*phi_vector.at(idx_n)));
							Yml_imag.push_back(NPml*sin(double(idx_m)*phi_vector.at(idx_n)));
						} else {
							NPml = pow (-1.0,-idx_m)*gsl_sf_legendre_sphPlm(idx_l,-idx_m,cos(theta_vector.at(idx_n)));
							Yml_real.push_back(NPml*cos(double(-idx_m)*phi_vector.at(idx_n)));
							Yml_imag.push_back(NPml*sin(double(-idx_m)*phi_vector.at(idx_n)));
						}
						idx_n++;
					}//while(dist_vector.at(idx_n)<idx_r)
					double 	a_ml_real,a_ml_imag;
					a_ml_real=accumulate(Yml_real.begin(),Yml_real.end(),0.0);
					a_ml_imag=accumulate(Yml_imag.begin(),Yml_imag.end(),0.0);
					a_ml_pow.push_back(pow(a_ml_real,2.0)+pow(a_ml_imag,2.0));
					if(idx_l!=max_l-1) idx_n-=Yml_real.size();
				}//finish traversing frequency range
				SH_descriptor[(idx_r-1)*max_l+idx_l] = sqrt(accumulate(a_ml_pow.begin(),a_ml_pow.end(),0.0));
			}	//end of for(int idx_l = 0;idx_l<max_l;idx_l++)
		}// for(unsigned int idx_r = 0;idx_r<RADIUS;idx_r++)

		//TEST 1 rotational invariant
		//double phi = 0.0;
		//double theta = 0.0;
		//for(int idx_r=1;idx_r<=max_r;idx_r++)
		//{
		//	theta+=M_PI/32.0;
		//	//for each frequency
		//	for(int idx_l = 0;idx_l<max_l;idx_l++)
		//	{
		//		//initial Y_ml = N*P(m,l,cos(theta))*exp(i*m*phi)
		//		//function gsl_sf_legendre_sphPlm returns value of N*P(m,l,cos(theta))
		//		//cos(x) input should be a radian
		//		//theta_vector and phi_vector are radian

		//		double Yml_real = 0.0;
		//		double Yml_img = 0.0;
		//		double NPml = 0.0;
		//		//traverse frequency range
		//		for(int idx_m = -idx_l;idx_m<=idx_l;idx_m++)
		//		{
		//			if(idx_m>=0)
		//			{
		//				NPml = gsl_sf_legendre_sphPlm(idx_l,idx_m,cos(theta));
		//				Yml_real += NPml*cos(double(idx_m)*phi);
		//				Yml_img  +=	NPml*sin(double(idx_m)*phi);
		//			}
		//			else
		//			{
		//				NPml = pow(-1.0,idx_m)*gsl_sf_legendre_sphPlm(idx_l,-idx_m,cos(theta));
		//				Yml_real += NPml*cos(double(-idx_m)*phi);
		//				Yml_img  -=	NPml*sin(double(-idx_m)*phi);
		//			}
		//		}//finish traversing frequency range
		//		SH_descriptor[(idx_r-1)*max_l+idx_l] += pow(Yml_real,2.0) + pow(Yml_img,2.0);
		//	}//end of for(int idx_l = 0;idx_l<max_l;idx_l++)
		//}

		////TEST 2 reconstruct
		//int re_l = 16;
		//double * YML_real;
		//YML_real = new double [2*re_l+1]();
		//double * YML_imag;
		//YML_imag = new double [2*re_l+1]();
		//double Yml_real = 0.0;
		//double Yml_img = 0.0;
		//double NPml = 0.0;
		//theta = M_PI/4.0;
		//phi = 0.0;
		//for(int idx_m=-re_l;idx_m<=re_l;idx_m++)
		//{
		//	if(idx_m>=0)
		//	{
		//		NPml = gsl_sf_legendre_sphPlm(re_l,idx_m,cos(theta));
		//		YML_real[idx_m+16] = NPml*cos(double(idx_m)*phi);
		//		YML_imag[idx_m+16] = NPml*sin(double(idx_m)*phi);
		//	}
		//	else
		//	{
		//		NPml = pow(-1.0,idx_m)*gsl_sf_legendre_sphPlm(re_l,-idx_m,cos(theta));
		//		YML_real[idx_m+16] = NPml*cos(double(-idx_m)*phi);
		//		YML_imag[idx_m+16] = -NPml*sin(double(-idx_m)*phi);
		//	}
		//}

		/*write out the data*/
		WriteSH(write_filename,max_r,max_l,SH_descriptor);

		delete [] SH_descriptor;
	}//end of if(get_polar&&get_sorted)

	SPHARM_CONTROL = false;
}

/*batch transform*/
void BatchTrans(void)
{
	int file_num = 5;

	for(unsigned int file_id=3;file_id<=file_num;file_id++)
	{
		string id = static_cast<ostringstream*>( &(ostringstream() << file_id) )->str();
		string read_filename = "./MeshData/data ("+id+").stl";
		string write_filename = "./MeshSHData/SH"+id+".txt";

		MyMesh mesh_to_transform;
		OpenMesh::IO::read_mesh(mesh_to_transform, read_filename);

		vector<Point> grid_points_to_transform;
		vector<double> dist_vector_to_transform;

		NormalizeMesh(mesh_to_transform);
		RasterizeMesh(mesh_to_transform,grid_points_to_transform,dist_vector_to_transform);
		ComputeSpharm(grid_points_to_transform,dist_vector_to_transform,write_filename);

	}

	//ifstream fin("hello.txt");
	//if (!fin)
	//{
	//	std::cout << "can not open this file" << endl;
	//}

	BATCH_CONTROL = false;
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