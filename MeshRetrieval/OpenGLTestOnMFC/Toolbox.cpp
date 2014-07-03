#include "stdafx.h"
#include "OpenGLControl.h"
#include ".\openglcontrol.h"
#include "MeshOperation.h"

#include <math.h>
#include <random>

/*load histogram data*/
void loadHistogram(string filname,double *histogram)
{
	ifstream myfile (filname);
	string line;
	int count = 0;
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			*(histogram+count) = atof(line.c_str());
			count++;
		}
		myfile.close();
	}
}

/*Calculate similarity*/
//similarity_halfcircle = sum(hist_test.*hist)/(norm(hist_test)*norm(hist))
double similarity(double *histogram_test,double *histogram_sketch)
{
	double norm_test  = 0.0,norm_sketch = 0.0;
	double similarity = 0.0;
	for (int i = 0; i < 143; i++)
	{
		norm_test   += pow(*(histogram_test+i),2);
		norm_sketch += pow(*(histogram_sketch+i),2);		
	}
	norm_test   = sqrt(norm_test);
	norm_sketch = sqrt(norm_sketch);

	for (int i = 0; i < 143; i++)
	{
		similarity += (*(histogram_test+i)) * (*(histogram_sketch+i))/(norm_test*norm_sketch);
	}
	return similarity;
}

//swap
void swap(double &x,double &y)
{
	double temp = x;
	x = y;
	y = temp;
}

//Quick Sort of array and get the index of original order of array
void qsort_getid(double array[],double id_array[], int left_id, int right_id)
{
	if(left_id >= right_id)
		return;
	int flag = left_id;
	for(int i = left_id+1; i<=right_id; i++)
	{
		if(array[left_id]>array[i])
		{
			flag++;
			swap(array[flag],array[i]);
			swap(id_array[flag],id_array[i]);
		}
	}
	//swap(array,flag,left_id);
	swap(array[flag],array[left_id]);
	swap(id_array[flag],id_array[left_id]);
	qsort_getid(array,id_array,left_id,flag-1);
	qsort_getid(array,id_array,flag+1,right_id);
}

//according to distance vector sort phi and theta vectors
bool qsortPolarCoordinate(int left_id, int right_id,
						  vector<double> &dist_vector,vector<double> &phi_vector, vector<double> &theta_vector)
{
	if(left_id >= right_id)	return true;

	int flag = left_id;
	for(int i = left_id+1; i<=right_id; i++)
	{
		  if(dist_vector.at(left_id)>dist_vector.at(i))
		{
			flag++;
			swap(dist_vector.at(flag),dist_vector.at(i));
			swap(phi_vector.at(flag),phi_vector.at(i));
			swap(theta_vector.at(flag),theta_vector.at(i));
		}
	}

	swap(dist_vector.at(flag),dist_vector.at(left_id));
	swap(phi_vector.at(flag),phi_vector.at(left_id));
	swap(theta_vector.at(flag),theta_vector.at(left_id));

	qsortPolarCoordinate(left_id,flag-1,dist_vector,phi_vector,theta_vector);
	qsortPolarCoordinate(flag+1,right_id,dist_vector,phi_vector,theta_vector);
}

//get polar coordinate 
bool  GetPolarCoordinate(vector<double> &grid_id_x, vector<double> &grid_id_y,vector<double> &grid_id_z,vector<double> &dist_vector,
						 vector<double> &phi_vector, vector<double> &theta_vector)
{
	for (unsigned int i = 0; i<grid_id_x.size();i++)
	{
		//phi   = atan(y/x);
		//theta = acos(z/radious);
		//atan(sqrt(3.0))=1.047=pi/3
		//acos(0.5)      =1.047=pi/3
		double phi = atan(double(grid_id_y.at(i)/grid_id_x.at(i)));
		double theta = acos(double(grid_id_z.at(i)/dist_vector.at(i)));
		phi_vector.push_back(phi);
		theta_vector.push_back(theta);
	}

	if(phi_vector.size()==theta_vector.size() && phi_vector.size()>0)	 
		return true;
	else     															 
		return false;
}

//get the sum of a double vector
double getVectorSum(vector<double> input_vector) 
{
	double sum = 0.0;
	for(unsigned int i = 0;i<input_vector.size();i++)
	{
		sum += input_vector.at(i);
	}
	return sum;
}

/*Round*/
double round(double number)
{
	double temp = number;
	//ceil(1.7)-0.5=1.5<1.7;  ceil(1.3)-0.5=1.5>1.3;
	if((ceil(number)-0.5)>temp) number = floor(number);
	else				        number = ceil(number);
	return number;
}

/*find the max distance of the model*/
void FindMaxMin(MyMesh &mesh, double &x_max, double &y_max, double &z_max, double &x_min, double &y_min, double &z_min)
{
	//initial 
	MyMesh::VertexIter v_it1 = mesh.vertices_begin();
	x_max = mesh.point(v_it1).data()[0];
	y_max = mesh.point(v_it1).data()[1];
	z_max = mesh.point(v_it1).data()[2];
	x_min = mesh.point(v_it1).data()[0];
	y_min = mesh.point(v_it1).data()[1];
	z_min = mesh.point(v_it1).data()[2];

	for (MyMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end(); ++v_it)
	{
		//max
		if(mesh.point(v_it).data()[0]>x_max) x_max = mesh.point(v_it).data()[0];
		if(mesh.point(v_it).data()[1]>y_max) y_max = mesh.point(v_it).data()[1];
		if(mesh.point(v_it).data()[2]>z_max) z_max = mesh.point(v_it).data()[2];

		//min
		if(mesh.point(v_it).data()[0]<x_min) x_min = mesh.point(v_it).data()[0];
		if(mesh.point(v_it).data()[1]<y_min) y_min = mesh.point(v_it).data()[1];
		if(mesh.point(v_it).data()[2]<z_min) z_min = mesh.point(v_it).data()[2];
	}
}