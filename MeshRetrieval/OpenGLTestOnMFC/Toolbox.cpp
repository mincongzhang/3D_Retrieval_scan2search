#include "stdafx.h"
#include "OpenGLControl.h"
#include ".\openglcontrol.h"
#include "MeshOperation.h"

/*Quick Sort of array and get the index of original order of array*/
//test: can be replaced by map?
void getSortedID(double array[],double id_array[], int left_id, int right_id)
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
	getSortedID(array,id_array,left_id,flag-1);
	getSortedID(array,id_array,flag+1,right_id);
}

void  GetPolarCoordinate(vector<Point> &grid_points,vector<PolarPoint> &grid_polar_points)
{
	for (unsigned int i = 0; i<grid_points.size();i++)
	{
		//1. get distance
		double current_dist = pow(grid_points.at(i).x(),2.0)+pow(grid_points.at(i).y(),2.0)
			                 +pow(grid_points.at(i).z(),2.0);
		current_dist = sqrt(current_dist);

		//2. get phi and theta
		//phi   = atan(y/x);
		//theta = acos(z/radious);
		//atan(sqrt(3.0))=1.047=pi/3
		//acos(0.5)      =1.047=pi/3
		double phi = atan(double(grid_points.at(i).y()/grid_points.at(i).x()));
		double theta = acos(double(grid_points.at(i).z()/current_dist));

		//3. assign to polar point
		PolarPoint tmp_point(phi,theta,current_dist);
		grid_polar_points.push_back(tmp_point);
	}
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

/*get the sum of a double vector */
/*replaced by build in function accumulate*/
//double getVectorSum(vector<double> input_vector) 
//{
//	double sum = 0.0;
//	for(unsigned int i = 0;i<input_vector.size();i++)
//	{
//		sum += input_vector.at(i);
//	}
//	return sum;
//}

/*according to distance vector sort phi and theta vectors*/
/*replaced by PolarPoint class and STL sorting function*/
//bool qsortPolarCoordinate(int left_id, int right_id,
//						  vector<double> &dist_vector,vector<double> &phi_vector, vector<double> &theta_vector)
//{
//	if(left_id >= right_id)	return true;
//
//	int flag = left_id;
//	for(int i = left_id+1; i<=right_id; i++)
//	{
//		  if(dist_vector.at(left_id)>dist_vector.at(i))
//		{
//			flag++;
//			swap(dist_vector.at(flag),dist_vector.at(i));
//			swap(phi_vector.at(flag),phi_vector.at(i));
//			swap(theta_vector.at(flag),theta_vector.at(i));
//		}
//	}
//
//	swap(dist_vector.at(flag),dist_vector.at(left_id));
//	swap(phi_vector.at(flag),phi_vector.at(left_id));
//	swap(theta_vector.at(flag),theta_vector.at(left_id));
//
//	qsortPolarCoordinate(left_id,flag-1,dist_vector,phi_vector,theta_vector);
//	qsortPolarCoordinate(flag+1,right_id,dist_vector,phi_vector,theta_vector);
//}