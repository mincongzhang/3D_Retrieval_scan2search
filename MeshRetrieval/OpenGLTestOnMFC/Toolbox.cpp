
#include "stdafx.h"
#include "OpenGLControl.h"
#include ".\openglcontrol.h"
#include "MeshOperation.h"

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

//Swap ith and jth elements in arrat[]
void swap(double array[], int i, int j)
{
	double temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}

//Quick Sort of arry and get the index of original order of array
void qsort_getid(double array[],double id_array[], int left_id, int right_id)
{
	if(left_id >= right_id)
		return;
	int flag = left_id;
	for(int i = left_id+1; i<=right_id; i++)
	{
		if(array[left_id]>array[i])
		{
			flag = flag+1;
			swap(array,flag, i);
			swap(id_array,flag, i);
		}
	}
	swap(array,flag,left_id);
	swap(id_array,flag,left_id);
	qsort_getid(array,id_array,left_id,flag-1);
	qsort_getid(array,id_array,flag+1,right_id);
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
void FindMaxMin(MyMesh &mesh, float &x_max, float &y_max, float &z_max, float &x_min, float &y_min, float &z_min)
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