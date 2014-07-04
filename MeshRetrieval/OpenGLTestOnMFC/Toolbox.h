#pragma once
#include "afxwin.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#undef min
#undef max
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/IO/exporter/ExporterT.hh>

using namespace std; // make std:: accessible
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

void loadHistogram(string filname,double *histogram);
double similarity(double *histogram_test,double *histogram_sketch);
void qsort_getid(double array[],double id_array[], int left_id, int right_id);
bool  GetPolarCoordinate(vector<double> &grid_id_x, vector<double> &grid_id_y,vector<double> &grid_id_z,vector<double> &dist_vector,
						 vector<double> &phi_vector, vector<double> &theta_vector);
bool qsortPolarCoordinate(int left_id, int right_id,
						  vector<double> &dist_vector,vector<double> &phi_vector, vector<double> &theta_vector);
double getVectorSum(vector<double> input_vector);
double round(double number);
void FindMaxMin(MyMesh &mesh, double &x_max, double &y_max, double &z_max, double &x_min, double &y_min, double &z_min);