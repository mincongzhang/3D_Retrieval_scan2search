#include "stdafx.h"
#include "PolarPoint.h"

PolarPoint::PolarPoint(double phi,double theta,double distance){
	//constructor
	_phi = phi;
	_theta = theta;
	_distance = distance;
}

PolarPoint::PolarPoint(const PolarPoint & p){
	//initialisation from another point
	_phi = p._phi;
	_theta = p._theta;
	_distance = p._distance;
}

PolarPoint & PolarPoint::operator=(const PolarPoint & p){
	if(this == & p)	 return *this;
	_phi = p._phi;
	_theta = p._theta;
	_distance = p._distance;

	return *this;
}

bool PolarPoint::operator>(const PolarPoint & p){
	if(this->_distance > p._distance) 
		return true;
	else
		return false;
}

bool PolarPoint::operator<(const PolarPoint & p){
	if(this->_distance < p._distance) 
		return true;
	else
		return false;
}