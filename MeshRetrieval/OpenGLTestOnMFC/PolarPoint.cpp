#include "stdafx.h"
#include "PolarPoint.h"

PolarPoint::PolarPoint(double phi,double theta,double distance){
	//constructor
	phi_ = phi;
	theta_ = theta;
	distance_ = distance;
}

PolarPoint::PolarPoint(const PolarPoint & p){
	//initialisation from another point
	phi_ = p.phi_;
	theta_ = p.theta_;
	distance_ = p.distance_;
}

PolarPoint & PolarPoint::operator=(const PolarPoint & p){
	if(this == & p)	 return *this;
	phi_ = p.phi_;
	theta_ = p.theta_;
	distance_ = p.distance_;

	return *this;
}

bool PolarPoint::operator>(const PolarPoint & p){
	if(this->distance_ > p.distance_) 
		return true;
	else
		return false;
}

bool PolarPoint::operator<(const PolarPoint & p){
	if(this->distance_ < p.distance_) 
		return true;
	else
		return false;
}