#ifndef _POLARPOINT  
#define _POLARPOINT  

class PolarPoint{
	double phi_;
	double theta_;
	double distance_;

public:
	PolarPoint(){phi_=0.0;theta_=0.0;distance_=0.0;}; //constructor
	PolarPoint(double,double,double); //constructor
	PolarPoint(const PolarPoint&); //initialisation from another point

	double& phi() {return phi_;}
	double& theta() {return theta_;}
	double& distance() {return distance_;}

	PolarPoint & operator=(const PolarPoint&);//assignment
	bool operator>(const PolarPoint &);//>
	bool operator<(const PolarPoint &);//<
};

#endif