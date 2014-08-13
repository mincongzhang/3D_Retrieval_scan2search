#ifndef _POLARPOINT  
#define _POLARPOINT  

class PolarPoint{
	double _phi;
	double _theta;
	double _distance;

public:
	PolarPoint(){_phi=0.0;_theta=0.0;_distance=0.0;}; //constructor
	PolarPoint(double,double,double); //constructor
	PolarPoint(const PolarPoint&); //initialisation from another point

	double& phi() {return _phi;}
	double& theta() {return _theta;}
	double& distance() {return _distance;}

	PolarPoint & operator=(const PolarPoint&);//assignment
	bool operator>(const PolarPoint &);//>
	bool operator<(const PolarPoint &);//<
};

#endif