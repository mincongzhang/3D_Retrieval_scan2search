#ifndef _VECTOR
#define _VECTOR

#include <iostream>

using namespace std;

#include "point.h"

class Vector{
	float X;
	float Y;
	float Z;
	
public:
	Vector() {X=Y=Z=0.0;}
	Vector(float,float,float);
	
	Vector(const Vector&);	  		//initialisation from another
	Vector& operator=(const Vector&); 	//assignment
	
	float squarednorm();			//returns square of norm
	float norm();				//returns norm
	void normalise();			//normalises 
	Vector& normalised();			//returns normalised vector
	friend Vector vadd(Vector&,Vector&);	       	//add two vectors
	friend Vector vminus(Vector&,Vector&);	       	//subtract
	friend float operator^(Vector&,Vector&);   //dot product
	friend Vector operator*(Vector&,Vector&); //cross product
	friend Vector operator*(Vector&, float);//Vector = vector * scalar
	
	float &x() {return X;}
	float &y() {return Y;}
	float &z() {return Z;}
		
	friend ostream& operator<<(ostream&,Vector);  //writing
	friend istream& operator>>(istream&,Vector&); //reading
};

#endif