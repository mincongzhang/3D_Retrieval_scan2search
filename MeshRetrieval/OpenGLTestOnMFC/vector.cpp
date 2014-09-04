#include "stdafx.h"
#include "Vector.h"
#include "math.h"


Vector::Vector(float x, float y, float z)
//constructor
{
	X = x;
	Y = y;
	Z = z;
}

Vector::Vector(const Vector& p)
//initialisation from another Vector
{
	X = p.X;
	Y = p.Y;
	Z = p.Z;
}

Vector& Vector::operator=(const Vector& p)
//assignment
{
	
	if(this == &p) return *this;
	
	X = p.X;
	Y = p.Y;
	Z = p.Z;
	
	return *this;
}

float Vector::squarednorm()
{
	return X*X+Y*Y+Z*Z;
}


float Vector::norm()
//returns norm of vector
{
	return (float)sqrt(squarednorm());
}

void Vector::normalise()
//normalises the vector
{
	float d = this->norm();
	
	X /= d;
	Y /= d;
	Z /= d;
}

Vector& Vector::normalised()
//returns the normalised vector
{
	float d = this->norm();
	
	X /= d;
	Y /= d;
	Z /= d;
	
	return *this;
}
	

ostream& operator<<(ostream& s,Vector p)
//writing
{	
	s << "[" << p.X << "," << p.Y << "," << p.Z << "]";
	return s;
}


istream& operator>>(istream& s,Vector& p)
//reading
{
	s >> p.X;
	s >> p.Y;
	s >> p.Z;
	return s;
}

Vector operator+(Vector& v1, Vector& v2)
//returns v1 + v2
{
	Vector v;
	
	v.x() = v1.x()+v2.x();
	v.y() = v1.y()+v2.y();
	v.z() = v1.z()+v2.z();
	
	return v;
}

Vector vadd(Vector& p1, Vector& p2)
//returns p1 + p2
{
	Vector v;
	
	v.x() = p1.x()+p2.x();
	v.y() = p1.y()+p2.y();
	v.z() = p1.z()+p2.z();
	
	return v;
}

Vector operator-(Vector& v1, Vector& v2)
//returns v1 - v2
{
	Vector v;
	
	v.x() = v1.x()-v2.x();
	v.y() = v1.y()-v2.y();
	v.z() = v1.z()-v2.z();
	
	return v;
}

Vector vminus(Vector& p1, Vector& p2)
//returns p1 + p2
{
	Vector v;
	
	v.x() = p1.x()-p2.x();
	v.y() = p1.y()-p2.y();
	v.z() = p1.z()-p2.z();
	
	return v;
}

float operator^(Vector& u, Vector& v)
//returns dot product
{
	return u.X*v.X + u.Y*v.Y + u.Z*v.Z;
}

Vector operator*(Vector& v1, Vector& v2)
//returns cross product
{
	Vector vout;
	
	vout.x() = v1.Y*v2.Z - v1.Z*v2.Y;
	vout.y() = v1.Z*v2.X - v1.X*v2.Z;
	vout.z() = v1.X*v2.Y - v1.Y*v2.X;
	
	return vout;
}
	
	
Vector operator*(Vector& v1, float a)
//vector by scalar
{
	Vector vout;
	
	vout.X = v1.X*a;
	vout.Y = v1.Y*a;
	vout.Z = v1.Z*a;
	
	return vout;
}
