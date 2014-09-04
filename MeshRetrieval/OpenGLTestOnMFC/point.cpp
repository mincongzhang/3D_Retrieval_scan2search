#include "stdafx.h"
#include "Point.h"
#include "Vector.h"

Point::Point(float x, float y, float z)
//constructor
{
	X = x;
	Y = y;
	Z = z;
}

Point::Point(const Point& p)
//initialisation from another point
{
	X = p.X;
	Y = p.Y;
	Z = p.Z;
}

Point& Point::operator=(const Point& p)
//assignment
{
	
	if(this == &p) return *this;
	
	X = p.X;
	Y = p.Y;
	Z = p.Z;
	
	return *this;
}

Vector operator-(Point& p1, Point& p2)
//returns p1 - p2
{
	Vector v(
		p1.x()-p2.x(),
		p1.y()-p2.y(),
		p1.z()-p2.z()
		);
	
	return v;
}

Vector operator+(Point& p1, Point& p2)
//returns p1 + p2
{
	Vector v(
		p1.x()+p2.x(),
		p1.y()+p2.y(),
		p1.z()+p2.z()
		);
	
	return v;
}

ostream& operator<<(ostream& s,Point p)
//writing
{	
	s << "(" << p.X << "," << p.Y << "," << p.Z << ")";
	return s;
}



istream& operator>>(istream& s,Point& p)
//reading
{
	s >> p.X;
	s >> p.Y;
	s >> p.Z;
	return s;
}
