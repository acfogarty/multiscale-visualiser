#include <iostream>
#include <math.h>
#include "structs.h"

// Vector3 methods
Vector3 Vector3::cross(Vector3 b) //return this cross b
{
    Vector3 c(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
    return c;
}

float Vector3 :: dot(Vector3 b) // return this dotted with b
{
    return x * b.x + y * b.y + z * b.z;
}

//******************************************************************
void Vector3::normalize()//adjust this vector to unit length
{		
    double sizeSq = x * x + y * y + z * z;
    if(sizeSq < 0.0000001)
    {
	std::cerr << "\nnormalize() sees vector (0,0,0)!";
 	return; // does nothing to zero vectors;
    }
    float scaleFactor = 1.0/(float)sqrt(sizeSq);
    x *= scaleFactor; y *= scaleFactor; z *= scaleFactor;
}

float Vector3::length()
{
    double sizeSqr = x * x + y * y + z * z;
    return sqrt(sizeSqr);
}
