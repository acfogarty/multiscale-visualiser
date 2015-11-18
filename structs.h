#ifndef STRUCT_H
#define STRUCT_H

//############################### class Point2 #######################
class Point2{
	public:float x,y;
	Point2(float xx, float yy){x = xx; y = yy;}
	Point2(){ x = y = 0;}
	void set(float xx, float yy){ x = xx; y = yy;}
};

//@@@@@@@@@@@@@@@@@@ Point3 class @@@@@@@@@@@@@@@@
class Point3{ 
public: 
	float x,y,z;
	void set(float dx, float dy, float dz){x = dx; y = dy; z = dz;}
	void set(Point3& p){x = p.x; y = p.y; z = p.z;} 	
        Point3(float xx, float yy, float zz){x = xx; y = yy; z = zz;}
	Point3(){x = y = z = 0;}
	void build4tuple(float v[])
	{// load 4-tuple with this color: v[3] = 1 for homogeneous
		v[0] = x; v[1] = y; v[2] = z; v[3] = 1.0f;
	}	
};
//@@@@@@@@@@@@@@@@@@ Vector3 class @@@@@@@@@@@@@@@@
class Vector3{ 
public: 
	float x,y,z;
	void set(float dx, float dy, float dz){ x = dx; y = dy; z = dz;} 
	void set(Vector3& v){ x = v.x; y = v.y; z = v.z;}
	void flip(){x = -x; y = -y; z = -z;} // reverse this vector
	void setDiff(Point3& a, Point3& b)//set to difference a - b
	    { x = a.x - b.x; y = a.y - b.y; z = a.z - b.z;}
	void normalize();//adjust this vector to unit length
	//Vector3(float xx, float yy, float zz){x = xx; y = yy; z = zz;}
	Vector3(float xx, float yy, float zz):x(xx),y(yy),z(zz) {}
	Vector3(const Vector3& v):x(v.x),y(v.y),z(v.z) {}
	Vector3(){x = y = z = 0.0;} //default constructor
	Vector3 cross(Vector3 b); //return this cross b
	float dot(Vector3 b); // return this dotted with b
};

#endif
