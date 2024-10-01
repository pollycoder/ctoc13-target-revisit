#ifndef _MATRIX3D_H_
#define _MATRIX3D_H_


#include<iostream>
#include<assert.h>
#include<math.h>
#include"MyVector.h"
//////////////////////////////////////////////////////////////////////
//CLASS Mat3D
 
class Mat3D
{	
public:
	Mat3D();//按9个元素的一维数组构造.
	Mat3D(double angle, int axis);//构造正交矩阵,angle为逆时针转解,axis=1,2,3表示x,y,z轴
	Mat3D(const Mat3D& TheMat3D);
	Mat3D(const double& m11, const double& m12, const double& m13, const double& m21, const double& m22, const double& m23,
		   const double& m31, const double& m32, const double& m33);
	~Mat3D();
	
	inline double& operator()(int n, int m);//n=1,2,3;m=1,2,3.
	inline const double& operator()(int n, int m) const;
	inline Mat3D& operator=(const Mat3D& Mat);

	void SetArray(int n, const double* Array);
	void SetArray(const MyVec& Array);
	void GetArray(int n, double* Array) const;
	MyVec GetArray() const;

	double Det() const;//求行列式值
	double trace() const;
	Mat3D Inver() const;//求逆
	Mat3D Trans() const;//求转置



	inline bool operator==(const Mat3D& Mat) const;
	inline bool operator!=(const Mat3D& Mat) const;



	const Mat3D	operator +() const;	
	const Mat3D&	operator+=(const Mat3D& Mat);
	const Mat3D	operator +(const Mat3D& Mat) const;
	const Mat3D&	operator+=(const double& Num);
	const Mat3D	operator +(const double& Num) const;
	friend Mat3D operator+(const double& Num, const Mat3D& Mat);

	const Mat3D	operator -() const;
	const Mat3D&	operator-=(const Mat3D& Mat);
	const Mat3D	operator -(const Mat3D& Mat) const;
	const Mat3D&	operator-=(const double& Num);
	const Mat3D	operator -(const double& Num) const;
	friend Mat3D operator-(const double& Num, const Mat3D& Mat);

	const Mat3D&	operator*=(const Mat3D& Mat);
	const Mat3D	operator*(const Mat3D& Mat) const;
	const MyVec operator*(const MyVec& Vec) const;
	const Mat3D& operator*=(const double& Num);
	const Mat3D	operator*(const double& Num) const;	
	friend MyVec operator*(const MyVec& Vec, const Mat3D& Mat);
	friend Mat3D operator*(const double& Num, const Mat3D& Mat);

	const Mat3D&	operator/=(const double& Num);
	const Mat3D	operator/(const double& Num) const;

	const Mat3D Absol() const;
	const double Max() const;
	const double Min() const;
	

//全局运算符重载函数
	//inline ostream & operator<<(ostream &output, Mat3D& Mat);
	//inline istream & operator>>(istream &input, Mat3D& Mat);

private:
	double* array;
};

inline double& Mat3D::operator()(int n, int m)
{
	assert(n>=1&&n<=3&&m>=1&&m<=3);
	int i=(n-1)*3+m-1;
	return array[i];
}


inline const double&	Mat3D::operator()(int n, int m) const
{
	assert(n>=1&&n<=3&&m>=1&&m<=3);
	int i=(n-1)*3+m-1;
	return array[i];
}

inline Mat3D& Mat3D::operator=(const Mat3D& Mat)
{
	for(int i=0;i<9;i++)
		array[i]=Mat.array[i];
	return *this;
}

inline bool Mat3D::operator==(const Mat3D& Mat) const
{	
	for(int i=0;i<9;i++)
	{
		if (array[i]!=Mat.array[i])
			return false;
	}
	return true;
}

inline bool Mat3D::operator!=(const Mat3D& Mat) const
{
	return(!(*this==Mat));
}


std::istream & operator>>(std::istream &input, Mat3D& Mat);

std::ostream & operator<<(std::ostream &output, Mat3D& Mat);


#endif