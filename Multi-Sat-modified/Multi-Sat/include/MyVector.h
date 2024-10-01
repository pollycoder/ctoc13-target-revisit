#ifndef _MYVECTOR_H_
#define _MYVECTOR_H_


#include<iostream>
#include<assert.h>
#include<math.h>

//////////////////////////////////////////////////////////////////////
//TEMPLATE CLASS MyVec

class MyVec
{	
public:
	explicit MyVec(int arraysize=3);
	MyVec(const MyVec& TheVector);
	MyVec(const double& x, const double& y, const double& z);
	MyVec(const double& x, const double& y, const double& z, const double& vx, const double& vy, const double& vz);
	MyVec(const double& x, const double& y, const double& z, const double& vx, const double& vy, const double& vz, const double& m);
	~MyVec();
	 int GetSize() const;
	void SetValue(int n, const double* Value);
	void GetValue(int n, double* Value) const;

	 double& operator[](int i);
	 const double&	operator[](int i) const;
	 double& operator()(int i);
	 const double&	operator()(int i) const;

	 bool operator==(const MyVec& Vec) const;
	 bool operator!=(const MyVec& Vec) const;

	 MyVec& operator=(const MyVec& Vec);

	 const MyVec	operator +() const;	
	 const MyVec&	operator+=(const MyVec& Vec);
	 const MyVec	operator +(const MyVec& Vec) const;
	 const MyVec&	operator+=(const double& Num);
	 const MyVec	operator +(const double& Num) const;
	friend  MyVec operator+(const double& Num, const MyVec& Vec);

	 const MyVec	operator -() const;
	 const MyVec&	operator-=(const MyVec& Vec);
	 const MyVec	operator -(const MyVec& Vec) const;
	 const MyVec&	operator-=(const double& Num);
	 const MyVec	operator -(const double& Num) const;
	friend  MyVec operator-(const double& Num, const MyVec& Vec);

	 const MyVec&	operator*=(const MyVec& Vec);
	 const MyVec	operator*(const MyVec& Vec) const;
	 const MyVec	operator*(const double& Num) const;
	 const MyVec&	operator*=(const double& Num);
	friend  MyVec operator*(const double& Num, const MyVec& Vec);

	 const MyVec&	operator/=(const MyVec& Vec);
	 const MyVec	operator/(const MyVec& Vec) const;
	 const MyVec&	operator/=(const double& Num);
	 const MyVec	operator/(const double& Num) const;
	friend  MyVec operator/(const double& Num, const MyVec& Vec);

	 const double Dot(const MyVec& Vec) const;
	 const MyVec Cross(const MyVec& Vec) const;
	 const MyVec Absol() const;
	 const double Norm1() const;
	 const double NormInf() const;
	 const double Norm2() const;
	 const double Max() const;
	 const double Min() const;
	

//全局运算符重载函数
	// ostream & operator<<(ostream &output, MyVec& Vec);
	// istream & operator>>(istream &input, MyVec& Vec);

private:
	double* array;
	int size;

};

std::istream & operator>>(std::istream &input, MyVec& Vec);

std::ostream & operator<<(std::ostream &output, MyVec& Vec);

#endif