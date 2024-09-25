#include"Mat3D.h"
//#include "StdAfx.h"

//构造函数
Mat3D::Mat3D()
{
	array=new double[9];
	assert(array!=NULL);
	for(int i=0;i<9;i++)
		array[i]=0;
}

Mat3D::Mat3D(double angle, int axis)
{
	assert(axis==1||axis==2||axis==3);
	array=new double[9];
	assert(array!=NULL);
	for(int i=0;i<9;i++)
		array[i]=0;
	if(axis==1)
	{
		array[0]=1.0;
		array[4]=array[8]=cos(angle);
		array[5]=sin(angle);
		array[7]=-array[5];
	}
	if(axis==2)
	{
		array[4]=1.0;
		array[0]=array[8]=cos(angle);
		array[6]=sin(angle);
		array[2]=-array[6];
	}
	if(axis==3)
	{
		array[8]=1.0;
		array[0]=array[4]=cos(angle);
		array[1]=sin(angle);
		array[3]=-array[1];
	}
}

//副本构造函数
Mat3D::Mat3D(const Mat3D& TheMat3D)
{
	int n=9;
	array=new double[n];
	double* srcptr=TheMat3D.array;
	double* destptr=array;
	while(n--)
		*destptr++=*srcptr++;
//	for(int j=0;j<9;j++)
//		array[j]=TheMat3D.array[j];
}


//3X3矩阵副本构造函数
Mat3D::Mat3D(const double& m11, const double& m12, const double& m13, const double& m21, const double& m22, const double& m23,
		   const double& m31, const double& m32, const double& m33)
{
	array=new double[9];
	assert(array!=NULL);
	array[0]=m11;
	array[1]=m12;
	array[2]=m13;
	array[3]=m21;
	array[4]=m22;
	array[5]=m23;
	array[6]=m31;
	array[7]=m32;
	array[8]=m33;
}

//析构函数
Mat3D::~Mat3D()
{
	delete[] array;
}


void Mat3D::SetArray(int n, const double* Array)
{
	assert(n==9);
	for(int i=0;i<9;i++)
		array[i]=Array[i];
}


void Mat3D::SetArray(const MyVec& Array)
{
	assert(Array.GetSize()==9);
	for(int i=0;i<9;i++)
		array[i]=Array[i];
}


void Mat3D::GetArray(int n, double* Array) const
{
	assert(n==9);
	for(int i=0;i<9;i++)
		Array[i]=array[i];
}


MyVec Mat3D::GetArray() const
{
	MyVec temp(9);
	for(int i=0;i<9;i++)
		temp[i]=array[i];
	return temp;
}


double Mat3D::Det() const
{
	double Temp1=array[0]*(array[4]*array[8]-array[5]*array[7]);
	double Temp2=array[3]*(array[2]*array[7]-array[1]*array[8]);
	double Temp3=array[6]*(array[1]*array[5]-array[2]*array[4]);
	return Temp1+Temp2+Temp3;
}


double Mat3D::trace() const
{
	return array[0] + array[4] + array[8];
}


Mat3D Mat3D::Inver() const
{
	double detvalue=this->Det();
	if(fabs(detvalue)<1.0E-14)
		std::cout<<"矩阵奇异,不能求逆"<<std::endl;
	assert(detvalue!=0);
	double na[9]=
	{
		array[4]*array[8]-array[5]*array[7], array[2]*array[7]-array[1]*array[8], array[1]*array[5]-array[2]*array[4],
		array[5]*array[6]-array[3]*array[8], array[0]*array[8]-array[2]*array[6], array[2]*array[3]-array[0]*array[5],
		array[3]*array[7]-array[4]*array[6], array[1]*array[6]-array[0]*array[7], array[0]*array[4]-array[1]*array[3]
	};
	Mat3D temp;
	temp.SetArray(9, na);
	return temp/detvalue;
}


Mat3D Mat3D::Trans() const
{
	Mat3D Tran;
	for(int i=1;i<=3;i++)
		for(int j=1;j<=3;j++)
			Tran(i,j)=(*this)(j,i);
	return Tran;
}

const Mat3D	Mat3D::operator +() const
{
	return *this;
}
	
const Mat3D& Mat3D::operator+=(const Mat3D& Mat)
{
	for(int i=0;i<9;i++)
		array[i]+=Mat.array[i];
	return *this;
}

const Mat3D	Mat3D::operator +(const Mat3D& Mat) const
{
	Mat3D temp=*this;
	temp+=Mat;
	return temp;
}

const Mat3D& Mat3D::operator+=(const double& Num)
{
	for(int i=0;i<9;i++)
		array[i]+=Num;
	return *this;
}

const Mat3D	Mat3D::operator +(const double& Num) const
{
	Mat3D temp=*this;
	temp+=Num;
	return temp;
}

Mat3D operator+(const double& Num, const Mat3D& Mat)
{
	Mat3D temp;
	for(int i=0;i<9;i++)
		temp.array[i]=Num+Mat.array[i];
	return temp;
}

const Mat3D	Mat3D::operator -() const
{
	Mat3D Mat=*this;
	for(int i=0;i<9;i++)
		Mat.array[i]=-array[i];
	return Mat;
}

const Mat3D& Mat3D::operator-=(const Mat3D& Mat)
{
	for(int i=0;i<9;i++)
		array[i]-=Mat.array[i];
	return *this;
}

const Mat3D	Mat3D::operator -(const Mat3D& Mat) const
{
	Mat3D temp=*this;
	temp-=Mat;
	return temp;
}

const Mat3D& Mat3D::operator-=(const double& Num)
{
	for(int i=0;i<9;i++)
		array[i]-=Num;
	return *this;
}

const Mat3D	Mat3D::operator -(const double& Num) const
{
	Mat3D temp=*this;
	temp-=Num;
	return temp;
}

Mat3D operator-(const double& Num, const Mat3D& Mat)
{
	Mat3D temp;
	for(int i=0;i<9;i++)
		temp.array[i]=Num-Mat.array[i];
	return temp;
}

const Mat3D& Mat3D::operator*=(const Mat3D& Mat)
{
	Mat3D temp=*this;
	for(int i=1;i<=3;i++)
		for(int j=1;j<=3;j++)
			(*this)(i,j)=temp(i,1)*Mat(1,j)+temp(i,2)*Mat(2,j)+temp(i,3)*Mat(3,j);
	return *this;
}

const Mat3D	Mat3D::operator*(const Mat3D& Mat) const
{
	Mat3D temp=*this;
	temp*=Mat;
	return temp;
}

const MyVec Mat3D::operator*(const MyVec& Vec) const
{
	MyVec temp;
	for(int i=1;i<=3;i++)
		for(int j=1;j<=3;j++)
		temp(i)+=(*this)(i,j)*Vec(j);
	return temp;
}

const Mat3D& Mat3D::operator*=(const double& Num)
{
	for(int i=0;i<9;i++)
		array[i]*=Num;
	return *this;
}
const Mat3D	Mat3D::operator*(const double& Num) const
{
	Mat3D temp=*this;
	temp*=Num;
	return temp;
}

MyVec operator*(const MyVec& Vec, const Mat3D& Mat)
{
	MyVec temp;
	for(int i=1;i<=3;i++)
		for(int j=1;j<=3;j++)
		temp(i)+=Mat(j,i)*Vec(j);
	return temp;
}

Mat3D operator*(const double& Num, const Mat3D& Mat)
{
	return Mat*Num;
}

const Mat3D& Mat3D::operator/=(const double& Num)
{
	assert(Num!=0);
	for(int i=0;i<9;i++)
		array[i]/=Num;
	return *this;
}


const Mat3D	Mat3D::operator/(const double& Num) const
{
	Mat3D temp=*this;
	temp/=Num;
	return temp;
}

const Mat3D Mat3D::Absol() const
{
	Mat3D temp=*this;
	for(int i=0;i<9;i++)
		temp.array[i]=fabs(temp.array[i]);
	return temp;
}

const double Mat3D::Max() const
{
	double temp=array[0];
	for(int i=1;i<9;i++)
		temp=(temp<array[i])?array[i]:temp;
	return temp;
}

const double Mat3D::Min() const
{
	double temp=array[0];
	for(int i=1;i<9;i++)
		temp=(temp>array[i])?array[i]:temp;
	return temp;
}

std::istream & operator>>(std::istream &input, Mat3D& Mat)
{
  for (int i=1;i<=3;i++)
		for(int j=1;j<=3;j++)
			input>>Mat(i,j);
    return input;
}


std::ostream & operator<<(std::ostream &output, Mat3D& Mat)
{
    for (int i=1;i<=3;i++)
	{
		for(int j=1;j<=3;j++)
			output<<Mat(i,j)<<"  ";
		output<<std::endl;
	}
    return output;
}
