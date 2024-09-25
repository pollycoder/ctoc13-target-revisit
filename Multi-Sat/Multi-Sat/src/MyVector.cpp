#include"MyVector.h"
//#include "StdAfx.h"

//���캯��

MyVec::MyVec(int arraysize)
{
	assert(arraysize>=1);
	size=arraysize;
	array=new double[size];
	assert(array!=NULL);
	for(int i=0;i<size;i++)
		array[i]=0;
}

//�������캯��

MyVec::MyVec(const MyVec& TheVector)
{

	int n=TheVector.size;
	size=n;
	array=new double[n];
	assert(array!=NULL);
	double* srcptr=TheVector.array;
	double* destptr=array;
	while(n--)
		*destptr++=*srcptr++;
}

//��άʸ���������캯��

MyVec::MyVec(const double& x, const double& y, const double& z)
{
	size=3;
	array=new double[size];
	assert(array!=NULL);
	array[0]=x;
	array[1]=y;
	array[2]=z;
}

//��άʸ���������캯��

MyVec::MyVec(const double& x, const double& y, const double& z, const double& vx, const double& vy, const double& vz)
{
	size=6;
	array=new double[size];
	assert(array!=NULL);
	array[0]=x;
	array[1]=y;
	array[2]=z;
	array[3]=vx;
	array[4]=vy;
	array[5]=vz;
}

MyVec::MyVec(const double& x, const double& y, const double& z, const double& vx, const double& vy, 
			 const double& vz, const double& m)
{
	size=7;
	array=new double[size];
	assert(array!=NULL);
	array[0]=x;
	array[1]=y;
	array[2]=z;
	array[3]=vx;
	array[4]=vy;
	array[5]=vz;
	array[6]=m;
}

//��������
MyVec::~MyVec()
{
	delete[] array;
}

//��ȡʸ����ά��

 int MyVec::GetSize() const
{
	return size;
}

//��һά�����ʸ����ֵ

void MyVec::SetValue(int n, const double* Value)
{
/*	assert(size==n);
	for(int i=0;i<size;i++)
		array[i]=Value[i];*/
	delete[] array;
	array=new double[n];
	size=n;
	assert(array!=NULL);
	double* destptr=array;
	while(n--)
		*destptr++=*Value++;
}

//��ʸ��ֵ����һά����

void MyVec::GetValue(int n, double* Value) const
{
	assert(size==n);
	for(int i=0;i<size;i++)
		Value[i]=array[i];
}

//�±���������غ���

 double& MyVec::operator[](int i)
{
	assert((i>=0)&&(i<size));
	return array[i];
}

//�±���������غ���

 const double& MyVec::operator[](int i) const
{
	assert((i>=0) && (i<size));
	return array[i];
}

//�±���������غ���

 double& MyVec::operator()(int i)
{
	assert((i>=1) && (i<=size));
	return array[i-1];
}

//�±���������غ���

 const double& MyVec::operator()(int i) const
{
	assert((i>=1) && (i<=size));
	return array[i-1];
}

//==��������غ���

 bool MyVec::operator==(const MyVec& Vec) const
{
	assert(Vec.size==size);
	for(int i=0;i<size;i++)
	{
		if (array[i]!=Vec.array[i])
			return false;
	}
	return true;
}

//!=��������غ���

 bool MyVec::operator !=(const MyVec& Vec) const
{
	return(!(*this==Vec));
}

//��ֵ��������غ���

MyVec& MyVec::operator=(const MyVec& Vec)
{
	delete[] array;
	size=Vec.size;
	array=new double[size];
	for(int i=0;i<size;i++)
		array[i]=Vec.array[i];
	return *this;
}

//ȡ����������غ���

 const MyVec MyVec::operator+() const
{
	return *this;
}

//ȡ����������غ���

 const MyVec MyVec::operator-() const
{
	MyVec Vec=*this;
	for(int i=0;i<size;i++)
		Vec.array[i]=-array[i];
	return Vec;
}

//ʸ��+=��������غ���

 const MyVec& MyVec::operator+=(const MyVec& Vec)
{
	assert(Vec.size==size);
	for(int i=0;i<size;i++)
		array[i]+=Vec.array[i];
	return *this;
}

//ʸ���ӷ���������غ���

 const MyVec MyVec::operator+(const MyVec& Vec) const
{
	assert(size==Vec.size);
	MyVec temp=*this;
	temp+=Vec;
	return temp;
}

//ʸ�������+=��������غ���

 const MyVec& MyVec::operator+=(const double& Num)
{
	for(int i=0;i<size;i++)
		array[i]+=Num;
	return *this;
}

//ʸ�������+��������غ���

 const MyVec MyVec::operator+(const double& Num) const
{
	MyVec temp=*this;
	temp+=Num;
	return temp;
}

//������ʸ��+��������غ���

 MyVec operator+(const double& Num, const MyVec& Vec)
{
	MyVec temp(Vec.size);
	for(int i=0;i<Vec.size;i++)
		temp.array[i]=Num+Vec.array[i];
	return temp;
}

//ʸ��-=��������غ���

 const MyVec& MyVec::operator-=(const MyVec& Vec)
{
	assert(Vec.size==size);
	for(int i=0;i<size;i++)
		array[i]-=Vec.array[i];
	return *this;
}

//ʸ��-��������غ���	

 const MyVec MyVec::operator-(const MyVec& Vec) const
{
	assert(size==Vec.size);
	MyVec temp=*this;
	temp-=Vec;
	return temp;
}	

//ʸ�������-=��������غ���

 const MyVec& MyVec::operator-=(const double& Num)
{
	for(int i=0;i<size;i++)
		array[i]-=Num;
	return *this;
}

//ʸ�������-��������غ���

 const MyVec MyVec::operator-(const double& Num) const
{
	MyVec temp=*this;
	temp-=Num;
	return temp;
}

//������ʸ��-��������غ���

 MyVec operator-(const double& Num, const MyVec& Vec)
{
	MyVec temp(Vec.size);
	for(int i=0;i<Vec.size;i++)
		temp.array[i]=Num-Vec.array[i];
	return temp;
}

//ʸ��*=��������غ���

 const MyVec& MyVec::operator*=(const MyVec& Vec)
{
	assert(size==Vec.size);
	MyVec temp=*this;
	for(int i=0;i<size;i++)
		array[i]=temp.array[i]*Vec.array[i];
	return *this;
}	

//ʸ��*��������غ���

 const MyVec MyVec::operator*(const MyVec& Vec) const
{
	assert(size==Vec.size);
	MyVec temp=*this;
	temp*=Vec;
	return temp;
}

//ʸ�������*��������غ���

 const MyVec MyVec::operator*(const double& Num) const
{
	MyVec temp=*this;
	temp*=Num;
	return temp;
}

//ʸ�������*=��������غ���	

 const MyVec& MyVec::operator*=(const double& Num)
{
	for(int i=0;i<size;i++)
		array[i]*=Num;
	return *this;
}

//������ʸ��*=��������غ���

 MyVec operator*(const double& Num, const MyVec& Vec)
{
	MyVec temp(Vec.size);
	for(int i=0;i<Vec.size;i++)
		temp.array[i]=Num*Vec.array[i];
	return temp;
}

//ʸ��/=��������غ���

 const MyVec& MyVec::operator/=(const MyVec& Vec)
{
	MyVec absvec=Vec.Absol();
	double absmin=absvec.Min();
	assert((size==Vec.size)&&(absmin>1.0E-14));
	MyVec temp=*this;
	for(int i=0;i<size;i++)
		array[i]=temp.array[i]/Vec.array[i];
	return *this;
}	

//ʸ��/��������غ���

 const MyVec MyVec::operator/(const MyVec& Vec) const
{
	assert(size==Vec.size);
	MyVec temp=*this;
	temp/=Vec;
	return temp;
}

//ʸ�������/��������غ���

 const MyVec MyVec::operator/(const double& Num) const
{
	MyVec temp=*this;
	temp/=Num;
	return temp;
}

//ʸ�������/=��������غ���

 const MyVec& MyVec::operator/=(const double& Num)
{
	assert(Num!=0);
	for(int i=0;i<size;i++)
		array[i]/=Num;
	return *this;
}

//������ʸ��/=��������غ���

 MyVec operator/(const double& Num, const MyVec& Vec)
{
	MyVec absvec=Vec.Absol();
	double absmin=absvec.Min();
	assert(absmin!=0);
	MyVec temp(Vec.size);
	for(int i=0;i<Vec.size;i++)
		temp.array[i]=Num/Vec.array[i];
	return temp;
}

//ʸ���ڻ�

 const double MyVec::Dot(const MyVec& Vec) const
{
	assert(size==Vec.size);
	double temp=0;
	for(int i=0;i<size;i++)
		temp+=array[i]*Vec.array[i];
	return temp;
}

//��άʸ�����

 const MyVec MyVec::Cross(const MyVec& Vec) const
{
	assert((size==3)&&(3==Vec.size));
	MyVec temp(3);
	temp[0]=array[1]*Vec.array[2]-array[2]*Vec.array[1];
	temp[1]=array[2]*Vec.array[0]-array[0]*Vec.array[2];
	temp[2]=array[0]*Vec.array[1]-array[1]*Vec.array[0];
	return temp;
}

//ʸ������ֵ

 const MyVec MyVec::Absol() const
{
	MyVec temp=*this;
	for(int i=0;i<size;i++)
		temp.array[i]=(temp.array[i]<0)?(-temp.array[i]):temp.array[i];
	return temp;
}

//ʸ��1����

 const double MyVec::Norm1() const
{
	double value=0;
	MyVec temp=this->Absol();
	for(int i=0;i<size;i++)
		value+=temp[i];
	return value;
}

//ʸ�������

 const double MyVec::NormInf() const
{
	MyVec temp=this->Absol();
	return temp.Max();
}

//ʸ��2����,ģ

 const double MyVec::Norm2() const
{
	double temp=Dot(*this);
	return sqrt(temp);
}

//ʸ���������ֵ

 const double MyVec::Max() const
{
	double temp=array[0];
	for(int i=0;i<size;i++)
		temp=(temp<array[i])?array[i]:temp;
	return temp;
}

//ʸ��������Сֵ

 const double MyVec::Min() const
{
	double temp=array[0];
	for(int i=0;i<size;i++)
		temp=(temp>array[i])?array[i]:temp;
	return temp;
}

//������>>��������غ���

std::istream & operator>>(std::istream &input, MyVec& Vec)
{
    int N=Vec.GetSize();
    for (int i=0;i<N;i++)
		input>>Vec[i];
    return input;
}

//�����<<��������غ���

std::ostream & operator<<(std::ostream &output, MyVec& Vec)
{
    int N=Vec.GetSize();
//	output<<endl;
    for (int i=0;i<N;i++)
		output<<Vec[i]<<std::endl;
    return output;
}
