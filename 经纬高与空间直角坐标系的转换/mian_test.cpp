#include<iostream>
#include"Deployment.h"
using namespace std;

int main()
{
	//原点经纬高
	double Lg0 = 124;
	double Bg0 = 28;
	double hg0 = 0;
	double theta_azg0 = 0;    

	Deployment deployment(Lg0, Bg0, hg0, theta_azg0);


	//经纬高    --->   xyz
	cout << "xyz  ---->  经纬高" << endl;
	const int n = 4;
	double LBH0[n][3] = { {124,28,0},{125,28,0}, {123,28,0}, {124,27,0} };   //测试输入数据    顺序是 经纬高
	double Txyz[n][3];       //测试输出数据    顺序是 xyz    z为高度

	for (int i = 0; i < n; i++) {
		//转换 
		
		deployment.GetXYZFromLBH(Txyz[i][0], Txyz[i][1], Txyz[i][2]/*高度*/, LBH0[i][0] * d2r/*单位:弧度*/, LBH0[i][1] * d2r/*单位：弧度*/, LBH0[i][2]);


		//输出
		cout <<"LBH: "<< LBH0[i][0] << "     " << LBH0[i][1] << "     " << LBH0[i][2] << "  " << endl;
		cout << "xyz: " << Txyz[i][0] << "  "<<Txyz[i][1] << "  "<<Txyz[i][2] << "  " << endl;
		cout << endl;
	}





	//xyz  ---->  经纬高
	cout << "xyz  ---->  经纬高" << endl;
	double xyz0[n][3] = { {0,0,0}, {1500,0,2400},{10000,0,2400},{10000,10000,2400} };   //测试输入数据    顺序是 xyz    z为高度
	double TLBH[n][3];      //测试输出数据    顺序是  经纬高

	for (int i = 0; i < n; i++) {
		//转换 
		deployment.GetLBHFromXYZ(TLBH[i][0], TLBH[i][1], TLBH[i][2], xyz0[i][0], xyz0[i][1], xyz0[i][2]/*高度*/);
		TLBH[i][0] *= r2d;   //转换出来的数据为弧度      度=弧度*r2d
		TLBH[i][1] *= r2d;


		//输出
		cout << "xyz: " << xyz0[i][0] << "     " << xyz0[i][1] << "     " << xyz0[i][2] << "  " << endl;
		cout << "LBH: " << TLBH[i][0] << "  " << TLBH[i][1] << "  " << TLBH[i][2] << "  " << endl;
		cout << endl;
	}



	//system("pause");
	return 0;
}