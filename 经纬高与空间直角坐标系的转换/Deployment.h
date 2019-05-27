// Deployment.h

#ifndef DEPLOYMENT_H
#define DEPLOYMENT_H

// ==================================================
// 常用变量

// 圆周率
static const double PI = 3.141592653589793;
// 将角度转为弧度
static const double d2r = PI / 180;
// 将弧度转为角度
static const double r2d = 1 / d2r;
// 自然对数的底
static const double e = 2.71828182845905;


class Deployment
{
	// 地心赤道旋转坐标系 Se: ze 轴指向北极, xe 轴沿赤道平面与 Greenwich 子午面的相交线, ye 轴由右手法则确定.
	// 当地铅垂坐标系 Sv: zv 轴当地铅垂向下, xv 与 yv 在当地水平面内, xv 指向北, yv 指向东.
	// 物体测量坐标系 Sg: 由 Sv 绕 zv 轴旋转角度 theta_azg 使 xv 轴指向物体测量轴线即得 Sg.
public:
// Constructor
	// Lg: 物体部署点经度 [deg]
	// Bg: 物体部署点大地纬度 [deg]
	// hg: 物体部署点高程 [m]
	// theta_azg: 物体测量坐标系 x 轴与当地铅垂坐标系 x 轴的夹角, 由当地铅垂坐标系 x 轴转向物体坐标系 x 轴为正 [deg]
	Deployment(double Lg = 116.38, double Bg = 39.9, double hg = 0, double theta_azg = 0);
	~Deployment();

// Attributes
	// 设置地球椭球模型长半轴和偏心率
	// RE: 地球椭球模型长半轴 [m]
	// f: 地球椭球模型偏心率
	void SetREandf(double RE, double f);

	// 获取地球椭球模型长半轴和偏心率
	// RE: 地球椭球模型长半轴 [m]
	// f: 地球椭球模型偏心率
	void GetREandf(double &RE, double &f) const;

private:
	double m_Lg;	// 物体部署点经度 [rad]
	double m_Bg;	// 物体部署点大地纬度 [rad]
	double m_hg;	// 物体部署点高程 [m]
	double m_theta_azg;	// 物体测量坐标系 x 轴与当地铅垂坐标系 x 轴的夹角, 由当地铅垂坐标系 x 轴转向物体测量坐标系 x 轴为正 [rad]
	double m_RE;	// 地球椭球模型长半轴 [m]
	double m_f;		// 地球椭球模型偏心率

	// 坐标转换矩阵 Cge 中的各元素
	double m_c11;
	double m_c12;
	double m_c13;
	double m_c21;
	double m_c22;
	double m_c23;
	double m_c31;
	double m_c32;
	double m_c33;

	// 地心指向物体测量坐标系原点的矢径在地心直角坐标系中的坐标
	double m_xRoge;
	double m_yRoge;
	double m_zRoge;

// Implementation
public:
	// 根据矢径端点的经度,大地纬度和高程计算此矢径在地心坐标系中的坐标向量
	// vector_x, vector_y, vector_z: 矢径在地心坐标系中的坐标向量 [m]
	// L: 经度[rad]
	// B: 大地纬度[rad]
	// h: 高程 [m]
	void GetEcFromLBH(double &vector_x, double &vector_y, double &vector_z, double L, double B, double h) const;

	// 根据矢径端点在地心坐标系中的坐标向量计算此矢径的高程,经度和大地纬度
	// L: 经度[rad]
	// B: 大地纬度[rad]
	// h: 高程 [m]
	// vector_x, vector_y, vector_z: 矢径在地心坐标系中的坐标向量 [m]
	void GetLBHFromEc(double &L, double &B, double &h, double vector_x, double vector_y, double vector_z) const;

	// 根据矢径端点的经度, 大地纬度和高程计算矢径端点在物体测量坐标系中的坐标
	// xrg, yrg, zrg: 矢径端点在物体测量坐标系中的坐标向量 [m]
	// L: 矢径端点的经度[rad]
	// B: 矢径端点的大地纬度[rad]
	// h: 矢径端点的高程 [m]
	void GetXYZFromLBH(double &xrg, double &yrg, double &zrg, double L, double B, double h) const;

	// 根据矢径端点在物体测量坐标系中的坐标计算矢径端点的经度, 大地纬度和高程
	// L: 矢径端点的经度[rad]
	// B: 矢径端点的大地纬度[rad]
	// h: 矢径端点的高程 [m]
	// xrg, yrg, zrg: 矢径端点在物体测量坐标系中的坐标向量 [m]
	void GetLBHFromXYZ(double &L, double &B, double &h, double xrg, double yrg, double zrg) const;

	// 将矢量在物体测量坐标系中的坐标分量转换为矢量在地心赤道旋转坐标系中的坐标分量
	// xrhoe, yrhoe, zrhoe: 矢量在地心赤道旋转坐标系中的坐标分量
	// xrhog, yrhog, zrhog: 矢量在物体测量坐标系中的坐标分量
	void GetEcVFromGcV(double &xrhoe, double &yrhoe, double &zrhoe,	double xrhog, double yrhog, double zrhog) const;

	// 将矢量在地心赤道旋转坐标系中的坐标分量转换为矢量在物体测量坐标系中的坐标分量
	// xrhog, yrhog, zrhog: 矢量在物体测量坐标系中的坐标分量
	// xrhoe, yrhoe, zrhoe: 矢量在地心赤道旋转坐标系中的坐标分量
	void GetGcVFromEcV(double &xrhog, double &yrhog, double &zrhog,	double xrhoe, double yrhoe, double zrhoe) const;

	// 将矢量在当地铅垂坐标系中的坐标分量转换为矢量在地心赤道旋转坐标系中的坐标分量
	// xrhoe, yrhoe, zrhoe: 矢量在地心赤道旋转坐标系中的坐标分量
	// xlocal, ylocal, zlocal: 矢量在当地铅垂坐标系中的坐标分量
	// L: 当地铅垂坐标系原点经度[rad]
	// B: 当地铅垂坐标系原点大地纬度[rad]
	void GetEcVFromLcV(double &xrhoe, double &yrhoe, double &zrhoe,
		double xlocal, double ylocal, double zlocal, double L, double B) const;

	// 将矢量在地心赤道旋转坐标系中的坐标分量转换为矢量在当地铅垂坐标系中的坐标分量
	// xlocal, ylocal, zlocal: 矢量在当地铅垂坐标系中的坐标分量
	// xrhoe, yrhoe, zrhoe: 矢量在地心赤道旋转坐标系中的坐标分量
	// L: 当地铅垂坐标系原点经度[rad]
	// B: 当地铅垂坐标系原点大地纬度[rad]
	void GetLcVFromEcV(double &xlocal, double &ylocal, double &zlocal,
		double xrhoe, double yrhoe, double zrhoe, double L, double B) const;

	// 将矢量在当地铅垂坐标系中的坐标分量转换为矢量在物体测量坐标系中的坐标分量
	// xrhog, yrhog, zrhog: 矢量在物体测量坐标系中的坐标分量
	// xlocal, ylocal, zlocal: 矢量在当地铅垂坐标系中的坐标分量
	// L: 当地铅垂坐标系原点经度[rad]
	// B: 当地铅垂坐标系原点大地纬度[rad]
	void GetGcVFromLcV(double &xrhog, double &yrhog, double &zrhog,
		double xlocal, double ylocal, double zlocal, double L, double B) const;

	// 将矢量在物体测量坐标系中的坐标分量转换为矢量在当地铅垂坐标系中的坐标分量
	// xlocal, ylocal, zlocal: 矢量在当地铅垂坐标系中的坐标分量
	// xrhog, yrhog, zrhog: 矢量在物体测量坐标系中的坐标分量
	// L: 当地铅垂坐标系原点经度[rad]
	// B: 当地铅垂坐标系原点大地纬度[rad]
	void GetLcVFromGcV(double &xlocal, double &ylocal, double &zlocal,
		double xrhog, double yrhog, double zrhog, double L, double B) const;
};

#endif