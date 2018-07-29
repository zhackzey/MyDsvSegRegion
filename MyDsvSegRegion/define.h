#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <io.h>
#include <memory.h>
#include <conio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <opencv2/opencv.hpp>
/*
#if _DEBUG
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_highgui249d.lib")
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_imgproc249d.lib")
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_core249d.lib")
#else
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_highgui249.lib")
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_imgproc249.lib")
#pragma comment(lib, "C:/opencv-2.4.9/opencv/build/x64/vc10/lib/opencv_core249.lib")
#endif
*/
#pragma comment(lib, "D:/opencv3/opencv/build/x64/vc14/lib/opencv_world341d.lib")


using namespace std;
using namespace cv;

const int maxn = 2000000000;
const double topi = acos(-1.0) / 180.0;	// pi/180
#define BOUND(x,min,max) ((x) < (min) ? (min) : ((x) > (max) ? (max) : (x)))
#define	nint(x)			(int)((x>0)?(x+0.5):(x-0.5))
#define	sqr(x)			((x)*(x))

struct point2d
{
	double x;
	double y;
};

struct point3d
{
	double x;
	double y;
	double z;
};

typedef double  MATRIX[3][3];

/*
 * 3D 的车辆轨迹信息
 * ang 角度
 * shv 位移
 * rot 用于坐标系变换
 */
typedef struct {
	point3d			ang;
	point3d			shv;
	MATRIX			rot;
} TRANSINFO;

typedef double  MAT2D[2][2];

/*
 * 2D 的车辆轨迹信息
 * ang 角度
 * rot 用于坐标变换
 */
typedef struct {
	double			ang;
	point2d			shv;
	MAT2D			rot;
} TRANS2D;

#define	PNTS_PER_LINE		32
#define	LINES_PER_BLK		12
#define	PTNUM_PER_BLK		(32*12)
#define	BKNUM_PER_FRM		580
#define	SCANDATASIZE		(BKNUM_PER_FRM*LINES_PER_BLK/2)

//for vel64
//HORIERRFACTOR=tan（水平角分辨率=0.1度）*（放大系数=2.0）=0.0018*5
#define	HORIERRFACTOR	0.02	//0.006
//VERTERRFACTOR=tan（垂直角分辨率=0.38度）*（放大系数=1.5）=0.0067*5
#define	VERTERRFACTOR	0.035	//0.035
#define	BASEERROR		0.3
#define	MAXSMOOTHERR	1.0
#define	MAXDISTHRE		2.0


#define	M_PI		3.1415926536

#define	INVALIDDOUBLE		99999999.9


typedef struct {
	float			x, y, z;
	u_char			i;		//反射率
} point3fi;

typedef struct {
	int x, y;
} point2i;

/*
 * 一个VDN数据里有32*12个激光点
 * 每个激光点包含的三维坐标x,y,z和一个i
 */
typedef struct {
	long			millisec;
	point3fi		points[PTNUM_PER_BLK];
} ONEVDNDATA;

/*
 * 一个DSV数据结构体里有32*12个激光点
 * 其实是在VDN数据的基础上加入了车辆位置信息
 * ONEDSVDATA = TRANSINFO + ONEVDNDATA
 */
typedef struct {
	point3d			ang;
	point3d			shv;
	long			millisec;
	point3fi		points[PTNUM_PER_BLK];
	MATRIX			rot;
} ONEDSVDATA;

/*
 * 一帧DSV数据 包含580 个DSV数据块
 * 直接是一个ONEDSVDATA 数据块的数组
 */
typedef struct {
	ONEDSVDATA		dsv[BKNUM_PER_FRM];
} ONEDSVFRAME;

typedef	struct {
	unsigned short	lab;
	point2i		dmin;
	point2i		dmax;
	point3fi	maxxp, maxyp, maxzp;
	point3fi	minxp, minyp, minzp;
	point3d		cp;
	int			ptnum;
	point3d		norm;
	double		var;
} SEGBUF;

//vel64
#define	VMINANG		(-21.627*M_PI/180.0)
#define	VMAXANG		(2.432*M_PI/180.0)

typedef struct {
	int				wid;	//距离图像的长，实际长LINES_PER_BLK*BKNUM_PER_FRM/2;	
							//为减少计算量，实际长/2
	int				len;	//距离图像的宽，64
	double			h0;		//水平零位角
	double			v0;		//垂直零位角
	double			hres;	//水平角解像度
	double			vres;	//垂直角解像度
	point3fi		*pts;	//每个像素对应的激光点（CorrectPoints后）
	point2i			*idx;	//每个像素对应的激光点扫描序号
	BYTE			*di;		//for data alignment only
	int				*regionID;	//每个像素点对应到区域编号
	int				regnum;		//总区域数
	SEGBUF			*segbuf;	//每个区域数据块
	IplImage		*rMap;		//可视化距离图像用
	IplImage		*lMap;		//可视化分割结果用
} RMAP;

#define	WIDSIZ		25 //120.0		//DEM宽、单位：米
#define	LENSIZ		25 //200.0		//DEM长、单位：米
#define	PIXSIZ		0.25
#define	POSOBSMINHEIGHT	0.6		//0.6m
#define	VEHICLEHEIGHT	3.0		//3.0m
#define	NEARVEHICLEDIS	6.0		//5.0m


typedef struct {
	int			x0, x1;		//DEM中的地面点开始和结束像素序号[0,dm.wid)
							//	int			y;			//DEM中的像素序号[0,dm.len)
	double		h;			//中心点位置（(x0+x1)/2))的地面高度
	double		dl;			//纵向距离y处的地面扫描线与前一条扫描线间的正常水平距离
							//（前一条扫描线为与车体更近的那条，两扫描线间角度d_ang=(VMAXANG-VMINANG)/63)
} CENTERLN;
/*以下部分处理新加入的features*/
#define MAX_PTS_PER_GRID 500
typedef struct {
	int				wid;			//DEM宽、像素数
	int				len;			//DEM长、像素数
	double			*demg;			//ground
	int				*demgnum;
	double			*demhmin;		//非地面DEM最低高度
	double			*demhmax;		//非地面DEM最高高度
	int				*demhnum;		//落入该栅格中的非地面点数目
	BYTE			*lab;			//每个像素点对应的标签、包括TRAVESABLE、NONTRAVESABLE
	double			*groll;			//该像素邻域平面拟合的横滚角
	double			*gpitch;		//该像素邻域平面拟合的俯仰角
	BYTE			*sublab;		//路面附近像素的对应的子标签
									//TRAVESABLE包括：FLATGROUND、DOWNSLOPE、UPSLOPE、
									//LEFTSIDESLOPE、RIGHTSIDESLOPE、EDGEPOINTS
									// NONTRAVESABLE包括：POSSIOBSTA、NEGATOBSTA

	double			*lpr;			//probability of the lab 为该属性的概率
	double			*WX, *WY, *WZ;	//用于地面拟合的临时buffer
	CENTERLN		*centerln;		//道路中心线信息、辅助用，不准确
	IplImage		*lmap;			//可视化用
	IplImage		*smap;			//可视化用
	TRANS2D			trans;			//该DEM的车体位姿，对应于每一帧第0个block
	bool			dataon;

	/*
	 * 一下部分是新加入的feature
	 */

	double			*meanHeight;	// 每一个栅格的平均高度
	
	double			**ptsHeight;	// 用于存储打在每一个栅格中的激光点的高度序列
	double			*heightVariance;// 每一个栅格的方差
	int				*demnum;		// 每一个栅格的激光点数
	bool			*visible;		// 可见性（待定）
} DMAP;

#define UNKNOWN			0
#define NONVALID		-9999
#define EDGEPT			-9

#define TRAVESABLE		1
#define NONTRAVESABLE	2
#define POSSIOBSTA		3
#define	NEGATOBSTA		4
#define HANGDOWNTR		5
#define HANGDOWNUN		6
#define	FLATGROUND		10
#define DOWNSLOPE		11
#define UPSLOPE			12
#define	LEFTSIDESLOPE	13
#define	RIGHTSIDESLOPE	14
#define	EDGEPOINTS		15

#define	istravesable(x)			(x==TRAVESABLE||x==DOWNSLOPE||x==UPSLOPE||x==SIDESLOPE)

extern RMAP	rm;
extern TRANSINFO calibInfo;
extern ONEDSVFRAME	*onefrm;

void rMatrixInit(MATRIX &rt);
void rMatrixmulti(MATRIX &r, MATRIX &rt);
void createRotMatrix_ZYX(MATRIX &rt, double rotateX, double rotateY, double rotateZ);
void createRotMatrix_XYZ(MATRIX &rt, double rotateX, double rotateY, double rotateZ);
void createRotMatrix_ZXY(MATRIX &rt, double rotateX, double rotateY, double rotateZ);
void shiftPoint3d(point3d &pt, point3d &sh);
void rotatePoint3d(point3d &pt, MATRIX &a);
double normalize2d(point2d *p);
double ppDistance2d(point2d *p1, point2d *p2);
double innerProduct2d(point2d *v1, point2d *v2);
double ppDistance3fi(point3fi *pt1, point3fi *pt2);
double p2r(point3fi *pt1);
void rotatePoint3fi(point3fi &pt, MATRIX &a);

BOOL ContourSegger();
void SmoothingData();
void Region2Seg();
void EstimateSeg();
void ContourExtraction();
UINT RegionGrow();
void EdgeGrow();
void ClassiSeg();
void OutputLog(char *filename, char *str);

void INVshiftPoint3d(point3d &pt, point3d &sh);
void INVrotatePoint3d(point3d &pt, MATRIX &a);
void shiftPoint3fi(point3fi &pt, point3d &sh);
void rotatePoint3fi(point3fi &pt, MATRIX &a);

void Calculate_Plane(int Points_Total, double *X_Coord, double *Y_Coord, double *Z_Coord,
	int Origin_Flag, double Plane_Eq[4]);
void Calculate_Residuals(double *X, double *Y, double *Z, double Equation[4],
	double *Error, int PointsTotal);
void shiftPoint2d(point2d &pt, point2d &sh);
void rotatePoint2d(point2d &pt, MAT2D &a);

void DrawRangeView();
void GenerateRangeView();
void InitRmap(RMAP *rm);
void ReleaseRmap(RMAP *rm);

void DrawDem(DMAP &m);
void CopyGloDem(DMAP *tar, DMAP *src);
void ZeroGloDem(DMAP *m);
void InitDmap(DMAP *dm);
void ReleaseDmap(DMAP *dm);
void PredictGloDem(DMAP &gmtar, DMAP &gmtmp);
void UpdateGloDem(DMAP &glo, DMAP &loc);
void GenerateLocDem(DMAP &loc);
void CallbackLocDem(int event, int x, int y, int flags, void *ustc);
void LabelRoadSurface(DMAP &glo);
void LabelObstacle(DMAP &glo);
void ExtractRoadCenterline(DMAP &glo);
/*以下部分处理新加入的features*/
void SaveDEM(DMAP & dm, int a, int b);