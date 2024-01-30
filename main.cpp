#include <graphics.h>		// 引用 EasyX 绘图库头文件
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <map>
#include <vector>
#include <sstream>
#include <random>
#include <numbers>
#include <list>
#include <chrono>

using namespace std;
//定义全局的向量，方便后续坐标变换的计算
struct vec
{
	double xyz[3];
	double eps = 1e-8;
	vec(double a = 0,double b = 0,double c = 0)
	{
		xyz[0] = a; xyz[1] = b; xyz[2] = c;
	}
	vec(double a[])
	{
		xyz[0] = a[0]; xyz[1] = a[1]; xyz[2] = a[2];
	}

	//默认拷贝构造函数

	//向量相等当且仅当各坐标相等
	bool operator == (const vec& b)
	{
		for (int i = 0; i < 3; i++)
		{
			if (fabs(xyz[i] - b.xyz[i]) > eps)
				return false;
		}
		return true;
	}

	vec operator +(const vec& b)
	{
		return vec(xyz[0] + b.xyz[0],xyz[1] + b.xyz[1],xyz[2] + b.xyz[2]);
	}

	vec operator -(const vec& b)
	{
		return vec(xyz[0] - b.xyz[0], xyz[1] - b.xyz[1], xyz[2] - b.xyz[2]);
	}

	//向量数乘c
	vec operator *(double c)const
	{
		return vec(xyz[0] * c, xyz[1] * c, xyz[2] * c);
	}

	//向量长度，2范数
	double norm2()const
	{
		return sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	}

	//向量叉积，直接用3阶行列式推出来即可。
	vec cross_product(const vec& b)const
	{
		return vec(xyz[1] * b.xyz[2] - xyz[2] * b.xyz[1], xyz[2] * b.xyz[0] - xyz[0] * b.xyz[2], xyz[0] * b.xyz[1] - xyz[1] * b.xyz[0]);
	}
	//向量点积
	double dot_product(const vec& b)const
	{
		double ans = 0;
		for (int i = 0; i < 3; i++)
		{
			ans += xyz[i] * b.xyz[i];
		}
		return ans;
	}
	//三阶行列式，等价于先做叉积，再做点积。
	static double determinant(const vec& a,const vec&b, const vec&c)
	{
		return (a.cross_product(b)).dot_product(c);
	}

	static vec get_zero_vector()
	{
		return vec(0, 0, 0);
	}
};

//vfromobj,来自obj文件的顶点v
/*struct vfromobj
{
	vec v;
	vfromobj(double x = 0, double y = 0, double z = 0)
	{
		v = vec(x,y,z);
	}
};

//vnfromobj,来自obj文件的顶点法向vn
struct vnfromobj
{
	vec vn;
	vnfromobj(double x = 0, double y = 0, double z = 0)
	{
		vn = vec(x, y, z);
	}
};
*/

//index pair,用于辅助存储facet中描述polygon顶点编号index的数据结构;
//facetfromobj,来自obj文件的facet，可能是polygon

struct indpair
{
	//v index, vt(顶点纹理坐标) index, vn index
	int vi, vti, vni;
	indpair(int a = 0,int b = 0,int c = 0)
	{
		vi = a, vti = b, vni = c;
	}
	//默认拷贝构造函数
};
struct facetfromobj
{
	//可能不止3个顶点
	vector<indpair> pointlist;
	//int colorind;
};

struct geometryfromobj
{
	vector<vec> vlist;
	vector<vec> vnlist;
	list<facetfromobj> flist;

	void read_from_file(string filename = ".\\soccerball.obj")
	{
		//map<int,int>mp;
		fstream fs;
		fs.open(filename, std::fstream::in);
		int linn = 0;
		if (fs.is_open())
		{
			//读取obj文件每一行
			while (!fs.eof())
			{
				linn++;
				string a;
				getline(fs, a);
				string b;
				int pos = a.find_first_of(' ');
				if (pos == string::npos)
					continue;
				b.assign(a.substr(0, pos));
				a = a.substr(pos + 1);

				//开始解析 a为指令，b为后续数据
				if (b == "v")
				{
					double x, y, z;
					if(sscanf_s(a.c_str(), "%lf%lf%lf", &x,&y,&z)!=3)
						break;
					vec v(x,y,z);
					vlist.push_back(v);
				}
				else if (b == "vn")
				{
					double x, y, z;
					if(sscanf_s(a.c_str(), "%lf%lf%lf", &x, &y, &z)!=3)
						break;
					vec vn(x, y, z);
					vnlist.push_back(vn);
				}
				else if (b == "f")
				{
					facetfromobj f;
					//此处未使用纹理坐标，所以只使用x,和z来输入v和vn
					int x, y, z;
					stringstream ss;
					ss << a;
					//facet不知道有几个顶点，故采用字符串流循环读入
					while (!ss.eof())
					{
						string tmp="";

						ss >> x;
						tmp.push_back(ss.get());
						tmp.push_back(ss.get());
						ss >> z;
						//若读取出错
						if (ss.fail()||tmp != "//")
							break;
						//读取成功
						f.pointlist.push_back(indpair(x, 0, z));
					}
					//facet面片存储
					flist.push_back(f);
				}
				
			}
			fs.close();
		}
		print("read_from_file");
	}
	//多边形mesh转为triangle mesh
	void triangle_mesh()
	{
		//用新list来装每个facet新产生的额外三角形面片
		list<facetfromobj> nfacetlist;
		for (list<facetfromobj>::iterator it = flist.begin(); it != flist.end(); it++)
		{			

			//当前facet的端点列表
			vector<indpair>& plst = it->pointlist;
			int pcnt = plst.size();
			if (pcnt > 3)
			{
				//除第一个三角形以外，新增三角形面片
				for (int j = 2; j + 1 < pcnt; j++)
				{
					facetfromobj ntriangle;
					ntriangle.pointlist.push_back(plst[0]);
					ntriangle.pointlist.push_back(plst[j]);
					ntriangle.pointlist.push_back(plst[j + 1]);
					nfacetlist.push_back(ntriangle);
				}
				//只保留第一个三角形的三个点
				plst.resize(3);
			}
		}
		flist.splice(flist.end(), nfacetlist);
		print("triangle_mesh");
	}
	//打印输出当前面数，点数
	void print(string funcname)
	{
		printf("-------------------------------\n");
		printf((funcname + "已完成\n").c_str());
		printf("当前顶点v数量%d 法向vn数量%d 面片f数量%d\n", vlist.size(), vnlist.size(), flist.size());
		map<int, int>mp;
		printf("总面片数量=%d\n", flist.size());
		for (list<facetfromobj>::iterator it = flist.begin(); it != flist.end(); it++)
		{
			int cnt = it->pointlist.size();
			if (mp.find(cnt) == mp.end())
			{
				mp[cnt] = 0;
			}
			mp[cnt]++;
		}
		for (map<int, int>::iterator a = mp.begin(); a != mp.end(); a++)
		{
			printf("%d点面片数为%d\n", a->first, a->second);
		}
	}
};

//齐次坐标
struct hvec
{
	double xyzw[4];
	double eps = 1e-8;
	hvec(double a = 0, double b = 0, double c = 0,double d = 1)
	{
		xyzw[0] = a; xyzw[1] = b; xyzw[2] = c, xyzw[3] = d;
	}
	//从三维vec向量解析过来，默认w=1
	hvec(const vec& b)
	{
		for (int i = 0; i < 3; i++)
		{
			xyzw[i] = b.xyz[i];
		}
		xyzw[3] = 1;
	}
	//拷贝构造
	hvec(const hvec& b)
	{
		memcpy(xyzw, b.xyzw, sizeof(xyzw));
	}
	//向量相等当且仅当回到三维后相等
	bool operator == (const hvec& b)
	{
		for (int i = 0; i < 4; i++)
		{
			if (fabs(xyzw[i] / xyzw[3] - b.xyzw[i] / b.xyzw[3]) >= eps)
				return false;
		}
		return true;
	}
	hvec operator +(const hvec& b)
	{
		return hvec(xyzw[0] + b.xyzw[0], xyzw[1] + b.xyzw[1], xyzw[2] + b.xyzw[2], xyzw[3] + b.xyzw[3]);
	}
	hvec operator -(const hvec& b)
	{
		return hvec(xyzw[0] - b.xyzw[0], xyzw[1] - b.xyzw[1], xyzw[2] - b.xyzw[2], xyzw[3] - b.xyzw[3]);
	}
	//向量数乘c
	hvec mulc(double c)
	{
		return hvec(xyzw[0] * c, xyzw[1] * c, xyzw[2] * c, xyzw[3] * c);
	}
	//令w=1
	void normalize()
	{
		if (fabs(xyzw[3]) < eps)
		{
			//printf("w为0，无法归一化\n");
			for (int i = 0; i < 3; i++)
			{
				xyzw[i] = DBL_MAX;
			}
			return;
		}
		for (int i = 0; i < 3; i++)
		{
			xyzw[i] = xyzw[i] / xyzw[3];
		}
		xyzw[3] = 1;
	}

	//仅在x,y分量上做归一化，令w=1，用于保留投影变换后的z值。
	void normalize_x_y()
	{
		if (fabs(xyzw[3]) < eps)
		{
			//printf("w为0，无法xy归一化\n");
			for (int i = 0; i < 2; i++)
			{
				xyzw[i] = DBL_MAX;
			}
			return;
		}
		for (int i = 0; i < 2; i++)
		{
			xyzw[i] = xyzw[i] / xyzw[3];
		}
		xyzw[3] = 1;
	}

	//甩掉w坐标，转换为3维坐标
	vec convert_to_vec()const
	{
		return vec(xyzw[0], xyzw[1], xyzw[2]);
	}
};
struct hmat
{
	//4*4齐次坐标系的变换。
	double A[4][4];
	hmat()
	{
		memset(A, 0, sizeof(A));
	}
	hmat(const hmat& b)
	{
		memcpy(A, b.A, sizeof(A));
	}
	//矩阵乘向量
	hvec operator*(const hvec& x)const
	{
		hvec b(0, 0, 0, 0);
		for (int i = 0; i < 4; i++)
		{
			for (int k = 0; k < 4; k++)
			{
				b.xyzw[i] += A[i][k] * x.xyzw[k];
			}
		}
		return b;
	}
	//矩阵乘矩阵
	hmat operator*(const hmat& b)const
	{
		const double(*B)[4] = b.A;
		hmat c;
		double(*C)[4] = c.A;
		for (int k = 0; k < 4; k++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		return c;
	}
	//平移坐标变换的矩阵
	static hmat get_translation_hmat(const vec& b)
	{
		hmat T;
		double t[][4] = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
		for (int i = 0; i < 3; i++)
			t[i][3] = b.xyz[i];
		memcpy(T.A, t, sizeof(T.A));
		return T;
	}
	//绕x轴逆时针旋转theta角度,单位为弧度制。
	static hmat get_x_axis_rot_hmat(double theta)
	{
		hmat T;
		double t[][4] = { {1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1} };
		t[1][1] = t[2][2] = cos(theta);
		t[2][1] = sin(theta);
		t[1][2] = -t[2][1];
		memcpy(T.A, t, sizeof(T.A));
		return T;
	}
	//绕y轴逆时针旋转theta角度,单位为弧度制。
	static hmat get_y_axis_rot_hmat(double theta)
	{
		hmat T;
		double t[][4] = { {0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,1} };
		t[0][0] = t[2][2] = cos(theta);
		t[0][2] = sin(theta);
		t[2][0] = -t[0][2];
		memcpy(T.A, t, sizeof(T.A));
		return T;
	}
	//构造单位矩阵
	//static hmat I;
	//static bool constructed;
	static hmat get_unity()
	{
		double I[][4] = { {1,0,0,0},{0,1,0,0 },{0,0,1,0},{0,0,0,1} };
		hmat T;
		memcpy(T.A, I, sizeof(T.A));
		return T;
	}

	void print()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				printf("%f ", A[i][j]);
			}
			printf("\n");
		}
	}
};

//三角形面片
struct triangle
{
	//三个端点坐标
	indpair ends[3];
	int color;

	//从facetfromobj拷贝
	triangle(const facetfromobj& b)
	{
		if (b.pointlist.size() >= 3)
		{
			for (int i = 0; i < 3; i++)
				ends[i] = b.pointlist[i];
		}
	}
	//默认拷贝构造函数
	/*triangle(const triangle& b)
	{
		memcpy(ends,b.ends,sizeof(ends));

	}*/
};

//
// 
struct boundingbox
{
	//AABB boundingbox,b[i][0]——第i维坐标最小值，b[i][1],第i维坐标最大值
	double b[3][2];
	boundingbox()
	{
		memset(b, 0, sizeof(b));
	}
	boundingbox(double x, double y, double z, double w, double m, double n)
	{
		b[0][0] = x; b[0][1] = y; b[1][0] = z; b[1][1] = w; b[2][0] = m; b[2][1] = n;
	}
	boundingbox(boundingbox& c)
	{
		memcpy(b, c.b, sizeof(b));
	}
	void print()const
	{
		for (int i = 0; i < 3; i++)
		{
			printf("%f %f|", b[i][0],b[i][1]);
		}
		printf("\n");
	}
	//静态的计算任意顶点集合的boundingbox,这个boundingbox边界上可以有顶点。
	static boundingbox get_bounding_box(const vector<hvec>& vlist)
	{
		boundingbox a;
		memset(a.b, 0, sizeof(a.b));

		for (int i = 0; i < vlist.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				//取xyz坐标时要注意除以w
				double tt = vlist[i].xyzw[j] / vlist[i].xyzw[3];
				if (tt < a.b[j][0])
					a.b[j][0] = tt;
				if (tt > a.b[j][1])
					a.b[j][1] = tt;
			}
		}
		return a;
	}

	//用于未来实现视域四棱体精准裁剪所有现存面片。
	bool is_intersect_with_triangle(const triangle &tri, const vector<hvec>& vlist)const
	{
		return false;
	}

	bool is_contains_triangle(const triangle& tri, const vector<hvec>& vlist)const
	{
		for (int i = 0; i < 3; i++)
		{
			const hvec& p = vlist[tri.ends[i].vi - 1];
			
			//三角形每个点p都完全在包围盒里才返回真，不能在边界上

			for (int j = 0; j < 3; j++)
			{
				if (p.xyzw[j] <= b[j][0] || p.xyzw[j] >= b[j][1])
				{
					return false;
				}
			}
		}
		return true;
	}
};



//用于单个模型的核心数据结构。
struct geometry
{
	vector<hvec> vlist;
	vector<hvec> vnlist;
	list<triangle> flist;

	//将geometryfromobj转化为geometry
	void assign(const geometryfromobj& b)
	{
		vlist.assign(b.vlist.begin(), b.vlist.end());
		vnlist.assign(b.vnlist.begin(), b.vnlist.end());
		flist.assign(b.flist.begin(), b.flist.end());
		print("assign");
	}
	void print(string funcname)
	{
		printf("-------------------------------\n");
		printf((funcname + "已完成\n").c_str());
		printf("当前顶点v数量%d 顶点法向vn数量%d 面片f数量%d\n", vlist.size(), vnlist.size(), flist.size());
	}
	//将模型居中于坐标原点
	void centered()
	{
		//计算当前模型的包围盒
		boundingbox bb = boundingbox::get_bounding_box(vlist);
		//得到包围盒的中点
		vec vs(-(bb.b[0][0] + bb.b[0][1])/2, -(bb.b[1][0] + bb.b[1][1]) / 2, -(bb.b[2][0] + bb.b[2][1]) / 2);
		//得到仿射变换矩阵
		hmat T = hmat::get_translation_hmat(vs);

		//开始变换
		for (int i = 0; i < vlist.size(); i++)
		{
			vlist[i] = T * vlist[i];
		}
		for (int i = 0; i < vnlist.size(); i++)
		{
			vnlist[i] = T * vnlist[i];
		}
		print("centered");
	}
	/*
	boundingbox apply_transform_and_get_bounding_box(hmat t)
	{
		boundingbox a;
		memset(mm, 0, sizeof(mm));
		for (int i = 0; i < vlist.size(); i++)
		{
			t.
			vlist[i];
			for (int j = 0; j < 3; j++)
			{
				double tt = vlist[i].xyzw[j];
				if (tt > mm[0][j])
					mm[0][j] = tt;
				if (tt < mm[1][j])
					mm[1][j] = tt;
			}
		}
		return
	}*/
};



//完整场景
struct scene {

	//存储场景中所有顶点的齐次坐标
	vector<hvec> vlist;
	//存储场景中所有顶点的法向坐标
	vector<hvec> vnlist;
	//存储场景中所有三角面片的顶点编号
	list<triangle> flist;
	
	void print(string funcname)
	{
		printf("-------------------------------\n");
		printf((funcname + " 已完成!\n").c_str());
		printf("当前顶点v数量%d 法向vn数量%d 面片f数量%d\n", vlist.size(), vnlist.size(), flist.size());
	}

	//场景造型，随机在世界坐标系中复制生成单个模型
	//模型变换，使用仿射变换来复制基础模型,base->指向基础模型的引用,n->复制n次
	void generate_from_base(const geometry& base,int n)
	{
		boundingbox bbox = boundingbox::get_bounding_box(base.vlist);


		double range[3];
		//计算xyz坐标的随机范围range 
		for (int i = 0; i < 3; i++)
		{
			range[i] = 2 * pow(n, 1.0 / 3)*bbox.b[i][1];
		}

		unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();

		std::default_random_engine generator;
		//std::default_random_engine generator2;

		//旋转角度theta的分布
		std::uniform_real_distribution<double> thetapdf(0, 2 * std::numbers::pi);
		
		//coordinate（坐标）的分布
		std::uniform_real_distribution<double> coordpdf(-1.0, 1.0);
	
		for (int i = 0; i < n; i++)
		{
			double xtheta = thetapdf(generator);
			double ytheta = thetapdf(generator);
			//旋转矩阵
			hmat R = hmat::get_y_axis_rot_hmat(ytheta) * hmat::get_x_axis_rot_hmat(xtheta);
			
			double coord[3];
			for (int j = 0; j < 3; j++)
			{
				coord[j] = coordpdf(generator)*range[j];
			}

			//先旋转，再平移。
			hmat T = hmat::get_translation_hmat(vec(coord)) * R;

			//对每个模型duplicate开始变换,顶点v需要做T变换，法向vn需要做R变换，纹理vt(如果有)不需要做变换。
			//triangle中顶点的编号从1开始。因此每个新的duplicate的顶点编号需要加上场景中已有的顶点数量。
			//pvcnt->present v number目前场景scene中的v的数量;pvncnt->present vn 的数量
			int pvcnt = vlist.size();
			int pvncnt = vnlist.size();

			for (int j = 0; j < base.vlist.size(); j++)
			{
				hvec nv = T * base.vlist[j];
				vlist.push_back(nv);
			}

			for (int j = 0; j < base.vnlist.size(); j++)
			{
				hvec nv = R * base.vnlist[j];
				vnlist.push_back(nv);
			}

			for (list<triangle>::const_iterator it = base.flist.cbegin(); it != base.flist.cend(); it++)
			{
				//k < 3,因为每个triangle面片都只有三个顶点。
				triangle ntriangle = *it;
				for (int k = 0; k < 3; k++)
				{
					//新的顶点编号
					ntriangle.ends[k].vi += pvcnt;
					ntriangle.ends[k].vni += pvncnt;
				}
				ntriangle.color = i + 1;
				//新的三角形面片汇总到本场景中
				flist.push_back(ntriangle);
			}
		}
		print("generate_from_base");
	}

};


//2维整数坐标
struct vec2d
{
	int xy[2];
	vec2d(int a = 0, int b = 0)
	{
		xy[0] = a;
		xy[1] = b;
	}
	bool operator < (const vec2d &b)const
	{
		if (xy[0] != b.xy[0])
			return xy[0] < b.xy[0];
		return xy[1] < b.xy[1];
	}
	vec2d operator - (const vec2d& b)const
	{
		vec2d nv;
		nv.xy[0] = xy[0] - b.xy[0];
		nv.xy[1] = xy[1] - b.xy[1];
		return nv;
	}
	vec2d operator + (const vec2d& b)const
	{
		vec2d nv;
		nv.xy[0] = xy[0] + b.xy[0];
		nv.xy[1] = xy[1] + b.xy[1];
		return nv;
	}
	vec2d operator / (int c)const
	{
		vec2d nv;
		nv.xy[0] = xy[0]/c;
		nv.xy[1] = xy[1]/c;
		return nv;
	}
	bool operator == (const vec2d& b)const
	{
		return (xy[0] == b.xy[0]) && (xy[1] == b.xy[1]);
	}

	//vec构造函数重载中无法使用定义在后面的struct的类型。除非使用class，并在vec结构体前面加类声明。
	//会破坏后面的代码，所以在本结构体中使用convert方法做替代。
	static vec convert_to_vec(const vec2d& b)
	{
		//二维向量对应的三维齐次坐标
		return vec(b.xy[0], b.xy[1], 1);
	}
};

//2维正整数盒子
struct box2d
{
	//ends[0]最靠近原点的向量，ends[1]最远离原点的向量。闭区间。
	vec2d ends[2];

	//判断2维点是否在2d盒子内,边界也认为包含
	bool is_point_within(const vec2d& b)const
	{
		for (int i = 0; i < 2; i++)
		{
			if (b.xy[i] < ends[0].xy[i] || b.xy[i] > ends[1].xy[i])
			{
				return false;
			}
		}
		return true;
	}
	//判断浮点2维坐标是否在本盒子内部。盒子的边界像素取像素中心坐标。
	bool is_point_within(double xy[2])const
	{
		for (int i = 0; i < 2; i++)
		{
			//边界不算Within
			if (xy[i] <= ends[0].xy[i] || xy[i] >= ends[1].xy[i])
				return false;
		}
		return true;
	}
	//判断2d盒子是否在本盒子内部，边界也认为包含
	bool is_box2d_within(const box2d &b)const
	{
		vec2d points[2][2];
		points[1][1] = b.ends[1];
		points[0][0] = points[0][1] = points[1][0] = b.ends[0];
		points[0][1].xy[1] = points[1][1].xy[1];
		points[1][0].xy[0] = points[1][1].xy[0];
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (!is_point_within(points[i][j]))
				{
					return false;
				}
			}
		}
		return true;
	}

	//判断三角形是否完全在本盒子内部，边界也认为包含
	bool is_triangle_within(const triangle& b, const vector<hvec>& v_in_window_list,int xmax,int ymax)const
	{
		for (int i = 0; i < 3; i++)
		{
			int x = min(xmax - 1, max(0, (int)floor(v_in_window_list[b.ends[i].vi - 1].xyzw[0]+0.5)));
			int y = min(ymax - 1, max(0, (int)floor(v_in_window_list[b.ends[i].vi - 1].xyzw[1]+0.5)));

			if (!is_point_within(vec2d(x,y)))
			{
				return false;
			}
		}
		return true;
	}

	//判断本2d盒子是否与三角形有公共点。使用重心坐标来判断。在double顶点坐标空间内，有公共点当且仅当，存在盒子顶点的像素中心在三角形内部，或者存在三角形顶点在盒子内部。
	//具体计算参考维基的重心坐标章节。
	bool is_intersect_with_triangle(const triangle& b, const vector<hvec>& v_in_window_list, int xmax, int ymax)const
	{

		vec points[3];
		for (int i = 0; i < 3; i++)
		{
			points[i] = vec(v_in_window_list[b.ends[i].vi - 1].xyzw[0], v_in_window_list[b.ends[i].vi - 1].xyzw[1], 1);
			
			
			double xy[2];
			xy[0] = points[i].xyz[0];
			xy[1] = points[i].xyz[1];
			/*错误
			//如果三角形角点在整数2d盒子内，则一定相交
			int x = min(xmax - 1, max(0, (int)floor(points[i].xyz[0] + 0.5)));
			int y = min(ymax - 1, max(0, (int)floor(points[i].xyz[1] + 0.5)));*/
			
			//三角形的角点要完全在像素中心的网格内部，但此处还是会多引入不存在的点，需要再想。
			if (is_point_within(xy))
			{
				return true;
			}
		}

		//算重心坐标。
		//三角形abc的面积
		double abc = vec::determinant(points[0],points[1],points[2]);
		
		/*错误
		//四个顶点的3维齐次坐标,为防止漏掉一些像素，取整数坐标外移0.5作为角点坐标。
		vec corner[2][2];
		corner[0][0] = corner[0][1] = corner[1][0] = (vec2d::convert_to_vec(ends[0])+vec(-0.5,-0.5,0));
		corner[0][1] = corner[0][1] + vec(0, (ends[1] - ends[0]).xy[1] + 1, 0);
		corner[1][0] = corner[1][0] + vec((ends[1] - ends[0]).xy[0] + 1, 0 , 0);
		corner[1][1] = (vec2d::convert_to_vec(ends[1]) + vec(0.5,0.5));
		*/
		
		vec corner[2][2];
		corner[0][0] = corner[0][1] = corner[1][0] = vec2d::convert_to_vec(ends[0]);
		corner[0][1] = corner[0][1] + vec(0, (ends[1] - ends[0]).xy[1] + 1, 0);
		corner[1][0] = corner[1][0] + vec((ends[1] - ends[0]).xy[0] + 1, 0, 0);
		corner[1][1] = vec2d::convert_to_vec(ends[1]);

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				//系数
				double lmd[3];
				//pbc的面积/abc的面积，得到系数。apc,abp同理
				lmd[0]  = vec::determinant(corner[i][j], points[1], points[2])/abc;
				lmd[1] = vec::determinant(points[0], corner[i][j], points[2])/abc;
				lmd[2] = vec::determinant(points[0], points[1], corner[i][j])/abc;
				bool in = true;

				for (int k = 0; k < 3; k++)
				{
					//只要不在三角形外部，则一定相交。在三角形边上也认为相交，见wiki重心坐标
					if (lmd[k] < 0)
					{
						in = false;
						break;
					}
				}
				if (in)
				{
					return true;
				}
			}
		}
		return false;
	}
};

//层次z_buffer
struct node
{
	//进入节点已知的数据。
	node* parent;
	//左上角坐标coord[0], 右下角坐标coord[1],
	box2d box;

	//进入节点未知待求解的数据。kids递归向下求解，z回溯求解。
	node* kids[2][2];
	double z;
};
struct zpyramid
{
	node* root;
	//ends0为zpyramid最靠近坐标原点的整数下标向量，ends1为最原理原点的整数下标向量。
	zpyramid(const vec2d &ends0,const vec2d &ends1)
	{
		root = new node();
		//左上角
		root->box.ends[0] = ends0;
		//右下角
		root->box.ends[1] = ends1;
		root->z = -DBL_MAX;
	}
	//单个像素点位置映射到堆中的叶节点。
	map< vec2d , node* > mp;
	//递归向下建立pyramid结点，回溯构建z值小顶堆。返回值为最小z值。
	double make_up_heap(node* p)
	{
		//printf("%d %d %d %d\n", p->box.ends[0].xy[0], p->box.ends[0].xy[1], p->box.ends[1].xy[0], p->box.ends[1].xy[1]);
		double minz = DBL_MAX;
		//距离坐标原点的左下指向右上的向量。
		vec2d corner[2][2];
		corner[1][1] = p->box.ends[1];
		corner[0][0] = corner[0][1] = corner[1][0] = p->box.ends[0];
		corner[0][1].xy[1] = corner[1][1].xy[1];
		corner[1][0].xy[0] = corner[1][1].xy[0];

		//坐标左下角
		if (!(corner[0][0] == corner[1][1]))
		{
			p->kids[0][0] = new node();
			p->kids[0][0]->parent = p;
			p->kids[0][0]->box.ends[0] = corner[0][0];
			p->kids[0][0]->box.ends[1] = (corner[0][0]+corner[1][1])/2;
			minz = fmin(minz, make_up_heap(p->kids[0][0]));
		}
		else
		{
			//此时当前节点为单个像素，不可进一步拆分。
			minz = -DBL_MAX;
			mp[corner[0][0]] = p;
		}

		//坐标右下角
		if (!(corner[0][0] == corner[1][0]))
		{
			p->kids[1][0] = new node();
			p->kids[1][0]->parent = p;
			p->kids[1][0]->box.ends[0] = (corner[0][0] + corner[1][0])/2 + vec2d(1,0);
			p->kids[1][0]->box.ends[1] = (corner[1][0] + corner[1][1]) / 2;
			minz = fmin(minz, make_up_heap(p->kids[1][0]));
		}

		//坐标左上角

		if (!(corner[0][0] == corner[0][1]))
		{
			p->kids[0][1] = new node();
			p->kids[0][1]->parent = p;
			p->kids[0][1]->box.ends[0] = (corner[0][0] + corner[0][1]) / 2 + vec2d(0,1);
			p->kids[0][1]->box.ends[1] = (corner[0][1] + corner[1][1]) / 2;
			minz = fmin(minz, make_up_heap(p->kids[0][1]));
		}

		//坐标右上角
		if (corner[0][0] != corner[0][1] && corner[0][0] != corner[1][0])
		{
			p->kids[1][1] = new node();
			p->kids[1][1]->parent = p;
			p->kids[1][1]->box.ends[0] = (corner[0][0] + corner[1][1]) / 2 + vec2d(1, 1);
			p->kids[1][1]->box.ends[1] = corner[1][1];
			minz = fmin(minz, make_up_heap(p->kids[1][1]));
		}

		return p->z = minz;
	}

	//得到单个像素的z值
	double get_pixel_z(const vec2d& pixel)
	{
		if (mp.find(pixel) == mp.end())
		{
			printf("z_pyramid中找不到该像素对应的节点 # %d # %d #\n", pixel.xy[0], pixel.xy[1]);
			return 0;
		}
		return mp[pixel]->z;
	}

	//更新单个像素的z_buffer值
	void update(const vec2d& pixel, double nz)
	{
		if (mp.find(pixel) == mp.end())
		{
			printf("update_z_pyramid中找不到该像素对应的节点 # %d # %d #\n",pixel.xy[0],pixel.xy[1]);
			return;
		}
		node* p = mp[pixel];
		p->z = nz;
		while (p->parent != NULL)
		{
			node* par = p->parent;
			//重新求一遍最小值父节点z_buffer的最小值。
			double nmin = p->z;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (par->kids[i][j] != NULL)
					{
						if (par->kids[i][j]->z < nmin)
						{
							nmin = par->kids[i][j]->z;
						}
					}
				}
			}
			//如果更新成功，则继续向上更新，否则退出。
			if (nmin > par->z)
			{
				par->z = nmin;
				p = par;
			}
			else
			{
				break;
			}
		}
	}

	//判断覆盖任意子矩阵的最小pyramid节点,并返回该节点z值,a为左上角，b为右下角。
	double min_enclosing_z(const box2d &a, node *p)
	{
		//若子矩阵完全在p的任意子节点范围内，则往下递归。否则，返回当前p->z值；
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (p->kids[i][j] == NULL)
					continue;
				if (p->kids[i][j]->box.is_box2d_within(a))
				{
					return min_enclosing_z(a,p->kids[i][j]);
				}
			}
		}
		return p->z;
	}
	//判断覆盖三角面片的最小pyramid节点
	double min_enclosing_z(const triangle& a, const vector<hvec> & v_in_window_list, node* p,int xmax,int ymax)
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (p->kids[i][j] == NULL)
					continue;
				if (p->kids[i][j]->box.is_triangle_within(a,v_in_window_list,xmax,ymax))
				{
					return min_enclosing_z(a,v_in_window_list,p->kids[i][j],xmax,ymax);
				}
			}
		}
		return p->z;
	}

	//a与节点p相交的部分是否真的有可见点。
	//假设进入的时候，三角形一定与当前节点的2d盒子有交集。
	//xmax,ymax->2d范围最大的整数坐标+1。triz->三角形a的最大z值
	bool is_visible(const triangle& a, const vector<hvec>& v_in_window_list, node* p, double triz ,int xmax, int ymax)
	{
		//如果相交的部分的最小z都更大，那么肯定不可见。
		if (p->z >= triz)
			return false;
		//如果已经是单个像素了，并且前一条不满足，则一定可见。
		if (p->box.ends[0] == p->box.ends[1])
		{
			return true;
		}

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (p->kids[i][j] == NULL)
					continue;
				//如果这个盒子与a相交
				if (p->kids[i][j]->box.is_intersect_with_triangle(a,v_in_window_list,xmax,ymax))
				{
					//递归向下，如果相交的部分还有可见点，那么当前节点也有可见点
					if (is_visible(a, v_in_window_list, p->kids[i][j], triz,xmax,ymax))
					{
						return true;
					}
				}
			}
		}
		//相交的所有部分都不可见。
		return false;
	}
};

struct pipeline
{
	//存储当前场景中所有顶点的齐次坐标
	vector<hvec> vlist;
	//存储变换到屏幕坐标系后顶点的齐次坐标
	vector<hvec> v_in_window_list;

	//存储当前场景中所有顶点的法向坐标
	vector<hvec> vnlist;
	//存储当前场景中所有三角面片
	list<triangle> flist;

	//Transform存储取景变换，投影变换,视窗变换等后续所有变换的复合变换。

	//cuvn视点坐标向量
	hvec c, u, v, n;

	//投影平面到视点的距离，应当为正数。
	double d;

	//ouvn坐标中，视域四棱台的大小。umax,vmax分别为视窗u,v坐标最大值，front,back是z值范围，以下四个坐标都为正数。正负轴对称。
	double umax;
	double vmax;
	double front;
	double back;

	//oxyz坐标中的四棱体范围，都为正数。正负轴对称。
	double xmax;
	double ymax;
	double zmax;

	//屏幕分辨率
	int xpixelnum;
	int ypixelnum;

	//场景8叉停止细分的三角面片数量
	const int stop_triangle_num = 10;


	//显示帧缓冲器
	int** framebuffer;

	void print(string funcname)
	{
		printf("-------------------------------\n");
		printf((funcname+"已完成\n").c_str());

	}

	//取景变换，世界坐标变换到ouvn坐标,viewpoint->视点坐标，n->视线方向,up->向上的方向
	void view_transform(const vec& viewpoint, vec n, const vec& up, hmat& transform)
	{
		if (n == vec::get_zero_vector() || up == vec::get_zero_vector())
		{
			printf("视线方向和up方向不能为0\n");
			return;
		}
		n = n * (1.0 / n.norm2());
		vec v = n.cross_product(up);
		if (v == vec::get_zero_vector())
		{
			printf("视线方向与up方向不能在一条直线上\n");
			return;
		}
		v = v * (1.0 / v.norm2());
		vec u = v.cross_product(n);
		u = u * (1.0 / u.norm2());

		this->c = viewpoint;
		this->u = u;
		this->v = v;
		this->n = n;

		//构造坐标变换矩阵
		hmat T;

		for (int j = 0; j < 3; j++)
			T.A[0][j] = u.xyz[j];
		for (int j = 0; j < 3; j++)
			T.A[1][j] = v.xyz[j];
		for (int j = 0; j < 3; j++)
			T.A[2][j] = n.xyz[j];
		T.A[3][3] = 1;

		hvec c = T * hvec(viewpoint);
		for (int i = 0; i < 3; i++)
			T.A[i][3] = -c.xyzw[i];

		//与现有变换复合
		transform = T * transform;

		print("viewtransform");
	}

	//视域四棱台,仍然在ouvn坐标系中,u为up方向,theta为相应维度视线与水平线的最大夹角,theta<pi/2。
	//front->视域四棱台前面，back->视域四棱台后面
	void viewing_frustum(double utheta, double vtheta, double front, double back)
	{
		if (front < 0 || back < 0 || front > back)
		{
			printf("front或back有误\n");
			return;
		}
		if (utheta > 0.5 * std::numbers::pi || vtheta > 0.5 * std::numbers::pi || utheta < 0 || vtheta < 0)
		{
			printf("视野夹角utheta或vtheta有误\n");
			return;
		}
		double d = (front + back) / 2;
		umax = d * tan(utheta);
		vmax = d * tan(vtheta);
		this->front = front;
		this->back = back;
		this->d = d;
		print("viewing_frustum");
	}

	//透视投影变换,d为投影平面，d>0
	void projection(hmat& transform)
	{
		if (d < 0)
		{
			printf("投影平面应当与视点为正距离");
			return;
		}

		hmat T = hmat::get_unity();
		T.A[3][2] = 1 / d;
		T.A[3][3] = 0;
		transform = T * transform;

		this->d = d;

		print("projection");
	}


	//真正通过矩阵计算应用变换transform到场景中，并将结果填充到vlist，vnlist，flist
	//transfrom直接应用至v,transform去掉仿射偏移后应用至vn。注意透视投影会影响法向量的方向。
	//分块矩阵可以证明正确性：Ax+b/(c'x+d)，仿射偏移部分为b/c'x+d，因此只需要令b=0即可。

	//第一次应用transform变换。
	void apply_transform_from_scene(const scene& s, const hmat& transform)
	{
		for (int i = 0; i < s.vlist.size(); i++)
		{
			hvec p = transform * s.vlist[i];
			//重要：在投影变换之后，把真实的z值求出来。可以证明，包含投影变换的复合transform应用于初始数据的齐次坐标后，得到的x,y,z分量仍然等于仅对其做仿射坐标变换的分量。
			//因此，复合投影变换之后将transform应用至初始数据，得到真实的z值，为视域四棱台clipping做准备。
			p.normalize_x_y();
			vlist.push_back(p);
		}

		//bremoved -> 去掉b的变换矩阵
		hmat bremoved = transform;
		for (int i = 0; i < 3; i++)
			bremoved.A[i][3] = 0;

		for (int i = 0; i < s.vnlist.size(); i++)
		{
			hvec p = bremoved * s.vnlist[i];
			p.normalize();
			vnlist.push_back(p);
		}

		//深度拷贝triangle面片
		flist.assign(s.flist.cbegin(), s.flist.cend());

		print("apply_transform_from_scene");
	}

	//在本场景现有的坐标上应用T变换，非第一次应用变换。
	void apply_transform(const hmat& transform)
	{
		for (int i = 0; i < vlist.size(); i++)
		{
			hvec p = transform * vlist[i];
			vlist[i] = p;
		}

		//bremoved -> 去掉b的变换矩阵
		hmat bremoved = transform;
		for (int i = 0; i < 3; i++)
			bremoved.A[i][3] = 0;

		for (int i = 0; i < vnlist.size(); i++)
		{
			hvec p = bremoved * vnlist[i];
			vnlist[i] = p;
		}
	}

	//裁剪掉视域四棱台以外的部分。
	void clipping()
	{

		boundingbox frustum(-umax, umax, -vmax, vmax, front, back);

		for (list<triangle>::iterator it = flist.begin(); it != flist.end(); )
		{

			//此处考虑到裁剪算法的复杂性，只保留完全在视域四棱体内部的三角面片。未来可以将真正的三维裁剪在此处应用。
			if (!frustum.is_contains_triangle(*it, vlist))
			{
				//注意erase的返回值指向被删除的下一个元素。
				it = (flist.erase(it));
				continue;
			}
			it++;
		}

		//根据剩余的triangle，删除相应的顶点坐标和顶点法向，同时修改triangle中顶点的编号。
		vector<hvec> nvlist;
		vector<hvec> nvnlist;
		//map建立旧点编号到新的点编号的映射关系。
		map<int, int> vmp;
		map<int, int> vnmp;
		for (list<triangle>::iterator it = flist.begin(); it != flist.end(); it++)
		{
			for (int i = 0; i < 3; i++)
			{
				int vi = it->ends[i].vi - 1;
				if (vmp.find(vi) == vmp.end())
				{
					nvlist.push_back(vlist[vi]);
					vmp[vi] = nvlist.size();
				}
				int nvi = vmp[vi];
				it->ends[i].vi = nvi;

				int vni = it->ends[i].vni - 1;
				if (vnmp.find(vni) == vnmp.end())
				{
					nvnlist.push_back(vnlist[vni]);
					vnmp[vni] = nvnlist.size();
				}
				int nvni = vnmp[vni];
				it->ends[i].vni = nvni;
			}
		}

		vlist.assign(nvlist.begin(), nvlist.end());
		vnlist.assign(nvnlist.begin(), nvnlist.end());
		print("clipping");
	}

	//ouvn变换到oxyz；oxyz的坐标原点放到投影平面上，与视点距离d。
	void ouvn_to_oxyz(hmat& transform)
	{
		hmat T;
		T.A[0][1] = T.A[1][0] = T.A[3][3] = 1;
		T.A[2][2] = -1;
		T.A[2][3] = d;
		transform = T * transform;
		xmax = vmax;
		ymax = umax;
		//小心精度误差。
		zmax = fmax(d - front, back - d);
		print("ouvn_to_oxyz");
	}

	//视窗变换->变换到像素为单位的屏幕坐标系,屏幕坐标系原点在左上角。x向下为正，y向右为正。
	//屏幕的像素为xp*yp,例如1024*768
	void window_transform(int xp, int yp, hmat& transform)
	{
		if (xp % 2 && yp % 2)
		{
			printf("屏幕像素需要为偶数!\n");
			return;
		}

		//一半像素
		xp /= 2;
		yp /= 2;
		hmat T;
		//景物空间坐标直接与像素中心对齐，注意此处的坐标变换。
		/*T.A[0][1] = -yp / ymax;
		T.A[1][0] = xp / xmax;
		T.A[0][3] = yp - 0.5;
		T.A[1][3] = xp - 0.5;
		T.A[2][2] = T.A[3][3] = 1;
		*/
		T.A[0][0] = -yp / xmax;
		T.A[1][1] = -xp / ymax;
		T.A[0][3] = yp - 0.5;
		T.A[1][3] = xp - 0.5;
		T.A[2][2] = T.A[3][3] = 1;

		transform = T * transform;
		//print_transform("window_transform");
	}

	//消隐
	map< vec2d, node* >* surface_visibility(hmat& transform, int xpixelnum, int ypixelnum)
	{
		this->xpixelnum = xpixelnum; this->ypixelnum = ypixelnum;

		//帧缓冲区构建
		framebuffer = new int* [ypixelnum];
		for (int i = 0; i < ypixelnum; i++)
			framebuffer[i] = new int[xpixelnum]();


		//将顶点列表预先备份变换到图像坐标系，供scan_convert函数使用。
		hmat T = transform;
		window_transform(xpixelnum, ypixelnum, T);

		for (int i = 0; i < this->vlist.size(); i++)
		{
			hvec tmp = T * this->vlist[i];

			//变换后的坐标不能超过0-n-1的范围。
			//tmp.xyzw[0] = min(double(ypixelnum - 1),max(0,floor(tmp.xyzw[0] + 0.5)));
			//tmp.xyzw[1] = min(double(xpixelnum - 1), max(0, floor(tmp.xyzw[1] + 0.5)));

			this->v_in_window_list.push_back(tmp);
		}

		//建立zpyramid，注意屏幕坐标系与分辨率参数的对应关系。
		zpyramid zpy (vec2d(0, 0), vec2d(ypixelnum - 1, xpixelnum - 1));
		zpy.make_up_heap(zpy.root);

		boundingbox bb(-xmax, xmax, -ymax, ymax, -zmax, zmax);

		//八叉树过程中会不断将三角面片list 8拆分，故使用原有的面片的备份进行运算。
		list<triangle> mylist(flist.begin(), flist.end());
		octree(transform,bb, mylist, stop_triangle_num, zpy);
		map< vec2d, node* >* msp = new map< vec2d, node* >(zpy.mp);
		print("surface_visibility");
		return msp;
	}


	static bool mycmp(const hvec& a, const hvec& b)
	{
		const double eps = 1e-6;
		if (a.xyzw[1] != b.xyzw[1])
			return a.xyzw[1] > b.xyzw[1];
		return a.xyzw[0] < b.xyzw[0];
	}
	//四舍五入找最近的整数，并保证在0——maxn-1的范围。
	inline static int double_to_legal_int(double a, int maxn)
	{
		return max(0, min(maxn - 1, (int)floor(a + 0.5)));
	}

	//计算z_buffer跳过了多少面片
	int cnta = 0;
	int cntb = 0;


	void scan_convert(const list<triangle>& mylist, zpyramid& zpy)
	{
		for (list<triangle>::const_iterator it = mylist.cbegin(); it != mylist.cend(); it++)
		{
			hvec p[3];
			double zmax = -DBL_MAX;
			for (int i = 0; i < 3; i++)
			{
				p[i] = this->v_in_window_list[it->ends[i].vi - 1];
				zmax = fmax(zmax, p[i].xyzw[2]);
			}
			//判断是否真的可见。
			if (!zpy.is_visible(*it, v_in_window_list, zpy.root, zmax, xmax, ymax))
			{
				cntb++;
				continue;
			}

			//只画顶端的顶点
			for (int i = 0; i < 1; i++)
			{
				int x = double_to_legal_int(p[i].xyzw[0], ypixelnum);
				int y = double_to_legal_int(p[i].xyzw[1], xpixelnum);
				if (p[i].xyzw[2] > zpy.get_pixel_z(vec2d(x,y)))
				{
					//当前像素的值需要被更新。
					framebuffer[x][y] = it->color;
					zpy.update(vec2d(x, y), p[i].xyzw[2]);
				}
			}


			//y值大的往前放，y值相同，x小的往前放。 
			sort(p, p + 3, mycmp);

			//normal 当前三角形面片的法向
			vec l = (p[1] - p[0]).convert_to_vec();
			vec r = (p[2] - p[0]).convert_to_vec();
			//l向量在x方向前进的终止位置。
			double endx = p[1].xyzw[0];

			vec normal = l.cross_product(r);
			if (normal.xyz[2] == 0)
			{
				//投影到xy平面上后三点共线。
				continue;
			}
			double dzx = -normal.xyz[0] / normal.xyz[2];
			double dzy = -normal.xyz[1] / normal.xyz[2];

			//除了三点共线的情况以外，只可能l水平或者第二段l水平。具体画图讨论。
			if (l.xyz[0] == 0 )
			{
				//水平线
				l = r;
				r = (p[2] - p[1]).convert_to_vec();
				endx = p[2].xyzw[0];
			}

			double dxl = l.xyz[0] / l.xyz[1];
			double dxr = r.xyz[0] / r.xyz[1];

			double xl = p[0].xyzw[0];
			double xr = xl;
			double zl = p[0].xyzw[2];
			double zr = zl;
			
			//topminus1,顶点的下面一个整数坐标。
			int topminus1 =  double_to_legal_int(p[0].xyzw[1],xpixelnum)-1;
			if (topminus1 < 0)
			{
				//顶点已在边缘。
				continue;
			}
			double ydelta = topminus1 - p[0].xyzw[1];
			
			//当前y的整数坐标
			int yp = topminus1;
			
			while (true)
			{
				double xldelta = ydelta * dxl;
				double xrdelta = ydelta * dxr;
				xl += xldelta;
				xr += xrdelta;
				zl += dzx * xldelta + dzy * ydelta;
				zr += dzx * xrdelta + dzy * ydelta;
				
				//也可以在此处用yp<p[2].xyzw[1]来作为结束条件
				//此处我用endx来进入判断结束和换l向量的情况。

				if ((xl - endx) * l.xyz[0] > 0)
				{
					//结束或者需要换l向量了
					//结束
					if (yp < p[2].xyzw[1])
					{
						break;
					}
					//换l向量
					l = (p[2] - p[1]).convert_to_vec();
					dxl = l.xyz[0] / l.xyz[1];
					ydelta = yp - p[1].xyzw[1];
					xldelta = ydelta * dxl;
					
					xl = p[1].xyzw[0] + xldelta;
					zl = p[1].xyzw[2] + dzx * xldelta + dzy * ydelta;
					endx = p[2].xyzw[0];
				}

				pair<double, double> lr[2];
				if (xl < xr)
				{
					lr[0] = make_pair(xl, zl);
					lr[1] = make_pair(xr, zr);
				}
				else
				{
					lr[0] = make_pair(xr, zr);
					lr[1] = make_pair(xl, zl);
				}
				//int roundxl = max(0, min(ypixelnum - 1, (int)ceil(lr[0].first - 1)));
				//int roundxr = max(0, min(ypixelnum - 1, (int)floor(lr[1].first + 1)));
				
				//int roundxl = double_to_legal_int(lr[0].first,ypixelnum);
				//int roundxr = double_to_legal_int(lr[1].first, ypixelnum);

				//left向上取整，right向下取整。
				int roundxl = (int)ceil(lr[0].first);
				int roundxr = (int)floor(lr[1].first);

				//当前像素点的z值。
				double zz = lr[0].second;
				//每行的初始横向x增量
				double deltax = roundxl - lr[0].first;
				for (int i = roundxl; i <= roundxr; i++)
				{
					zz += deltax * dzx;
					//当前像素的x坐标,y坐标,z值。
					i, yp, zz;
					vec2d points2d(i,yp);
					if (zz > zpy.get_pixel_z(points2d))
					{
						//当前像素的值需要被更新。
						framebuffer[i][yp] = it->color;
						zpy.update(points2d, zz);
					}
					deltax = 1;
				}

				ydelta = -1;
				//列坐标减一
				yp--;
			}

		}
	}

	//场景8叉树 + 层次z_buffer，开始！递归向下，一边建立8叉树，一边消隐。
	//在oxyz空间做8叉树,当box内三角面片数量少于threshold时，停止向下划分空间。transform是为了ouvn->oxyz
	void octree(hmat &transform, const boundingbox &box,list<triangle> &myflist, int threshold, zpyramid &zpy)
	{
		//不要octree，直接扫描转换。
		//scan_convert(myflist, zpy);
		//return;

		//剪枝
		if (myflist.size() == 0)
		{
			return;
		}
		if (myflist.size() < threshold)
		{
			//扫描转换myflist
			scan_convert(myflist, zpy);
			return;
		}

		//检查当前节点是否被完全遮挡。
		//当前box左上角和右下角在oxyz空间的齐次坐标,其z值一样，是当前box靠向视点这一面的z值。
		hvec lefttop(box.b[0][0],box.b[1][1],box.b[2][1]);
		hvec rightdown(box.b[0][1],box.b[1][0],box.b[2][1]);


		hmat T = transform;
		window_transform(xpixelnum, ypixelnum, T);

		lefttop = T * lefttop;
		rightdown = T * rightdown;
		
		//boxinwindow->当前8叉树节点对应的cube离视点最近的面c转换到屏幕坐标系所得到的box投影到2维xy平面上,注意此处较小坐标向下取整，较大坐标向上取整，保证变换后的投影不小于面c的投影。
		box2d boxinwindow;
		
		double minx = min(lefttop.xyzw[0], rightdown.xyzw[0]);
		double maxx = max(lefttop.xyzw[0], rightdown.xyzw[0]);
		double miny = min(lefttop.xyzw[1], rightdown.xyzw[1]);
		double maxy = max(lefttop.xyzw[1], rightdown.xyzw[1]);
		
		int xmin = (int)floor(minx + 0.5);
		if (xmin <= minx)
			xmin += 1;
		int xmax = (int)floor(maxx + 0.5);
		if (xmax >= maxx)
			xmax -= 1;

		int ymin = (int)floor(miny + 0.5);
		if (ymin <= miny)
			ymin += 1;
		int ymax = (int)floor(maxy + 0.5);
		if (ymax >= maxy)
			ymax -= 1;
		/*
		int xmin = double_to_legal_int(floor(minx), ypixelnum);
		int ymin = double_to_legal_int(floor(miny), xpixelnum);
		int xmax = double_to_legal_int(ceil(maxx), ypixelnum);
		int ymax = double_to_legal_int(ceil(maxy), xpixelnum);
		*/
	/*	if (xmin < 0 || ymin < 0 || xmax >= ypixelnum || ymax >= xpixelnum)
		{
			printf("no!!!!!\n");
			system("pause");
		}*/

		boxinwindow.ends[0] = vec2d(xmin, ymin);
		boxinwindow.ends[1] = vec2d(xmax, ymax);

		//若当前盒子完全不可见，则退出。
		double z =  zpy.min_enclosing_z(boxinwindow, zpy.root);

		if (z > lefttop.xyzw[2])
		{
			cnta += myflist.size();
			//printf("###%d\n", cnta);
			//system("pause");
			return;
		}
		

		//构造8个子盒子
		//3个维度的最小、中间、最大3个值。
		double points[3][3];
		for (int i = 0; i < 3; i++)
		{
			points[i][0] = box.b[i][0];
			points[i][2] = box.b[i][1];
			points[i][1] = (points[i][0] + points[i][2]) / 2;
		}
		boundingbox subbox[2][2][2];
		//z值从大的开始，也就是离视点更近的box
		for (int k = 1; k >= 0; k--)
		{
			//x值从小到大
			for (int i = 0; i < 2; i++)
			{
				//y值从小到大
				for (int j = 0; j < 2; j++)
				{
					//boundingbox[i][j]
					subbox[i][j][k] = boundingbox(points[0][i], points[0][i + 1], points[1][j], points[1][j + 1], points[2][k], points[2][k + 1]);
					//subbox[i][j][k].print();
					//system("pause");
				}
			}
		}

		list<triangle> sublist[2][2][2];
		for (list<triangle>::iterator it = myflist.begin(); it != myflist.end();)
		{ 
			bool used = false;
			for (int k = 1; k >= 0; k--)
			{
				for (int i = 0; i < 2; i++)
				{
					for (int j = 0; j < 2; j++)
					{
						if (subbox[i][j][k].is_contains_triangle(*it,this->vlist))
						{
							sublist[i][j][k].push_back(*it);
							it = myflist.erase(it);
							used = true;
							break;
						}
					}
					if (used)
						break;
				}
				if (used)
					break;
			}
			if (!used)
			{
				it++;
			}
		}

		//扫描转换myflist
		scan_convert(myflist, zpy);

		//接着处理子节点
		for (int k = 1; k >= 0; k--)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					octree(transform,subbox[i][j][k], sublist[i][j][k], threshold, zpy);
				}
			}
		}
	}

	void print_framebuffer()
	{
		for (int i = 0; i < ypixelnum; i++)
		{
			for (int j = 0; j < xpixelnum; j++)
			{
				printf("%d ", framebuffer[i][j]);
			}
			printf("\n");
		}
	}
};

vec get_a_coord(string name)
{
	while (1)
	{
		printf( ("请输入" + name).c_str() );
		double x[3];
		scanf_s("%lf%lf%lf", x, x + 1, x + 2);
		if (x[0] > 100 || x[0] < -100 || x[1] > 100 || x[1] < -100 || x[2] > 100 || x[2] < -100)
		{
			printf("坐标不合法\n");
			continue;
		}
		else
		{
			return vec(x);
		}
	}
}
void userinterface(vec &o,vec &n, vec&up)
{
	printf("-----------------------------\n");
	printf("是否需要更改默认视点坐标系？1-需要，0-不需要\n");
	int c;
	while (1)
	{
		scanf_s("%d", &c);
		if (c == 0 || c == 1)
		{
			if (c == 1)
			{
				o = get_a_coord("原点坐标（示例：0 0 0）\n");
				n = get_a_coord("视点方向（示例：1 0 0）\n");
				up = get_a_coord("照相机up方向（示例：0 0 1）\n");
			}
			break;
		}
		printf("请重新输入\n");
	}
}

geometryfromobj soccerobj;
geometry soccer;

scene soccerfield;
pipeline mypipeline;
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	//控制台
	FILE* fp;
	AllocConsole();
	freopen_s(&fp, "CONIN$", "r", stdin);
	freopen_s(&fp, "CONOUT$", "w", stdout);
	freopen_s(&fp, "CONOUT$", "w", stderr);

	//用于计时
	using namespace chrono;
	auto start = system_clock::now();

	//读取基准模型
	soccerobj.read_from_file();
	//多边形面片全部转换为三角面片
	soccerobj.triangle_mesh();
	soccer.assign(soccerobj);

	//模型居中到坐标原点
	soccer.centered();

	//模型变换：在1000个随机位置复制基准模型，构建场景。
	soccerfield.generate_from_base(soccer,1000);
	
	auto check1 = system_clock::now();
	printf("模型变换用时%f秒\n", double(duration_cast<microseconds>(check1 - start).count() * microseconds::period::num / microseconds::period::den));


	//视点ouvn坐标系的oo->o,nn->n,up->up三个向量
	vec oo(0, 0, 0);
	vec nn(1, 0, 0);
	vec up(0, 0, 1);
	
	userinterface(oo, nn, up);



	//取景变换
	hmat I = hmat::get_unity();
	mypipeline.view_transform(oo,nn,up,I);
	//视域四棱台参数设置
	mypipeline.viewing_frustum(1.0 / 3 * std::numbers::pi, 1.0 / 3 * std::numbers::pi, 1, 1000);
	//投影变换
	mypipeline.projection(I);
	//将前面的变换真正应用到场景的点上,得到ouvn下的坐标。
	mypipeline.apply_transform_from_scene(soccerfield,I);
	//根据视域四棱台做裁剪
	mypipeline.clipping();
	//将ouvn下的点转换到oxyz坐标中，方便做消隐。
	I = hmat::get_unity();
	mypipeline.ouvn_to_oxyz(I);

	printf("目前场景中待扫描转换的三角面片总数%d个\n", mypipeline.flist.size());
	
	//在oxyz坐标系中做消隐,其中包含场景八叉树的调用，层次z_buffer的调用，视窗变换的调用。输出为帧缓冲器中的消隐结果。
	map< vec2d, node* >* msp = mypipeline.surface_visibility(I,1024,768);
	
	printf("\n场景八叉树跳过的三角面片数量 = %d 层次z_buffer跳过的三角面片数量 = %d\n", mypipeline.cnta, mypipeline.cntb);
	
	auto check2 = system_clock::now();
	printf("取景变换、投影变换、裁剪、消隐总用时%f秒\n", double(duration_cast<microseconds>(check2 - start).count() * microseconds::period::num / microseconds::period::den));


	system("pause");



	// 初始化绘图窗口
	initgraph(1024, 768);

	// 获取指向显示缓冲区的指针
	DWORD* pMem = GetImageBuffer();
	
	// 直接对显示缓冲区赋值
	for (int i = 0; i < 768; i++)
	{
		for (int j = 0; j < 1024; j++)
		{
			if (mypipeline.framebuffer[i][j] != 0)
			{
				pMem[i * 1024 + j] = BGR(RGB(mypipeline.framebuffer[i][j]/16*10, mypipeline.framebuffer[i][j]/4%16*10, mypipeline.framebuffer[i][j] %4*10));
				/*if (mypipeline.framebuffer[i][j] == -1)
				{
					//如果是顶点，则换红色单独显示。
					pMem[i * 1024 + j] = BGR(RGB(255, 0, 0));
				}
				else
				{
					pMem[i * 1024 + j] = BGR(RGB(255, 255, 255));
				}*/
			}
		}
	}

	/*
	将每个像素的z_buffer值拿出来显示
	for(int i = 0; i < 768; i++)
		for (int j = 0; j < 1024; j++)
		{
			int c = int(msp->at(vec2d(i, j))->z+mypipeline.zmax)%200+50;
			pMem[i * 1024 + j] = BGR(RGB(c,c,c));
		}
	*/
	// 使显示缓冲区生效
	FlushBatchDraw();

	// 保存绘制的图像
	//saveimage(_T("C:\\Users\\luotong\\Desktop\\test.bmp"));
	saveimage(_T(".\\test.bmp"));
	system("pause");
	closegraph();
	return 0;
}
