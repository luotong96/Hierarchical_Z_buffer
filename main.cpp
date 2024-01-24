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
using namespace std;
//定义全局的向量，方便后续坐标变换的计算
struct vec
{
	double xyz[3];
	double eps = 1e-6;
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

	void read_from_file(string filename = "..\\soccerball.obj")
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
		printf((funcname + "已完成\n").c_str());
		printf("当前v数量%d vn数量%d f数量%d\n", vlist.size(), vnlist.size(), flist.size());
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
	double eps = 1e-6;
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
			printf("w为0，无法归一化\n");
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
			printf("w为0，无法归一化\n");
			return;
		}
		for (int i = 0; i < 2; i++)
		{
			xyzw[i] = xyzw[i] / xyzw[3];
		}
		xyzw[3] = 1;
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
	//int color

	//从facetfromobj拷贝
	triangle(const facetfromobj& b)
	{
		if (b.pointlist.size() >= 3)
		{
			for (int i = 0; i < 3; i++)
				ends[i] = b.pointlist[i];
		}
	}
	//拷贝构造函数
	triangle(const triangle & b )
	{
		memcpy(ends,b.ends,sizeof(ends));
	}
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

	//静态的计算任意顶点集合的boundingbox
	static boundingbox get_bounding_box(const vector<hvec>& vlist)
	{
		boundingbox a;
		double(*mm)[2] = a.b;
		memset(mm, 0, sizeof(a.b));
		for (int i = 0; i < vlist.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				//取xyz坐标时要注意除以w
				double tt = vlist[i].xyzw[j] / vlist[i].xyzw[3];
				if (tt < mm[j][0])
					mm[j][0] = tt;
				if (tt > mm[j][1])
					mm[j][1] = tt;
			}
		}
		return a;
	}

	//用于未来实现视域四棱体精准裁剪所有现存面片。
	bool is_intersec_with_triangle(const triangle &tri, const vector<hvec>& vlist)const
	{
		return false;
	}

	bool is_contains_triangle(const triangle& tri, const vector<hvec>& vlist)const
	{
		for (int i = 0; i < 3; i++)
		{
			const hvec& p = vlist[tri.ends[i].vi];
			
			//三角形每个点p都在包围盒里才返回真。

			for (int j = 0; j < 3; j++)
			{
				if (p.xyzw[j] < b[j][0] || p.xyzw[j] > b[j][1])
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
		printf((funcname + "已完成\n").c_str());
		printf("当前v数量%d vn数量%d f数量%d\n", vlist.size(), vnlist.size(), flist.size());
	}
	//将模型居中于坐标原点
	void centered()
	{
		//计算当前模型的包围盒
		boundingbox bb = boundingbox::get_bounding_box(vlist);
		//得到包围盒的中点
		vec vs((bb.b[0][0] - bb.b[0][1])/2, (bb.b[1][0] - bb.b[1][1]) / 2, (bb.b[2][0] - bb.b[2][1]) / 2);
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
		printf((funcname + " 已完成!\n").c_str());
		printf("当前v数量%d vn数量%d f数量%d\n", vlist.size(), vnlist.size(), flist.size());
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

		std::default_random_engine generator;
		std::default_random_engine generator2;

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
				coord[j] = coordpdf(generator2)*range[j];
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
};

//层次z_buffer
struct node
{
	//进入节点已知的数据。
	node* parent;
	//左上角坐标coord[0], 右下角坐标coord[1],
	vec2d coord[2];

	//进入节点未知待求解的数据。kids递归向下求解，z回溯求解。
	node* kids[2][2];
	double z;
};
struct zpyramid
{
	node* root;
	zpyramid(int a,int b,int c,int d)
	{
		root = new node();
		//左上角
		root->coord[0] = vec2d(a,b);
		//右下角
		root->coord[1] = vec2d(c,d);
	}
	//单个像素点位置映射到堆中的叶节点。
	map< vec2d , node* > mp;
	//递归向下建立pyramid结点，回溯构建z值小顶堆。返回值为最小z值。
	double make_up_heap(node* p)
	{
		double minz = DBL_MAX;

		int x[2][2];
		x[0][0] = (p->coord[0].xy[0] + p->coord[1].xy[0]) / 2;
		x[0][1] = x[0][0] + 1;
		x[1][0] = (p->coord[0].xy[1] + p->coord[1].xy[1]) / 2;
		x[1][1] = x[1][0] + 1;

		//左上角
		if (!(x[0][1] > p->coord[1].xy[0] && x[1][1] > p->coord[1].xy[1]))
		{
			p->kids[0][0] = new node();
			p->kids[0][0]->parent = p;
			p->kids[0][0]->coord[0] = p->coord[0];
			p->kids[0][0]->coord[1] = vec2d(x[0][0],x[1][0]);
			minz = fmin(minz, make_up_heap(p->kids[0][0]));
		}
		else
		{
			//此时当前节点为单个像素，不可进一步拆分。
			minz = DBL_MIN;
			mp[p->coord[0]] = p;
		}
		//右上角
		if (x[1][1] <= p->coord[1].xy[1])
		{
			p->kids[0][1] = new node();
			p->kids[0][1]->parent = p;
			p->kids[0][1]->coord[0] = vec2d(p->coord[0].xy[0],x[1][1]);
			p->kids[0][1]->coord[1] = vec2d(x[0][0], p->coord[1].xy[1]);
			minz = fmin(minz, make_up_heap(p->kids[0][1]));
		}
		//左下角
		if (x[0][1] <= p->coord[1].xy[0])
		{
			p->kids[1][0] = new node();
			p->kids[1][0]->parent = p;
			p->kids[1][0]->coord[0] = vec2d(x[0][1],p->coord[0].xy[1]);
			p->kids[1][0]->coord[1] = vec2d(p->coord[1].xy[0], x[1][0]);
			minz = fmin(minz, make_up_heap(p->kids[1][0]));
		}
		//右下角
		if (x[0][1] <= p->coord[1].xy[0] && x[1][1] <= p->coord[1].xy[1])
		{
			p->kids[1][1] = new node();
			p->kids[1][1]->parent = p;
			p->kids[1][1]->coord[0] = vec2d(x[0][1], x[1][1]);
			p->kids[1][1]->coord[1] = p->coord[1];
			minz = fmin(minz, make_up_heap(p->kids[1][1]));
		}
		return p->z = minz;
	}

	//更新单个像素的z_buffer值
	void update(const vec2d& pixel, double nz)
	{
		if (mp.find(pixel) == mp.end())
		{
			printf("z_pyramid中找不到该像素对应的节点\n");
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
			if (nmin < par->z)
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


};

struct pipeline
{
	//存储当前场景中所有顶点的齐次坐标
	vector<hvec> vlist;
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

	void print_transform(string funcname)
	{
		printf(("已完成"+funcname).c_str());
		//printf(" transform当前为：\n");
		//transform.print();

	}

	//取景变换，世界坐标变换到ouvn坐标,viewpoint->视点坐标，n->视点方向,up->向上的方向
	void view_transform(const vec& viewpoint, vec n, const vec& up,hmat& transform)
	{
		if (n == vec::get_zero_vector() || up == vec::get_zero_vector())
		{
			printf("视点方向和up方向不能为0\n");
			return;
		}
		n = n * (1.0 / n.norm2());
		vec v = n.cross_product(up);
		if (v == vec::get_zero_vector())
		{
			printf("视线方向与up方向不能在一条直线上\n");
			return;
		}
		v = v * (1.0/v.norm2());
		vec u = v.cross_product(n);

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

		for (int i = 0; i < 3; i++)
			T.A[i][3] = -viewpoint.xyz[i];

		//与现有变换复合
		transform = T * transform;

		print_transform("viewtransform");
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

	}

	//透视投影变换,d为投影平面，d>0
	void projection(hmat &transform)
	{
		if (d < 0)
		{
			printf("投影平面应当与视点为正距离");
			return;
		}

		hmat T = hmat::get_unity();
		T.A[3][2] = 1/d;
		T.A[3][3] = 0;
		transform = T * transform;

		this->d = d;

		print_transform("projection");
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
			p.normalize_x_y();
			vnlist.push_back(p);
		}

		//深度拷贝triangle面片
		flist.assign(s.flist.cbegin(), s.flist.cend());
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

		boundingbox frustum(-umax,umax,-vmax,vmax,front,back);

		for (list<triangle>::iterator it = flist.begin(); it != flist.end(); )
		{

			//此处考虑到裁剪算法的复杂性，只保留完全在视域四棱体内部的三角面片。未来可以将真正的三维裁剪在此处应用。
			if (!frustum.is_contains_triangle(*it,vlist))
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
				int vi = it->ends[i].vi;
				if (vmp.find(vi) == vmp.end())
				{
					nvlist.push_back(vlist[vi]);
					vmp[vi] = nvlist.size();
				}
				int nvi = vmp[vi];
				it->ends[i].vi = nvi;

				int vni = it->ends[i].vni;
				if (vnmp.find(vni) == vnmp.end())
				{
					nvnlist.push_back(vnlist[vni]);
					vnmp[vni] = nvnlist.size();
				}
				int nvni = vnmp[vni];
				it->ends[i].vni = nvni;
			}
		}

		vlist.assign(nvlist.begin(),nvlist.end());
		vnlist.assign(nvnlist.begin(), nvnlist.end());
	}   

	//ouvn变换到oxyz；oxyz的坐标原点放到投影平面上，与视点距离d。
	void ouvn_to_oxyz(hmat &transform)
	{
		hmat T;
		T.A[0][1] = T.A[1][0] = T.A[3][3] = 1;
		T.A[2][2] = -1;
		T.A[2][3] = d;
		transform = T * transform;
		xmax = vmax;
		ymax = umax;
		//小心精度误差。
		zmax = fmax(d - front,back - d);
		print_transform("ouvn_to_oxyz");
	}

	//视窗变换->变换到像素为单位的屏幕坐标系,屏幕坐标系原点在左上角。x向下为正，y向右为正。
	//屏幕的像素为xp*yp,例如1024*768
	void window_transform(int xp,int yp,hmat &transform)
	{
		if (xp % 2 && yp % 2)
		{
			printf("屏幕像素需要为偶数!\n");
			return;
		}
		this->xpixelnum = xp;
		this->ypixelnum = yp;
		//一半像素
		xp /= 2;
		yp /= 2;
		hmat T;
		//景物空间坐标直接与像素中心对齐，注意此处的坐标变换。
		T.A[0][1] = -yp / ymax;
		T.A[1][0] = xp / xmax;
		T.A[0][3] = yp - 0.5;
		T.A[1][3] = xp - 0.5;
		T.A[2][2] = T.A[3][3] = 1;

		transform = T * transform;
		//print_transform("window_transform");
	}

	//场景8叉树 + 层次z_buffer，开始！递归向下，一边建立8叉树，一边消隐。
	//在oxyz空间做8叉树
	void octree(const boundingbox &box,int threshold)
	{
		//检查当前节点是否被完全遮挡。
		//if()
	}
};

geometryfromobj soccerobj;
geometry soccer;

scene soccerfield;
pipeline mypipeline;
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	FILE* fp;

	AllocConsole();
	freopen_s(&fp, "CONIN$", "r", stdin);
	freopen_s(&fp, "CONOUT$", "w", stdout);
	freopen_s(&fp, "CONOUT$", "w", stderr);
	/*initgraph(640, 480);	// 创建绘图窗口，分辨率 640x480
	circle(320, 240, 100);	// 画圆，圆心 (320, 240)，半径 100
	Sleep(5000);			// 延时 5000 毫秒
	closegraph();			// 关闭图形窗口
	*/
	soccerobj.read_from_file();
	soccerobj.triangle_mesh();
	soccer.assign(soccerobj);
	soccer.centered();
	//soccer.print("");

	soccerfield.generate_from_base(soccer,1000);
	
	vec oo(1, 1, 1);
	vec nn(1, 1, 1);
	vec uu(1, 2, 1);
	//mypipeline.viewtransform(oo,nn,uu);
	//mypipeline.projection(5);

	//soccerfield.print("");

	system("pause");
	return 0;
}
