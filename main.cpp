#include <graphics.h>		// 引用 EasyX 绘图库头文件
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <map>
#include <vector>
#include <sstream>
#include <cassert>
using namespace std;
//定义全局的向量，方便后续坐标变换的计算
struct vec
{
	double x, y, z;
	double eps = 1e-6;
	vec(double a = 0,double b = 0,double c = 0)
	{
		x = a; y = b; z = c;
	}
	//向量相等当且仅当各坐标相等
	bool operator == (const vec& b)
	{
		if (fabs(x - b.x) < eps && fabs(y - b.y) < eps && fabs(z - b.z) < eps)
			return true;
		return false;
	}
	vec operator +(const vec& b)
	{
		return vec(x + b.x,y + b.y,z + b.z);
	}
	vec operator -(const vec& b)
	{
		return vec(x - b.x, y - b.y, z - b.z);
	}
	//向量数乘c
	vec mulc(double c)
	{
		return vec(x * c, y * c, z * c);
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
	vector<facetfromobj> flist;

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
		int n = flist.size();
		for (int i = 0; i < n; i++)
		{
			//用新vector来装每个facet新产生的额外三角形面片
			vector<facetfromobj> nfacetlist;
			
			//当前facet的端点列表
			vector<indpair>& plst = flist[i].pointlist;
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

			//新产生的三角形放回原flist中
			for (int j = 0; j < nfacetlist.size(); j++)
			{
				flist.push_back(nfacetlist[j]);
			}
		}
		print("triangle_mesh");
	}
	//打印输出当前面数，点数
	void print(string funcname)
	{
		printf((funcname + "已完成\n").c_str());
		printf("当前v数量%d vn数量%d f数量%d\n", vlist.size(), vnlist.size(), flist.size());
		map<int, int>mp;
		printf("总面片数量=%d\n", flist.size());
		for (int i = 0; i < flist.size(); i++)
		{
			int cnt = flist[i].pointlist.size();
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
		xyzw[0] = b.x; xyzw[1] = b.y; xyzw[2] = b.z;
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
};

//
// 
struct boundingbox
{
	//AABB boundingbox,b[i][0]——第i维坐标最小值，b[i][1],第i维坐标最大值
	double b[3][2];
};
//用于单个模型的核心数据结构。
struct geometry
{
	vector<hvec> vlist;
	vector<hvec> vnlist;
	vector<triangle> flist;

	void assign(const geometryfromobj& b)
	{
		vlist.assign(b.vlist.begin(), b.vlist.end());
		vnlist.assign(b.vnlist.begin(), b.vnlist.end());
		flist.assign(b.flist.begin(), b.flist.end());
	}
	void print(string funcname)
	{
		printf((funcname + "已完成\n").c_str());
		printf("当前v数量%d vn数量%d f数量%d\n", vlist.size(), vnlist.size(), flist.size());
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
		memcpy(A,b.A,sizeof(A));
	}
	//矩阵乘向量
	hvec operator*(const hvec& x)
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
	hmat operator*(const hmat& b)
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
};

//场景造型，随机在世界坐标系中复制生成单个模型
struct scene {
	//指向基础模型的引用
	const geometry& base;
	
	vector<hvec> vlist;
	vector<hvec> vnlist;
	vector<triangle> flist;
	//使用仿射变换来复制模型
	void make_duplicates(int n)
	{
		//计算bounding 
		pow(n, 1.0 / 3);
	}
};



geometryfromobj soccerobj;
geometry soccer;
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
	soccer.print("");

	double mm[2][3];
	memset(mm, 0, sizeof(mm));
	for (int i = 0; i < soccer.vlist.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double tt = soccer.vlist[i].xyzw[j];
			if (tt > mm[0][j])
				mm[0][j] = tt;
			if (tt < mm[1][j])
				mm[1][j] = tt;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
			printf(" %f", mm[i][j]);
		printf("\n");
	}
	system("pause");
	return 0;
}
