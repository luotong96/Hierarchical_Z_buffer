#include <graphics.h>		// 引用 EasyX 绘图库头文件
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <map>
#include <vector>
#include <sstream>
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
struct vfromobj
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
};

struct geometry
{
	vector<vfromobj> vlist;
	vector<vnfromobj> vnlist;
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
					vfromobj v(x,y,z);
					vlist.push_back(v);
				}
				else if (b == "vn")
				{
					double x, y, z;
					if(sscanf_s(a.c_str(), "%lf%lf%lf", &x, &y, &z)!=3)
						break;
					vnfromobj vn(x, y, z);
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

					/*int tp = f.pointlist.size();
					if (tp == 5)
					{
						printf("linn = %d\n", linn);
					}
					if (mp.find(tp) == mp.end())
					{

						printf("size = %d\n",tp);
						mp[tp] = 0;
					}
					mp[tp]++;
					*/
					flist.push_back(f);
				}
				
			}
			//printf("facet = %d %d %d\n", mp[3], mp[4], mp[5]);
			fs.close();
		}
		printf("###%d %d %d########\n", vlist.size(),vnlist.size(),flist.size());
	}
};


/*void read_obj()
{

	map<string, bool> mp;
	mp.clear();
	fstream fs;
	//fs.open("..\\t3qqxibgic-CenterCitySciFi\\Center city Sci-Fi\\Center City Sci-Fi.obj", std::fstream::in);
	fs.open("..\\soccerball.obj", std::fstream::in);
	while (!fs.eof())
	{
		string a;
		getline(fs, a);
		string b;
		int pos = a.find_first_of(' ');
		b.assign(a.substr(0, pos));
		if (b != "g")
			continue;
		if (mp.find(a) == mp.end())
		{
			cout << a << endl;
			mp[a] = true;
		}
	}
	fs.close();
}*/
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
	//read_obj();
	soccer.read_from_file();

	system("pause");
	return 0;
}
