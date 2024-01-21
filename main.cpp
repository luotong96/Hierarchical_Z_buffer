#include <graphics.h>		// 引用 EasyX 绘图库头文件

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	initgraph(640, 480);	// 创建绘图窗口，分辨率 640x480
	circle(320, 240, 100);	// 画圆，圆心 (320, 240)，半径 100
	Sleep(5000);			// 延时 5000 毫秒
	closegraph();			// 关闭图形窗口
	return 0;
}


