层次zbuffer算法实现总说明及实验报告

一、总体说明

本实验采用c++（ISOC++20标准）作为编程语言，除c++STL标准库外，仅使用easyx将帧缓存数据输出到屏幕上。包括读入文件、向量计算在内的所有基础操作均由本人独立完成代码实现，保留著作权。编码环境：visual studio2019，版本管理：git，代码托管：https://github.com/luotong96/Hierarchical_Z_Buffer

二、数据结构说明

本实验数据结构繁多，主要采用了结构体，堆，数组等数据结构。数据传递大量使用引用类型。具体参见代码及其注释。

三、实验流程

a)   将soccerball.obj文件中的数据读入程序，得到1760个顶点和1992个多边形面片。

b)   对每个多边形，以其中一个顶点为基准，将多边形面片直接拆分成多个三角形。总共得到3516个三角形面片。

c)   将模型局部坐标系的原点移动到模型中心，使其三角形面片在8个象限中尽量均匀分布。

d)   在1000个随机位置（修改源代码647行，将时间作为随机生成器种子，可真正实现随机效果）复制基准模型，构建完整场景soccerfield。总顶点数量：1760000，总三角面片数量：3516000。

e)   定义矩阵transform = 单位矩阵。

f)    根据视点坐标基ouvn计算用于坐标变换的矩阵T，并将其复合至transform中。

g)   设置视域四棱台的参数（x,y,方向视角大小，棱台”前面”和”后面”的坐标），取

d=（”前面”+“后面“）/2来作为投影平面的位置。（即保证在透视投影中一半的内容缩小，一半的内容放大）。并计算视窗大小。

h)   根据投影平面的位置z=d，构造投影矩阵T，并将其复合至transform中。令transform=单位矩阵

i)    将transform应用至场景soccerfield中的所有顶点。完成取景变换和投影变换，并保留z值。所有新的顶点和面片转移到mypipeline数据结构中。

j)    裁剪掉视域四棱台以外的所有面片。考虑到裁剪算法的复杂性，本实验只保留完全在视域四棱体内部的三角面片，舍弃与四棱台相交的三角面片。最终得到**1020782**个三角面片需要处理。

k)   将坐标原点放至投影平面，z轴正方向为视线方向的反方向：构建ouvn到oxyz坐标的变换矩阵T，并将其复合至transform。

l)    视窗变换：将oxyz坐标系变换到屏幕坐标系。直接将最终坐标与屏幕像素中心点对齐。构建该变换的矩阵T，并将其复合至transform。此处变换的矩阵比较复杂易错。

m)  面消隐算法：

​         i.     设置图像分辨率1024*768

​        ii.     将trasnform应用至当前所有顶点。完成ouvn->oxyz和视窗变换的复合变换。新的顶点坐标存于v_in_window_list。

​        iii.     构建层次z-buffer数据结构zpyramid。其本质是一个小顶堆，建构于屏幕空间的像素点坐标上，在octree和scan_convert方法中都有使用。具体构建方法参考源代码915行起。

​        iv.     场景八叉树划分octree

1. 若当前box中三角面片小于threshold，直接扫描转换其所有面片。
2. 当前box投影到屏幕空间坐标系后，询问zpyramid是否被遮挡。若遮挡，则直接返回。
3. 构建8个子box，将当前box待处理的三角面片依次划分到其归属的子box中。要求三角面片完全在一个box内部，才认为其属于这个box。（本实验视域四棱台采用的裁剪方法也满足该条性质）
4. 扫描转换当前box中与任一分割面有交点的三角面片。
5. 按从离视点近到离视点远的顺序对子box递归处理。

​        v.     扫描转换scan_convert

1. 计算三角面片三个顶点的z值最大值。
2. 向zpyramid询问当前三角面片是否有像素可见——zpyramid.isvisible()函数

a)   判断三角面片与当前box的每个子box是否相交。对相交的box递归向下判断。相交判断采用重心坐标系，判断是否有box顶点在三角形内，或者三角形顶点在box内。

b)   若递归回溯时发现三角形与当前box相交的部分至少一个像素可见。则直接返回可见。

c)   否则，返回不可见。

3. 若三角面片不可见，跳过当前三角面片。
4. 否则，扫描转换当前三角形。

a)   先画三角形的上端点。

b)   扫描线从上往下扫描，得到xl,xr。只画完全在[xl,xr]内部的像素中心整点。（这个地方可能累积了不少浮点计算误差）

c)   三角形有两段，中间要换一条边，扫描线接着往下扫。具体见代码1500行起。

d)   三角形面片颜色取决于构建soccerfield时其所属模型的编号。

​        vi.     场景八叉树和层次zbuffer的实现原理参考于：

Greene N, Kass M, Miller G. Hierarchical Z-buffer visibility[C]//Proceedings of the 20th annual conference on Computer graphics and interactive techniques. 1993: 231-238.

四、实验结果

使用动态链接版本

1、 场景8叉树+层次Z-buffer

运行情况：场景8叉树跳过面片47953个，层次z_buffer跳过面片600239个

取景变换、投影变换、裁剪、消隐总用时49秒

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image002.png)

2、 不要场景8叉只要层次zbuffer

运行情况：场景8叉树跳过面片0个，层次z_buffer跳过面片694667个

取景变换、投影变换、裁剪、消隐总用时42秒

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image004.png)

3、 只要8叉序不要判断已被遮挡

运行情况：场景8叉树跳过面片0个，层次z_buffer跳过面片647327个

取景变换、投影变换、裁剪、消隐总用时48秒

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image006.png)

4、 绘制结果

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image008.jpg)

五、实验分析

由实验结果可知，层次z_buffer能够快速有效拒绝被遮挡的三角面片。但是本次实验的场景8叉树效果不佳。对比上面三项不同数据可以发现，使用场景8叉的计算时间花费已经超过他所带来的收益。仅使用层次z_buffer，虽然最终渲染的面片数量更多，但是总时间花费更小。可能的原因是本次实现的扫描转换不够精准，导致图像上总是有大量的散点被判断为镂空。所以一个8叉节点总是无法有效拒绝掉大量被遮挡的三角面片。因为该8叉节点的最小z值总是离视点无穷远，这是这些未被填充的镂空散点造成的。仅采用场景8叉序对场景进行扫描转换，可以拒绝约5万个面片，但时间花费仍然较大，高于收益。

​    仔细检查扫描转换的代码后，暂未发现明显的错误。这些镂空散点可能源于浮点运算误差。尝试在计算中加入eps偏移量，发现能够减少部分镂空点，但不能完解决。也可能是扫描转换代码本身不够完备，或者我未捕捉到基准模型本身的三角片特性，这还需要进一步的调试和优化。故本次实验的完整模式与简单模式的加速比<1。这一部分还望老师指点。

六、链接说明

使用静态链接的程序运行速度极快，4秒即可完成绘制。

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image010.png)

七、用户界面及使用说明：

源代码依赖于easyx库，可从EasyX_2023大暑版.exe安装。

运行时project.exe需和soccerball.obj放在同一个目录下，运行完毕后会在同一目录下生成绘制结果test.bmp。

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image012.png)

默认视点坐标系为：o(0, 0, 0)，n(1, 0, 0)，up(0, 0, 1)

![img](file:///C:/Users/luotong/AppData/Local/Temp/msohtmlclip1/01/clip_image014.png)

 

 
