#include "mex.h" 
#include <math.h>
#include <iostream>

#define sqrt3 1.7321
#define zero 1e-7
#define threshold 0.45

int MCon=0, NCon=0,PCon=0;
int count=0;


bool location(double *depth,double x,double y,double z,double &Val)
{
	//基于质心的二维和三维插值算法

	//MCon----y, NCon----x,PCon----z
	int xmax=0,ymax=0,zmax=0;
	//double dist[8],ratio[8];
	int coordinate[24];
	//int tempDepth[8];

	const double deltax=1,deltay=1,deltaz=1;
	double x1=0,x2=0,y1=0,y2=0,z1=0,z2=0;
	double weight[8];	
	
	Val=0;

	if (x<1||x>=NCon||y<1||y>=MCon||z<1||z>=PCon)
	{
		Val=-1000;
		return false;
	} 
	else
	{
		//找到点（x，y，z）最近的八个整数点
		xmax=ceil(x);
		if (xmax==x)
			xmax++;
		ymax=ceil(y);
		if (ymax==y)
			ymax++;
		zmax=ceil(z);
		if (zmax==z)
			zmax++;
		
		
		/****************************基于质心的二维和三维插值算法*******************************/
		//坐标
		coordinate[0] =xmax-1;	coordinate[1] =ymax-1;	coordinate[2] =zmax-1;
		coordinate[3] =xmax;	coordinate[4] =ymax-1;	coordinate[5] =zmax-1;
		coordinate[6] =xmax;	coordinate[7] =ymax;	coordinate[8] =zmax-1;
		coordinate[9] =xmax-1;	coordinate[10]=ymax;	coordinate[11]=zmax-1;
		coordinate[12]=xmax-1;	coordinate[13]=ymax-1;	coordinate[14]=zmax;
		coordinate[15]=xmax;	coordinate[16]=ymax-1;	coordinate[17]=zmax;
		coordinate[18]=xmax;	coordinate[19]=ymax;	coordinate[20]=zmax;
		coordinate[21]=xmax-1;	coordinate[22]=ymax;	coordinate[23]=zmax;
		//前向差值
		x2=xmax-x;	x1=x-(xmax-1);
		y2=ymax-y;	y1=y-(ymax-1);
		z2=zmax-z;	z1=z-(zmax-1);
		weight[0]=(x2*y2*z2)/(deltax*deltay*deltaz);
		weight[1]=(x1*y2*z2)/(deltax*deltay*deltaz);
		weight[2]=(x1*y1*z2)/(deltax*deltay*deltaz);
		weight[3]=(x2*y1*z2)/(deltax*deltay*deltaz);
		weight[4]=(x2*y2*z1)/(deltax*deltay*deltaz);
		weight[5]=(x1*y2*z1)/(deltax*deltay*deltaz);
		weight[6]=(x1*y1*z1)/(deltax*deltay*deltaz);
		weight[7]=(x2*y1*z1)/(deltax*deltay*deltaz);
		
		for (int i=0;i<8;i++)
		{

			Val=Val+weight[i]*
				depth[(coordinate[1+i*3]-1)+
					  (coordinate[0+i*3]-1)*MCon+
					  (coordinate[2+i*3]-1)*MCon*NCon];
		}
		return true;


		/****************************计算距离推导权重*******************************/
/*
		coordinate[0] =xmax-1;	coordinate[1] =ymax-1;	coordinate[2] =zmax-1;
		coordinate[3] =xmax;	coordinate[4] =ymax-1;	coordinate[5] =zmax-1;
		coordinate[6] =xmax;	coordinate[7] =ymax;	coordinate[8] =zmax-1;
		coordinate[9] =xmax-1;	coordinate[10]=ymax;	coordinate[11]=zmax-1;
		coordinate[12]=xmax-1;	coordinate[13]=ymax-1;	coordinate[14]=zmax;
		coordinate[15]=xmax;	coordinate[16]=ymax-1;	coordinate[17]=zmax;
		coordinate[18]=xmax;	coordinate[19]=ymax;	coordinate[20]=zmax;
		coordinate[21]=xmax-1;	coordinate[22]=ymax;	coordinate[23]=zmax;
		//计算（x，y，z）到八个整数点的距离  和 八个点的值
		for (int i=0;i<8;i++)
		{		
			//八个整数点的距离
			dist[i]=(coordinate[0+i*3]-x)*(coordinate[0+i*3]-x)+
					(coordinate[1+i*3]-y)*(coordinate[1+i*3]-y)+
					(coordinate[2+i*3]-z)*(coordinate[2+i*3]-z);
			dist[i]=sqrt(dist[i]);
			if (dist[i]<threshold)
			{
				Val=depth[(coordinate[1+i*3]-1)+
						  (coordinate[0+i*3]-1)*MCon+
						  (coordinate[2+i*3]-1)*MCon*NCon];
				//mexPrintf("the test count: %d\n",++count);
				return true;
			}
			//八个点的值
			tempDepth[i]=depth[(coordinate[1+i*3]-1)+
							   (coordinate[0+i*3]-1)*MCon+
							   (coordinate[2+i*3]-1)*MCon*NCon];
		}
		//计算（x，y，z）到八个整数点的距离占总和的比例
		for (int i=0;i<8;i++)
		{
			ratio[i]=(sqrt3-dist[i])/sqrt3;
			Val=Val+tempDepth[i]*ratio[i];
		}
*/
	}

	

}
void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
	//    const EPS=1.0E-10;
	//    const double EPS=1.0E-6;    
	double *depth,*point,*size;
	double *valide,*sizeVal,Val=-1000; 
	int mCount=0,nCount=0,pCount=0;
	int Md,Nd,Mp,Np,M2,N2; 
	double x,y,z;
	count=0;
	//异常处理 
	//异常处理 
	if(nrhs!=3)
		mexErrMsgTxt("至少输入三个参数\n"); 
	if(!mxIsDouble(prhs[0])) 
		mexErrMsgTxt("the Input Matrix must be double!\n"); 
	depth=mxGetPr(prhs[0]);
	point=mxGetPr(prhs[1]);
	size=mxGetPr(prhs[2]); 

	Md=mxGetM(prhs[0]); 
	Nd=mxGetN(prhs[0]);	
	Mp=mxGetM(prhs[1]); 
	Np=mxGetN(prhs[1]);
	M2=mxGetM(prhs[2]); 
	N2=mxGetN(prhs[2]);
	mCount=(int)size[0];
	nCount=(int)size[1];
	pCount=(int)size[2];
	if ((Md*Nd!=mCount*nCount*pCount)&&(Md!=mCount))
	{
		mexErrMsgTxt("输入的depth的尺寸大小与size不一致！\n"); 
	}
	//mexPrintf("depth 的大小为 %4d*%4d*%4d \n",mCount,nCount,pCount);
	MCon=mCount;
	NCon=nCount;
	PCon=pCount;
	//mexPrintf("点的个数为 %4d \n",Mp);

	//plhs[0]=mxCreateDoubleMatrix(Mp,1,mxREAL); 
	//valide=mxGetPr(plhs[0]); 
	plhs[0]=mxCreateDoubleMatrix(Mp,1,mxREAL); 
	valide=mxGetPr(plhs[0]);
	plhs[1]=mxCreateDoubleMatrix(Mp,1,mxREAL); 
	sizeVal=mxGetPr(plhs[1]);

	/*for(int i=0;i<M2;i++)
	{
		x=point[i];
		y=point[i+M2];
		z=point[i+M2*2];
		location(scale,x,y,z,loc);
		val[i]=depth[(loc[2]-1)*M1*M1+(loc[0]-1)*M1+loc[1]-1];
	}*/
	for (int i=0;i<Mp;i++)
	{
		x=point[i+Mp*0];
		y=point[i+Mp*1];
		z=point[i+Mp*2];
		if (!location(depth,x,y,z,Val))
		{
			valide[i]=1;//标记出超出范围的点
		}		
		else valide[i]=0;
		sizeVal[i]=Val;
	}
	
} 
