#include "mex.h" 
#include <math.h>
#include <iostream>

int scaleCount=0;

void location(double *scale,double x,double y,double z,int *loc)
{
	double deltax=abs(scale[0]-x),deltax0=0;
	double deltay=abs(scale[0]-y),deltay0=0;
	double deltaz=abs(scale[0]-z),deltaz0=0;
	//X 计算在depth中的所在：列
	for(int i=1;i<scaleCount;i++)
	{	
		deltax0=abs(scale[i]-x);
		if (deltax>=deltax0)
		{
			deltax=deltax0;
		} 
		else
		{
			loc[0]=i;
			break;
		}

		if (i==scaleCount-1)
		{
			loc[0]=scaleCount;
			break;
		}
		
	}
	//Y 计算在depth中的所在：行
	for(int i=1;i<scaleCount;i++)
	{	
		deltay0=abs(scale[i]-y);
		if (deltay>=deltay0)
		{
			deltay=deltay0;
		} 
		else
		{
			loc[1]=i;
			break;
		}

		if (i==scaleCount-1)
		{
			loc[1]=scaleCount;
			break;
		}

	}
	//Z
	for(int i=1;i<scaleCount;i++)
	{	
		deltaz0=abs(scale[i]-z);
		if (deltaz>=deltaz0)
		{
			deltaz=deltaz0;
		} 
		else
		{
			loc[2]=i;
			break;
		}

		if (i==scaleCount-1)
		{
			loc[2]=scaleCount;
			break;
		}

	}	
}
void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
	//    const EPS=1.0E-10;
	//    const double EPS=1.0E-6;    
	double *scale,*depth,*point; 
	double *val; 
	double x=0,y=0,z=0;
	int M0,N0,M1,N1,M2,N2; 
	int loc[3]={0,0,0};
	//异常处理 
	//异常处理 
	if(nrhs!=3)
		mexErrMsgTxt("至少输入三个参数\n"); 
	if(!mxIsDouble(prhs[0])) 
		mexErrMsgTxt("the Input Matrix must be double!\n"); 
	scale=mxGetPr(prhs[0]); 
	depth=mxGetPr(prhs[1]);
	point=mxGetPr(prhs[2]);

	M0=mxGetM(prhs[0]); 
	N0=mxGetN(prhs[0]); 
	scaleCount=M0;
	M1=mxGetM(prhs[1]); 
	N1=mxGetN(prhs[1]);
	M2=mxGetM(prhs[2]); 
	N2=mxGetN(prhs[2]);

	if(M0!=M1)
		mexErrMsgTxt("输入参数大小不一致\n"); 
	mexPrintf("depth 的大小为 %4d*%4d*%4d \n",M0,M0,M0);
	mexPrintf("点的个数为 %4d \n",M2);

	plhs[0]=mxCreateDoubleMatrix(M2,1,mxREAL); 
	val=mxGetPr(plhs[0]); 

	for(int i=0;i<M2;i++)
	{
		x=point[i];
		y=point[i+M2];
		z=point[i+M2*2];
		location(scale,x,y,z,loc);
		val[i]=depth[(loc[2]-1)*M1*M1+(loc[0]-1)*M1+loc[1]-1];
	}

	
	
} 
