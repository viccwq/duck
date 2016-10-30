#include "mex.h" 
#include <math.h>
#include <iostream>

#define sqrt3 1.7321
#define zero 1e-7
#define threshold 0.45

int MCon=0,NCon=0,PCon=0;
double *Out_Nodes,*Norms;
int SizeNorms;


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
	}
}

bool findNormalIndex2(double *out_nodes,double *norms,int size,int x,int y,int z,double normalxyz[][3],int *enable)
{
	int size2=size*2;
	int size1=size*1;
	int size0=0;
	int xyz[8][3]={	{x,	 y,  z  },
					{x,  y,  z+1},
					{x,  y+1,z  },
					{x,  y+1,z+1},
					{x+1,y,  z  },
					{x+1,y,  z+1},
					{x+1,y+1,z  },
					{x+1,y+1,z+1}};
	/*XYZ[0][0]=x;		XYZ[0][1]=y;		XYZ[0][2]=z;
	XYZ[1][0]=x;		XYZ[1][1]=y;		XYZ[1][2]=z+1;
	XYZ[2][0]=x;		XYZ[2][1]=y+1;		XYZ[2][2]=z;
	XYZ[3][0]=x;		XYZ[3][1]=y+1;		XYZ[3][2]=z+1;
	XYZ[4][0]=x+1;		XYZ[4][1]=y;		XYZ[4][2]=z;
	XYZ[5][0]=x+1;		XYZ[5][1]=y;		XYZ[5][2]=z+1;
	XYZ[6][0]=x+1;		XYZ[6][1]=y+1;		XYZ[6][2]=z;
	XYZ[7][0]=x+1;		XYZ[7][1]=y+1;		XYZ[7][2]=z+1;*/

	int j=0;
	for (int i=0;i<size;i++)
	{
		if ((out_nodes[i+size0]==xyz[j][0])&&
			(out_nodes[i+size1]==xyz[j][1])&&
			(out_nodes[i+size2]==xyz[j][2]))
		{
			normalxyz[j][0]=norms[i+size0];
			normalxyz[j][1]=norms[i+size1];
			normalxyz[j][2]=norms[i+size2];
			enable[j]=1;
			j++;
			if (j==8)
			{
				return true;
			}
		} 
	}
	if (j<8)
	{
		mexErrMsgTxt("未能找到该点的周围八个点的法向量\n"); 
	}
	return false;
}
bool findNormalIndex(double *out_nodes,double *norms,int size,int *node,double *normal)
{
	//double *out_nodes,*norms,查询的点坐标和法向量库
	//int size,第一二参数的大小
	//int *node,待查询的点
	//double *normal返回的指针
	int i=0;
	normal[0]=0;
	normal[1]=0;
	normal[2]=0;

	for (i=0;i<size;i++)
	{
		if ((out_nodes[i]==node[0])&&
			(out_nodes[i+1*size]==node[1])&&
			(out_nodes[i+2*size]==node[2]))
		{
			normal[0]=norms[i];
			normal[1]=norms[i+1*size];
			normal[2]=norms[i+2*size];
			return true;
		} 
	}
	if (i==size)
	{
		return false;
	}
}
//Out_Nodes变成了顺序排列
bool averageNormal2(double end_x,double end_y,double end_z,double *normal)
{
	normal[0]=0;
	normal[1]=0;
	normal[2]=0;
	double normalxyz[8][3];//x1y1z1x2y2z2x3y3z3x4y4z4...
	for (int i=0;i<8;i++)
	{
		for (int j=0;j<3;j++)
		{
			normalxyz[i][j]=0;
		}
	}

	int enable[8]={0,0,0,0,0,0,0,0},sume=0;

	int x=floor(end_x);
	int y=floor(end_y);
	int z=floor(end_z);

	findNormalIndex2(Out_Nodes,Norms,SizeNorms,x,y,z,normalxyz,enable);
	for (int i=0;i<8;i++)
	{
		sume=enable[i]+sume;
	}


	if (sume==8)//采用插值法
	{
		const double deltax=1,deltay=1,deltaz=1;
		double x1=0,x2=0,y1=0,y2=0,z1=0,z2=0;
		double weight[8],weig=0;
		x1=end_x-x;	x2=x+1-end_x;
		y1=end_y-y;	y2=y+1-end_y;
		z1=end_z-z;	z2=z+1-end_z;
		//weight中点的顺序按照：
		//{x,  y,  z  },
		//{x,  y,  z+1},
		//{x,  y+1,z  },
		//{x,  y+1,z+1},
		//{x+1,y,  z  },
		//{x+1,y,  z+1},
		//{x+1,y+1,z  },
		//{x+1,y+1,z+1}
		weight[0]=(x2*y2*z2)/(deltax*deltay*deltaz);
		weight[1]=(x2*y2*z1)/(deltax*deltay*deltaz);
		weight[2]=(x2*y1*z2)/(deltax*deltay*deltaz);
		weight[3]=(x2*y1*z1)/(deltax*deltay*deltaz);
		weight[4]=(x1*y2*z2)/(deltax*deltay*deltaz);
		weight[5]=(x1*y2*z1)/(deltax*deltay*deltaz);
		weight[6]=(x1*y1*z2)/(deltax*deltay*deltaz);
		weight[7]=(x1*y1*z1)/(deltax*deltay*deltaz);
		
		for (int i=0;i<8;i++)
		{
			normal[0]=normal[0]+weight[i]*normalxyz[i][0];
			normal[1]=normal[1]+weight[i]*normalxyz[i][1];
			normal[2]=normal[2]+weight[i]*normalxyz[i][2];		
		}
		return true;
	} 
	else//采用平均法
	{
		mexErrMsgTxt("XXXXXXXXXXXXX\n"); 
		return false;
	}
}

bool averageNormal(int x,int y,int z,double *normal)
{
	double normalx=0,normaly=0,normalz=0;
	double normalxyz[3]={0,0,0};
	int node[3],count=0;
	//***1
	node[0]=x;		node[1]=y;		node[2]=z;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***2
	node[0]=x+1;	node[1]=y;		node[2]=z;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***3
	node[0]=x+1;	node[1]=y+1;	node[2]=z;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***4
	node[0]=x;		node[1]=y+1;	node[2]=z;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***5
	node[0]=x;		node[1]=y;		node[2]=z+1;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***6
	node[0]=x+1;	node[1]=y;		node[2]=z+1;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***7
	node[0]=x+1;	node[1]=y+1;	node[2]=z+1;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	//***8
	node[0]=x;		node[1]=y+1;	node[2]=z+1;
	if (findNormalIndex(Out_Nodes,Norms,SizeNorms,node,normalxyz)
		)
	{
		count++;
		normalx=normalx+normalxyz[0];
		normaly=normaly+normalxyz[1];
		normalz=normalz+normalxyz[2];
	}
	if (count!=0)
	{
		normal[0]=normalx/count;
		normal[1]=normaly/count;
		normal[2]=normalz/count;
		return true;
	} 
	else
	{
		return false;
	}
}


void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
	//    const EPS=1.0E-10;
	//    const double EPS=1.0E-6;    
	double *depth,*sizeDep,*end_node,*vec,*out_nodes,*norms;
	double *final_node;
	int mCount=0,nCount=0,pCount=0;

	int Md,Nd,Men,Nen,Mv,Nv,Mnor,Nnor; 
	double x,y,z;
	int CountRep=0;
	//异常处理 
	//异常处理 
	if(nrhs!=6)
		mexErrMsgTxt("至少输入六个参数\n"); 
	if(!mxIsDouble(prhs[0])) 
		mexErrMsgTxt("the Input Matrix must be double!\n"); 
	depth=mxGetPr(prhs[0]);
	sizeDep=mxGetPr(prhs[1]); 
	end_node=mxGetPr(prhs[2]);
	vec=mxGetPr(prhs[3]); 
	out_nodes=mxGetPr(prhs[4]); 
	norms=mxGetPr(prhs[5]); 
 

	Md=mxGetM(prhs[0]); 
	Nd=mxGetN(prhs[0]);	
	mCount=(int)sizeDep[0];
	nCount=(int)sizeDep[1];
	pCount=(int)sizeDep[2];
	if ((Md*Nd!=mCount*nCount*pCount)&&(Md!=mCount))
	{
		mexErrMsgTxt("输入的depth的尺寸大小与size不一致！\n"); 
	}
	//mexPrintf("depth 的大小为 %4d*%4d*%4d \n",mCount,nCount,pCount);


	Men=mxGetM(prhs[2]); 
	Nen=mxGetN(prhs[2]);
	Mv=mxGetM(prhs[3]); 
	Nv=mxGetN(prhs[3]);
	if ((Men!=Mv)&&(Nen!=Nv))
	{
		mexErrMsgTxt("输入的end_node和vec的大小不一致！\n"); 
	}
	int Men0=0;
	int Men1=Men*1;
	int Men2=Men*2;

	mexPrintf("共有%d个需要往回平移的点.\n",Men);
	Mnor=mxGetM(prhs[4]); 
	Nnor=mxGetN(prhs[4]);
	if ((Mnor!=mxGetM(prhs[5]))&&(Nnor!=mxGetN(prhs[5])))
	{
		mexErrMsgTxt("输入的out_nodes和norms的大小不一致！\n"); 
	}

	//初始化全局变量
	MCon=mCount;
	NCon=nCount;
	PCon=pCount;

	Out_Nodes=out_nodes;
	Norms=norms;
	SizeNorms=Mnor;
	//mexPrintf("点的个数为 %4d \n",Mp);


	plhs[0]=mxCreateDoubleMatrix(Men,3,mxREAL); 
	final_node=mxGetPr(plhs[0]);

//////////////////////////////////////////////////////////////////////////old
	//将节点整数化
/*
	int *nodex,*nodey,*nodez;
	nodex=new int[Men];
	nodey=new int[Men];
	nodez=new int[Men];
	for (int i=0;i<Men;i++)
	{
		nodex[i]=floor(end_node[i+0*Men]);
		nodey[i]=floor(end_node[i+1*Men]);
		nodez[i]=floor(end_node[i+2*Men]);
	}
	

	//获取每个点的平均法向量	
	double *normalx,*normaly,*normalz;
	normalx=new double[Men];
	normaly=new double[Men];
	normalz=new double[Men];
	for (int i=0;i<Men;i++)
	{		
		double normal[3];
		if (averageNormal(nodex[i],nodey[i],nodez[i],normal))
		{
			normalx[i]=normal[0];
			normaly[i]=normal[1];
			normalz[i]=normal[2];
		}
		else
		{
			mexErrMsgTxt("计算法向量错误\n"); 
		}
	}
*/
//////////////////////////////////////////////////////////////////////////new
	//获取每个点的平均法向量	
	double *normalx,*normaly,*normalz;
	normalx=new double[Men];
	normaly=new double[Men];
	normalz=new double[Men];
	for (int i=0;i<Men;i++)
	{		
		double normal[3];
		if (averageNormal2(end_node[i+Men0],end_node[i+Men1],end_node[i+Men2],normal))
		{
			normalx[i]=normal[0];
			normaly[i]=normal[1];
			normalz[i]=normal[2];
		}
		else
		{
			mexErrMsgTxt("计算法向量错误\n"); 
		}
	}

//////////////////////////////////////////////////////////////////////////


	

	//法向量方向纠正
	double *length=new double[Men];//近似平移回去的长度
	for (int i=0;i<Men;i++)
	{
		//求向量的数量积
		double cosSita,sinSita,modVec,modNorm,dot;//此时输入法向量的模为一
		dot=normalx[i]*vec[i+Men0]+
			normaly[i]*vec[i+Men1]+
			normalz[i]*vec[i+Men2];
		//modVec=sqrt(vec[i+Men0]*vec[i+Men0]+
		//			vec[i+Men1]*vec[i+Men1]+
		//			vec[i+Men2]*vec[i+Men2]);
		//if (modVec==0)
		//{
		//	mexPrintf("%d",i);
		//	mexErrMsgTxt("XXXXXXXXXXXXX222\n"); 
		//}
		//cosSita=dot/modVec;


		//modNorm=sqrt(normalx[i]*normalx[i]+
		//			normaly[i]*normaly[i]+
		//			normalz[i]*normalz[i]);
		//cosSita=dot/modNorm;

		cosSita=dot;
		if (cosSita>1)
		{
			mexPrintf("%d",i);
			mexErrMsgTxt("XXXXXXXXXXXXX333\n"); 
		}
//		if (cosSita>0)	
//		{
			normalx[i]=normalx[i]*(-1);
			normaly[i]=normaly[i]*(-1);
			normalz[i]=normalz[i]*(-1);
//		}
//		else
//		{
//			CountRep++;
//		}
		sinSita=sqrt(1-cosSita*cosSita);
		length[i]=sinSita*modVec;
	}
	mexPrintf("共有%d个点的法向量没有纠正.\n",CountRep);
	CountRep=0;
	

	//找到返回后的点
	for (int i=0;i<Men;i++)
	{
		//法向量的沿着,应该是递增的
		const int Le=6;
		double sizeVal=0,valTest,step;
		int flagTest=0;
		double x,y,z;
		x=end_node[i+Men0];
		y=end_node[i+Men1];
		z=end_node[i+Men2];
		step=length[i]/Le;


		location(depth,x,y,z,sizeVal);
		valTest=sizeVal;
		for (int j=0;j<Le;j++)//判断单调性
		{
			x=x+(j+1)*step*normalx[i];
			y=y+(j+1)*step*normaly[i];
			z=z+(j+1)*step*normalz[i];
			location(depth,x,y,z,sizeVal);
			if (sizeVal>=valTest)
			{
				flagTest++;
			}
			valTest=sizeVal;
		}
		if (flagTest<Le-2)
		{
			normalx[i]=normalx[i]*(-1);
			normaly[i]=normaly[i]*(-1);
			normalz[i]=normalz[i]*(-1);
			mexPrintf("第%d个点的法向量进行了反修正.\n",i+1);
		}

		//找到深度值为0的点
		double k=1;
		x=end_node[i+Men0]+k*normalx[i];
		y=end_node[i+Men1]+k*normaly[i];
		z=end_node[i+Men2]+k*normalz[i];		
		location(depth,x,y,z,sizeVal);
		while (abs(sizeVal)>0.05)
		{
			if (sizeVal>0)
			{
				k=k/2;
				x=x-k*normalx[i];
				y=y-k*normaly[i];
				z=z-k*normalz[i];		
			} 
			else
			{
				if (k!=1)
					k=k/2;
				x=x+k*normalx[i];
				y=y+k*normaly[i];
				z=z+k*normalz[i];	
			}
			location(depth,x,y,z,sizeVal);
            CountRep++;
            if (CountRep==100)
            {                
                mexPrintf("第%d个点无法找回其返回点.\n",i+1);
                mexErrMsgTxt("程序终断!!!!!!!!!!!!!!!!!!\n");
           }
		}
		//返回值
		final_node[i+Men0]=x;
		final_node[i+Men1]=y;
		final_node[i+Men2]=z;
		CountRep=0;
	}	
	
	//delete [] nodex;
	//delete [] nodey;
	//delete [] nodez;

	delete [] normalx;
	delete [] normaly;
	delete [] normalz;

	delete [] length;
	
} 