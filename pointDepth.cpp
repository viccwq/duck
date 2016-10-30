#include "mex.h" 
#include <math.h>
#include <iostream>
void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
//    const EPS=1.0E-10;
//    const double EPS=1.0E-6;    
    double *tree,*point; 
	double *depth; 
    int M0,N0,M1,N1; 
    //异常处理 
    //异常处理 
    if(nrhs!=2)
        mexErrMsgTxt("USAGE: b=reverse(a)/n"); 
    if(!mxIsDouble(prhs[0])) 
        mexErrMsgTxt("the Input Matrix must be double!/n"); 
    tree=mxGetPr(prhs[0]); 
    point=mxGetPr(prhs[1]);

    M0=mxGetM(prhs[0]); 
    N0=mxGetN(prhs[0]); 
    M1=mxGetM(prhs[1]); 
    N1=mxGetN(prhs[1]);
    mexPrintf("There are %4d leaves\n",M0);
    mexPrintf("There are %4d points\n",M1);

    plhs[0]=mxCreateDoubleMatrix(M1,1,mxREAL); 
    depth=mxGetPr(plhs[0]); 

    for(int j=0;j<M1;j++)	//points 
	{
		for(int i=0;i<M0;i++)		//leaves
        {
			if ((point[j]>=tree[i])&&
				(point[j]<tree[i+M0])&&
				(point[j+M1]>=tree[i+2*M0])&&
				(point[j+M1]<tree[i+3*M0])&&
				(point[j+2*M1]>=tree[i+4*M0])&&
				(point[j+2*M1]<tree[i+5*M0]))
			{
				depth[j]=tree[i+6*M0];
				break;                   
			}
			else
				depth[j]=0;
        }
	}

                       
} 
