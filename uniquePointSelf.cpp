#include "mex.h" 
#include <math.h>
#include <iostream>
void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
//    const EPS=1.0E-10;
//    const double EPS=1.0E-6;    
    double *Point,*threshold; 
	double *outData; 
    int M0,N0,M1,N1; 
    //异常处理 
    //异常处理 
    if(nrhs!=2)
        mexErrMsgTxt("USAGE: b=reverse(a)/n"); 
    if(!mxIsDouble(prhs[0])) 
        mexErrMsgTxt("the Input Matrix must be double!/n"); 
    Point=mxGetPr(prhs[0]); 
    threshold=mxGetPr(prhs[1]);

    M0=mxGetM(prhs[0]); 
    N0=mxGetN(prhs[0]); 
    M1=mxGetM(prhs[1]); 
    N1=mxGetN(prhs[1]); 
    mexPrintf("Points:%d * %d\n",M0,N0);

    plhs[0]=mxCreateDoubleMatrix(M0,1,mxREAL); 
    //plhs[0]=mxCreateLogicalMatrix(M1,1);
    outData=mxGetPr(plhs[0]); 

	for (int i=0;i<M0;i++)
	{
		outData[i]=0;
	}

    for(int i=0;i<M0-1;i++) 
	{
		if(outData[i]==1)
		{
			continue;
		} 
		else
		{
			for(int j=i+1;j<M0;j++) 
			{
				if(outData[j]==1)
					continue;
				else
				{
					if(((Point[i]-Point[j])*(Point[i]-Point[j])+
						(Point[i+M0]-Point[j+M0])*(Point[i+M0]-Point[j+M0])+
						(Point[i+2*M0]-Point[j+2*M0])*(Point[i+2*M0]-Point[j+2*M0]))<threshold[0]*threshold[0])
						outData[j]=1;            

				}

			}
		}
	}

	int temp=0;
	for (int i=0;i<M0;i++)
	{
		temp=temp+outData[i];
	}
	mexPrintf("Have removed %d Points\n",temp);
		
/*     
    for(int i=0;i<M0;i++)
        for(int j=0;j<N0;j++)
        mexPrintf("fixPoint is :%f\n",fixPoint[i+M0*j]);
    for(int i=0;i<M1;i++)
        for(int j=0;j<N1;j++)
        mexPrintf("disPoint is :%f\n",disPoint[i+M1*j]);
*/
                       
} 
