#include "mex.h" 
#include <math.h>
#include <iostream>
void mexFunction(int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[]) 
{ 
//    const EPS=1.0E-10;
//    const double EPS=1.0E-6;    
    double *fixPoint,*disPoint,*threshold; 
	double *outData; 
    int M0,N0,M1,N1; 
    //异常处理 
    //异常处理 
    if(nrhs!=3)
        mexErrMsgTxt("USAGE: b=reverse(a)/n"); 
    if(!mxIsDouble(prhs[0])) 
        mexErrMsgTxt("the Input Matrix must be double!/n"); 
    fixPoint=mxGetPr(prhs[0]); 
    disPoint=mxGetPr(prhs[1]);
    threshold=mxGetPr(prhs[2]);

    M0=mxGetM(prhs[0]); 
    N0=mxGetN(prhs[0]); 
    M1=mxGetM(prhs[1]); 
    N1=mxGetN(prhs[1]); 
    mexPrintf("fixPoints:%d * %d\n",M0,N0);
    mexPrintf("distributePoints:%d * %d\n",M1,N1);
    plhs[0]=mxCreateDoubleMatrix(M1,1,mxREAL); 
    //plhs[0]=mxCreateLogicalMatrix(M1,1);
    outData=mxGetPr(plhs[0]); 
	for (int i=0;i<M1;i++)
	{
		outData[i]=0;
	}
    for(int i=0;i<M0;i++) 
	{
        for(int j=0;j<M1;j++) 
        {
            if(outData[j]==1)
                continue;
            else
            {
                if(((fixPoint[i]-disPoint[j])*(fixPoint[i]-disPoint[j])+
					(fixPoint[i+M0]-disPoint[j+M1])*(fixPoint[i+M0]-disPoint[j+M1])+
					(fixPoint[i+2*M0]-disPoint[j+2*M1])*(fixPoint[i+2*M0]-disPoint[j+2*M1])<threshold[0]*threshold[0]))
                    outData[j]=1;            
			}
                    
        }
	}
	int temp=0;
	for (int i=0;i<M1;i++)
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
