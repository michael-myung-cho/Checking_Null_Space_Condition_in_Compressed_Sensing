/* function [ portion ] = find_portion_3( InList_Portion, indices )
% find portion according to index1 and index2*/

#include "mex.h"
void find_portion_1(double *InList, double *Indices, double *alpha, int m, int n)
{
    int i;     
    for(i=0; i < m; i++)
    {
        if(InList[i] == Indices[0])
        {
            alpha[0] = InList[i+m];
        }
    }
}
void find_portion_2(double *InList, double *Indices, double *alpha, int m, int n)
{
    int i;     
    for(i=0; i < m; i++)
    {
        if((InList[i] == Indices[0]) && (InList[i+m] == Indices[1]))
        {
            alpha[0] = InList[i+m*2];
        }
    }
}
void find_portion_3(double *InList, double *Indices, double *alpha, int m, int n)
{
    int i;     
    for(i=0; i < m; i++)
    {
//        mexPrintf("[%d]%f, %f, %f - %f, %f, %f\n",i,*(InList + i),*(InList+i+m),*(InList+ i+2*m),*Indices,*(Indices+1),*(Indices+2));
        if((InList[i] == Indices[0]) && (InList[i+m] == Indices[1]) && (InList[i+2*m] == Indices[2]))
        {
            alpha[0] = InList[i+m*3];
        }
    }
}
void find_portion_4(double *InList, double *Indices, double *alpha, int m, int n)
{
    int i;     
    for(i=0; i < m; i++)
    {
//        mexPrintf("[%d]%f, %f, %f - %f, %f, %f\n",i,*(InList + i),*(InList+i+m),*(InList+ i+2*m),*Indices,*(Indices+1),*(Indices+2));
        if((InList[i] == Indices[0]) && (InList[i+m] == Indices[1]) && (InList[i+2*m] == Indices[2])&& (InList[i+3*m] == Indices[3]))
        {
            alpha[0] = InList[i+m*4];
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *InList;
    double *Indices;
    double *alpha;
    int L;
    int m;
    int n;

    if(nrhs != 3)
    {   
        mexErrMsgTxt("three input required");
    }else if ( nlhs > 1 )
    {
        mexErrMsgTxt("Too many output arguments");
    }
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    InList = mxGetPr(prhs[0]);
    Indices = mxGetPr(prhs[1]);
    L = mxGetScalar(prhs[2]);
    alpha = mxGetPr(plhs[0]);
    
    switch(L)
    {
        case 1:
            find_portion_1( InList, Indices, alpha, m, n); break;
        case 2:
            find_portion_2( InList, Indices, alpha, m, n); break;
        case 3:
            find_portion_3( InList, Indices, alpha, m, n); break;
        default:
            find_portion_4( InList, Indices, alpha, m, n); break;
    }
}

