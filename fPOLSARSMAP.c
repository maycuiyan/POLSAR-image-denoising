#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#define sqr(x) (x)*(x)
#define EPS 1e-16

void ENLCal(double T11[], double T12r[], double T12i[], double T13r[], double T13i[], double T22[], double T23r[], double T23i[], double T33[], double nu[], int mrows, int ncols, int r)
{
    int i,j;
    int ii,jj;
    int N;
    double *s11, *s12r, *s12i, *s13r, *s13i, *s22, *s23r, *s23i, *s33;
    double *ssqr11, *ssqr12r, *ssqr12i, *ssqr13r, *ssqr13i, *ssqr22, *ssqr23r, *ssqr23i, *ssqr33;
    double S11, S12r, S12i, S13r, S13i, S22, S23r, S23i, S33;
    double SSQR11, SSQR12r, SSQR12i, SSQR13r, SSQR13i, SSQR22, SSQR23r, SSQR23i, SSQR33;
    double tmp;
    
    
    s11 = (double*)calloc(mrows,sizeof(double));
    s12r = (double*)calloc(mrows,sizeof(double));
    s12i = (double*)calloc(mrows,sizeof(double));
    s13r = (double*)calloc(mrows,sizeof(double));
    s13i = (double*)calloc(mrows,sizeof(double));
    s22 = (double*)calloc(mrows,sizeof(double));
    s23r = (double*)calloc(mrows,sizeof(double));
    s23i = (double*)calloc(mrows,sizeof(double));
    s33 = (double*)calloc(mrows,sizeof(double));
    ssqr11 = (double*)calloc(mrows,sizeof(double));
    ssqr12r = (double*)calloc(mrows,sizeof(double));
    ssqr12i = (double*)calloc(mrows,sizeof(double));
    ssqr13r = (double*)calloc(mrows,sizeof(double));
    ssqr13i = (double*)calloc(mrows,sizeof(double));
    ssqr22 = (double*)calloc(mrows,sizeof(double));
    ssqr23r = (double*)calloc(mrows,sizeof(double));
    ssqr23i = (double*)calloc(mrows,sizeof(double));
    ssqr33 = (double*)calloc(mrows,sizeof(double));
    
    N = (2*r+1)*(2*r+1);
    
    
    
    for(j=0;j<=r;j++)
    {
        for(i=0;i<=r-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
 
            for(jj=0;jj<=j+r;jj++)
            {
                for(ii=0;ii<=i+r;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(j+r+1)*(i+r+1)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);
        }
        for(i=r;i<=mrows-r-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=0;jj<=j+r;jj++)
            {
                for(ii=i-r;ii<=i+r;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(j+r+1)*(2*r+1)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);
        }
        for(i=mrows-r;i<=mrows-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=0;jj<=j+r;jj++)
            {
                for(ii=i-r;ii<=mrows-1;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(j+r+1)*(mrows-i+r)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        } 
    }
    
    for(j=r+1;j<=ncols-r-1;j++)
    {
        for(i=0;i<=r-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=j-r;jj<=j+r;jj++)
            {
                for(ii=0;ii<=i+r;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(i+r+1)*(2*r+1)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);
        }
        
        i=r;
        S11 = s11[i];
        S12r = s12r[i];
        S12i = s12i[i];
        S13r = s13r[i];
        S13i = s13i[i];
        S22 = s22[i];
        S23r = s23r[i];
        S23i = s23i[i];
        S33 = s33[i];
        SSQR11 = ssqr11[i];
        SSQR12r = ssqr12r[i];
        SSQR12i = ssqr12i[i];
        SSQR13r = ssqr13r[i];
        SSQR13i = ssqr13i[i];
        SSQR22 = ssqr22[i];
        SSQR23r = ssqr23r[i];
        SSQR23i = ssqr23i[i];
        SSQR33 = ssqr33[i];
        
        
        s11[i] = 0;
        s12r[i] = 0;
        s12i[i] = 0;
        s13r[i] = 0;
        s13i[i] = 0;
        s22[i] = 0;
        s23r[i] = 0;
        s23i[i] = 0;
        s33[i] = 0;
        ssqr11[i] = 0;
        ssqr12r[i] = 0;
        ssqr12i[i] = 0;
        ssqr13r[i] = 0;
        ssqr13i[i] = 0;
        ssqr22[i] = 0;
        ssqr23r[i] = 0;
        ssqr23i[i] = 0;
        ssqr33[i] = 0;
        for(jj=j-r;jj<=j+r;jj++)
        {
            for(ii=0;ii<=i+r;ii++)
            {
                s11[i]=s11[i]+T11[ii+jj*mrows];
                s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                s22[i]=s22[i]+T22[ii+jj*mrows];
                s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                s33[i]=s33[i]+T33[ii+jj*mrows];
                ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
            }
        }
        nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(i+r+1)*(2*r+1)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        
        for(i=r+1;i<=mrows-r-1;i++)
        {
            
            tmp = S11; S11 = s11[i];
            s11[i]=s11[i-1]+S11-tmp-T11[(j-1-r)*mrows+(i+r)]+T11[(j-1-r)*mrows+(i-1-r)]-T11[(j+r)*mrows+(i-1-r)]+T11[(j+r)*mrows+(i+r)];
            tmp = S12r; S12r = s12r[i];
            s12r[i]=s12r[i-1]+S12r-tmp-T12r[(j-1-r)*mrows+(i+r)]+T12r[(j-1-r)*mrows+(i-1-r)]-T12r[(j+r)*mrows+(i-1-r)]+T12r[(j+r)*mrows+(i+r)];
            tmp = S12i; S12i = s12i[i];
            s12i[i]=s12i[i-1]+S12i-tmp-T12i[(j-1-r)*mrows+(i+r)]+T12i[(j-1-r)*mrows+(i-1-r)]-T12i[(j+r)*mrows+(i-1-r)]+T12i[(j+r)*mrows+(i+r)];
            tmp = S13r; S13r = s13r[i];
            s13r[i]=s13r[i-1]+S13r-tmp-T13r[(j-1-r)*mrows+(i+r)]+T13r[(j-1-r)*mrows+(i-1-r)]-T13r[(j+r)*mrows+(i-1-r)]+T13r[(j+r)*mrows+(i+r)];
            tmp = S13i; S13i = s13i[i];
            s13i[i]=s13i[i-1]+S13i-tmp-T13i[(j-1-r)*mrows+(i+r)]+T13i[(j-1-r)*mrows+(i-1-r)]-T13i[(j+r)*mrows+(i-1-r)]+T13i[(j+r)*mrows+(i+r)];
            tmp = S22; S22 = s22[i];
            s22[i]=s22[i-1]+S22-tmp-T22[(j-1-r)*mrows+(i+r)]+T22[(j-1-r)*mrows+(i-1-r)]-T22[(j+r)*mrows+(i-1-r)]+T22[(j+r)*mrows+(i+r)];
            tmp = S23r; S23r = s23r[i];
            s23r[i]=s23r[i-1]+S23r-tmp-T23r[(j-1-r)*mrows+(i+r)]+T23r[(j-1-r)*mrows+(i-1-r)]-T23r[(j+r)*mrows+(i-1-r)]+T23r[(j+r)*mrows+(i+r)];
            tmp = S23i; S23i = s23i[i];
            s23i[i]=s23i[i-1]+S23i-tmp-T23i[(j-1-r)*mrows+(i+r)]+T23i[(j-1-r)*mrows+(i-1-r)]-T23i[(j+r)*mrows+(i-1-r)]+T23i[(j+r)*mrows+(i+r)];
            tmp = S33; S33 = s33[i];
            s33[i]=s33[i-1]+S33-tmp-T33[(j-1-r)*mrows+(i+r)]+T33[(j-1-r)*mrows+(i-1-r)]-T33[(j+r)*mrows+(i-1-r)]+T33[(j+r)*mrows+(i+r)];
            
            tmp = SSQR11; SSQR11 = ssqr11[i];
            ssqr11[i]=ssqr11[i-1]+SSQR11-tmp-sqr(T11[(j-1-r)*mrows+(i+r)])+sqr(T11[(j-1-r)*mrows+(i-1-r)])-sqr(T11[(j+r)*mrows+(i-1-r)])+sqr(T11[(j+r)*mrows+(i+r)]);
            tmp = SSQR12r; SSQR12r = ssqr12r[i];
            ssqr12r[i]=ssqr12r[i-1]+SSQR12r-tmp-sqr(T12r[(j-1-r)*mrows+(i+r)])+sqr(T12r[(j-1-r)*mrows+(i-1-r)])-sqr(T12r[(j+r)*mrows+(i-1-r)])+sqr(T12r[(j+r)*mrows+(i+r)]);
            tmp = SSQR12i; SSQR12i = ssqr12i[i];
            ssqr12i[i]=ssqr12i[i-1]+SSQR12i-tmp-sqr(T12i[(j-1-r)*mrows+(i+r)])+sqr(T12i[(j-1-r)*mrows+(i-1-r)])-sqr(T12i[(j+r)*mrows+(i-1-r)])+sqr(T12i[(j+r)*mrows+(i+r)]);
            tmp = SSQR13r; SSQR13r = ssqr13r[i];
            ssqr13r[i]=ssqr13r[i-1]+SSQR13r-tmp-sqr(T13r[(j-1-r)*mrows+(i+r)])+sqr(T13r[(j-1-r)*mrows+(i-1-r)])-sqr(T13r[(j+r)*mrows+(i-1-r)])+sqr(T13r[(j+r)*mrows+(i+r)]);
            tmp = SSQR13i; SSQR13i = ssqr13i[i];
            ssqr13i[i]=ssqr13i[i-1]+SSQR13i-tmp-sqr(T13i[(j-1-r)*mrows+(i+r)])+sqr(T13i[(j-1-r)*mrows+(i-1-r)])-sqr(T13i[(j+r)*mrows+(i-1-r)])+sqr(T13i[(j+r)*mrows+(i+r)]);
            tmp = SSQR22; SSQR22 = ssqr22[i];
            ssqr22[i]=ssqr22[i-1]+SSQR22-tmp-sqr(T22[(j-1-r)*mrows+(i+r)])+sqr(T22[(j-1-r)*mrows+(i-1-r)])-sqr(T22[(j+r)*mrows+(i-1-r)])+sqr(T22[(j+r)*mrows+(i+r)]);
            tmp = SSQR23r; SSQR23r = ssqr23r[i];
            ssqr23r[i]=ssqr23r[i-1]+SSQR23r-tmp-sqr(T23r[(j-1-r)*mrows+(i+r)])+sqr(T23r[(j-1-r)*mrows+(i-1-r)])-sqr(T23r[(j+r)*mrows+(i-1-r)])+sqr(T23r[(j+r)*mrows+(i+r)]);
            tmp = SSQR23i; SSQR23i = ssqr23i[i];
            ssqr23i[i]=ssqr23i[i-1]+SSQR23i-tmp-sqr(T23i[(j-1-r)*mrows+(i+r)])+sqr(T23i[(j-1-r)*mrows+(i-1-r)])-sqr(T23i[(j+r)*mrows+(i-1-r)])+sqr(T23i[(j+r)*mrows+(i+r)]);
            tmp = SSQR33; SSQR33 = ssqr33[i];
            ssqr33[i]=ssqr33[i-1]+SSQR33-tmp-sqr(T33[(j-1-r)*mrows+(i+r)])+sqr(T33[(j-1-r)*mrows+(i-1-r)])-sqr(T33[(j+r)*mrows+(i-1-r)])+sqr(T33[(j+r)*mrows+(i+r)]);
            
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*N-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        }
        
        for(i=mrows-r;i<=mrows-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=j-r;jj<=j+r;jj++)
            {
                for(ii=i-r;ii<=mrows-1;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(mrows-i+r)*(2*r+1)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        }
    }
    
    for(j=ncols-r;j<=ncols-1;j++)
    {
        for(i=0;i<=r-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=j-r;jj<=ncols-1;jj++)
            {
                for(ii=0;ii<=i+r;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(i+r+1)*(ncols-j+r)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        }
        for(i=r;i<=mrows-r-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=j-r;jj<=ncols-1;jj++)
            {
                for(ii=i-r;ii<=i+r;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(2*r+1)*(ncols-j+r)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        }
        for(i=mrows-r;i<=mrows-1;i++)
        {
            s11[i] = 0;
            s12r[i] = 0;
            s12i[i] = 0;
            s13r[i] = 0;
            s13i[i] = 0;
            s22[i] = 0;
            s23r[i] = 0;
            s23i[i] = 0;
            s33[i] = 0;
            ssqr11[i] = 0;
            ssqr12r[i] = 0;
            ssqr12i[i] = 0;
            ssqr13r[i] = 0;
            ssqr13i[i] = 0;
            ssqr22[i] = 0;
            ssqr23r[i] = 0;
            ssqr23i[i] = 0;
            ssqr33[i] = 0;
            for(jj=j-r;jj<=ncols-1;jj++)
            {
                for(ii=i-r;ii<=mrows-1;ii++)
                {
                    s11[i]=s11[i]+T11[ii+jj*mrows];
                    s12r[i]=s12r[i]+T12r[ii+jj*mrows];
                    s12i[i]=s12i[i]+T12i[ii+jj*mrows];
                    s13r[i]=s13r[i]+T13r[ii+jj*mrows];
                    s13i[i]=s13i[i]+T13i[ii+jj*mrows];
                    s22[i]=s22[i]+T22[ii+jj*mrows];
                    s23r[i]=s23r[i]+T23r[ii+jj*mrows];
                    s23i[i]=s23i[i]+T23i[ii+jj*mrows];
                    s33[i]=s33[i]+T33[ii+jj*mrows];
                    ssqr11[i]=ssqr11[i]+sqr(T11[ii+jj*mrows]);
                    ssqr12r[i]=ssqr12r[i]+sqr(T12r[ii+jj*mrows]);
                    ssqr12i[i]=ssqr12i[i]+sqr(T12i[ii+jj*mrows]);
                    ssqr13r[i]=ssqr13r[i]+sqr(T13r[ii+jj*mrows]);
                    ssqr13i[i]=ssqr13i[i]+sqr(T13i[ii+jj*mrows]);
                    ssqr22[i]=ssqr22[i]+sqr(T22[ii+jj*mrows]);
                    ssqr23r[i]=ssqr23r[i]+sqr(T23r[ii+jj*mrows]);
                    ssqr23i[i]=ssqr23i[i]+sqr(T23i[ii+jj*mrows]);
                    ssqr33[i]=ssqr33[i]+sqr(T33[ii+jj*mrows]);
                }
            }
            nu[i+j*mrows] = sqr(s11[i]+s22[i]+s33[i])/((ssqr11[i]+ssqr22[i]+ssqr33[i]+2*(ssqr12r[i]+ssqr12i[i]+ssqr13r[i]+ssqr13i[i]+ssqr23r[i]+ssqr23i[i]))*(mrows-i+r)*(ncols-j+r)-(sqr(s11[i])+sqr(s22[i])+sqr(s33[i])+2*(sqr(s12r[i])+sqr(s12i[i])+sqr(s13r[i])+sqr(s13i[i])+sqr(s23r[i])+sqr(s23i[i])))+EPS);

        }
    }   
    
    free(s11);
    free(s12r);
    free(s12i);
    free(s13r);
    free(s13i);
    free(s23r);
    free(s23i);
    free(s33);
    free(ssqr11);
    free(ssqr12r);
    free(ssqr12i);
    free(ssqr13r);
    free(ssqr13i);
    free(ssqr23r);
    free(ssqr23i);
    free(ssqr33);
}


void POLSARSMAP(double T11[],double T12r[], double T12i[], double T13r[], double T13i[], double T22[], double T23r[], double T23i[], double T33[], double l_enl[], double Tf11[],double Tf12r[], double Tf12i[], double Tf13r[], double Tf13i[], double Tf22[], double Tf23r[], double Tf23i[], double Tf33[],double L, int mrows,int ncols, double kappa, double alpha)
{
    double *buf_t11, *buf_t12r, *buf_t12i, *buf_t13r, *buf_t13i, *buf_t22, *buf_t23r, *buf_t23i, *buf_t33;
    double x0_t11, x0_t12r, x0_t12i, x0_t13r, x0_t13i, x0_t22, x0_t23r, x0_t23i, x0_t33;
    double x1_t11, x1_t12r, x1_t12i, x1_t13r, x1_t13i, x1_t22, x1_t23r, x1_t23i, x1_t33;
    double nu;
    double wt;
    int i,j;
    buf_t11 = (double*)calloc(mrows,sizeof(double));
    buf_t12r = (double*)calloc(mrows,sizeof(double));
    buf_t12i = (double*)calloc(mrows,sizeof(double));
    buf_t13r = (double*)calloc(mrows,sizeof(double));
    buf_t13i = (double*)calloc(mrows,sizeof(double));
    buf_t22 = (double*)calloc(mrows,sizeof(double));
    buf_t23r = (double*)calloc(mrows,sizeof(double));
    buf_t23i = (double*)calloc(mrows,sizeof(double));
    buf_t33 = (double*)calloc(mrows,sizeof(double));
    // forward filtering
    j = 0;
    x0_t11 = T11[0];
    x0_t12r = T12r[0];
    x0_t12i = T12i[0];
    x0_t13r = T13r[0];
    x0_t13i = T13i[0];
    x0_t22 = T22[0];
    x0_t23r = T23r[0];
    x0_t23i = T23i[0];
    x0_t33 = T33[0];
    for (i=0;i<=mrows-1;i++)
    {
        nu = kappa*sqr(l_enl[j*mrows+i]/L);
        x1_t11 = L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
        Tf11[j*mrows+i] = x1_t11;
        x0_t11=x1_t11;
        //-------------
        x1_t12r = L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
        Tf12r[j*mrows+i] = x1_t12r;
        x0_t12r=x1_t12r;
        //-------------
        x1_t12i = L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
        Tf12i[j*mrows+i] = x1_t12i;
        x0_t12i=x1_t12i;
        //-------------
        x1_t13r = L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
        Tf13r[j*mrows+i] = x1_t13r;
        x0_t13r=x1_t13r;
        //-------------
        x1_t13i = L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
        Tf13i[j*mrows+i] = x1_t13i;
        x0_t13i=x1_t13i;
        //-------------
        x1_t22 = L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
        Tf22[j*mrows+i] = x1_t22;
        x0_t22=x1_t22;
        //-------------
        x1_t23r = L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
        Tf23r[j*mrows+i] = x1_t23r;
        x0_t13r=x1_t13r;
        //-------------
        x1_t23i = L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
        Tf23i[j*mrows+i] = x1_t23i;
        x0_t23i=x1_t23i;
        //-------------
        x1_t33 = L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
        Tf33[j*mrows+i] = x1_t33;
        x0_t33=x1_t33;
    }
    for (j=1;j<=ncols-1;j++)
    {
        if (j%2==0)
        {
            i=0;
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x0_t11=Tf11[(j-1)*mrows+i];
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            Tf11[j*mrows+i] = x1_t11;
            //--------------------------
            x0_t12r=Tf12r[(j-1)*mrows+i];
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            Tf12r[j*mrows+i] = x1_t12r;
            //--------------------------
            x0_t12i=Tf12i[(j-1)*mrows+i];
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            Tf12i[j*mrows+i] = x1_t12i;
            //--------------------------
            x0_t13r=Tf13r[(j-1)*mrows+i];
            x1_t13r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            Tf13r[j*mrows+i] = x1_t13r;
            //--------------------------
            x0_t13i=Tf13i[(j-1)*mrows+i];
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            Tf13i[j*mrows+i] = x1_t13i;
            //--------------------------
            x0_t22=Tf22[(j-1)*mrows+i];
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            Tf22[j*mrows+i] = x1_t22;
            //--------------------------
            x0_t23r=Tf23r[(j-1)*mrows+i];
            x1_t23r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            Tf23r[j*mrows+i] = x1_t23r;
            //--------------------------
            x0_t23i=Tf23i[(j-1)*mrows+i];
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            Tf23i[j*mrows+i] = x1_t23i;
            //--------------------------
            x0_t33=Tf33[(j-1)*mrows+i];
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            Tf33[j*mrows+i] = x1_t33; 
            
            for(i=1;i<=mrows-1;i++)
            {
                nu = kappa*sqr(l_enl[j*mrows+i]/L);
                x0_t11=alpha*Tf11[(j-1)*mrows+i]+(1-alpha)*Tf11[j*mrows+i-1];
                x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
                Tf11[j*mrows+i] = x1_t11;
                //----------------------------
                x0_t12r=alpha*Tf12r[(j-1)*mrows+i]+(1-alpha)*Tf12r[j*mrows+i-1];
                x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
                Tf12r[j*mrows+i] = x1_t12r;
                //----------------------------
                x0_t12i=alpha*Tf12i[(j-1)*mrows+i]+(1-alpha)*Tf12i[j*mrows+i-1];
                x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
                Tf12i[j*mrows+i] = x1_t12i;
                //----------------------------
                x0_t13r=alpha*Tf13r[(j-1)*mrows+i]+(1-alpha)*Tf13r[j*mrows+i-1];
                x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
                Tf13r[j*mrows+i] = x1_t13r;
                //----------------------------
                x0_t13i=alpha*Tf13i[(j-1)*mrows+i]+(1-alpha)*Tf13i[j*mrows+i-1];
                x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
                Tf13i[j*mrows+i] = x1_t13i;
                //----------------------------
                x0_t22=alpha*Tf22[(j-1)*mrows+i]+(1-alpha)*Tf22[j*mrows+i-1];
                x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
                Tf22[j*mrows+i] = x1_t22;
                //----------------------------
                x0_t23r=alpha*Tf23r[(j-1)*mrows+i]+(1-alpha)*Tf23r[j*mrows+i-1];
                x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
                Tf23r[j*mrows+i] = x1_t23r;
                //----------------------------
                x0_t23i=alpha*Tf23i[(j-1)*mrows+i]+(1-alpha)*Tf23i[j*mrows+i-1];
                x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
                Tf23i[j*mrows+i] = x1_t23i;
                //----------------------------
                x0_t33=alpha*Tf33[(j-1)*mrows+i]+(1-alpha)*Tf33[j*mrows+i-1];
                x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
                Tf33[j*mrows+i] = x1_t33;
            }
        }
        else
        {
            i=mrows-1;
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x0_t11=Tf11[(j-1)*mrows+i];
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            Tf11[j*mrows+i] = x1_t11;
            //------------------------------
            x0_t12r=Tf12r[(j-1)*mrows+i];
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            Tf12r[j*mrows+i] = x1_t12r;
            //------------------------------
            x0_t12i=Tf12i[(j-1)*mrows+i];
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            Tf12i[j*mrows+i] = x1_t12i;
            //------------------------------
            x0_t13r=Tf13r[(j-1)*mrows+i];
            x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            Tf13r[j*mrows+i] = x1_t13r;
            //------------------------------
            x0_t13i=Tf13i[(j-1)*mrows+i];
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            Tf13i[j*mrows+i] = x1_t13i;
            //-----------------------------
            x0_t22=Tf22[(j-1)*mrows+i];
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            Tf22[j*mrows+i] = x1_t22;
            //-----------------------------
            x0_t23r=Tf23r[(j-1)*mrows+i];
            x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            Tf23r[j*mrows+i] = x1_t23r;
            //------------------------------
            x0_t23i=Tf23i[(j-1)*mrows+i];
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            Tf23i[j*mrows+i] = x1_t23i;
            //------------------------------
            x0_t33=Tf33[(j-1)*mrows+i];
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            Tf33[j*mrows+i] = x1_t33;
            
            for(i=mrows-2;i>=0;i--)
            {
                nu = kappa*sqr(l_enl[j*mrows+i]/L);
                x0_t11 = alpha*Tf11[(j-1)*mrows+i]+(1-alpha)*Tf11[j*mrows+i+1];
                x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
                Tf11[j*mrows+i] = x1_t11;
                //----------------------------
                x0_t12r = alpha*Tf12r[(j-1)*mrows+i]+(1-alpha)*Tf12r[j*mrows+i+1];
                x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
                Tf12r[j*mrows+i] = x1_t12r;
                //----------------------------
                x0_t12i = alpha*Tf12i[(j-1)*mrows+i]+(1-alpha)*Tf12i[j*mrows+i+1];
                x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
                Tf12i[j*mrows+i] = x1_t12i;
                //----------------------------
                x0_t13r = alpha*Tf13r[(j-1)*mrows+i]+(1-alpha)*Tf13r[j*mrows+i+1];
                x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
                Tf13r[j*mrows+i] = x1_t13r;
                //----------------------------
                x0_t13i = alpha*Tf13i[(j-1)*mrows+i]+(1-alpha)*Tf13i[j*mrows+i+1];
                x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
                Tf13i[j*mrows+i] = x1_t13i;
                //----------------------------
                x0_t22 = alpha*Tf22[(j-1)*mrows+i]+(1-alpha)*Tf22[j*mrows+i+1];
                x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
                Tf22[j*mrows+i] = x1_t22;
                //----------------------------
                x0_t23r = alpha*Tf23r[(j-1)*mrows+i]+(1-alpha)*Tf23r[j*mrows+i+1];
                x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
                Tf23r[j*mrows+i] = x1_t23r;
                //----------------------------
                x0_t23i = alpha*Tf23i[(j-1)*mrows+i]+(1-alpha)*Tf23i[j*mrows+i+1];
                x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
                Tf23i[j*mrows+i] = x1_t23i;
                //----------------------------
                x0_t33 = alpha*Tf33[(j-1)*mrows+i]+(1-alpha)*Tf33[j*mrows+i+1];
                x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
                Tf33[j*mrows+i] = x1_t33;
            }
        }
    }
    
    
    // backward filtering
    j=ncols-1;
    if (ncols%2==0)
    {
        x0_t11 = Tf11[j*mrows];
        x0_t12r = Tf12r[j*mrows];
        x0_t12i = Tf12i[j*mrows];
        x0_t13r = Tf13r[j*mrows];
        x0_t13i = Tf13i[j*mrows];
        x0_t22 = Tf22[j*mrows];
        x0_t23r = Tf23r[j*mrows];
        x0_t23i = Tf23i[j*mrows];
        x0_t33 = Tf33[j*mrows];
        
        for (i=0;i<=mrows-1;i++)
        {
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            buf_t11[i] = x1_t11;
            x0_t11=x1_t11;
            //---------------------------------
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            buf_t12r[i] = x1_t12r;
            x0_t12r=x1_t12r;
            //---------------------------------
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            buf_t12i[i] = x1_t12i;
            x0_t12i=x1_t12i;
            //---------------------------------
            x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            buf_t13r[i] = x1_t13r;
            x0_t13r=x1_t13r;
            //---------------------------------
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            buf_t13i[i] = x1_t13i;
            x0_t13i=x1_t13i;
            //---------------------------------
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            buf_t22[i] = x1_t22;
            x0_t22=x1_t22;
            //---------------------------------
            x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            buf_t23r[i] = x1_t23r;
            x0_t23r=x1_t23r;
            //---------------------------------
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            buf_t23i[i] = x1_t23i;
            x0_t23i=x1_t23i;
            //---------------------------------
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            buf_t33[i] = x1_t33;
            x0_t33=x1_t33;
            //================================
            wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
            Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
            Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
            Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
            Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
            Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
            Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
            Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
            Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
            Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
        }
    }
    else
    {

        
        x0_t11 = Tf11[j*mrows+mrows-1];
        x0_t12r = Tf12r[j*mrows+mrows-1];
        x0_t12i = Tf12i[j*mrows+mrows-1];
        x0_t13r = Tf13r[j*mrows+mrows-1];
        x0_t13i = Tf13i[j*mrows+mrows-1];
        x0_t22 = Tf22[j*mrows+mrows-1];
        x0_t23r = Tf23r[j*mrows+mrows-1];
        x0_t23i = Tf23i[j*mrows+mrows-1];
        x0_t33 = Tf33[j*mrows+mrows-1];
        
        for (i=mrows-1;i>=0;i--)
        {
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            buf_t11[i] = x1_t11;
            x0_t11=x1_t11;
            //---------------------------------
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            buf_t12r[i] = x1_t12r;
            x0_t12r=x1_t12r;
            //---------------------------------
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            buf_t12i[i] = x1_t12i;
            x0_t12i=x1_t12i;
            //---------------------------------
            x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            buf_t13r[i] = x1_t13r;
            x0_t13r=x1_t13r;
            //---------------------------------
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            buf_t13i[i] = x1_t13i;
            x0_t13i=x1_t13i;
            //---------------------------------
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            buf_t22[i] = x1_t22;
            x0_t22=x1_t22;
            //---------------------------------
            x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            buf_t23r[i] = x1_t23r;
            x0_t23r=x1_t23r;
            //---------------------------------
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            buf_t23i[i] = x1_t23i;
            x0_t23i=x1_t23i;
            //---------------------------------
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            buf_t33[i] = x1_t33;
            x0_t33=x1_t33;
            //================================
            wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
            Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
            Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
            Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
            Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
            Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
            Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
            Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
            Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
            Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
        }
    }
    
    for(j=ncols-2;j>=0;j--)
    {
        if (j%2==1)
        {
            i=0;
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x0_t11=buf_t11[i];
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            buf_t11[i] = x1_t11;
            //---------------------------------
            x0_t12r=buf_t12r[i];
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            buf_t12r[i] = x1_t12r;
            //---------------------------------
            x0_t12i=buf_t12i[i];
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            buf_t12i[i] = x1_t12i;
            //---------------------------------
            x0_t13r=buf_t13r[i];
            x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            buf_t13r[i] = x1_t13r;
            //---------------------------------
            x0_t13i=buf_t13i[i];
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            buf_t13i[i] = x1_t13i;
            //---------------------------------
            x0_t22=buf_t22[i];
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            buf_t22[i] = x1_t22;
            //---------------------------------
            x0_t23r=buf_t23r[i];
            x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            buf_t23r[i] = x1_t23r;
            //---------------------------------
            x0_t23i=buf_t23i[i];
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            buf_t23i[i] = x1_t23i;
            //---------------------------------
            x0_t33=buf_t33[i];
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            buf_t33[i] = x1_t33;
            //================================
            wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
            Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
            Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
            Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
            Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
            Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
            Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
            Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
            Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
            Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
            
            for (i=1;i<mrows-1;i++)
            {
                nu = kappa*sqr(l_enl[j*mrows+i]/L);
                x0_t11=alpha*buf_t11[i]+(1-alpha)*buf_t11[i-1];
                x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
                buf_t11[i] = x1_t11;
                //-----------------------------
                x0_t12r=alpha*buf_t12r[i]+(1-alpha)*buf_t12r[i-1];
                x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
                buf_t12r[i] = x1_t12r;
                //-----------------------------
                x0_t12i=alpha*buf_t12i[i]+(1-alpha)*buf_t12i[i-1];
                x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
                buf_t12i[i] = x1_t12i;
                //-----------------------------
                x0_t13r=alpha*buf_t13r[i]+(1-alpha)*buf_t13r[i-1];
                x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
                buf_t13r[i] = x1_t13r;
                //-----------------------------
                x0_t13i=alpha*buf_t13i[i]+(1-alpha)*buf_t13i[i-1];
                x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
                buf_t13i[i] = x1_t13i;
                //-----------------------------
                x0_t22=alpha*buf_t22[i]+(1-alpha)*buf_t22[i-1];
                x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
                buf_t22[i] = x1_t22;
                //-----------------------------
                x0_t23r=alpha*buf_t23r[i]+(1-alpha)*buf_t23r[i-1];
                x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
                buf_t23r[i] = x1_t23r;
                //-----------------------------
                x0_t23i=alpha*buf_t23i[i]+(1-alpha)*buf_t23i[i-1];
                x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
                buf_t23i[i] = x1_t23i;
                //-----------------------------
                x0_t33=alpha*buf_t33[i]+(1-alpha)*buf_t33[i-1];
                x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
                buf_t33[i] = x1_t33;
                //================================
                wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
                Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
                Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
                Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
                Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
                Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
                Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
                Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
                Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
                Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
            }
        }
        else
        {
            i=mrows-1;
            nu = kappa*sqr(l_enl[j*mrows+i]/L);
            x0_t11=buf_t11[i];
            x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
            buf_t11[i] = x1_t11;
            //----------------------------
            x0_t12r=buf_t12r[i];
            x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
            buf_t12r[i] = x1_t12r;
            //----------------------------
            x0_t12i=buf_t12i[i];
            x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
            buf_t12i[i] = x1_t12i;
            //----------------------------
            x0_t13r=buf_t13r[i];
            x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
            buf_t13r[i] = x1_t13r;
            //----------------------------
            x0_t13i=buf_t13i[i];
            x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
            buf_t13i[i] = x1_t13i;
            //----------------------------
            x0_t22=buf_t22[i];
            x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
            buf_t22[i] = x1_t22;
            //----------------------------
            x0_t23r=buf_t23r[i];
            x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
            buf_t23r[i] = x1_t23r;
            //----------------------------
            x0_t23i=buf_t23i[i];
            x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
            buf_t23i[i] = x1_t23i;
            //----------------------------
            x0_t33=buf_t33[i];
            x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
            buf_t33[i] = x1_t33;
            //================================
            wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
            Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
            Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
            Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
            Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
            Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
            Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
            Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
            Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
            Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
            
            
            for (i=mrows-2;i>=0;i--)
            {
                nu = kappa*sqr(l_enl[j*mrows+i]/L);
                x0_t11=alpha*buf_t11[i]+(1-alpha)*buf_t11[i+1];
                x1_t11=L/(nu+L+2)*T11[j*mrows+i]+nu/(nu+L+2)*x0_t11;
                buf_t11[i] = x1_t11;
                //------------------------------
                x0_t12r=alpha*buf_t12r[i]+(1-alpha)*buf_t12r[i+1];
                x1_t12r=L/(nu+L+2)*T12r[j*mrows+i]+nu/(nu+L+2)*x0_t12r;
                buf_t12r[i] = x1_t12r;
                //------------------------------
                x0_t12i=alpha*buf_t12i[i]+(1-alpha)*buf_t12i[i+1];
                x1_t12i=L/(nu+L+2)*T12i[j*mrows+i]+nu/(nu+L+2)*x0_t12i;
                buf_t12i[i] = x1_t12i;
                //------------------------------
                x0_t13r=alpha*buf_t13r[i]+(1-alpha)*buf_t13r[i+1];
                x1_t13r=L/(nu+L+2)*T13r[j*mrows+i]+nu/(nu+L+2)*x0_t13r;
                buf_t13r[i] = x1_t13r;
                //------------------------------
                x0_t13i=alpha*buf_t13i[i]+(1-alpha)*buf_t13i[i+1];
                x1_t13i=L/(nu+L+2)*T13i[j*mrows+i]+nu/(nu+L+2)*x0_t13i;
                buf_t13i[i] = x1_t13i;
                //------------------------------
                x0_t22=alpha*buf_t22[i]+(1-alpha)*buf_t22[i+1];
                x1_t22=L/(nu+L+2)*T22[j*mrows+i]+nu/(nu+L+2)*x0_t22;
                buf_t22[i] = x1_t22;
                //------------------------------
                x0_t23r=alpha*buf_t23r[i]+(1-alpha)*buf_t23r[i+1];
                x1_t23r=L/(nu+L+2)*T23r[j*mrows+i]+nu/(nu+L+2)*x0_t23r;
                buf_t23r[i] = x1_t23r;
                //------------------------------
                x0_t23i=alpha*buf_t23i[i]+(1-alpha)*buf_t23i[i+1];
                x1_t23i=L/(nu+L+2)*T23i[j*mrows+i]+nu/(nu+L+2)*x0_t23i;
                buf_t23i[i] = x1_t23i;
                //------------------------------
                x0_t33=alpha*buf_t33[i]+(1-alpha)*buf_t33[i+1];
                x1_t33=L/(nu+L+2)*T33[j*mrows+i]+nu/(nu+L+2)*x0_t33;
                buf_t33[i] = x1_t33;
                //================================
                wt = (x1_t11+x1_t22+x1_t33)/(x1_t11+x1_t22+x1_t33+Tf11[j*mrows+i]+Tf22[j*mrows+i]+Tf33[j*mrows+i]+EPS);
                Tf11[j*mrows+i] = Tf11[j*mrows+i]*wt + x1_t11*(1-wt);
                Tf12r[j*mrows+i] = Tf12r[j*mrows+i]*wt + x1_t12r*(1-wt);
                Tf12i[j*mrows+i] = Tf12i[j*mrows+i]*wt + x1_t12i*(1-wt);
                Tf13r[j*mrows+i] = Tf13r[j*mrows+i]*wt + x1_t13r*(1-wt);
                Tf13i[j*mrows+i] = Tf13i[j*mrows+i]*wt + x1_t13i*(1-wt);
                Tf22[j*mrows+i] = Tf22[j*mrows+i]*wt + x1_t22*(1-wt);
                Tf23r[j*mrows+i] = Tf23r[j*mrows+i]*wt + x1_t23r*(1-wt);
                Tf23i[j*mrows+i] = Tf23i[j*mrows+i]*wt + x1_t23i*(1-wt);
                Tf33[j*mrows+i] = Tf33[j*mrows+i]*wt + x1_t33*(1-wt);
            }
        }
    }
    free(buf_t11);
    free(buf_t12r);
    free(buf_t12i);
    free(buf_t13r);
    free(buf_t13i);
    free(buf_t22);
    free(buf_t23r);
    free(buf_t23i);
    free(buf_t33);
}






// [Tf11,Tf12,Tf13,Tf22,Tf23,Tf33] = fPOLSARSMAP(T11,T12,T13,T22,T23,T33,L,r,kappa)
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    int mrows,ncols,r;
    double L,kappa,alpha;
    double *T11,*T12r,*T12i,*T13r,*T13i,*T22,*T23r,*T23i,*T33;
    double *l_enl;
    double *Tf11,*Tf12r,*Tf12i,*Tf13r,*Tf13i,*Tf22,*Tf23r,*Tf23i,*Tf33;
    double duration;
    clock_t start, finish;
    
    mrows=mxGetM(prhs[0]);
    ncols=mxGetN(prhs[0]);
    L=mxGetScalar(prhs[6]);
    r = mxGetScalar(prhs[7]);
    kappa = mxGetScalar(prhs[8]);
    alpha = 0.38;
    
    plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(mrows,ncols, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(mrows,ncols, mxCOMPLEX);
    plhs[3] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(mrows,ncols, mxCOMPLEX);
    plhs[5] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    
    
    T11=mxGetPr(prhs[0]);
    T12r=mxGetPr(prhs[1]);
    T12i=mxGetPi(prhs[1]);
    T13r=mxGetPr(prhs[2]);
    T13i=mxGetPi(prhs[2]);
    T22=mxGetPr(prhs[3]);
    T23r=mxGetPr(prhs[4]);
    T23i=mxGetPi(prhs[4]);
    T33=mxGetPr(prhs[5]);
    
    Tf11=mxGetPr(plhs[0]);
    Tf12r=mxGetPr(plhs[1]);
    Tf12i=mxGetPi(plhs[1]);
    Tf13r=mxGetPr(plhs[2]);
    Tf13i=mxGetPi(plhs[2]);
    Tf22=mxGetPr(plhs[3]);
    Tf23r=mxGetPr(plhs[4]);
    Tf23i=mxGetPi(plhs[4]);
    Tf33=mxGetPr(plhs[5]);
    
    
    start = clock();
    l_enl = (double*)calloc(mrows*ncols,sizeof(double));
    ENLCal(T11,T12r,T12i,T13r,T13i,T22,T23r,T23i,T33,l_enl,mrows,ncols,r);
    POLSARSMAP(T11,T12r,T12i,T13r,T13i,T22,T23r,T23i,T33,l_enl,Tf11,Tf12r,Tf12i,Tf13r,Tf13i,Tf22,Tf23r,Tf23i,Tf33,L, mrows,ncols,kappa,alpha);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );
    free(l_enl);         
}

