#include <math.h>
#ifndef MAXDOUBLE
#define MAXDOUBLE 1e300
#endif

double gamma(double x)
{
    if(x<=0)return MAXDOUBLE;
    if(x>171)return MAXDOUBLE;
    double f=1.,z;
    while(x<7.0)
    {
        f*=x;
        x+=1.0;
    }
    z=1.0/(x*x);
    return 2.506628274631*exp(x*(log(x)-1.)+
        (((-0.000595238095238*z+0.000793650793651)*z-
        0.002777777777778)*z+0.083333333333333)/x)/
        (f*sqrt(x));
}

double polygamma(int n,double x)
{
    if((n<1)||(x<=0))return MAXDOUBLE;
    double f=0.,t,x2;
    int sfct=1,i;
    for(i=1;i<n;i++)sfct*=-i;
    t= n<4 ? 6.:9.*n;
    if(x<t){
        do{f+=pow(x,-n-1);x++;}while(x<t);
        f*=n;
    }
    x2=x*x;
    return sfct*(f+(((((( (2.08767569878681e-08-
        5.28419013868749e-10*(n+10.)*(n+11.)/x2)*
        (n+8.)*(n+9.)/x2-0.000000826719576719577)*
        (n+6.)*(n+7.)/x2+0.0000330687830687831)*
        (n+4.)*(n+5.)/x2-0.00138888888888889)*
        (n+2.)*(n+3.)/x2+0.0833333333333333)*
        (n+1.)/x+0.5)*n/x+1.)/pow(x,n));
}

double trigamma(double x){return polygamma(1,x);}
double pentagamma(double x){return polygamma(3,x);}

double x2aux(double x,int upper,double v,double vfac,int* ifail)
{
    const double acc=1e-15, big=1.3e154;
    double _x2aux;
    if(x==0.0)
    {
        _x2aux= upper ? 1.0 : 0.0;
        *ifail=0;
    }
    else
    {
        double gin,factor,a,b,term,rn,dif,pn[6];
        int i,undone;
        factor=v*log(x)-x-vfac;
        if(x>1.0&&x>=v)
        {
            a=2.0-v;b=a+x;
            term=pn[0]=1.0;pn[1]=x;pn[2]=x+1.0;pn[3]=x*b;
            gin=pn[2]/pn[3];undone=1;
            b+=2.0;
            while(undone)
            {
                pn[4]=b*pn[2]-a*term*pn[0];
                pn[5]=b*pn[3]-a*term*pn[1];
                if(pn[5]!=0)
                {
                    rn=pn[4]/pn[5];dif=fabs(gin-rn);
                    if(dif<=acc&&dif<=acc*rn)undone=0;
                    else gin=rn;
                }
                if(undone)
                {
                    for(i=0;i<4;i++)pn[i]=pn[i+2];
                    if(fabs(pn[4])>=big)
                        for(i=0;i<4;i++)pn[i]=pn[i]/big;
                    a+=1.0;term+=1.0;
                }
                b+=2.0;
            }
            gin=exp(factor+log(gin));
            *ifail= !upper&&gin>0.9 ? 5 : 0;
            _x2aux= upper ? gin : 1.0-gin;
        }
        else
        {
            gin=term=1.0;
            rn=v+1.0;
            do
            {
                term*=x/rn;gin+=term;
                rn+=1.0;
            }while(term>acc);
            gin=exp(factor)*gin/v;
            *ifail= upper&&gin>0.9 ? 5 : 0;
            _x2aux= upper ? 1.0-gin : gin;
        }
    }
    return _x2aux;
}

double igamma(double a,double x)
{
    if((a<=0)||(x<0))return MAXDOUBLE;
    int ifail;
    double gam=gamma(a);
    return x2aux(x,0,a,log(gam),&ifail)*gam;
}
