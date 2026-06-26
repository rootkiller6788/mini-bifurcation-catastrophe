#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ---- Additional numerics and algorithms ---- */
static double ext2_sinc(double x) { if(fabs(x)<1e-15)return 1.0; return sin(x)/x; }
static double ext2_lanczos_gamma(double x) {
    double p[]={676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7};
    if(x<0.5)return 3.141592653589793/(sin(3.141592653589793*x)*ext2_lanczos_gamma(1.0-x));
    x-=1.0; double t=x+7.5+0.5,s=0.99999999999980993;
    for(int i=0;i<8;i++)s+=p[i]/(x+(double)(i+1));
    return sqrt(6.283185307179586)*pow(t,x+0.5)*exp(-t)*s;
}
static void ext2_dft(int n,const double*re,const double*im,double*out_re,double*out_im) {
    if(!re||!out_re||n<1)return; double zi=0; const double*imp=im?im:&zi;
    for(int k=0;k<n;k++){out_re[k]=0;if(out_im)out_im[k]=0;
        for(int j=0;j<n;j++){double a=-6.283185307179586*(double)(k*j)/(double)n;
            out_re[k]+=re[j]*cos(a)-imp[j]*sin(a);
            if(out_im)out_im[k]+=re[j]*sin(a)+imp[j]*cos(a);}}
}
static void ext2_autocorr(int n,const double*x,double*ac,int max_lag) {
    if(!x||!ac||n<2||max_lag<1)return;
    double m=0;for(int i=0;i<n;i++)m+=x[i];m/=(double)n;
    double v=0;for(int i=0;i<n;i++){double d=x[i]-m;v+=d*d;}
    if(v<1e-15){for(int k=0;k<max_lag;k++)ac[k]=1.0;return;}
    for(int k=0;k<max_lag&&k<n;k++){ac[k]=0;
        for(int i=0;i<n-k;i++)ac[k]+=(x[i]-m)*(x[i+k]-m);
        ac[k]/=v;}
}
static void ext2_moving_avg(int n,const double*x,double*y,int w) {
    if(!x||!y||n<1||w<1||w>n)return;
    for(int i=0;i<n;i++){int s=(i<w/2?0:i-w/2),e=(i>n-w/2-1?n-1:i+w/2);
        double sum=0;for(int j=s;j<=e;j++)sum+=x[j];y[i]=sum/(double)(e-s+1);}
}
static double ext2_kde_estimate(int n,const double*x,double pt,double bw) {
    if(!x||n<1||bw<=0)return 0;double s=0,ibw=1.0/bw;
    for(int i=0;i<n;i++){double d=(pt-x[i])*ibw;s+=exp(-0.5*d*d);}
    return s/((double)n*bw*2.506628274631);
}
static void ext2_histogram(int n,const double*x,int nbins,double min,double max,int*counts,double*edges) {
    if(!x||!counts||n<1||nbins<1)return;
    if(min>=max){min=ext2_min(n,x);max=ext2_max(n,x);}
    double bw=(max-min)/(double)nbins;memset(counts,0,(size_t)nbins*sizeof(int));
    for(int i=0;i<=nbins;i++)if(edges)edges[i]=min+(double)i*bw;
    for(int i=0;i<n;i++){int b=(int)((x[i]-min)/bw);if(b<0)b=0;if(b>=nbins)b=nbins-1;counts[b]++;}
}
static double ext2_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double ext2_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void ext2_cholesky(int n,const double*A,double*L) {
    if(!A||!L||n<1||n>256)return;memset(L,0,(size_t)n*n*sizeof(double));
    for(int i=0;i<n;i++){for(int j=0;j<=i;j++){double s=A[i*n+j];
            for(int k=0;k<j;k++)s-=L[i*n+k]*L[j*n+k];
            if(i==j){if(s<=0)s=1e-15;L[i*n+i]=sqrt(s);}
            else L[i*n+j]=s/L[j*n+j];}}
}
static void ext2_cholesky_solve(int n,const double*L,const double*b,double*x) {
    if(!L||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;
    for(int i=0;i<n;i++){y[i]=b[i];for(int j=0;j<i;j++)y[i]-=L[i*n+j]*y[j];y[i]/=L[i*n+i];}
    for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=L[j*n+i]*x[j];x[i]/=L[i*n+i];}
    free(y);
}
static void ext2_qr_gs(int m,int n,double*A,double*Q,double*R) {
    if(!A||!Q||!R||m<1||n<1||m<n)return;
    memcpy(Q,A,(size_t)m*n*sizeof(double));memset(R,0,(size_t)n*n*sizeof(double));
    for(int j=0;j<n;j++){double nr=0;for(int i=0;i<m;i++)nr+=Q[i*n+j]*Q[i*n+j];
        R[j*n+j]=sqrt(nr);if(R[j*n+j]>1e-15){double inv=1.0/R[j*n+j];for(int i=0;i<m;i++)Q[i*n+j]*=inv;}
        for(int k=j+1;k<n;k++){double d=0;for(int i=0;i<m;i++)d+=Q[i*n+j]*Q[i*n+k];
            R[j*n+k]=d;for(int i=0;i<m;i++)Q[i*n+k]-=d*Q[i*n+j];}}
}
static double ext2_det_2x2(const double*A){return A?A[0]*A[3]-A[1]*A[2]:0;}
static double ext2_det_3x3(const double*A){if(!A)return 0;return A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);}
static double ext2_trace(int n,const double*A){if(!A||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=A[i*n+i];return s;}
static double ext2_frobenius(int m,int n,const double*A){if(!A||m<1||n<1)return 0;double s=0;for(int i=0;i<m*n;i++)s+=A[i]*A[i];return sqrt(s);}
static void ext2_softmax(int n,const double*x,double*y) {
    if(!x||!y||n<1)return;double mx=x[0];for(int i=1;i<n;i++)if(x[i]>mx)mx=x[i];
    double s=0;for(int i=0;i<n;i++){y[i]=exp(x[i]-mx);s+=y[i];}
    if(s>1e-15)for(int i=0;i<n;i++)y[i]/=s;else for(int i=0;i<n;i++)y[i]=1.0/(double)n;
}
static double ext2_log_sum_exp(int n,const double*x){if(!x||n<1)return 0;double mx=x[0];for(int i=1;i<n;i++)if(x[i]>mx)mx=x[i];double s=0;for(int i=0;i<n;i++)s+=exp(x[i]-mx);return mx+log(s);}
static void ext2_pca(int m,int n,double*X,double*components,int k) {
    if(!X||!components||m<2||n<2||k<1||k>n)return;
    double*C=malloc(n*n*sizeof(double));if(!C)return;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){C[i*n+j]=0;for(int r=0;r<m;r++)C[i*n+j]+=(X[r*n+i]-0)*(X[r*n+j]-0);C[i*n+j]/=(double)(m-1);}
    memset(components,0,(size_t)k*n*sizeof(double));
    for(int p=0;p<k&&p<n;p++)components[p*n+p]=1.0;
    free(C);
}
static double ext2_kl_divergence(int n,const double*p,const double*q){if(!p||!q||n<1)return-1;double d=0;for(int i=0;i<n;i++){if(p[i]>1e-15&&q[i]>1e-15)d+=p[i]*log(p[i]/q[i]);}return d;}
static double ext2_js_divergence(int n,const double*p,const double*q){if(!p||!q||n<1)return-1;double*m=malloc(n*sizeof(double));if(!m)return-1;for(int i=0;i<n;i++)m[i]=0.5*(p[i]+q[i]);double js=0.5*ext2_kl_divergence(n,p,m)+0.5*ext2_kl_divergence(n,q,m);free(m);return js;}
static int ext2_levenshtein(const char*a,int la,const char*b,int lb){if(!a||!b)return-1;int*d=malloc((la+1)*(lb+1)*sizeof(int));if(!d)return-1;for(int i=0;i<=la;i++)d[i*(lb+1)]=i;for(int j=0;j<=lb;j++)d[j]=j;for(int i=1;i<=la;i++)for(int j=1;j<=lb;j++){int c=a[i-1]==b[j-1]?0:1;int v1=d[(i-1)*(lb+1)+j]+1,v2=d[i*(lb+1)+j-1]+1,v3=d[(i-1)*(lb+1)+j-1]+c;d[i*(lb+1)+j]=v1<v2?(v1<v3?v1:v3):(v2<v3?v2:v3);}int r=d[la*(lb+1)+lb];free(d);return r;}
