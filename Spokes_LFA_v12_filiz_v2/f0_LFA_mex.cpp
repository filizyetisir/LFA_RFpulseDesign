


#include "mex.h"
#include "omp.h"

#include <complex>
#include <cmath>

using namespace std;

#define TWO_PI_GAMMA (2.675128976384781e+08);
#define NUM_THREADS 8


void f0_LFA_MZ_mex(double *f,double *df,double *b1s_re,double *b1s_im,double *grad,double *time_to_sinc_time,double *time_to_spoke,double *subpulse,
        double *b1maps_re,double *b1maps_im,double *target_LFA_re,double *target_LFA_im,
        int nnzp,int ntimes,int ncoils,int nspokes,double deltat,
        int comp_grad,int do_mls,int refoc_pulse);

void multiply_2spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> *a3,complex<double> *b3);

void multiply_3spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> a3,complex<double> b3,complex<double> *a4,complex<double> *b4);



void f0_LFA_MZ_mex(double *f,double *df,double *b1s_re,double *b1s_im,double *grad,double *time_to_sinc_time,double *time_to_spoke,double *subpulse,
        double *b1maps_re,double *b1maps_im,double *target_LFA_re,double *target_LFA_im,
        int nnzp,int ntimes,int ncoils,int nspokes,double deltat,
        int comp_grad,int do_mls,int refoc_pulse)
{
    
    // initialize objective function
    *f=0.0;

    complex<double> *a_forw,*b_forw,*as,*bs,*da_dre,*da_dim,*db_dre,*db_dim;
    complex<double> *datot_dre,*dbtot_dre,*datot_dim,*dbtot_dim;
            
    if(comp_grad==1){  // allocate some structures only if gradient is requested
        
        a_forw=new complex<double>[nnzp*ntimes];
        b_forw=new complex<double>[nnzp*ntimes];
        for(int i=0;i<nnzp;i++){  // only the first time point needs to be initialized
            a_forw[i]=complex<double>(1.0,0.0);
            b_forw[i]=complex<double>(0.0,0.0);
        }
        
        // Q matrices
        as=new complex<double>[nnzp*ntimes];
        bs=new complex<double>[nnzp*ntimes];                    
        
        // derivatives of Q matrices
        da_dre=new complex<double>[nnzp*ntimes*ncoils];
        da_dim=new complex<double>[nnzp*ntimes*ncoils];
        db_dre=new complex<double>[nnzp*ntimes*ncoils];
        db_dim=new complex<double>[nnzp*ntimes*ncoils];
        
        // derivatives of the total Q
        datot_dre=new complex<double>[nnzp*ncoils*nspokes];
        datot_dim=new complex<double>[nnzp*ncoils*nspokes];
        dbtot_dre=new complex<double>[nnzp*ncoils*nspokes];
        dbtot_dim=new complex<double>[nnzp*ncoils*nspokes];
        
        for(int i=0;i<nnzp*ncoils*nspokes;i++)
            datot_dre[i]=dbtot_dre[i]=datot_dim[i]=dbtot_dim[i]=complex<double>(0.0,0.0);
    }
    
    complex<double> *a=new complex<double>[nnzp];
    complex<double> *b=new complex<double>[nnzp];
    for(int i=0;i<nnzp;i++){
        a[i]=complex<double>(1.0,0.0);
        b[i]=complex<double>(0.0,0.0);
    }

    // FORWARD BLOCH SIMULATION    
    complex<double> pure_imag=complex<double>(0.0,1.0);
    complex<double> *y=new complex<double>[nnzp];

#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i=0;i<nnzp;i++){
        
        for(int time=0;time<ntimes;time++){
            
            int sinctime=(int)( time_to_sinc_time[time] );      
            
            double norm=pow(b1s_re[time*nnzp+i],2) + pow(b1s_im[time*nnzp+i],2) + pow(abs(grad[time*nnzp+i]),2);
            norm=sqrt(norm);
            
            complex<double> a2,b2;
            double n1,n2,n3,phi,cosphi,sinphi;
            if(norm==0){
                a2=1.0;
                b2=0.0;
            }else{
                phi=-deltat*norm*TWO_PI_GAMMA;
                
                n1=b1s_re[time*nnzp+i] / norm;
                n2=b1s_im[time*nnzp+i] / norm;
                n3=grad[time*nnzp+i] / norm;    
                
                // new CK parameters
                cosphi=cos(phi/2.0);
                sinphi=sin(phi/2.0);
                a2=cosphi - pure_imag*n3*sinphi;
                b2=-pure_imag*( n1+pure_imag*n2 )*sinphi;              
                
                // forward Q product
                multiply_2spinors(a[i],b[i],a2,b2,a+i,b+i);
            }
                       
            if(comp_grad==1){
                
                // store Q matrix
                as[i+time*nnzp]=a2;
                bs[i+time*nnzp]=b2;        
                
                // store one time late forward Q matrix product
                if(time<ntimes-1){                    
                    a_forw[i+(time+1)*nnzp]=a[i];
                    b_forw[i+(time+1)*nnzp]=b[i];
                }
                
                if(sinctime>0){
                    
                    complex<double> b1s_=complex<double>( b1s_re[time*nnzp+i],b1s_im[time*nnzp+i] );
                    for(int j=0;j<ncoils;j++){
                        
                        if(norm>0){                    
                            
                            complex<double> b1maps_=complex<double>( b1maps_re[i+j*nnzp],b1maps_im[i+j*nnzp] );                          
                            
                            // derivative of phi wrt spokes amplitudes
                            complex<double> dphi=-deltat * conj(b1s_)/norm * subpulse[sinctime-1] * b1maps_ * TWO_PI_GAMMA;
                            double dphi_dre=real(dphi);
                            double dphi_dim=-imag(dphi);

                            // derivative of n (rotation axis) wrt spokes amplitudes
                            complex<double> tmp=-conj(b1s_)/pow(norm,3) * subpulse[sinctime-1] * b1maps_;
                            double tmp3=subpulse[sinctime-1]/norm;            
                            
                            double dnx_dre=real(tmp)*real(b1s_)  + tmp3*real(b1maps_);
                            double dnx_dim=-imag(tmp)*real(b1s_) - tmp3*imag(b1maps_);
                            
                            double dny_dre=real(tmp)*imag(b1s_) + tmp3*imag(b1maps_);
                            double dny_dim=-imag(tmp)*imag(b1s_) + tmp3*real(b1maps_);
                            
                            complex<double> dnz=tmp * grad[time*nnzp+i];
                            double dnz_dre=real(dnz);
                            double dnz_dim=-imag(dnz);                                                   
                            
                            // derivative of Q 
                            da_dre[i + time*nnzp + j*ntimes*nnzp]=dphi_dre/2*( -sinphi - pure_imag*n3*cosphi ) - pure_imag*sinphi*dnz_dre;
                            da_dim[i + time*nnzp + j*ntimes*nnzp]=dphi_dim/2*( -sinphi - pure_imag*n3*cosphi ) - pure_imag*sinphi*dnz_dim;            
                            
                            db_dre[i + time*nnzp + j*ntimes*nnzp]=-pure_imag*( dnx_dre + pure_imag*dny_dre ) * sinphi - pure_imag/2.0*( n1+pure_imag*n2 )*cosphi*dphi_dre;                                
                            db_dim[i + time*nnzp + j*ntimes*nnzp]=-pure_imag*( dnx_dim + pure_imag*dny_dim ) * sinphi - pure_imag/2.0*( n1+pure_imag*n2 )*cosphi*dphi_dim;
                        }
                        
                    }  // for(j=0;j<ncoils;j++)
                    
                }  // if(sinctime>0)
                
            }  // if(comp_grad==1)            

        }  // for(time=0;times<ntimes;time++)   
        
    }  // for(i=0;i<nnzp;i++) 
    
    
    
    // objective function
    complex<double> opt_metric;
    for(int i=0;i<nnzp;i++){
        
        if(refoc_pulse==0){
            if(do_mls==1)
                opt_metric=pow(abs(a[i]),2) - pow(abs(b[i]),2);  // mz=|alpha|^2-|beta|^2
            else
                opt_metric=2.0*conj(a[i])*b[i];  // 2*conj(alpha)*beta
        }else{
            if(do_mls==1)
                opt_metric=pow(abs(b[i]),2);  // |beta|^2
            else
                opt_metric=pow(b[i],2);  // beta^2
        }
        
        y[i]=opt_metric - complex<double>(target_LFA_re[i],target_LFA_im[i]);
        *f += pow( abs(y[i]),2 );
    }
    

    // BACKWARD BLOCH SIMULATION        
    if(comp_grad==1){
        
        for(int i=0;i<nnzp;i++){
            a[i]=complex<double>(1.0,0.0);
            b[i]=complex<double>(0.0,0.0);
        }
        
#pragma omp parallel for num_threads(NUM_THREADS)
        for(int i=0;i<nnzp;i++){
            
            for(int time=ntimes-1;time>=0;time--){            
                
                int sp_num=(int)( time_to_spoke[time] );
                            
                complex<double> tmp3,tmp4;                
                if(sp_num>0){               
                    
                    for(int j=0;j<ncoils;j++){
                    
                        // update derivatives of total Q matrices (sandwitch product)
                        multiply_3spinors(a_forw[i+time*nnzp],b_forw[i+time*nnzp],da_dre[i + time*nnzp + j*ntimes*nnzp],db_dre[i + time*nnzp + j*ntimes*nnzp],a[i],b[i],&tmp3,&tmp4);                           
                        datot_dre[i + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dre[i + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3;
                        dbtot_dre[i + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dre[i + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4;            
                        
                        multiply_3spinors(a_forw[i+time*nnzp],b_forw[i+time*nnzp],da_dim[i + time*nnzp + j*ntimes*nnzp],db_dim[i + time*nnzp + j*ntimes*nnzp],a[i],b[i],&tmp3,&tmp4);
                        datot_dim[i + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dim[i + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3;
                        dbtot_dim[i + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dim[i + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4;      
                        
                    }  // for(j=0;j<ncoils;j++)
                    
                }  // if(sp_num>0)

                // keep track of the backward Q matrices product
                multiply_2spinors( as[i+time*nnzp],bs[i+time*nnzp],a[i],b[i],a+i,b+i );
                
            }  // for(time=ntimes-1;time>=0;time--)
            
        }  // for(i=0;i<nnzp;i++)
        
       
        // form gradient out of the previously computed structures
         for(int j=0;j<nspokes;j++){
             
            for(int k=0;k<ncoils;k++){                
                df[k+j*ncoils]=df[k+j*ncoils + ncoils*nspokes]=0.0;                
                
                for(int i=0;i<nnzp;i++){                                
                    
                    if(refoc_pulse==0){
                        if(do_mls==1){
                            df[k+j*ncoils] += 4.0 * real( datot_dre[i+k*nnzp+j*nnzp*ncoils]*conj(a[i]) - dbtot_dre[i+k*nnzp+j*nnzp*ncoils]*conj(b[i]) ) * real(y[i]);
                            df[k+j*ncoils + ncoils*nspokes] += 4.0 * real( datot_dim[i+k*nnzp+j*nnzp*ncoils]*conj(a[i]) - dbtot_dim[i+k*nnzp+j*nnzp*ncoils]*conj(b[i]) ) * real(y[i]);
                        }else{
                            df[k+j*ncoils] += 4.0 * real( ( conj(datot_dre[i+k*nnzp+j*nnzp*ncoils])*b[i] + dbtot_dre[i+k*nnzp+j*nnzp*ncoils]*conj(a[i]) )*conj(y[i]) );
                            df[k+j*ncoils + ncoils*nspokes] += 4.0 * real( ( conj(datot_dim[i+k*nnzp+j*nnzp*ncoils])*b[i] + dbtot_dim[i+k*nnzp+j*nnzp*ncoils]*conj(a[i]) )*conj(y[i]) );
                        }
                    }else{                    
                        if(do_mls==1){
                            df[k+j*ncoils] += 4.0 * real( dbtot_dre[i+k*nnzp+j*nnzp*ncoils]*conj(b[i]) ) * real(y[i]);
                            df[k+j*ncoils + ncoils*nspokes] += 4.0 * real( dbtot_dim[i+k*nnzp+j*nnzp*ncoils]*conj(b[i]) ) * real(y[i]);
                        }else{
                            df[k+j*ncoils] += 4.0 * real( dbtot_dre[i+k*nnzp+j*nnzp*ncoils]*b[i]*conj(y[i]) );
                            df[k+j*ncoils + ncoils*nspokes] += 4.0 * real( dbtot_dim[i+k*nnzp+j*nnzp*ncoils]*b[i]*conj(y[i]) );
                        }                        
                    }
                    
                }
            }            
        }
        
    }  // if(comp_grad==1)
    
    
    delete [] a;
    delete [] b;    
    delete [] y;
    
    if(comp_grad==1){
        delete [] a_forw;
        delete [] b_forw;
        
        delete [] as;
        delete [] bs;
        
        delete [] da_dre;
        delete [] da_dim;
        delete [] db_dre;
        delete [] db_dim;
        
        delete [] datot_dre;
        delete [] datot_dim;
        delete [] dbtot_dre;
        delete [] dbtot_dim;
    }    
    
}





void multiply_2spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> *a3,complex<double> *b3)
{
    *a3=a2*a1 - conj(b2)*b1;
    *b3=b2*a1 + conj(a2)*b1;
}




void multiply_3spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> a3,complex<double> b3,complex<double> *a4,complex<double> *b4)
{
    complex<double> tmp1=a2*a1 - conj(b2)*b1;
    complex<double> tmp2=b2*a1 + conj(a2)*b1;
    
    *a4=a3*tmp1 - conj(b3)*tmp2;
    *b4=b3*tmp1 + conj(a3)*tmp2;
}



void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    if(nrhs!=12){
        printf("Usage:\n\t[f df]=function f0_LFA_MZ_mex(b1s,grad,time_to_sinc_time,time_to_spoke,subpulse,b1maps,target_LFA,nspokes,deltat,comp_grad,do_mls,refoc_pulse)\n");                
        mexErrMsgTxt("Wrong number of input parameters.");        
    }

    // STRUCTURES
    
    // b1s
    double *b1s_re=(double *)mxGetPr(prhs[0]);  
    double *b1s_im=(double *)mxGetPi(prhs[0]);  

    // grad
    double *grad=(double *)mxGetPr(prhs[1]);  
    
    // time_to_sinctime
    double *time_to_sinc_time=(double *)mxGetPr(prhs[2]);  

    // time_to_spoke
    double *time_to_spoke=(double *)mxGetPr(prhs[3]);
    
    // subpulse
    double *subpulse=(double *)mxGetPr(prhs[4]);
    
    // b1maps
    double *b1maps_re=(double *)mxGetPr(prhs[5]);  
    double *b1maps_im=(double *)mxGetPi(prhs[5]);  

    // target_LFA
    double *target_LFA_re=(double *)mxGetPr(prhs[6]);
    double *target_LFA_im=(double *)mxGetPi(prhs[6]);    
    
    // SCALARS
    
    // nspokes
    int nspokes=(int)mxGetScalar(prhs[7]);

    // deltat
    double deltat=(double)mxGetScalar(prhs[8]);

    // comp_grad
    int comp_grad=(int)mxGetScalar(prhs[9]);
   
    // do_mls
    int do_mls=(int)mxGetScalar(prhs[10]);

    // refoc_pulse
    int refoc_pulse=(int)mxGetScalar(prhs[11]);
    

    // CONSISTENCY CHECK
    int nnzp=(int)mxGetM(prhs[0]);
    if( (int)mxGetM(prhs[1])!=nnzp || (int)mxGetM(prhs[5])!=nnzp || (int)mxGetM(prhs[6])!=nnzp ){
        mexErrMsgTxt("Inconsistent dimension of input structures (nnzp).");
    }

    int ntimes=(int)mxGetN(prhs[0]);    
    if( (int)mxGetN(prhs[1])!=ntimes || (int)mxGetM(prhs[2])!=ntimes || (int)mxGetM(prhs[3])!=ntimes  ){
        mexErrMsgTxt("Inconsistent dimension of input structures (ntimes).");
    }

    int ncoils=mxGetN(prhs[5]);
    
    
    // ALLOCATE MEMORY
    double *f,*df;
    
    plhs[0]=mxCreateDoubleMatrix( 1,1,mxREAL);
    f=(double *)mxGetPr(plhs[0]);
    
    plhs[1]=mxCreateDoubleMatrix( 2*nspokes*ncoils,1,mxREAL);
    df=(double *)mxGetPr(plhs[1]);        
    
    if(nlhs>2){
        mexErrMsgTxt("Too many output variables.");
    }
    
    // PERFORM COMPUTATION        
    // printf("nnzp=%d  ntimes=%d  ncoils=%d  nspokes=%d  deltat=%e  comp_grad=%d  \n",nnzp,ntimes,ncoils,nspokes,deltat,comp_grad);
    
    f0_LFA_MZ_mex(f,df,b1s_re,b1s_im,grad,time_to_sinc_time,time_to_spoke,subpulse,b1maps_re,b1maps_im,target_LFA_re,target_LFA_im,
        nnzp,ntimes,ncoils,nspokes,deltat,comp_grad,do_mls,refoc_pulse);

}



























