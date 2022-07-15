


#include "mex.h"
#include "omp.h"

#include <complex>
#include <cmath>

using namespace std;


#define TWO_PI_GAMMA (2.675128976384781e+08);
#define NUM_THREADS 8  // number of CPUs to be used in parallel using OMP


void hessian_LFA_noB0_mex(double *btotspokes_re,double *btotspokes_im,double *q_gblips_re,double *q_gblips_im,double *b1maps_re,double *b1maps_im,double *target_LFA_re,double *target_LFA_im,
        int nnzp,int nspokes,int ncoils,double sumsinc,double deltat,int do_mls,int refoc_pulse,
        double *h);

void multiply_2spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> *a3,complex<double> *b3);

void multiply_3spinors(complex<double> a1,complex<double> b1,complex<double> a2,complex<double> b2,complex<double> a3,complex<double> b3,complex<double> *a4,complex<double> *b4);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[]);



void hessian_LFA_noB0_mex(double *btotspokes_re,double *btotspokes_im,double *q_gblips_re,double *q_gblips_im,double *b1maps_re,double *b1maps_im,double *target_LFA_re,double *target_LFA_im,
        int nnzp,int nspokes,int ncoils,double sumsinc,double deltat,int do_mls,int refoc_pulse,
        double *h)
{

    // double T0=omp_get_wtime();
    
    // MEMORY ALLOCATION total ~ 140 MB (nnzp=3702 ncoils=8 nspokes=2)
    complex<double> *as_gblip=new complex<double>[nnzp*nspokes];  // 0.11 MB
    complex<double> *bs_gblip=new complex<double>[nnzp*nspokes];  // 0.11 MB
    
    complex<double> *a_forw=new complex<double>[nnzp*nspokes];  // 0.11 MB
    complex<double> *b_forw=new complex<double>[nnzp*nspokes];  // 0.11 MB
    
    for(int i=0;i<nnzp;i++){  // only need to initialize the first spoke of the one step late forward products
        a_forw[i]=complex<double>(1.0,0.0);
        b_forw[i]=complex<double>(0.0,0.0);
    }
    
    double *dphi_dre=new double[nnzp*ncoils*nspokes];  // 1 MB
    double *dphi_dim=new double[nnzp*ncoils*nspokes];  // 1 MB

    complex<double> *dexppha_dre=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *dexppha_dim=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB

    double *h_phi_1=new double[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    double *h_phi_2=new double[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    double *h_phi_3=new double[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB

    complex<double> *h_exppha_1=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    complex<double> *h_exppha_2=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    complex<double> *h_exppha_3=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB

    complex<double> *h_a_1=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB 
    complex<double> *h_a_2=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    complex<double> *h_a_3=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB

    complex<double> *h_b_1=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    complex<double> *h_b_2=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB
    complex<double> *h_b_3=new complex<double>[nspokes*nnzp*ncoils*ncoils];  // 7.6 MB

    complex<double> *da_dre=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *da_dim=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB

    complex<double> *db_dre=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *db_dim=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
        
    complex<double> *a=new complex<double>[nnzp];  // 0.06 MB
    complex<double> *b=new complex<double>[nnzp];  // 0.06 MB
    for(int i=0;i<nnzp;i++){
        a[i]=complex<double>(1.0,0.0);
        b[i]=complex<double>(0.0,0.0);
    }
     
    int ncunk=nspokes*ncoils;
    complex<double> *ha=new complex<double>[nnzp*(2*ncunk)*(2*ncunk)];  // 15 MB
    complex<double> *hb=new complex<double>[nnzp*(2*ncunk)*(2*ncunk)];  // 15 MB
    
    complex<double> *ja_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *ja_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    
    complex<double> *jb_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *jb_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    
    complex<double> *a_tmp1_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *a_tmp1_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB

    complex<double> *a_tmp2_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *a_tmp2_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB

    complex<double> *b_tmp1_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *b_tmp1_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB

    complex<double> *b_tmp2_re=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    complex<double> *b_tmp2_im=new complex<double>[nnzp*ncoils*nspokes];  // 1 MB
    
    
    // double T1=omp_get_wtime();

    complex<double> pure_imag=complex<double>(0.0,1.0);
    
#pragma omp parallel for num_threads(NUM_THREADS)
    for(int i=0;i<nnzp;i++){  

        // FORWARD BLOCH SIMULATION    
        for(int sp=0;sp<nspokes;sp++){
            
            complex<double> btot=complex<double>( btotspokes_re[i+sp*nnzp],btotspokes_im[i+sp*nnzp] );
            double btot_mag=abs(btot);
            complex<double> exp_btot_pha=exp( pure_imag*arg(btot) );
            
            double phi=-deltat*btot_mag*TWO_PI_GAMMA;  // see Pauly et al.            
            double cosphi=cos(phi/2.0);  // no need to recompute these several times...
            double sinphi=sin(phi/2.0);  
            
            // CK parameters
            complex<double> a2=cosphi;
            complex<double> b2=-pure_imag*exp_btot_pha*sinphi;    
            
            // CK parameters incorporating gradient blips (a gradien blip is a a pure alpha with beta=0)
            complex<double> q_gblip=complex<double>( q_gblips_re[i+sp*nnzp],q_gblips_im[i+sp*nnzp] );
            as_gblip[i+sp*nnzp]=a2 * q_gblip;
            bs_gblip[i+sp*nnzp]=b2 * conj(q_gblip);
            
            // forward Q products (one step late)
            if(sp<nspokes-1)
                multiply_2spinors(a_forw[i+sp*nnzp],b_forw[i+sp*nnzp],as_gblip[i+sp*nnzp],bs_gblip[i+sp*nnzp],a_forw+i+(sp+1)*nnzp,b_forw+i+(sp+1)*nnzp);    
                
			// Jacobian terms (@ each spatial position)
            complex<double> b1map,tmp;
            int ind,ind2;
            for(int j=0;j<ncoils;j++){
                b1map=complex<double>( b1maps_re[i+j*nnzp],b1maps_im[i+j*nnzp] );
                tmp=b1map*conj(exp_btot_pha);

                ind=i + (j+sp*ncoils)*nnzp;
                
                dphi_dre[ind]=-deltat*abs(sumsinc)*real(tmp)*TWO_PI_GAMMA;
                dphi_dim[ind]=deltat*abs(sumsinc)*imag(tmp)*TWO_PI_GAMMA;
                                
                da_dre[ind]=-0.5*sinphi*dphi_dre[ind];
                da_dim[ind]=-0.5*sinphi*dphi_dim[ind];    
                
                if(btot_mag==0){
                    dexppha_dre[ind]=0.0;
                    dexppha_dim[ind]=0.0;                
                }else{                    
                    dexppha_dre[ind]=( sumsinc*b1map - exp_btot_pha*real( sumsinc*b1map*conj(exp_btot_pha) ) ) / btot_mag;  // derivative of the term exp(1j*angle(btot))
                    dexppha_dim[ind]=( pure_imag*sumsinc*b1map - exp_btot_pha*real( pure_imag*sumsinc*b1map*conj(exp_btot_pha) ) ) / btot_mag;
                }    
                
                db_dre[ind]=-pure_imag*sinphi*dexppha_dre[ind] - 0.5*pure_imag*cosphi*exp_btot_pha*dphi_dre[ind];
                db_dim[ind]=-pure_imag*sinphi*dexppha_dim[ind] - 0.5*pure_imag*cosphi*exp_btot_pha*dphi_dim[ind];                
            }            
            
            // Hessian of phi (@ each spatial position)
            double norm=0.0;
            if(btot_mag>0.0)
                norm=-deltat*sumsinc*sumsinc/abs(btot)*TWO_PI_GAMMA;
            
            complex<double> tmp_j,tmp_k;
            for(int j=0;j<ncoils;j++){
                tmp_j=complex<double>( b1maps_re[i+j*nnzp],b1maps_im[i+j*nnzp] )*conj(exp_btot_pha);                
                for(int k=0;k<ncoils;k++){
                    tmp_k=complex<double>( b1maps_re[i+k*nnzp],b1maps_im[i+k*nnzp] )*conj(exp_btot_pha);
                    tmp=complex<double>( b1maps_re[i+j*nnzp],b1maps_im[i+j*nnzp] ) * complex<double>( b1maps_re[i+k*nnzp],-b1maps_im[i+k*nnzp] );  // b1maps(:,j).*conj(b1maps(:,k))                    
                    ind=sp + i*nspokes + k*nspokes*nnzp + j*nspokes*nnzp*ncoils;
                    if(k>=j){
                        h_phi_1[ind]=norm*( real(tmp) - real(tmp_j)*real(tmp_k) );            
                        h_phi_3[ind]=norm*( real(tmp) - imag(tmp_j)*imag(tmp_k) );                                                        
                        if(k!=j){
                            ind2=sp + i*nspokes + j*nspokes*nnzp + k*nspokes*nnzp*ncoils;
                            h_phi_1[ind2]=h_phi_1[ind];  // symmetric term
                            h_phi_3[ind2]=h_phi_3[ind];  // symmetric term
                        }
                    }
                    h_phi_2[ind]=norm*( -imag(tmp) + imag(tmp_j)*real(tmp_k) );  // not symmetric
                }   
            }
            
            // Hessian of a (@ each spatial location)
            int ind_j,ind_k;
            for(int j=0;j<ncoils;j++){
                ind_j=i + (j+sp*ncoils)*nnzp;
                for(int k=0;k<ncoils;k++){
                    ind_k=i + (k+sp*ncoils)*nnzp;
                    ind=sp + i*nspokes + k*nspokes*nnzp + j*nspokes*nnzp*ncoils;
                    if(k>=j){
                        h_a_1[ind]=-0.25*cosphi*dphi_dre[ind_j]*dphi_dre[ind_k] - 0.5*sinphi*h_phi_1[ind];  // hessian term #1 (dre_dre)
                        h_a_3[ind]=-0.25*cosphi*dphi_dim[ind_j]*dphi_dim[ind_k] - 0.5*sinphi*h_phi_3[ind];  // hessian term #3 (dim_dim)
                    }
                    if(k!=j){
                        ind2=sp + i*nspokes + j*nspokes*nnzp + k*nspokes*nnzp*ncoils;
                        h_a_1[ind2]=h_a_1[ind];  // symmetric term
                        h_a_3[ind2]=h_a_3[ind];  // symmetric term
                    }                    
                    h_a_2[ind]=-0.25*cosphi*dphi_dim[ind_j]*dphi_dre[ind_k] - 0.5*sinphi*h_phi_2[ind];  // hessian term #2 (dre_dim)
                }
            }
            
            // Hessian of the term exp(1j*angle(btot)) (@ each spatial position)
            double V=btot_mag;
            if(V==0){
                for(int j=0;j<ncoils*ncoils;j++){
                    h_exppha_1[sp+i*nspokes+j]=h_exppha_2[sp+i*nspokes+j]=h_exppha_3[sp+i*nspokes+j]=0.0;
                }
            }else{            
                complex<double> U,dU,dV,b1map_j,b1map_k;
                for(int j=0;j<ncoils;j++){
                    ind_j=i + (j+sp*ncoils)*nnzp;
                    b1map_j=complex<double>( b1maps_re[i+j*nnzp],b1maps_im[i+j*nnzp] );
                    for(int k=0;k<ncoils;k++){
                        ind_k=i + (k+sp*ncoils)*nnzp;                        
                        b1map_k=complex<double>( b1maps_re[i+k*nnzp],b1maps_im[i+k*nnzp] );
                        
                        ind=sp + i*nspokes + k*nspokes*nnzp + j*nspokes*nnzp*ncoils;                        
                        if(k>=j){                            
                            // term #1
                            U=sumsinc*( b1map_j - exp_btot_pha*real( b1map_j*conj(exp_btot_pha) ) );
                            dU=-dexppha_dre[ind_k]*real( sumsinc*b1map_j*conj(exp_btot_pha) ) - exp_btot_pha*real( sumsinc*b1map_j*conj(dexppha_dre[ind_k] ) );
                            dV=real( sumsinc*b1map_k*conj(exp_btot_pha) );
                            h_exppha_1[ind]=( dU*V-dV*U )/(V*V);
                            
                            // term #3
                            U=sumsinc*( pure_imag*b1map_j - exp_btot_pha*real( pure_imag*b1map_j*conj(exp_btot_pha) ) );
                            dU=dexppha_dim[ind_k]*imag( sumsinc*b1map_j*conj(exp_btot_pha) ) + exp_btot_pha*imag( sumsinc*b1map_j*conj(dexppha_dim[ind_k] ) );
                            dV=-imag( sumsinc*b1map_k*conj(exp_btot_pha) );
                            h_exppha_3[ind]=( dU*V-dV*U )/(V*V);                    
                            
                            if(k!=j){
                                ind2=sp + i*nspokes + j*nspokes*nnzp + k*nspokes*nnzp*ncoils;
                                h_exppha_1[ind2]=h_exppha_1[ind];  // symmetric term
                                h_exppha_3[ind2]=h_exppha_3[ind];  // symmetric term
                            }
                        }            
                        
                        // term #2
                        U=sumsinc*( pure_imag*b1map_j - exp_btot_pha*real( pure_imag*b1map_j*conj(exp_btot_pha) ) );
                        dU=dexppha_dre[ind_k]*imag( sumsinc*b1map_j*conj(exp_btot_pha) ) + exp_btot_pha*imag( sumsinc*b1map_j*conj(dexppha_dre[ind_k] ) );
                        dV=real( sumsinc*b1map_k*conj(exp_btot_pha) );
                        h_exppha_2[ind]=( dU*V-dV*U )/(V*V);
                        
                    }
                }
            }
                        
            // Hessian of b (@ each spatial position)
            complex<double> tmp1,tmp2,tmp3,tmp4;
            for(int j=0;j<ncoils;j++){
                ind_j=i + (j+sp*ncoils)*nnzp;
                for(int k=0;k<ncoils;k++){
                    ind_k=i + (k+sp*ncoils)*nnzp;
                    
                    ind=sp + i*nspokes + k*nnzp*nspokes + j*nnzp*ncoils*nspokes;                    
                    if(k>=j){                        
                        // term #1
                        tmp1=-0.5*pure_imag*cosphi*( dphi_dre[ind_k]*dexppha_dre[ind_j] + dphi_dre[ind_j]*dexppha_dre[ind_k] );
                        tmp2=-pure_imag*sinphi*h_exppha_1[ind];
                        tmp3=0.25*pure_imag*sinphi*exp_btot_pha*dphi_dre[ind_k]*dphi_dre[ind_j];
                        tmp4=-0.5*pure_imag*cosphi*exp_btot_pha*h_phi_1[ind];
                        h_b_1[ind]=tmp1 + tmp2 + tmp3 + tmp4;
                        
                        // term #3
                        tmp1=-0.5*pure_imag*cosphi*( dphi_dim[ind_k]*dexppha_dim[ind_j] + dphi_dim[ind_j]*dexppha_dim[ind_k] );
                        tmp2=-pure_imag*sinphi*h_exppha_3[ind];
                        tmp3=0.25*pure_imag*sinphi*exp_btot_pha*dphi_dim[ind_k]*dphi_dim[ind_j];
                        tmp4=-0.5*pure_imag*cosphi*exp_btot_pha*h_phi_3[ind];
                        h_b_3[ind]=tmp1 + tmp2 + tmp3 + tmp4;
                        
                        if(k!=j){
                            ind2=sp + i*nspokes + j*nnzp*nspokes + k*nnzp*ncoils*nspokes;
                            h_b_1[ind2]=h_b_1[ind];  // symmetric term
                            h_b_3[ind2]=h_b_3[ind];  // symmetric term
                        }
                    }
                    
                    // term #2
                    tmp1=-0.5*pure_imag*cosphi*( dphi_dre[ind_k]*dexppha_dim[ind_j] + dphi_dim[ind_j]*dexppha_dre[ind_k] );
                    tmp2=-pure_imag*sinphi*h_exppha_2[ind];
                    tmp3=0.25*pure_imag*sinphi*exp_btot_pha*dphi_dre[ind_k]*dphi_dim[ind_j];
                    tmp4=-0.5*pure_imag*cosphi*exp_btot_pha*h_phi_2[ind];
                    h_b_2[ind]=tmp1 + tmp2 + tmp3 + tmp4;
                }
            }

            // incorporate gradient blips in Hessian and Jacobian terms
            for(int j=0;j<ncoils;j++){
                
                for(int k=0;k<ncoils;k++){
                    ind=sp + i*nspokes + k*nnzp*nspokes + j*nspokes*nnzp*ncoils;
                    h_a_1[ind]=h_a_1[ind] * q_gblip;
                    h_a_2[ind]=h_a_2[ind] * q_gblip;
                    h_a_3[ind]=h_a_3[ind] * q_gblip;
                    
                    h_b_1[ind]=h_b_1[ind] * conj(q_gblip);
                    h_b_2[ind]=h_b_2[ind] * conj(q_gblip);
                    h_b_3[ind]=h_b_3[ind] * conj(q_gblip);                    
                }
                
                ind=i + (j+sp*ncoils)*nnzp;
                da_dre[ind]=da_dre[ind] * q_gblip;
                da_dim[ind]=da_dim[ind] * q_gblip;
                
                db_dre[ind]=db_dre[ind] * conj(q_gblip);
                db_dim[ind]=db_dim[ind] * conj(q_gblip);
            }
            
        }  // for(int sp=0;sp<nspokes;sp++)

        
        
        // BACKWARD BLOCH SIMULATION         
        for(int sp=nspokes-1;sp>=0;sp--){
                 
            int ind=i + sp*nnzp;
            
            for(int j=0;j<ncoils;j++){
                int ind_j=j + sp*ncoils;
                int ind2=i + (j+sp*ncoils)*nnzp;
                        
                // Jacobian (sandwitch product of backward Q product, Jacobian of current Q matrix and forward Q products) 
                multiply_3spinors(a_forw[ind],b_forw[ind],da_dre[ind2],db_dre[ind2],a[i],b[i],ja_re+ind2,jb_re+ind2);
                multiply_3spinors(a_forw[ind],b_forw[ind],da_dim[ind2],db_dim[ind2],a[i],b[i],ja_im+ind2,jb_im+ind2);
                
                // Hessian term (derivatives taken within a single spoke)
                for(int k=0;k<ncoils;k++){
                    int ind_k=k + sp*ncoils;
                    int ind3=sp + i*nspokes + k*nspokes*nnzp + j*nspokes*nnzp*ncoils;
                    
                    complex<double> tmp_a_1,tmp_a_2,tmp_a_3,tmp_b_1,tmp_b_2,tmp_b_3;
                    multiply_3spinors(a_forw[ind],b_forw[ind],h_a_1[ind3],h_b_1[ind3],a[i],b[i],&tmp_a_1,&tmp_b_1);
                    multiply_3spinors(a_forw[ind],b_forw[ind],h_a_2[ind3],h_b_2[ind3],a[i],b[i],&tmp_a_2,&tmp_b_2);
                    multiply_3spinors(a_forw[ind],b_forw[ind],h_a_3[ind3],h_b_3[ind3],a[i],b[i],&tmp_a_3,&tmp_b_3);
                    
                    ha[i + (ind_k       + ind_j        *2*ncunk)*nnzp]=tmp_a_1;
                    ha[i + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp]=tmp_a_2;
                    ha[i + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp]=tmp_a_2;
                    ha[i + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp]=tmp_a_3;                    

                    hb[i + (ind_k       + ind_j        *2*ncunk)*nnzp]=tmp_b_1;
                    hb[i + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp]=tmp_b_2;
                    hb[i + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp]=tmp_b_2;
                    hb[i + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp]=tmp_b_3;
                }
                
            }                
            
            // Jacobian terms (derivatives taken for different spokes)
            if(sp>0){
                
                // product between left-hand jacobian term and backward Q products
                for(int j=0;j<ncoils;j++){
                    ind=i + (j+sp*ncoils)*nnzp;
                    multiply_2spinors(da_dre[ind],db_dre[ind],a[i],b[i],a_tmp1_re+ind,b_tmp1_re+ind);
                    multiply_2spinors(da_dim[ind],db_dim[ind],a[i],b[i],a_tmp1_im+ind,b_tmp1_im+ind);  
                } 
                
                for(int sp2=sp-1;sp2>=0;sp2--){
                    
                    // compute the product of Q matrices sandwitched between indices i and j
                    complex<double> a_sandw=complex<double>(1.0,0.0);
                    complex<double> b_sandw=complex<double>(0.0,0.0);                        
                    for(int sp3=sp-1;sp3>=sp2+1;sp3--) 
                        multiply_2spinors(as_gblip[i+sp3*nnzp],bs_gblip[i+sp3*nnzp],a_sandw,b_sandw,&a_sandw,&b_sandw);
                                        
                    // product between sandwictch term, right-hand jacobian term and forward Q products 
                    for(int j=0;j<ncoils;j++){
                        ind=i + (j+sp2*ncoils)*nnzp;
                        multiply_3spinors(a_forw[i+sp2*nnzp],b_forw[i+sp2*nnzp],da_dre[ind],db_dre[ind],a_sandw,b_sandw,a_tmp2_re+ind,b_tmp2_re+ind);
                        multiply_3spinors(a_forw[i+sp2*nnzp],b_forw[i+sp2*nnzp],da_dim[ind],db_dim[ind],a_sandw,b_sandw,a_tmp2_im+ind,b_tmp2_im+ind);
                    }

                    // outter products between the columns of the TMP1 and TMP2 Q terms
                    for(int j=0;j<ncoils;j++){
                        int ind_j=i + (j+sp*ncoils)*nnzp;
                        for(int k=0;k<ncoils;k++){
                            int ind_k=i + (k+sp2*ncoils)*nnzp;
                            
                            complex<double> tmp_a_1,tmp_a_2,tmp_a_3,tmp_a_4,tmp_b_1,tmp_b_2,tmp_b_3,tmp_b_4;
                            multiply_2spinors(a_tmp2_re[ind_k],b_tmp2_re[ind_k] , a_tmp1_re[ind_j],b_tmp1_re[ind_j] , &tmp_a_1,&tmp_b_1 );
                            multiply_2spinors(a_tmp2_re[ind_k],b_tmp2_re[ind_k] , a_tmp1_im[ind_j],b_tmp1_im[ind_j] , &tmp_a_2,&tmp_b_2 );
                            multiply_2spinors(a_tmp2_im[ind_k],b_tmp2_im[ind_k] , a_tmp1_re[ind_j],b_tmp1_re[ind_j] , &tmp_a_3,&tmp_b_3 );
                            multiply_2spinors(a_tmp2_im[ind_k],b_tmp2_im[ind_k] , a_tmp1_im[ind_j],b_tmp1_im[ind_j] , &tmp_a_4,&tmp_b_4 );

                            // half of the Hessian terms are filled by symmetry
                            int ind_j2= j + sp*ncoils;
                            int ind_k2= k + sp*ncoils;
                            int ind_j3=j + sp2*ncoils;
                            int ind_k3=k + sp2*ncoils;

                            ha[i + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = ha[i + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_a_1;
                            ha[i + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha[i + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_a_2;
                            ha[i + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = ha[i + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_3;
                            ha[i + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha[i + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_4;

							hb[i + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = hb[i + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_b_1;
                            hb[i + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb[i + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_b_2;
                            hb[i + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = hb[i + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_3;
                            hb[i + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb[i + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_4;
                        }                        
                    }
                    
                }  // for(int sp2=sp-1;sp>=0;sp2--)

			}  // if(sp>0)
        
            // keep track of the backward Q products
            multiply_2spinors(as_gblip[i+sp*nnzp],bs_gblip[i+sp*nnzp],a[i],b[i],a+i,b+i);

        }  // for(int sp=nspokes-1;sp>=0;sp++)
        
    }  // for(int i=0;i<nnzp;i++)


    // double T2=omp_get_wtime();
    
    // COMBINE ALL TERMS TO FORM THE TOTAL HESSIAN
    complex<double> *doptmetric_dre=new complex<double>[nnzp*ncoils*nspokes];
    complex<double> *doptmetric_dim=new complex<double>[nnzp*ncoils*nspokes];
    
    // initialize total Hessian
    for(int i=0;i<(2*ncunk)*(2*ncunk);i++)
        h[i]=0.0;

#pragma omp parallel num_threads(NUM_THREADS)
    {            
            double *h_omp=new double[(2*ncunk)*(2*ncunk)];
            for(int i=0;i<(2*ncunk)*(2*ncunk);i++)
                h_omp[i]=0.0;
             
            #pragma omp for
            for(int i=0;i<nnzp;i++){
                
                // y
                complex<double> y;
                complex<double> target_LFA=complex<double>(target_LFA_re[i],target_LFA_im[i]);
                if(refoc_pulse==0){
                    if(do_mls==1)
                        y=pow(abs(a[i]),2) - pow(abs(b[i]),2) - target_LFA;
                    else
                        y=2.0*conj(a[i])*b[i] - target_LFA;
                }else{
                    if(do_mls==1)
                        y=pow(abs(b[i]),2) - target_LFA;
                    else
                        y=pow(b[i],2) - target_LFA;
                }
                
                // Jacobian terms of optimization metric
                for(int j=0;j<ncunk;j++){
                    
                    if(refoc_pulse==0){                        
                        if(do_mls==1){
                            doptmetric_dre[i+j*nnzp]=2.0*real( ja_re[i+j*nnzp]*conj(a[i]) - jb_re[i+j*nnzp]*conj(b[i]) );                            
                            doptmetric_dim[i+j*nnzp]=2.0*real( ja_im[i+j*nnzp]*conj(a[i]) - jb_im[i+j*nnzp]*conj(b[i]) );                    
                        }else{
                            doptmetric_dre[i+j*nnzp]=2.0*( conj(ja_re[i+j*nnzp])*b[i] + jb_re[i+j*nnzp]*conj(a[i]) );    
                            doptmetric_dim[i+j*nnzp]=2.0*( conj(ja_im[i+j*nnzp])*b[i] + jb_im[i+j*nnzp]*conj(a[i]) );    
                        }                            
                    }else{                        
                        if(do_mls==1){
                            doptmetric_dre[i+j*nnzp]=2.0*real( jb_re[i+j*nnzp]*conj(b[i]) );                            
                            doptmetric_dim[i+j*nnzp]=2.0*real( jb_im[i+j*nnzp]*conj(b[i]) );                    
                        }else{
                            doptmetric_dre[i+j*nnzp]=2.0*jb_re[i+j*nnzp]*b[i];                            
                            doptmetric_dim[i+j*nnzp]=2.0*jb_im[i+j*nnzp]*b[i];
                        }
                    }
                }
                
                for(int j=0;j<ncunk;j++){
                    for(int k=0;k<ncunk;k++){
                        h_omp[ k       + j         *2*ncunk ]+=2.0*real( conj(doptmetric_dre[i+k*nnzp])*doptmetric_dre[i+j*nnzp] );
                        double tmp=2.0*real( conj(doptmetric_dre[i+k*nnzp])*doptmetric_dim[i+j*nnzp] );
                        h_omp[ k       + (j+ncunk)*2*ncunk ]+=tmp;
                        h_omp[ j+ncunk + k        *2*ncunk ]+=tmp;                
                        h_omp[ k+ncunk + (j+ncunk)*2*ncunk ]+=2.0*real( conj(doptmetric_dim[i+k*nnzp])*doptmetric_dim[i+j*nnzp] );
                    }        
                }
                
                // Hessian terms of optimization metric
                for(int j=0;j<ncunk;j++){
                    
                    for(int k=0;k<ncunk;k++){        
                        
                        if(refoc_pulse==0){                    
                            if(do_mls==1){
                                h_omp[k       + j        *2*ncunk] += 4.0*real( y*conj(a[i])*ha[i + (k       + j        *2*ncunk)*nnzp] - y*conj(b[i])*hb[i + (k       + j        *2*ncunk)*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( y*conj(a[i])*ha[i + (k+ncunk + j        *2*ncunk)*nnzp] - y*conj(b[i])*hb[i + (k+ncunk + j        *2*ncunk)*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(a[i])*ha[i + (k       + (j+ncunk)*2*ncunk)*nnzp] - y*conj(b[i])*hb[i + (k       + (j+ncunk)*2*ncunk)*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(a[i])*ha[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp] - y*conj(b[i])*hb[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp] );

                                h_omp[k       + j        *2*ncunk] += 4.0*real( y*conj(ja_re[i+k*nnzp])*ja_re[i+j*nnzp] - y*conj(jb_re[i+k*nnzp])*jb_re[i+j*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( y*conj(ja_im[i+k*nnzp])*ja_re[i+j*nnzp] - y*conj(jb_im[i+k*nnzp])*jb_re[i+j*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(ja_re[i+k*nnzp])*ja_im[i+j*nnzp] - y*conj(jb_re[i+k*nnzp])*jb_im[i+j*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(ja_im[i+k*nnzp])*ja_im[i+j*nnzp] - y*conj(jb_im[i+k*nnzp])*jb_im[i+j*nnzp] );                                    
                            }else{
                                h_omp[k       + j        *2*ncunk] += 4.0*real( conj(y)*( b[i]*conj(ha[i + (k       + j        *2*ncunk)*nnzp]) + conj(a[i])*hb[i + (k       + j        *2*ncunk)*nnzp] ) );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( conj(y)*( b[i]*conj(ha[i + (k+ncunk + j        *2*ncunk)*nnzp]) + conj(a[i])*hb[i + (k+ncunk + j        *2*ncunk)*nnzp] ) );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*( b[i]*conj(ha[i + (k       + (j+ncunk)*2*ncunk)*nnzp]) + conj(a[i])*hb[i + (k       + (j+ncunk)*2*ncunk)*nnzp] ) );                                
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*( b[i]*conj(ha[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp]) + conj(a[i])*hb[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp] ) );
                                
                                h_omp[k       + j        *2*ncunk] += 4.0*real( conj(y)*conj(ja_re[i+k*nnzp])*jb_re[i+j*nnzp] ) + 4.0*real( conj(y)*conj(ja_re[i+j*nnzp])*jb_re[i+k*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( conj(y)*conj(ja_im[i+k*nnzp])*jb_re[i+j*nnzp] ) + 4.0*real( conj(y)*conj(ja_re[i+j*nnzp])*jb_im[i+k*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*conj(ja_re[i+k*nnzp])*jb_im[i+j*nnzp] ) + 4.0*real( conj(y)*conj(ja_im[i+j*nnzp])*jb_re[i+k*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*conj(ja_im[i+k*nnzp])*jb_im[i+j*nnzp] ) + 4.0*real( conj(y)*conj(ja_im[i+j*nnzp])*jb_im[i+k*nnzp] );                    
                            }                            
                        }else{
                            if(do_mls==1){
                                h_omp[k       + j        *2*ncunk] += 4.0*real( y*conj(b[i])*hb[i + (k       + j        *2*ncunk)*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( y*conj(b[i])*hb[i + (k+ncunk + j        *2*ncunk)*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(b[i])*hb[i + (k       + (j+ncunk)*2*ncunk)*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(b[i])*hb[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp] );

                                h_omp[k       + j        *2*ncunk] += 4.0*real( y*conj(jb_re[i+k*nnzp])*jb_re[i+j*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( y*conj(jb_im[i+k*nnzp])*jb_re[i+j*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(jb_re[i+k*nnzp])*jb_im[i+j*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( y*conj(jb_im[i+k*nnzp])*jb_im[i+j*nnzp] );                    
                            }else{
                                h_omp[k       + j        *2*ncunk] += 4.0*real( conj(y)*b[i]*hb[i + (k       + j        *2*ncunk)*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( conj(y)*b[i]*hb[i + (k+ncunk + j        *2*ncunk)*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*b[i]*hb[i + (k       + (j+ncunk)*2*ncunk)*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*b[i]*hb[i + (k+ncunk + (j+ncunk)*2*ncunk)*nnzp] );

                                h_omp[k       + j        *2*ncunk] += 4.0*real( conj(y)*jb_re[i+k*nnzp]*jb_re[i+j*nnzp] );
                                h_omp[k+ncunk + j        *2*ncunk] += 4.0*real( conj(y)*jb_im[i+k*nnzp]*jb_re[i+j*nnzp] );
                                h_omp[k       + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*jb_re[i+k*nnzp]*jb_im[i+j*nnzp] );
                                h_omp[k+ncunk + (j+ncunk)*2*ncunk] += 4.0*real( conj(y)*jb_im[i+k*nnzp]*jb_im[i+j*nnzp] );                    
                            }
                        }        
                        
                    }
                }
                
            }  // for(int i=0;i<nnzp;i++) 
                        
            
            // fill the hessian using the partial hessian computed by each threads
            #pragma omp critical
            {
                for(int i=0;i<(2*ncunk)*(2*ncunk);i++)
                    h[i]+=h_omp[i];
            }
            delete [] h_omp;
    }        
        
        
    // double T3=omp_get_wtime();
    
            
    // clean up memory
    delete [] as_gblip;
    delete [] bs_gblip;
       
    delete [] a_forw;
    delete [] b_forw; 
        
    delete [] dphi_dre;
    delete [] dphi_dim;

    delete [] dexppha_dre;
    delete [] dexppha_dim;

    delete [] h_phi_1;
    delete [] h_phi_2;
    delete [] h_phi_3;

    delete [] h_exppha_1;
    delete [] h_exppha_2;
    delete [] h_exppha_3;

    delete [] h_a_1;
    delete [] h_a_2;
    delete [] h_a_3;

    delete [] h_b_1;
    delete [] h_b_2;
    delete [] h_b_3;

    delete [] da_dre;
    delete [] da_dim;

    delete [] db_dre;
    delete [] db_dim;
    
    delete [] a;
    delete [] b;
    
    delete [] ha;
    delete [] hb;
    
    delete [] ja_re;
    delete [] ja_im;
    
    delete [] jb_re;
    delete [] jb_im;
    
    delete [] a_tmp1_re;
    delete [] a_tmp1_im;

    delete [] a_tmp2_re;
    delete [] a_tmp2_im;

    delete [] b_tmp1_re;
    delete [] b_tmp1_im;

    delete [] b_tmp2_re;
    delete [] b_tmp2_im;
    
    delete [] doptmetric_dre; 
    delete [] doptmetric_dim;
    
    // double T4=omp_get_wtime();    
    // double T_tot=T4-T0;
    
    /*
    printf("T_mem\t=\t%f\t%3.3f%%\nT_bloch\t=\t%f\t%3.3f%%\nT_comb\t=\t%f\t%3.3f%%\nT_clean\t=\t%f\t%3.3f%%\n",
            T1-T0,(T1-T0)/T_tot*100.0,
            T2-T1,(T2-T1)/T_tot*100.0,
            T3-T2,(T3-T2)/T_tot*100.0,
            T4-T3,(T4-T3)/T_tot*100.0);
    printf("TOTAL\t=\t%f\t100%%\n",T_tot);
    */
    
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
    
    if(nrhs!=8){
        printf("Usage:\n\th=hessian_LFA_noB0_mex(btotspokes,q_gblips,b1maps,target_LFA,sumsinc,deltat,do_mls,refoc_pulse)\n");                
        mexErrMsgTxt("Wrong number of input parameters.");        
    }

    // STRUCTURES
    
    // btotspokes
    double *btotspokes_re=(double *)mxGetPr(prhs[0]);  
    double *btotspokes_im=(double *)mxGetPi(prhs[0]);  

    // q_gblips
    double *q_gblips_re=(double *)mxGetPr(prhs[1]);  
    double *q_gblips_im=(double *)mxGetPi(prhs[1]);  

    // b1maps
    double *b1maps_re=(double *)mxGetPr(prhs[2]);  
    double *b1maps_im=(double *)mxGetPi(prhs[2]);  

    // target_LFA
    double *target_LFA_re=(double *)mxGetPr(prhs[3]);  
    double *target_LFA_im=(double *)mxGetPi(prhs[3]); 
    
    
    // SCALARS
    
    // sumsinc
    double sumsinc=(double)mxGetScalar(prhs[4]);

    // deltat
    double deltat=(double)mxGetScalar(prhs[5]);

    // do_mls
    int do_mls=(int)mxGetScalar(prhs[6]);

    // refoc_pulse
    int refoc_pulse=(int)mxGetScalar(prhs[7]);

    // CONSISTENCY CHECK    
    int nnzp=(int)mxGetM(prhs[0]);
    if( (int)mxGetM(prhs[1])!=nnzp || (int)mxGetM(prhs[2])!=nnzp || (int)mxGetM(prhs[3])!=nnzp ){
        mexErrMsgTxt("Inconsistent dimension of input structures (nnzp).");
    }

    int nspokes=(int)mxGetN(prhs[0]);
    if( (int)mxGetN(prhs[1])!=nspokes ){
        mexErrMsgTxt("Inconsistent dimension of input structures (ntimes).");
    }

    int ncoils=mxGetN(prhs[2]);
    
    if( (int)mxGetN(prhs[3])!=1 ){
        mexErrMsgTxt("target_FA needs to be a complex column vector.");
    }
    

    // ALLOCATE MEMORY
    int ncunk=nspokes*ncoils;
    
    plhs[0]=mxCreateDoubleMatrix( 2*ncunk,2*ncunk,mxREAL);
    double *h=(double *)mxGetPr(plhs[0]);
    
    if(nlhs>1){
        mexErrMsgTxt("Too many output variables.");
    }
    
    
    // PERFORM COMPUTATION
    // printf("nnzp=%d  ncoils=%d  nspokes=%d  sumsinc=%e\n",nnzp,ncoils,nspokes,sumsinc);    
    hessian_LFA_noB0_mex(btotspokes_re,btotspokes_im,q_gblips_re,q_gblips_im,b1maps_re,b1maps_im,target_LFA_re,target_LFA_im,
            nnzp,nspokes,ncoils,sumsinc,deltat,do_mls,refoc_pulse,h);
    
}










