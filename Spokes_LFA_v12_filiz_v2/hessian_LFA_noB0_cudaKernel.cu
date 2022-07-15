



#include <cmath>

#include <stdio.h>


#define TWO_PI_GAMMA (2.675128976384781e+08)



__device__ void multiply_2spinors(double a1_re,double a1_im,double b1_re,double b1_im,
                                  double a2_re,double a2_im,double b2_re,double b2_im,
                                  double *a3_re,double *a3_im,double *b3_re,double *b3_im)
{
  *a3_re = a2_re*a1_re - a2_im*a1_im - b2_re*b1_re - b2_im*b1_im;
  *a3_im = a2_im*a1_re + a2_re*a1_im + b2_im*b1_re - b2_re*b1_im;
  
  *b3_re = b2_re*a1_re - b2_im*a1_im + a2_re*b1_re + a2_im*b1_im;
  *b3_im = b2_im*a1_re + b2_re*a1_im - a2_im*b1_re + a2_re*b1_im;
}




__device__ void multiply_3spinors(double a1_re,double a1_im,double b1_re,double b1_im,
                                  double a2_re,double a2_im,double b2_re,double b2_im,
                                  double a3_re,double a3_im,double b3_re,double b3_im,
                                  double *a4_re,double *a4_im,double *b4_re,double *b4_im)
{
  
  double tmp1_re = a2_re*a1_re - a2_im*a1_im - b2_re*b1_re - b2_im*b1_im;
  double tmp1_im = a2_im*a1_re + a2_re*a1_im + b2_im*b1_re - b2_re*b1_im;
  
  double tmp2_re = b2_re*a1_re - b2_im*a1_im + a2_re*b1_re + a2_im*b1_im;
  double tmp2_im = b2_im*a1_re + b2_re*a1_im - a2_im*b1_re + a2_re*b1_im;
  
  *a4_re = a3_re*tmp1_re - a3_im*tmp1_im - b3_re*tmp2_re - b3_im*tmp2_im;
  *a4_im = a3_im*tmp1_re + a3_re*tmp1_im + b3_im*tmp2_re - b3_re*tmp2_im;
  
  *b4_re = b3_re*tmp1_re - b3_im*tmp1_im + a3_re*tmp2_re + a3_im*tmp2_im;
  *b4_im = b3_im*tmp1_re + b3_re*tmp1_im - a3_im*tmp2_re + a3_re*tmp2_im;    
}



// set inc_gblip==1 for incorporation of the gradient blip in the CK parameters
__device__ void compute_caley_klein_params(double b1_re,double b1_im,double deltat,
					   double gblip_re,double gblip_im,int inc_gblip,
					   double *b1_mag,double *b1_pha,double *cos_b1_pha,double *sin_b1_pha,
					   double *phi,double *cosphi,double *sinphi,
					   double *a_re,double *a_im,double *b_re,double *b_im)
{

  *b1_mag=sqrt( b1_re*b1_re + b1_im*b1_im );
    
  *phi=-deltat*(*b1_mag)*TWO_PI_GAMMA;  // see Pauly et al.
  sincos( (*phi)/2.0,sinphi,cosphi );
  
  if(inc_gblip==1){
    *a_re = (*cosphi) * gblip_re;  // a2 = a2 * Qblip
    *a_im = (*cosphi) * gblip_im;
  }else{
    *a_re = (*cosphi);
    *a_im = 0.0;
  }


  *b1_pha = atan2( b1_im,b1_re );
  sincos( *b1_pha,sin_b1_pha,cos_b1_pha );
    
  double tmp_re= (*sin_b1_pha)*(*sinphi);
  double tmp_im=-(*cos_b1_pha)*(*sinphi);
  
  if(inc_gblip==1){
    *b_re = tmp_re*gblip_re + tmp_im*gblip_im;  // b2 = b2 * conj(Qblip)
    *b_im = tmp_im*gblip_re - tmp_re*gblip_im;
  }else{
    *b_re = tmp_re;
    *b_im = tmp_im;
  }

}




__device__  void compute_jacobian_Q_matrix( double b1maps_re,double b1maps_im,
					    double cos_btot_pha,double sin_btot_pha,
					    double cosphi,double sinphi,
					    double btot_mag,double deltat,double sumsinc,
					    double gblip_re,double gblip_im,
					    double *da_dre_re,double *da_dre_im,double *da_dim_re,double *da_dim_im,
					    double *db_dre_re,double *db_dre_im,double *db_dim_re,double *db_dim_im )
{

  double tmp_re = b1maps_re*cos_btot_pha + b1maps_im*sin_btot_pha;
  double tmp_im = b1maps_im*cos_btot_pha - b1maps_re*sin_btot_pha;
  
  double dphi_dre = -deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_re;
  double dphi_dim =  deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_im;

  // derivative of alpha
  double tmp_re_re = -0.5*sinphi*dphi_dre;
  double tmp_re_im = 0.0;

  double tmp_im_re = -0.5*sinphi*dphi_dim;
  double tmp_im_im = 0.0;
  
  // gradient blip: alpha = alpha * gblip
  *da_dre_re = tmp_re_re * gblip_re;
  *da_dre_im = tmp_re_re * gblip_im;

  *da_dim_re = tmp_im_re * gblip_re;
  *da_dim_im = tmp_im_re * gblip_im;
    
  // derivative of the term exp(1j*angle(btot))
  double dexppha_dre_re=0.0,dexppha_dre_im=0.0,dexppha_dim_re=0.0,dexppha_dim_im=0.0;
  if(btot_mag>0){
    dexppha_dre_re = sumsinc/btot_mag*( b1maps_re - cos_btot_pha*tmp_re );
    dexppha_dre_im = sumsinc/btot_mag*( b1maps_im - sin_btot_pha*tmp_re );

    dexppha_dim_re = sumsinc/btot_mag*( -b1maps_im + cos_btot_pha*tmp_im );
    dexppha_dim_im = sumsinc/btot_mag*(  b1maps_re + sin_btot_pha*tmp_im );
  }

  // derivative of beta
  tmp_re_re = sinphi*dexppha_dre_im + 0.5*cosphi*sin_btot_pha*dphi_dre;
  tmp_re_im = -sinphi*dexppha_dre_re - 0.5*cosphi*cos_btot_pha*dphi_dre;

  tmp_im_re = sinphi*dexppha_dim_im + 0.5*cosphi*sin_btot_pha*dphi_dim;
  tmp_im_im = -sinphi*dexppha_dim_re - 0.5*cosphi*cos_btot_pha*dphi_dim;

  // gradient blip: beta = beta * conj(gblip)
  *db_dre_re = tmp_re_re*gblip_re + tmp_re_im*gblip_im;
  *db_dre_im = tmp_re_im*gblip_re - tmp_re_re*gblip_im;

  *db_dim_re = tmp_im_re*gblip_re + tmp_im_im*gblip_im;
  *db_dim_im = tmp_im_im*gblip_re - tmp_im_re*gblip_im;
  
}





__device__ void compute_hessian_Q_matrix( double cosphi,double sinphi,
					  double b1maps_j_re,double b1maps_j_im,
					  double b1maps_k_re,double b1maps_k_im,
					  double cos_btot_pha,double sin_btot_pha, double btot_mag,
					  double deltat,double sumsinc,				
					  double gblip_re,double gblip_im,
					  double *h_a_1_re,double *h_a_1_im,double *h_b_1_re,double *h_b_1_im,
					  double *h_a_2_re,double *h_a_2_im,double *h_b_2_re,double *h_b_2_im,
					  double *h_a_3_re,double *h_a_3_im,double *h_b_3_re,double *h_b_3_im )
{

  double tmp_j_re = b1maps_j_re*cos_btot_pha + b1maps_j_im*sin_btot_pha;
  double tmp_j_im = b1maps_j_im*cos_btot_pha - b1maps_j_re*sin_btot_pha;

  double tmp_k_re = b1maps_k_re*cos_btot_pha + b1maps_k_im*sin_btot_pha;
  double tmp_k_im = b1maps_k_im*cos_btot_pha - b1maps_k_re*sin_btot_pha;

  double dphi_dre_j = -deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_j_re;
  double dphi_dim_j =  deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_j_im;

  double dphi_dre_k = -deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_k_re;
  double dphi_dim_k =  deltat*abs(sumsinc)*TWO_PI_GAMMA*tmp_k_im;


  // HESSIAN OF PHI
  double norm=0.0;
  if(btot_mag>0.0)
    norm=-deltat*sumsinc*sumsinc/btot_mag*TWO_PI_GAMMA;

  // b1maps(:,j).*conj(b1maps(:,k))
  double tmp_re = b1maps_j_re*b1maps_k_re + b1maps_j_im*b1maps_k_im;
  double tmp_im = b1maps_j_im*b1maps_k_re - b1maps_j_re*b1maps_k_im;

  double h_phi_1 = norm*( tmp_re - tmp_j_re*tmp_k_re );  
  double h_phi_3 = norm*( tmp_re - tmp_j_im*tmp_k_im );
  double h_phi_2 = norm*( -tmp_im + tmp_j_im*tmp_k_re );


  // HESSIAN OF ALPHA
  double tmp_h_1_re = -0.25*cosphi*dphi_dre_j*dphi_dre_k - 0.5*sinphi*h_phi_1;  // hessian term #1 (dre_dre)
  double tmp_h_3_re = -0.25*cosphi*dphi_dim_j*dphi_dim_k - 0.5*sinphi*h_phi_3;  // hessian term #3 (dim_dim) 
  double tmp_h_2_re = -0.25*cosphi*dphi_dim_j*dphi_dre_k - 0.5*sinphi*h_phi_2;  // hessian term #2 (dre_dim)

  double tmp_h_1_im = 0.0;  // alpha is real when B0=0
  double tmp_h_2_im = 0.0;
  double tmp_h_3_im = 0.0;

  *h_a_1_re = tmp_h_1_re * gblip_re;
  *h_a_1_im = tmp_h_1_re * gblip_im;

  *h_a_2_re = tmp_h_2_re * gblip_re;
  *h_a_2_im = tmp_h_2_re * gblip_im;

  *h_a_3_re = tmp_h_3_re * gblip_re;
  *h_a_3_im = tmp_h_3_re * gblip_im;
  

  // HESSIAN OF THE TERM exp(1j*angle(btot))
  double V=btot_mag;
  double V_re=V;
  double V_im=0.0;

  double h_exppha_1_re=0.0,h_exppha_1_im=0.0;
  double h_exppha_2_re=0.0,h_exppha_2_im=0.0;
  double h_exppha_3_re=0.0,h_exppha_3_im=0.0;

  double dexppha_dre_j_re=0.0, dexppha_dre_j_im=0.0;
  double dexppha_dim_j_re=0.0, dexppha_dim_j_im=0.0;

  double dexppha_dre_k_re=0.0, dexppha_dre_k_im=0.0;
  double dexppha_dim_k_re=0.0, dexppha_dim_k_im=0.0;

  if(V>0){

    dexppha_dre_j_re = sumsinc/btot_mag*( b1maps_j_re - cos_btot_pha*tmp_j_re );
    dexppha_dre_j_im = sumsinc/btot_mag*( b1maps_j_im - sin_btot_pha*tmp_j_re );  
    dexppha_dim_j_re = sumsinc/btot_mag*( -b1maps_j_im + cos_btot_pha*tmp_j_im );
    dexppha_dim_j_im = sumsinc/btot_mag*(  b1maps_j_re + sin_btot_pha*tmp_j_im );

    dexppha_dre_k_re = sumsinc/btot_mag*( b1maps_k_re - cos_btot_pha*tmp_k_re );
    dexppha_dre_k_im = sumsinc/btot_mag*( b1maps_k_im - sin_btot_pha*tmp_k_re );  
    dexppha_dim_k_re = sumsinc/btot_mag*( -b1maps_k_im + cos_btot_pha*tmp_k_im );
    dexppha_dim_k_im = sumsinc/btot_mag*(  b1maps_k_re + sin_btot_pha*tmp_k_im );

    // term #1
    double U_re = sumsinc*( b1maps_j_re - cos_btot_pha*tmp_j_re );
    double U_im = sumsinc*( b1maps_j_im - sin_btot_pha*tmp_j_re );

    double dU_re = ( -dexppha_dre_k_re*tmp_j_re - cos_btot_pha*(b1maps_j_re*dexppha_dre_k_re + b1maps_j_im*dexppha_dre_k_im) )*sumsinc;
    double dU_im = ( -dexppha_dre_k_im*tmp_j_re - sin_btot_pha*(b1maps_j_re*dexppha_dre_k_re + b1maps_j_im*dexppha_dre_k_im) )*sumsinc;

    double dV_re = sumsinc*tmp_k_re;
    double dV_im = 0.0;
    
    h_exppha_1_re = ( dU_re*V_re - dU_im*V_im - dV_re*U_re + dV_im*U_im )/(V*V);
    h_exppha_1_im = ( dU_im*V_re + dU_re*V_im - dV_im*U_re - dV_re*U_im )/(V*V);
    
    // term #3
    U_re = sumsinc*( -b1maps_j_im + cos_btot_pha*tmp_j_im );
    U_im = sumsinc*(  b1maps_j_re + sin_btot_pha*tmp_j_im );

    dU_re = sumsinc*( dexppha_dim_k_re*tmp_j_im + cos_btot_pha*(b1maps_j_im*dexppha_dim_k_re - b1maps_j_re*dexppha_dim_k_im) );
    dU_im = sumsinc*( dexppha_dim_k_im*tmp_j_im + sin_btot_pha*(b1maps_j_im*dexppha_dim_k_re - b1maps_j_re*dexppha_dim_k_im) );

    dV_re = -sumsinc*tmp_k_im;
    dV_im = 0.0;

    h_exppha_3_re = ( dU_re*V_re - dU_im*V_im - dV_re*U_re + dV_im*U_im )/(V*V);
    h_exppha_3_im = ( dU_im*V_re + dU_re*V_im - dV_im*U_re - dV_re*U_im )/(V*V);
    
    // term #2
    U_re = sumsinc*( -b1maps_j_im + cos_btot_pha*tmp_j_im );
    U_im = sumsinc*(  b1maps_j_re + sin_btot_pha*tmp_j_im );

    dU_re = sumsinc*( dexppha_dre_k_re*tmp_j_im + cos_btot_pha*(b1maps_j_im*dexppha_dre_k_re - b1maps_j_re*dexppha_dre_k_im) );
    dU_im = sumsinc*( dexppha_dre_k_im*tmp_j_im + sin_btot_pha*(b1maps_j_im*dexppha_dre_k_re - b1maps_j_re*dexppha_dre_k_im) );

    dV_re = sumsinc*tmp_k_re;
    dV_im = 0.0;

    h_exppha_2_re = ( dU_re*V_re - dU_im*V_im - dV_re*U_re + dV_im*U_im )/(V*V);
    h_exppha_2_im = ( dU_im*V_re + dU_re*V_im - dV_im*U_re - dV_re*U_im )/(V*V);
  }


  // HESSIAN OF BETA

  // term #1
  double tmp1_re =  0.5*cosphi*( dphi_dre_k*dexppha_dre_j_im + dphi_dre_j*dexppha_dre_k_im );
  double tmp1_im = -0.5*cosphi*( dphi_dre_k*dexppha_dre_j_re + dphi_dre_j*dexppha_dre_k_re );
                                         
  double tmp2_re =  sinphi*h_exppha_1_im;
  double tmp2_im = -sinphi*h_exppha_1_re;

  double tmp3_re = -0.25*sinphi*dphi_dre_k*dphi_dre_j*sin_btot_pha;
  double tmp3_im =  0.25*sinphi*dphi_dre_k*dphi_dre_j*cos_btot_pha;

  double tmp4_re =  0.5*cosphi*h_phi_1*sin_btot_pha;
  double tmp4_im = -0.5*cosphi*h_phi_1*cos_btot_pha;

  tmp_h_1_re=tmp1_re + tmp2_re + tmp3_re + tmp4_re;
  tmp_h_1_im=tmp1_im + tmp2_im + tmp3_im + tmp4_im;
  
  // term #3
  tmp1_re =  0.5*cosphi*( dphi_dim_k*dexppha_dim_j_im + dphi_dim_j*dexppha_dim_k_im );
  tmp1_im = -0.5*cosphi*( dphi_dim_k*dexppha_dim_j_re + dphi_dim_j*dexppha_dim_k_re );

  tmp2_re =  sinphi*h_exppha_3_im;
  tmp2_im = -sinphi*h_exppha_3_re;                                         

  tmp3_re = -0.25*sinphi*dphi_dim_k*dphi_dim_j*sin_btot_pha;
  tmp3_im =  0.25*sinphi*dphi_dim_k*dphi_dim_j*cos_btot_pha;

  tmp4_re =  0.5*cosphi*h_phi_3*sin_btot_pha;
  tmp4_im = -0.5*cosphi*h_phi_3*cos_btot_pha;

  tmp_h_3_re=tmp1_re + tmp2_re + tmp3_re + tmp4_re;
  tmp_h_3_im=tmp1_im + tmp2_im + tmp3_im + tmp4_im;
  
  // term #2
  tmp1_re =  0.5*cosphi*( dphi_dre_k*dexppha_dim_j_im + dphi_dim_j*dexppha_dre_k_im );
  tmp1_im = -0.5*cosphi*( dphi_dre_k*dexppha_dim_j_re + dphi_dim_j*dexppha_dre_k_re );

  tmp2_re =  sinphi*h_exppha_2_im;
  tmp2_im = -sinphi*h_exppha_2_re;                          

  tmp3_re = -0.25*sinphi*dphi_dre_k*dphi_dim_j*sin_btot_pha;
  tmp3_im =  0.25*sinphi*dphi_dre_k*dphi_dim_j*cos_btot_pha;               

  tmp4_re =  0.5*cosphi*h_phi_2*sin_btot_pha;
  tmp4_im = -0.5*cosphi*h_phi_2*cos_btot_pha;

  tmp_h_2_re=tmp1_re + tmp2_re + tmp3_re + tmp4_re;
  tmp_h_2_im=tmp1_im + tmp2_im + tmp3_im + tmp4_im;


  // incorporate gradient blip into the hessian of b
  *h_b_1_re = tmp_h_1_re*gblip_re + tmp_h_1_im*gblip_im;
  *h_b_1_im = tmp_h_1_im*gblip_re - tmp_h_1_re*gblip_im;

  *h_b_2_re = tmp_h_2_re*gblip_re + tmp_h_2_im*gblip_im;
  *h_b_2_im = tmp_h_2_im*gblip_re - tmp_h_2_re*gblip_im;

  *h_b_3_re = tmp_h_3_re*gblip_re + tmp_h_3_im*gblip_im;
  *h_b_3_im = tmp_h_3_im*gblip_re - tmp_h_3_re*gblip_im;  

}





__global__ void hessian_LFA_noB0_cudaKernel( double *a_re,double *a_im,double *b_re,double *b_im,
					     double *ha_re,double *ha_im,
					     double *hb_re,double *hb_im,
					     double *ja_re_re,double *ja_re_im,double *ja_im_re,double *ja_im_im,
					     double *jb_re_re,double *jb_re_im,double *jb_im_re,double *jb_im_im,					     
					     double *a_forw_re,double *a_forw_im,double *b_forw_re,double *b_forw_im,
					     double *btotspokes_re,double *btotspokes_im,
					     double *q_gblips_re,double *q_gblips_im,
					     double *b1maps_re,double *b1maps_im,
					     int nnzp,int nspokes,int ncoils,double sumsinc,double deltat )
{


  int index=threadIdx.x + blockIdx.x * blockDim.x;  //  global thread index == voxel index  


  if(index<nnzp){

    double atot_re=1.0;
    double atot_im=0.0;
    
    double btot_re=0.0;
    double btot_im=0.0;

    a_forw_re[ index ] = 1.0;
    a_forw_im[ index ] = 0.0;
    
    b_forw_re[ index ] = 0.0;
    b_forw_im[ index ] = 0.0;
    
     
    int ncunk=nspokes*ncoils;

    double a_1_re,a_1_im,b_1_re,b_1_im,b1_1_mag,b1_1_pha,cos_b1_1_pha,sin_b1_1_pha,phi_1,cos_phi_1,sin_phi_1;
    
    // FORWARD BLOCH SIMULATION    
    for(int sp=0;sp<nspokes;sp++){
      
      // CK parameters
      compute_caley_klein_params( btotspokes_re[index+sp*nnzp],btotspokes_im[index+sp*nnzp],deltat,
				  q_gblips_re[index+sp*nnzp],q_gblips_im[index+sp*nnzp],1,
				  &b1_1_mag,&b1_1_pha,&cos_b1_1_pha,&sin_b1_1_pha,
				  &phi_1,&cos_phi_1,&sin_phi_1,
				  &a_1_re,&a_1_im,&b_1_re,&b_1_im );
      
      // forward Q products (one step late)
      if(sp<nspokes-1){

	multiply_2spinors( atot_re ,atot_im ,btot_re ,btot_im  ,
			   a_1_re  ,a_1_im  ,b_1_re  ,b_1_im   ,
			   &atot_re,&atot_im,&btot_re,&btot_im );

	a_forw_re[ index+(sp+1)*nnzp ] = atot_re;
	a_forw_im[ index+(sp+1)*nnzp ] = atot_im;

	b_forw_re[ index+(sp+1)*nnzp ] = btot_re;
	b_forw_im[ index+(sp+1)*nnzp ] = btot_im;
	
      }
                  
    }  // for(int sp=0;sp<nspokes;sp++)


    // BACKWARD BLOCH SIMULATION         
    atot_re=1.0;
    atot_im=0.0;
    
    btot_re=0.0;
    btot_im=0.0;
   
    double a_2_re,a_2_im,b_2_re,b_2_im,b1_2_mag,b1_2_pha,cos_b1_2_pha,sin_b1_2_pha,phi_2,cos_phi_2,sin_phi_2;
    double a_3_re,a_3_im,b_3_re,b_3_im,b1_3_mag,b1_3_pha,cos_b1_3_pha,sin_b1_3_pha,phi_3,cos_phi_3,sin_phi_3;

    double da_1_dre_re,da_1_dre_im,da_1_dim_re,da_1_dim_im;
    double db_1_dre_re,db_1_dre_im,db_1_dim_re,db_1_dim_im;

    double A_TMP1_RE_re[8],A_TMP1_RE_im[8];
    double B_TMP1_RE_re[8],B_TMP1_RE_im[8];
    
    double A_TMP1_IM_re[8],A_TMP1_IM_im[8];
    double B_TMP1_IM_re[8],B_TMP1_IM_im[8];

    double A_TMP2_RE_re[8],A_TMP2_RE_im[8];
    double B_TMP2_RE_re[8],B_TMP2_RE_im[8];
    
    double A_TMP2_IM_re[8],A_TMP2_IM_im[8];
    double B_TMP2_IM_re[8],B_TMP2_IM_im[8];

    for(int sp=nspokes-1;sp>=0;sp--){

      compute_caley_klein_params( btotspokes_re[index+sp*nnzp],btotspokes_im[index+sp*nnzp],deltat,
				  q_gblips_re[index+sp*nnzp],q_gblips_im[index+sp*nnzp],1,
				  &b1_1_mag,&b1_1_pha,&cos_b1_1_pha,&sin_b1_1_pha,
				  &phi_1,&cos_phi_1,&sin_phi_1,
				  &a_1_re,&a_1_im,&b_1_re,&b_1_im );
            
      int ind=index + sp*nnzp;
      
      for(int j=0;j<ncoils;j++){
	
	int ind_j=j + sp*ncoils;
	int ind2=index + (j+sp*ncoils)*nnzp;

	// Jacobian (sandwitch product of backward Q product, Jacobian of current Q matrix and forward Q products) 
	compute_jacobian_Q_matrix( b1maps_re[index+j*nnzp],b1maps_im[index+j*nnzp],
				   cos_b1_1_pha,sin_b1_1_pha,
				   cos_phi_1,sin_phi_1,
				   b1_1_mag,deltat,sumsinc,
				   q_gblips_re[index+sp*nnzp],q_gblips_im[index+sp*nnzp],
				   &da_1_dre_re,&da_1_dre_im,&da_1_dim_re,&da_1_dim_im ,
				   &db_1_dre_re,&db_1_dre_im,&db_1_dim_re,&db_1_dim_im );

	
	multiply_3spinors( a_forw_re[ind] ,a_forw_im[ind] ,b_forw_re[ind] ,b_forw_im[ind] ,
			   da_1_dre_re    ,da_1_dre_im    ,db_1_dre_re    ,db_1_dre_im    ,
			   atot_re        ,atot_im        ,btot_re        ,btot_im        ,
			   &ja_re_re[ind2],&ja_re_im[ind2],&jb_re_re[ind2],&jb_re_im[ind2]);

	multiply_3spinors( a_forw_re[ind] ,a_forw_im[ind] ,b_forw_re[ind] ,b_forw_im[ind] ,
			   da_1_dim_re    ,da_1_dim_im    ,db_1_dim_re    ,db_1_dim_im    ,
			   atot_re        ,atot_im        ,btot_re        ,btot_im        ,
			   &ja_im_re[ind2],&ja_im_im[ind2],&jb_im_re[ind2],&jb_im_im[ind2]);


	// Hessian term (within a single spoke)
	for(int k=0;k<ncoils;k++){
	  int ind_k=k + sp*ncoils;

	  double h_a_1_re,h_a_1_im,h_b_1_re,h_b_1_im;
	  double h_a_2_re,h_a_2_im,h_b_2_re,h_b_2_im;
	  double h_a_3_re,h_a_3_im,h_b_3_re,h_b_3_im;

	  compute_hessian_Q_matrix(cos_phi_1,sin_phi_1,
				   b1maps_re[index+j*nnzp],b1maps_im[index+j*nnzp],
				   b1maps_re[index+k*nnzp],b1maps_im[index+k*nnzp],
				   cos_b1_1_pha,sin_b1_1_pha,b1_1_mag,
				   deltat,sumsinc,
				   q_gblips_re[index+sp*nnzp],q_gblips_im[index+sp*nnzp],
				   &h_a_1_re,&h_a_1_im,&h_b_1_re,&h_b_1_im,
				   &h_a_2_re,&h_a_2_im,&h_b_2_re,&h_b_2_im,
				   &h_a_3_re,&h_a_3_im,&h_b_3_re,&h_b_3_im );

	  
	  double tmp_a_1_re,tmp_a_1_im,tmp_b_1_re,tmp_b_1_im;
	  double tmp_a_2_re,tmp_a_2_im,tmp_b_2_re,tmp_b_2_im;
	  double tmp_a_3_re,tmp_a_3_im,tmp_b_3_re,tmp_b_3_im;

	  multiply_3spinors( a_forw_re[ind],a_forw_im[ind],b_forw_re[ind],b_forw_im[ind],
			     h_a_1_re      ,h_a_1_im      ,h_b_1_re      ,h_b_1_im      ,
			     atot_re       ,atot_im       ,btot_re       ,btot_im       ,
			     &tmp_a_1_re   ,&tmp_a_1_im   ,&tmp_b_1_re   ,&tmp_b_1_im   );

	  multiply_3spinors( a_forw_re[ind],a_forw_im[ind],b_forw_re[ind],b_forw_im[ind],
			     h_a_2_re      ,h_a_2_im      ,h_b_2_re      ,h_b_2_im      ,
			     atot_re       ,atot_im       ,btot_re       ,btot_im       ,
			     &tmp_a_2_re   ,&tmp_a_2_im   ,&tmp_b_2_re   ,&tmp_b_2_im   );

	  multiply_3spinors( a_forw_re[ind],a_forw_im[ind],b_forw_re[ind],b_forw_im[ind],
			     h_a_3_re      ,h_a_3_im      ,h_b_3_re      ,h_b_3_im      ,
			     atot_re       ,atot_im       ,btot_re       ,btot_im       ,
			     &tmp_a_3_re   ,&tmp_a_3_im   ,&tmp_b_3_re   ,&tmp_b_3_im   );
          
	  // fill ha
	  ha_re[index + (ind_k       + ind_j        *2*ncunk)*nnzp] = tmp_a_1_re;
	  ha_re[index + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp] = tmp_a_2_re;
	  ha_re[index + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_a_2_re;
	  ha_re[index + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_a_3_re;                    

	  ha_im[index + (ind_k       + ind_j        *2*ncunk)*nnzp] = tmp_a_1_im;
	  ha_im[index + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp] = tmp_a_2_im;
	  ha_im[index + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_a_2_im;
	  ha_im[index + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_a_3_im;                    
	   
	  // fill hb
	  hb_re[index + (ind_k       + ind_j        *2*ncunk)*nnzp] = tmp_b_1_re;
	  hb_re[index + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp] = tmp_b_2_re;
	  hb_re[index + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_b_2_re;
	  hb_re[index + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_b_3_re;

	  hb_im[index + (ind_k       + ind_j        *2*ncunk)*nnzp] = tmp_b_1_im;
	  hb_im[index + (ind_j+ncunk + ind_k        *2*ncunk)*nnzp] = tmp_b_2_im;
	  hb_im[index + (ind_k       + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_b_2_im;
	  hb_im[index + (ind_k+ncunk + (ind_j+ncunk)*2*ncunk)*nnzp] = tmp_b_3_im;

	}
        
	
	// Jacobian terms (derivatives taken for different spokes)
	if(nspokes>0){
	  
	  // product between left-hand jacobian term and backward Q products
	  multiply_2spinors( da_1_dre_re     ,da_1_dre_im     ,db_1_dre_re     ,db_1_dre_im      ,
			     atot_re         ,atot_im         ,btot_re         ,btot_im          ,
			     &A_TMP1_RE_re[j],&A_TMP1_RE_im[j],&B_TMP1_RE_re[j],&B_TMP1_RE_im[j] );
	  
	  multiply_2spinors( da_1_dim_re     ,da_1_dim_im     ,db_1_dim_re     ,db_1_dim_im      ,
			     atot_re         ,atot_im         ,btot_re         ,btot_im          ,
			     &A_TMP1_IM_re[j],&A_TMP1_IM_im[j],&B_TMP1_IM_re[j],&B_TMP1_IM_im[j] );	  	 

	}

      }  // if(int j=0;j<ncoils;j++)


      if(nspokes>0){
        
	for(int sp2=sp-1;sp2>=0;sp2--){
	  
	  // compute the product of Q matrices sandwitched between indices "sp" and "sp2"
	  double a_sandw_re=1.0;
	  double a_sandw_im=0.0;

	  double b_sandw_re=0.0;
	  double b_sandw_im=0.0;

	  for(int sp3=sp-1;sp3>=sp2+1;sp3--){

	    compute_caley_klein_params( btotspokes_re[index+sp3*nnzp],btotspokes_im[index+sp3*nnzp],deltat,
					q_gblips_re[index+sp3*nnzp],q_gblips_im[index+sp3*nnzp],1,
					&b1_3_mag,&b1_3_pha,&cos_b1_3_pha,&sin_b1_3_pha,
					&phi_3,&cos_phi_3,&sin_phi_3,
					&a_3_re,&a_3_im,&b_3_re,&b_3_im );

	    multiply_2spinors(a_3_re     ,a_3_im     ,b_3_re     ,b_3_im      ,
			      a_sandw_re ,a_sandw_im ,b_sandw_re ,b_sandw_im  ,
			      &a_sandw_re,&a_sandw_im,&b_sandw_re,&b_sandw_im );
	  }			      
	  

	  // product between sandwitch term, right-hand jacobian term and forward Q products 
	  for(int j=0;j<ncoils;j++){
	    ind=index + (j+sp2*ncoils)*nnzp;

	    compute_caley_klein_params( btotspokes_re[index+sp2*nnzp],btotspokes_im[index+sp2*nnzp],deltat,
					q_gblips_re[index+sp2*nnzp],q_gblips_im[index+sp2*nnzp],1,
					&b1_2_mag,&b1_2_pha,&cos_b1_2_pha,&sin_b1_2_pha,
					&phi_2,&cos_phi_2,&sin_phi_2,
					&a_2_re,&a_2_im,&b_2_re,&b_2_im );

	    double da_2_dre_re,da_2_dre_im,da_2_dim_re,da_2_dim_im;
	    double db_2_dre_re,db_2_dre_im,db_2_dim_re,db_2_dim_im;

	    compute_jacobian_Q_matrix( b1maps_re[index+j*nnzp],b1maps_im[index+j*nnzp],
				       cos_b1_2_pha,sin_b1_2_pha,
				       cos_phi_2,sin_phi_2,
				       b1_2_mag,deltat,sumsinc,
				       q_gblips_re[index+sp2*nnzp],q_gblips_im[index+sp2*nnzp],
				       &da_2_dre_re,&da_2_dre_im,&da_2_dim_re,&da_2_dim_im ,
				       &db_2_dre_re,&db_2_dre_im,&db_2_dim_re,&db_2_dim_im );

	    multiply_3spinors( a_forw_re[index+sp2*nnzp],a_forw_im[index+sp2*nnzp],b_forw_re[index+sp2*nnzp],b_forw_im[index+sp2*nnzp],
			       da_2_dre_re              ,da_2_dre_im              ,db_2_dre_re              ,db_2_dre_im              ,
			       a_sandw_re               ,a_sandw_im               ,b_sandw_re               ,b_sandw_im               ,
			       &A_TMP2_RE_re[j]         ,&A_TMP2_RE_im[j]         ,&B_TMP2_RE_re[j]         ,&B_TMP2_RE_im[j]         );

	    multiply_3spinors( a_forw_re[index+sp2*nnzp],a_forw_im[index+sp2*nnzp],b_forw_re[index+sp2*nnzp],b_forw_im[index+sp2*nnzp],
			       da_2_dim_re              ,da_2_dim_im              ,db_2_dim_re              ,db_2_dim_im              ,
			       a_sandw_re               ,a_sandw_im               ,b_sandw_re               ,b_sandw_im               ,
			       &A_TMP2_IM_re[j]         ,&A_TMP2_IM_im[j]         ,&B_TMP2_IM_re[j]         ,&B_TMP2_IM_im[j]         );
	  }
	  
	  // outter products between the columns of the TMP1 and TMP2 jacobian terms
	  double tmp_a_1_re,tmp_a_1_im,tmp_b_1_re,tmp_b_1_im;
	  double tmp_a_2_re,tmp_a_2_im,tmp_b_2_re,tmp_b_2_im;
	  double tmp_a_3_re,tmp_a_3_im,tmp_b_3_re,tmp_b_3_im;
	  double tmp_a_4_re,tmp_a_4_im,tmp_b_4_re,tmp_b_4_im;

	  for(int j=0;j<ncoils;j++){
	  
	    for(int k=0;k<ncoils;k++){
             	      
	      multiply_2spinors( A_TMP2_RE_re[k],A_TMP2_RE_im[k],B_TMP2_RE_re[k],B_TMP2_RE_im[k] ,
				 A_TMP1_RE_re[j],A_TMP1_RE_im[j],B_TMP1_RE_re[j],B_TMP1_RE_im[j] ,
				 &tmp_a_1_re    ,&tmp_a_1_im    ,&tmp_b_1_re    ,&tmp_b_1_im     );

	      multiply_2spinors( A_TMP2_RE_re[k],A_TMP2_RE_im[k],B_TMP2_RE_re[k],B_TMP2_RE_im[k] ,
				 A_TMP1_IM_re[j],A_TMP1_IM_im[j],B_TMP1_IM_re[j],B_TMP1_IM_im[j] ,
				 &tmp_a_2_re    ,&tmp_a_2_im    ,&tmp_b_2_re    ,&tmp_b_2_im     );

	      multiply_2spinors( A_TMP2_IM_re[k],A_TMP2_IM_im[k],B_TMP2_IM_re[k],B_TMP2_IM_im[k] ,
				 A_TMP1_RE_re[j],A_TMP1_RE_im[j],B_TMP1_RE_re[j],B_TMP1_RE_im[j] ,
				 &tmp_a_3_re    ,&tmp_a_3_im    ,&tmp_b_3_re    ,&tmp_b_3_im     );

	      multiply_2spinors( A_TMP2_IM_re[k],A_TMP2_IM_im[k],B_TMP2_IM_re[k],B_TMP2_IM_im[k] ,
				 A_TMP1_IM_re[j],A_TMP1_IM_im[j],B_TMP1_IM_re[j],B_TMP1_IM_im[j] ,
				 &tmp_a_4_re    ,&tmp_a_4_im    ,&tmp_b_4_re    ,&tmp_b_4_im     );


	      // half of the Hessian terms are filled by symmetry
	      int ind_j2 = j + sp*ncoils;
	      int ind_k3 = k + sp2*ncoils;

	      ha_re[index + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = ha_re[index + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_a_1_re;
	      ha_re[index + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha_re[index + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_a_2_re;
	      ha_re[index + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = ha_re[index + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_3_re;
	      ha_re[index + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha_re[index + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_4_re;

	      ha_im[index + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = ha_im[index + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_a_1_im;
	      ha_im[index + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha_im[index + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_a_2_im;
	      ha_im[index + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = ha_im[index + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_3_im;
	      ha_im[index + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = ha_im[index + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_a_4_im;


	      hb_re[index + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = hb_re[index + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_b_1_re;
	      hb_re[index + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb_re[index + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_b_2_re;
	      hb_re[index + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = hb_re[index + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_3_re;
	      hb_re[index + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb_re[index + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_4_re;

	      hb_im[index + (ind_k3       + ind_j2        *2*ncunk)*nnzp] = hb_im[index + (ind_j2       + ind_k3        *2*ncunk)*nnzp] = tmp_b_1_im;
	      hb_im[index + (ind_k3       + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb_im[index + (ind_j2+ncunk + ind_k3        *2*ncunk)*nnzp] = tmp_b_2_im;
	      hb_im[index + (ind_k3+ncunk + ind_j2        *2*ncunk)*nnzp] = hb_im[index + (ind_j2       + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_3_im;
	      hb_im[index + (ind_k3+ncunk + (ind_j2+ncunk)*2*ncunk)*nnzp] = hb_im[index + (ind_j2+ncunk + (ind_k3+ncunk)*2*ncunk)*nnzp] = tmp_b_4_im;

	      
	    }  // for(int k=0;k<ncoils;k++)
	  
	  }  // for(int j=0;j<ncoils;j++)

          
	}  // for(int sp2=sp-1;sp>=0;sp2--)
	
      }  // if(sp>0)
      

      // update backward Q products
      multiply_2spinors( a_1_re  ,a_1_im  ,b_1_re  ,b_1_im   , 
			 atot_re ,atot_im ,btot_re ,btot_im  ,
			 &atot_re,&atot_im,&btot_re,&btot_im );
      
    }  // for(int sp=nspokes-1;sp>=0;sp++)


    a_re[index] = atot_re;
    a_im[index] = atot_im;

    b_re[index] = btot_re;
    b_im[index] = btot_im;
    

  }  // if(index<nnzp)

    
}







