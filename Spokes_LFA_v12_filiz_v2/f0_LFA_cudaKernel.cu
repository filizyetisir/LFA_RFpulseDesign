


#include <cmath>


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




// CUDA kernel
__global__ void f0_LFA_cudaKernel(double *a_re,double *a_im,double *b_re,double *b_im,
				  double *datot_dre_re,double *datot_dre_im,double *datot_dim_re,double *datot_dim_im,
				  double *dbtot_dre_re,double *dbtot_dre_im,double *dbtot_dim_re,double *dbtot_dim_im,
				  double *da_dre_re,double *da_dre_im,double *da_dim_re,double *da_dim_im,
				  double *db_dre_re,double *db_dre_im,double *db_dim_re,double *db_dim_im,
				  double *a_forw_re,double *a_forw_im,double *b_forw_re,double *b_forw_im,
				  double *as_re,double *as_im,double *bs_re,double *bs_im,
				  double *b1s_re,double *b1s_im,double *grad,
				  double *time_to_sinc_time,double *time_to_spoke,double *subpulse,
				  double *b1maps_re,double *b1maps_im,
				  int nnzp,int ntimes,int ncoils,int nspokes,double deltat,int comp_grad)
{


  int index=threadIdx.x + blockIdx.x * blockDim.x;  //  global thread index == voxel index


  if(index<nnzp){

    
    // FORWARD BLOCH SIMULATION
    
    for(int time=0;time<ntimes;time++){
      
      int sinctime=(int)( time_to_sinc_time[time] );      
      
      double norm=pow(b1s_re[time*nnzp+index],2) + pow(b1s_im[time*nnzp+index],2) + pow(abs(grad[time*nnzp+index]),2);
      norm=sqrt(norm);
      
      double a2_re,a2_im,b2_re,b2_im;
      double n1,n2,n3,phi,cosphi,sinphi;
      if(norm==0){
	
	a2_re=1.0;
	a2_im=0.0;
	
	b2_re=0.0;
	b2_im=0.0;	
	
      }else{
	
	phi=-deltat*norm*TWO_PI_GAMMA;
	
	n1=b1s_re[time*nnzp+index] / norm;
	n2=b1s_im[time*nnzp+index] / norm;
	n3=grad[time*nnzp+index] / norm;    
	
	// new CK parameters
	cosphi=cos(phi/2.0);
	sinphi=sin(phi/2.0);
	
	a2_re=cosphi;
	a2_im=-n3*sinphi;
	
	b2_re=n2*sinphi;              
	b2_re=-n1*sinphi;              
	
	// forward Q product
	multiply_2spinors( a_re[index],a_im[index],b_re[index],b_im[index],
			   a2_re,a2_im,b2_re,b2_im,
			   &a_re[index],&a_im[index],&b_re[index],&b_im[index]);
      }
      
      
      if( comp_grad==1 ){
	
	// store Q matrix
	as_re[index + time*nnzp] = a2_re;
	as_im[index + time*nnzp] = a2_im;
	
	bs_re[index + time*nnzp] = b2_re;
	bs_im[index + time*nnzp] = b2_im;
	
	// store one time late forward Q matrix product
	if( time<ntimes-1 ){
	  a_forw_re[index + (time+1)*nnzp]=a_re[index];  // better way to handle memory?
	  a_forw_im[index + (time+1)*nnzp]=a_im[index];  
	  
	  b_forw_re[index + (time+1)*nnzp]=b_re[index];  
	  b_forw_im[index + (time+1)*nnzp]=b_im[index];  
	}
	
	if( sinctime>0 ){
	  
	  for(int j=0;j<ncoils;j++){
	    
	    if(norm>0){                    
	      
	      // derivative of phi wrt spokes amplitudes
	      double dphi_dre = -deltat/norm*subpulse[sinctime-1]*TWO_PI_GAMMA * ( b1s_re[index+time*nnzp]*b1maps_re[index+j*nnzp] + b1s_im[index+time*nnzp]*b1maps_im[index+j*nnzp] );
	      double dphi_dim = deltat/norm*subpulse[sinctime-1]*TWO_PI_GAMMA * ( b1s_re[index+time*nnzp]*b1maps_im[index+j*nnzp] - b1s_im[index+time*nnzp]*b1maps_re[index+j*nnzp] );
	      
	      // derivative of n (rotation axis) wrt spokes amplitudes
	      double tmp_re=-subpulse[sinctime-1]/pow(norm,3) * ( b1s_re[index+time*nnzp]*b1maps_re[index+j*nnzp] + b1s_im[index+time*nnzp]*b1maps_im[index+j*nnzp] );
	      double tmp_im=-subpulse[sinctime-1]/pow(norm,3) * ( b1s_re[index+time*nnzp]*b1maps_im[index+j*nnzp] - b1s_im[index+time*nnzp]*b1maps_re[index+j*nnzp] );
	      
	      double tmp3=subpulse[sinctime-1]/norm;            
	      
	      double dnx_dre=tmp_re*b1s_re[index+time*nnzp]  + tmp3*b1maps_re[index+j*nnzp];
	      double dnx_dim=-tmp_im*b1s_re[index+time*nnzp] - tmp3*b1maps_im[index+j*nnzp];
	      
	      double dny_dre=tmp_re*b1s_im[index+time*nnzp] + tmp3*b1maps_im[index+j*nnzp];
	      double dny_dim=-tmp_im*b1s_im[index+time*nnzp] + tmp3*b1maps_re[index+j*nnzp];
	      
	      double dnz_dre=tmp_re * grad[index+time*nnzp];
	      double dnz_dim=-tmp_im * grad[index+time*nnzp];
	      
	      // derivative of Q 
	      da_dre_re[index + time*nnzp + j*ntimes*nnzp]=-dphi_dre/2.0*sinphi;
	      da_dre_im[index + time*nnzp + j*ntimes*nnzp]=-dphi_dre/2.0*n3*cosphi - sinphi*dnz_dre;
	      
	      da_dim_re[index + time*nnzp + j*ntimes*nnzp]=-dphi_dim/2.0*sinphi;
	      da_dim_im[index + time*nnzp + j*ntimes*nnzp]=-dphi_dim/2.0*n3*cosphi - sinphi*dnz_dim;            
	      
	      db_dre_re[index + time*nnzp + j*ntimes*nnzp]= dny_dre*sinphi + 0.5*n2*cosphi*dphi_dre;
	      db_dre_im[index + time*nnzp + j*ntimes*nnzp]=-dnx_dre*sinphi - 0.5*n1*cosphi*dphi_dre;                                
	      
	      db_dim_re[index + time*nnzp + j*ntimes*nnzp]= dny_dim*sinphi + 0.5*n2*cosphi*dphi_dim;
	      db_dim_im[index + time*nnzp + j*ntimes*nnzp]=-dnx_dim*sinphi - 0.5*n1*cosphi*dphi_dim;
	      
	    }	  
	    
	  }  // for(j=0;j<ncoils;j++)
	  
	}  // if(sinctime>0)
	
      }  // if(comp_grad==1)            
      
    }  // for(time=0;times<ntimes;time++)     
    
    
    // BACKWARD BLOCH SIMULATION        
    if(comp_grad==1){
      
      a_re[index]=1.0;
      a_im[index]=0.0;
      
      b_re[index]=0.0;
      b_im[index]=0.0;
      
      for(int time=ntimes-1;time>=0;time--){            
	
	int sp_num=(int)( time_to_spoke[time] );

	double tmp3_re,tmp3_im,tmp4_re,tmp4_im;
	if(sp_num>0){               
	  
	  for(int j=0;j<ncoils;j++){
	    
	    // update derivatives of total Q matrices (sandwitch product)
	    multiply_3spinors( a_forw_re[index+time*nnzp],a_forw_im[index+time*nnzp],b_forw_re[index+time*nnzp],b_forw_im[index+time*nnzp],
			       da_dre_re[index + time*nnzp + j*ntimes*nnzp],da_dre_im[index + time*nnzp + j*ntimes*nnzp],db_dre_re[index + time*nnzp + j*ntimes*nnzp],db_dre_im[index + time*nnzp + j*ntimes*nnzp],
			       a_re[index],a_im[index],b_re[index],b_im[index],
			       &tmp3_re,&tmp3_im,&tmp4_re,&tmp4_im);                           
	    
	    datot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3_re;
	    datot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3_im;
	    
	    dbtot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4_re;
	    dbtot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4_im;
            
	    multiply_3spinors( a_forw_re[index+time*nnzp],a_forw_im[index+time*nnzp],b_forw_re[index+time*nnzp],b_forw_im[index+time*nnzp],
			       da_dim_re[index + time*nnzp + j*ntimes*nnzp],da_dim_im[index + time*nnzp + j*ntimes*nnzp],db_dim_re[index + time*nnzp + j*ntimes*nnzp],db_dim_im[index + time*nnzp + j*ntimes*nnzp],
			       a_re[index],a_im[index],b_re[index],b_im[index],
			       &tmp3_re,&tmp3_im,&tmp4_re,&tmp4_im);                           
	    
	    datot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3_re;
	    datot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=datot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp3_im;
	    
	    dbtot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4_re;            
	    dbtot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils]=dbtot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] + tmp4_im;
	    
	  }  // for(j=0;j<ncoils;j++)
          
	}  // if(sp_num>0)
	
	// keep track of the backward Q matrices product
	multiply_2spinors( as_re[index+time*nnzp],as_im[index+time*nnzp],bs_re[index+time*nnzp],bs_im[index+time*nnzp],
			   a_re[index],a_im[index],b_re[index],b_im[index],
			   &a_re[index],&a_im[index],&b_re[index],&b_im[index]);
        
      }  // for(time=ntimes-1;time>=0;time--)
      
    }  // if(comp_grad==1)

  }  //  if(index<nnzp)
  
}




















