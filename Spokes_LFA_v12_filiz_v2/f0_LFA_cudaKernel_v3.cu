


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
				  double *a_forw_re,double *a_forw_im,double *b_forw_re,double *b_forw_im,
				  double *b1s_re,double *b1s_im,double *grad,
				  int *time_to_sinc_time,int *time_to_spoke,double *subpulse,
				  double *b1maps_re,double *b1maps_im,
				  int nnzp,int ntimes,int ncoils,int nspokes,double deltat,int comp_grad)
{

  int index=threadIdx.x + blockIdx.x * blockDim.x;  //  global thread index == voxel index


  double a2_re,a2_im,b2_re,b2_im;
  double n1,n2,n3,phi,cosphi,sinphi;
      

  if(index<nnzp){

    double atot_re=1.0;
    double atot_im=0.0;

    double btot_re=0.0;
    double btot_im=0.0;

    if( comp_grad==1 ){
      a_forw_re[index]=1.0;
      a_forw_im[index]=0.0;
      
      b_forw_re[index]=0.0;
      b_forw_im[index]=0.0;
    }

    
    // FORWARD BLOCH SIMULATION    
    for(int time=0;time<ntimes;time++){
      
      double norm=b1s_re[time*nnzp+index]*b1s_re[time*nnzp+index] + b1s_im[time*nnzp+index]*b1s_im[time*nnzp+index] + grad[time*nnzp+index]*grad[time*nnzp+index];
      norm=sqrt(norm);
      
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
	sincos( phi/2.0,&sinphi,&cosphi );
	
	a2_re=cosphi;
	a2_im=-n3*sinphi;
	
	b2_re=n2*sinphi;              
	b2_im=-n1*sinphi;              
	
	// forward Q product
	multiply_2spinors( atot_re,atot_im,btot_re,btot_im,
			   a2_re,a2_im,b2_re,b2_im,
			   &atot_re,&atot_im,&btot_re,&btot_im );
      }


      // if gradient is required, store one time late forward Q matrix product      
      if( comp_grad==1 && time<ntimes-1 ){
	  a_forw_re[index + (time+1)*nnzp]=atot_re;  // better way to handle memory?
	  a_forw_im[index + (time+1)*nnzp]=atot_im;  
	  
	  b_forw_re[index + (time+1)*nnzp]=btot_re;  
	  b_forw_im[index + (time+1)*nnzp]=btot_im;  	
      }

      
    }  // for(time=0;times<ntimes;time++)     
    

    a_re[index] = atot_re;
    a_im[index] = atot_im;

    b_re[index] = btot_re;
    b_im[index] = btot_im;   
    
    // BACKWARD BLOCH SIMULATION        
    if(comp_grad==1){

      atot_re = 1.0;
      atot_im = 0.0;

      btot_re = 0.0;
      btot_im = 0.0;

      for(int j=0;j<nspokes*ncoils;j++){
	datot_dre_re[index + j*nnzp] = 0.0;
	datot_dre_im[index + j*nnzp] = 0.0;
	datot_dim_re[index + j*nnzp] = 0.0;
	datot_dim_im[index + j*nnzp] = 0.0;

	dbtot_dre_re[index + j*nnzp] = 0.0;
	dbtot_dre_im[index + j*nnzp] = 0.0;
	dbtot_dim_re[index + j*nnzp] = 0.0;
	dbtot_dim_im[index + j*nnzp] = 0.0;
      }
      
      for(int time=ntimes-1;time>=0;time--){            
	
	int sp_num=time_to_spoke[time];
	int sinctime=time_to_sinc_time[time];

	// compute CK parameters
	double norm=b1s_re[time*nnzp+index]*b1s_re[time*nnzp+index] + b1s_im[time*nnzp+index]*b1s_im[time*nnzp+index] + grad[time*nnzp+index]*grad[time*nnzp+index];
	norm = sqrt(norm);

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
	  
	  sincos( phi/2.0,&sinphi,&cosphi );
	  
	  a2_re=cosphi;
	  a2_im=-n3*sinphi;
	  
	  b2_re=n2*sinphi;              
	  b2_im=-n1*sinphi;              	  
	}

	
	if( sp_num>0 && sinctime>0 && norm>0 ){               
	  
	  for(int j=0;j<ncoils;j++){
	      
	    // derivative of phi wrt spokes amplitudes
	    double dphi_dre = -deltat/norm*subpulse[sinctime-1]*TWO_PI_GAMMA * ( b1s_re[index+time*nnzp]*b1maps_re[index+j*nnzp] + b1s_im[index+time*nnzp]*b1maps_im[index+j*nnzp] );
	    double dphi_dim = deltat/norm*subpulse[sinctime-1]*TWO_PI_GAMMA * ( b1s_re[index+time*nnzp]*b1maps_im[index+j*nnzp] - b1s_im[index+time*nnzp]*b1maps_re[index+j*nnzp] );
	    
	    // derivative of n (rotation axis) wrt spokes amplitudes
	    double tmp_re=-subpulse[sinctime-1]/(norm*norm*norm) * ( b1s_re[index+time*nnzp]*b1maps_re[index+j*nnzp] + b1s_im[index+time*nnzp]*b1maps_im[index+j*nnzp] );
	    double tmp_im=-subpulse[sinctime-1]/(norm*norm*norm) * ( b1s_re[index+time*nnzp]*b1maps_im[index+j*nnzp] - b1s_im[index+time*nnzp]*b1maps_re[index+j*nnzp] );
	    
	    double tmp3=subpulse[sinctime-1]/norm;            
	    
	    double dnx_dre=tmp_re*b1s_re[index+time*nnzp]  + tmp3*b1maps_re[index+j*nnzp];
	    double dnx_dim=-tmp_im*b1s_re[index+time*nnzp] - tmp3*b1maps_im[index+j*nnzp];
	    
	    double dny_dre=tmp_re*b1s_im[index+time*nnzp] + tmp3*b1maps_im[index+j*nnzp];
	    double dny_dim=-tmp_im*b1s_im[index+time*nnzp] + tmp3*b1maps_re[index+j*nnzp];
	    
	    double dnz_dre=tmp_re * grad[index+time*nnzp];
	    double dnz_dim=-tmp_im * grad[index+time*nnzp];
	    
	    // derivative of Q 
	    double da_dre_re = -dphi_dre/2.0*sinphi;
	    double da_dre_im = -dphi_dre/2.0*n3*cosphi - sinphi*dnz_dre;
	    
	    double da_dim_re = -dphi_dim/2.0*sinphi;
	    double da_dim_im = -dphi_dim/2.0*n3*cosphi - sinphi*dnz_dim;            
	    
	    double db_dre_re =  dny_dre*sinphi + 0.5*n2*cosphi*dphi_dre;
	    double db_dre_im = -dnx_dre*sinphi - 0.5*n1*cosphi*dphi_dre;                                
	    
	    double db_dim_re =  dny_dim*sinphi + 0.5*n2*cosphi*dphi_dim;
	    double db_dim_im = -dnx_dim*sinphi - 0.5*n1*cosphi*dphi_dim;
	    
	    
	    // update derivatives of total Q matrices (sandwitch product)
	    double tmp3_re,tmp3_im,tmp4_re,tmp4_im;
	    multiply_3spinors( a_forw_re[index+time*nnzp],a_forw_im[index+time*nnzp],b_forw_re[index+time*nnzp],b_forw_im[index+time*nnzp],
			       da_dre_re                 ,da_dre_im                 ,db_dre_re                 ,db_dre_im                 ,
			       atot_re                   ,atot_im                   ,btot_re                   ,btot_im                   ,
			       &tmp3_re                  ,&tmp3_im                  ,&tmp4_re                  ,&tmp4_im                  );                           
	    	    
	    datot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp3_re;
	    datot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp3_im;
	    
	    dbtot_dre_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp4_re;
	    dbtot_dre_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp4_im;
            

	    multiply_3spinors( a_forw_re[index+time*nnzp],a_forw_im[index+time*nnzp],b_forw_re[index+time*nnzp],b_forw_im[index+time*nnzp],
			       da_dim_re                 ,da_dim_im                 ,db_dim_re                 ,db_dim_im                 ,
			       atot_re                   ,atot_im                   ,btot_re                   ,btot_im                   ,
			       &tmp3_re                  ,&tmp3_im                  ,&tmp4_re                  ,&tmp4_im                  );                           
	    
	    datot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp3_re;
	    datot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp3_im;
	    
	    dbtot_dim_re[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp4_re;            
	    dbtot_dim_im[index + j*nnzp + (sp_num-1)*nnzp*ncoils] += tmp4_im;
	    
	  }  // for(j=0;j<ncoils;j++)
          
	}  // if(sp_num>0)
	
	// keep track of the backward Q matrices product
	multiply_2spinors( a2_re   ,a2_im   ,b2_re   ,b2_im   ,
			   atot_re ,atot_im ,btot_re ,btot_im ,
			   &atot_re,&atot_im,&btot_re,&btot_im);
        
      }  // for(time=ntimes-1;time>=0;time--)
      
    }  // if(comp_grad==1)


  }  //  if(index<nnzp)

  
}




















