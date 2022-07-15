


#include <cmath>


#define TWO_PI_GAMMA (2.675128976384781e+08)




__device__ void multiply_2spinors(float a1_re,float a1_im,float b1_re,float b1_im,
				  float a2_re,float a2_im,float b2_re,float b2_im,
				  float *a3_re,float *a3_im,float *b3_re,float *b3_im)
{
  *a3_re = a2_re*a1_re - a2_im*a1_im - b2_re*b1_re - b2_im*b1_im;
  *a3_im = a2_im*a1_re + a2_re*a1_im + b2_im*b1_re - b2_re*b1_im;
  
  *b3_re = b2_re*a1_re - b2_im*a1_im + a2_re*b1_re + a2_im*b1_im;
  *b3_im = b2_im*a1_re + b2_re*a1_im - a2_im*b1_re + a2_re*b1_im;
}




__device__ void multiply_3spinors(float a1_re,float a1_im,float b1_re,float b1_im,
				  float a2_re,float a2_im,float b2_re,float b2_im,
				  float a3_re,float a3_im,float b3_re,float b3_im,
				  float *a4_re,float *a4_im,float *b4_re,float *b4_im)
{
  
  float tmp1_re = a2_re*a1_re - a2_im*a1_im - b2_re*b1_re - b2_im*b1_im;
  float tmp1_im = a2_im*a1_re + a2_re*a1_im + b2_im*b1_re - b2_re*b1_im;
  
  float tmp2_re = b2_re*a1_re - b2_im*a1_im + a2_re*b1_re + a2_im*b1_im;
  float tmp2_im = b2_im*a1_re + b2_re*a1_im - a2_im*b1_re + a2_re*b1_im;
  
  *a4_re = a3_re*tmp1_re - a3_im*tmp1_im - b3_re*tmp2_re - b3_im*tmp2_im;
  *a4_im = a3_im*tmp1_re + a3_re*tmp1_im + b3_im*tmp2_re - b3_re*tmp2_im;
  
  *b4_re = b3_re*tmp1_re - b3_im*tmp1_im + a3_re*tmp2_re + a3_im*tmp2_im;
  *b4_im = b3_im*tmp1_re + b3_re*tmp1_im - a3_im*tmp2_re + a3_re*tmp2_im;    
}




// CUDA kernel
__global__ void f0_LFA_cudaKernel(float *a_re,float *a_im,float *b_re,float *b_im,
				  float *b1s_re,float *b1s_im,float *grad,
				  float *time_to_sinc_time,float *time_to_spoke,float *subpulse,
				  float *b1maps_re,float *b1maps_im,float *target_LFA_re,float *target_LFA_im,
				  int nnzp,int ntimes,int ncoils,int nspokes,float deltat,int comp_grad)
{


  int index=threadIdx.x + blockIdx.x * blockDim.x;  //  global thread index == voxel index


  if(index<nnzp){

    a_re[0]=1.0;
    a_im[0]=0.0;

    b_re[0]=0.0;
    b_im[0]=0.0;

    
    // FORWARD BLOCH SIMULATION
    for(int time=0;time<ntimes;time++){
      
      float norm=b1s_re[time*nnzp+index]*b1s_re[time*nnzp+index] + b1s_im[time*nnzp+index]*b1s_im[time*nnzp+index] + grad[time*nnzp+index]*grad[time*nnzp+index];
      norm=sqrt(norm);
      
      float a2_re,a2_im,b2_re,b2_im;
      float n1,n2,n3,phi,cosphi,sinphi;
      if(norm==0){
	
	a2_re=1.0;
	a2_im=0.0;
	
	b2_re=0.0;
	b2_im=0.0;	
	
      }else{
	
	phi=-deltat*norm*TWO_PI_GAMMA;
	
	n1=b1s_re[time*nnzp+index] / norm;  // not good to read these from global memory!! IMPROVE!!
	n2=b1s_im[time*nnzp+index] / norm;
	n3=grad[time*nnzp+index] / norm;    
	
	// new CK parameters
	sincos( phi/2.0,&sinphi,&cosphi );
	
	a2_re=cosphi;
	a2_im=-n3*sinphi;
	
	b2_re=n2*sinphi;              
	b2_im=-n1*sinphi;              
	
	// forward Q product
	multiply_2spinors( a_re[index],a_im[index],b_re[index],b_im[index],
			   a2_re,a2_im,b2_re,b2_im,
			   &a_re[index],&a_im[index],&b_re[index],&b_im[index]);
      }      
      
    }  // for(time=0;times<ntimes;time++)     
    
  }  //  if(index<nnzp)
  
}




















