
/*
Add description and legal texts
*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>


#define DT double
#define MAX(x,y)  ((x) >= (y) ? (x) : (y))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#define MEM_OPTIMIZED 1
#define THREADSX 32
#define THREADSY THREADSX

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600

#else
  static __inline__ __device__ double atomicAdd(double *address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    if (val==0.0)
      return __longlong_as_double(old);
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
  }


#endif

void Check_CUDA_Error(const char *message)
{
	cudaError_t error = cudaGetLastError();
	if(error!=cudaSuccess) {
		fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
		exit(-1);
	}
}


__device__ void computeVariogram(int i, int  j,const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    float *sh_np,float *sh_dis,float *sh_tm,float *sh_hm,float *sh_gam,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr, int sh_pos){

//    int half_nd = nd/2;
//    float dx,dy,dz;

    float dx,dy,dz,dxs,dys,dzs,hs;
    int id,ii,il,it,iv;
    int lagbeg,lagend,ilag;
    float band,dcazm,dcdec,dxy,vrh,vrhpr,vrt,vrtpr,h;
    int omni;

    dx  = d_x[i] - d_x[j];
    dy  = d_y[i] - d_y[j];
    dz  = d_z[i] - d_z[j];
    dxs = dx*dx;
    dys = dy*dy;
    dzs = dz*dz;
    hs  = dxs + dys + dzs;

    if(hs <= dismxs)
    {
        if(hs < 0.0) hs = 0.0;
        h   = sqrtf(hs);


    //
    // Determine which lag this is and skip if outside the defined distance
    // tolerance:
    //
        if(h<=EPSLON){
            lagbeg = 1;
            lagend = 1;
        }
        else{
            lagbeg = -1;
            lagend = -1;
            for(ilag=2;ilag<=nlag+2;ilag++){
                if(h>=(xlag*(float)(ilag-2)-xltol) && h<=(xlag*(float)(ilag-2)+xltol)){
                    if(lagbeg<0) lagbeg = ilag;
                    lagend = ilag;
                }
            }

        }
        if(lagend>=0)
        {
        //			printf("dx=%f dy=%f dz=%fh=%f lagbeg=%d lagend=%d\n",dx,dy,dz,h,lagbeg,lagend);


        //
        // Definition of the direction corresponding to the current pair. All
        // directions are considered (overlapping of direction tolerance cones
        // is allowed):
        //


            for(id=0;id<ndir;id++){
            //
            // Check for an acceptable azimuth angle:
            //
                dxy = sqrtf(MAX((dxs+dys),0.0));
                if(dxy<EPSLON){
                    dcazm = 1.0;
                }
                else{
                    dcazm = (dx*d_uvxazm[id]+dy*d_uvyazm[id])/dxy;
                }
                if(fabsf(dcazm)>=d_csatol[id])
                {
            //
            // Check the horizontal bandwidth criteria (maximum deviation
            // perpendicular to the specified direction azimuth):
            //
                    band = d_uvxazm[id]*dy - d_uvyazm[id]*dx;
                    if(fabsf(band)<d_bandwh[id])
                    {
                        //fprintf(stdout,"dxy=%f\tdcazm=%f\tband=%f\n",dxy,dcazm,band);


                //
                // Check for an acceptable dip angle:
                //
                        if(dcazm<0.0) dxy = -dxy;
                        if(lagbeg==1)
                            dcdec = 0.0;
                        else{
                            dcdec = (dxy*d_uvhdec[id]+dz*d_uvzdec[id])/h;

                        }
                        band = d_uvhdec[id]*dz - d_uvzdec[id]*dxy;
                        if(fabsf(dcdec)>=d_csdtol[id] && fabsf(band)<=d_bandwd[id])
                        {
                    //
                    // Check the vertical bandwidth criteria (maximum deviation perpendicular
                    // to the specified dip direction):
                    //

                        //
                        // Check whether or not an omni-directional variogram is being computed:
                        //
                                omni = 0;
                                if(d_atol[id]>=90.0) omni = 1;
                        //
                        // This direction is acceptable - go ahead and compute all variograms:
                        //

                            //printf("dxy=%f dcazm=%f uvxazm[0]=%f uvyazm[0]=%f band=%f dcdec=%f omni=%d csdtol[0]=%f\n",dxy,dcazm,uvxazm[0],uvyazm[0],band,dcdec,omni,csdtol[0]);

                        //				fprintf(stdout,"dcazm=%f\tdcdec=%f\n",dcazm,dcdec);

                            for(iv=0;iv<nvarg;iv++){
                    //
                    // For this variogram, sort out which is the tail and the head value:
                    //
                                it = d_ivtype[iv];
                                if(dcazm>=0.0 && dcdec>=0.0){
                                    ii = d_ivtail[iv]-1;
                                    vrh   = d_vr[i+ii*(maxdat)];
                                    ii = d_ivhead[iv]-1;
                                    vrt   = d_vr[j+ii*(maxdat)];
                                    if(omni || it==2){
                                        ii    = d_ivhead[iv]-1;
                                        vrtpr = d_vr[i+ii*(maxdat)];
                                        ii    = d_ivtail[iv]-1;
                                        vrhpr = d_vr[j+ii*(maxdat)];
                                    }
                                }
                                else{
                                    ii = d_ivtail[iv]-1;
                                    vrh   = d_vr[j+ii*(maxdat)];
                                    ii = d_ivhead[iv]-1;
                                    vrt   = d_vr[i+ii*(maxdat)];
                                    if(omni || it==2){
                                        ii    = d_ivhead[iv]-1;
                                        vrtpr = d_vr[j+ii*(maxdat)];
                                        ii    = d_ivtail[iv]-1;
                                        vrhpr = d_vr[i+ii*(maxdat)];
                                    }
                                }
                    //
                    // Reject this pair on the basis of missing values:
                    //
                                if(vrt>=tmin && vrh>=tmin && vrt<=tmax && vrh<=tmax && it!=2 || (vrtpr>=tmin && vrhpr>=tmin && vrtpr<=tmax && vrhpr<=tmax))
                                {
                                    if(it==1 || it==5 || it>=9){
                                        for(il=lagbeg;il<=lagend;il++){
                                            ii = (id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2)+il -1;



                                            atomicAdd(&sh_np[ii + mxdlv*sh_pos],1.0);
                                            atomicAdd(&sh_dis[ii + mxdlv*sh_pos],(h));
                                            atomicAdd(&sh_tm[ii + mxdlv*sh_pos],(vrt));
                                            atomicAdd(&sh_hm[ii + mxdlv*sh_pos],(vrh));
                                            atomicAdd(&sh_gam[ii + mxdlv*sh_pos],((vrh-vrt)*(vrh-vrt)));

                                            if(omni){
                                                if(vrtpr>=tmin && vrhpr>=tmin && vrtpr<tmax && vrhpr<tmax){

                                                    atomicAdd(&sh_np[ii + mxdlv*sh_pos],1.0);
                                                    atomicAdd(&sh_dis[ii + mxdlv*sh_pos],(h));
                                                    atomicAdd(&sh_tm[ii + mxdlv*sh_pos],(vrtpr));
                                                    atomicAdd(&sh_hm[ii + mxdlv*sh_pos],(vrhpr));
                                                    atomicAdd(&sh_gam[ii + mxdlv*sh_pos],((vrhpr-vrtpr)*(vrhpr-vrtpr)));

                                                }
                                            }
                                        }
                                    }

                                    // The Traditional Cross Semivariogram:
                //
                                    else if(it==2){
                                        for(il=lagbeg;il<=lagend;il++){
                                            ii = (id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2)+il -1;
                                            atomicAdd(&sh_np[ii + mxdlv*sh_pos],1.0);
                                            atomicAdd(&sh_dis[ii + mxdlv*sh_pos],(h));
                                            atomicAdd(&sh_tm[ii + mxdlv*sh_pos],(0.5*(vrt+vrtpr)));
                                            atomicAdd(&sh_hm[ii + mxdlv*sh_pos],(0.5*(vrh+vrhpr)));
                                            atomicAdd(&sh_gam[ii + mxdlv*sh_pos],((vrhpr-vrh)*(vrt-vrtpr)));

                                        }
                                    }
				/*

					Note: 
					If new spatial measure are requiered, they must be implemented here following the 
					previous examples, with it=1,2,5,9.

				*/
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}







__device__ void computePointsValues(int idx, int  idy,const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    float *sh_np,float *sh_dis,float *sh_tm,float *sh_hm,float *sh_gam,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr, int sh_pos,int half_nd){

    int i,j;
    j = idx + half_nd;
    i = idy;
    computeVariogram(i,j,nd,irepo,maxdat,MAXVAR,
        d_x,d_y,d_z,
        EPSLON,nlag,xlag,xltol,
        mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
        dismxs,tmax,tmin,ndir,nvarg,
        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
        d_csatol, d_csdtol, d_bandwh, d_bandwd,
        d_atol,
        d_ivtype, d_ivtail, d_ivhead,
        d_vr,sh_pos);

    if (idx > idy){
        i = idy;
        j = idx;

        computeVariogram(i,j,nd,irepo,maxdat,MAXVAR,
            d_x,d_y,d_z,
            EPSLON,nlag,xlag,xltol,
            mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
            dismxs,tmax,tmin,ndir,nvarg,
            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
            d_csatol, d_csdtol, d_bandwh, d_bandwd,
            d_atol,
            d_ivtype, d_ivtail, d_ivhead,
            d_vr,sh_pos);

    } else{
        if (idx == idy){
            i = idy;
            j = idy;

            computeVariogram(i,j,nd,irepo,maxdat,MAXVAR,
                d_x,d_y,d_z,
                EPSLON,nlag,xlag,xltol,
                mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
                dismxs,tmax,tmin,ndir,nvarg,
                d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
                d_csatol, d_csdtol, d_bandwh, d_bandwd,
                d_atol,
                d_ivtype, d_ivtail, d_ivhead,
                d_vr,sh_pos);
        }
        i = idx + half_nd;
        j = idy + half_nd;

        computeVariogram(i,j,nd,irepo,maxdat,MAXVAR,
            d_x,d_y,d_z,
            EPSLON,nlag,xlag,xltol,
            mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
            dismxs,tmax,tmin,ndir,nvarg,
            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
            d_csatol, d_csdtol, d_bandwh, d_bandwd,
            d_atol,
            d_ivtype, d_ivtail, d_ivhead,
            d_vr,sh_pos);
    }


}






__global__ void variogramKernelMemoryOptimized(const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    DT *d_np, DT *d_dis, DT *d_gam, DT *d_hm, DT *d_tm,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr,int chunks_sh_mem,int frac_nd){

    int tidx=threadIdx.x;
    int tidy=threadIdx.y;
    int bidx=blockIdx.x;
    int bidy=blockIdx.y;
    int bdimx=blockDim.x;
    int bdimy=blockDim.y;
    int idx = bidx*bdimx + tidx;
    int idy = bidy*bdimy + tidy;
    int threadId = tidx + bdimx*tidy;
    int half_nd = nd/2;
    int sh_pos = tidy%chunks_sh_mem;
    int i,j;
    int num_threads = bdimx*bdimy;
    extern __shared__ float buffer[];
    float *sh_np = &buffer[0];
    float *sh_dis = &buffer[chunks_sh_mem*mxdlv];
    float *sh_gam = &buffer[2*chunks_sh_mem*mxdlv];
    float *sh_hm = &buffer[3*chunks_sh_mem*mxdlv];
    float *sh_tm = &buffer[4*chunks_sh_mem*mxdlv];
    int init_sh_mem = threadId;

    while (init_sh_mem < chunks_sh_mem*mxdlv){
        sh_np[init_sh_mem] = 0;
        sh_dis[init_sh_mem] = 0.0;
        sh_gam[init_sh_mem] = 0.0;
        sh_hm[init_sh_mem] = 0.0;
        sh_tm[init_sh_mem] = 0.0;
        init_sh_mem += num_threads;
    }

    __syncthreads();

    if (idx < frac_nd && idy < frac_nd){
        for (i = idx; i < half_nd; i += frac_nd){
            for (j = idy; j < half_nd; j += frac_nd){
                computePointsValues(i,j,nd,irepo,maxdat,MAXVAR,
                    d_x,d_y,d_z,
                    EPSLON,nlag,xlag,xltol,
                    mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
                    dismxs,tmax,tmin,ndir,nvarg,
                    d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
                    d_csatol, d_csdtol, d_bandwh, d_bandwd,
                    d_atol,
                    d_ivtype, d_ivtail, d_ivhead,
                    d_vr,sh_pos,half_nd);
            }
        }
    }


    __syncthreads();

    if (threadId < mxdlv){
        sh_np[threadId] += sh_np[threadId + mxdlv] + sh_np[threadId + 2*mxdlv] + sh_np[threadId + 3*mxdlv];
        sh_dis[threadId] += sh_dis[threadId + mxdlv] + sh_dis[threadId + 2*mxdlv] + sh_dis[threadId + 3*mxdlv];
        sh_tm[threadId] += sh_tm[threadId + mxdlv] + sh_tm[threadId + 2*mxdlv] + sh_tm[threadId + 3*mxdlv];
        sh_hm[threadId] += sh_hm[threadId + mxdlv] + sh_hm[threadId + 2*mxdlv] + sh_hm[threadId + 3*mxdlv];
        sh_gam[threadId] += sh_gam[threadId + mxdlv] + sh_gam[threadId + 2*mxdlv] + sh_gam[threadId + 3*mxdlv];

        atomicAdd(&d_np[threadId],sh_np[threadId]);
        atomicAdd(&d_dis[threadId],sh_dis[threadId]);
        atomicAdd(&d_tm[threadId],sh_tm[threadId]);
        atomicAdd(&d_hm[threadId],sh_hm[threadId]);
        atomicAdd(&d_gam[threadId],sh_gam[threadId]);
    }

}







__global__ void variogramKernel(    const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    DT *d_np, DT *d_dis, DT *d_gam, DT *d_hm, DT *d_tm,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr,int frac_nd){

    int tidx=threadIdx.x;
    int tidy=threadIdx.y;
    int bidx=blockIdx.x;
    int bidy=blockIdx.y;
    int bdimx=blockDim.x;
    int bdimy=blockDim.y;
    int idx = bidx*bdimx + tidx;
    int idy = bidy*bdimy + tidy;
    int threadId = tidx + bdimx*tidy;
    int half_nd = nd/2;
    int i,j;
    int num_threads = bdimx*bdimy;
    extern __shared__ float buffer[];
    float *sh_np = &buffer[0];
    float *sh_dis = &buffer[mxdlv];
    float *sh_gam = &buffer[2*mxdlv];
    float *sh_hm = &buffer[3*mxdlv];
    float *sh_tm = &buffer[4*mxdlv];

    int init_sh_mem = threadId;
    while (init_sh_mem < mxdlv){
        sh_np[init_sh_mem] = 0;
        sh_dis[init_sh_mem] = 0.0;
        sh_gam[init_sh_mem] = 0.0;
        sh_hm[init_sh_mem] = 0.0;
        sh_tm[init_sh_mem] = 0.0;
        init_sh_mem += num_threads;
    }

    __syncthreads();

    if (idx < frac_nd && idy < frac_nd){
        for (i = idx; i < half_nd; i += frac_nd){
            for (j = idy; j < half_nd; j += frac_nd){
                computePointsValues(i,j,nd,irepo,maxdat,MAXVAR,
                    d_x,d_y,d_z,
                    EPSLON,nlag,xlag,xltol,
                    mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
                    dismxs,tmax,tmin,ndir,nvarg,
                    d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
                    d_csatol, d_csdtol, d_bandwh, d_bandwd,
                    d_atol,
                    d_ivtype, d_ivtail, d_ivhead,
                    d_vr,0,half_nd);
            }
        }
    }

    __syncthreads();

    if (threadId < mxdlv){

        atomicAdd(&d_np[threadId],sh_np[threadId]);
        atomicAdd(&d_dis[threadId],sh_dis[threadId]);
        atomicAdd(&d_tm[threadId],sh_tm[threadId]);
        atomicAdd(&d_hm[threadId],sh_hm[threadId]);
        atomicAdd(&d_gam[threadId],sh_gam[threadId]);
    }
}





extern "C" int extractstatisticscudawrapper_(
                            //      integer nd,irepo,maxdat,MAXVAR
                                int *nd, int *irepo, int *maxdat, int *MAXVAR,
                            //      real x(maxdat),y(maxdat),z(maxdat)
                                float *x, float *y, float *z,
                            //      real EPSLON
                                float *EPSLON,
                            //      integer nlag
                                int *nlag,
                            //      real xlag,xltol
                                float *xlag, float *xltol,
                            //      integer mxdlv
                                int *mxdlv,
                            //      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
                            //     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
                                double *np, double *dis, double *gam, double *hm, double *tm, double *hv, double *tv,
                            //      integer numThreads
                                int *numThreads,
                            //      real*8 reducedVariables(7,mxdlv,numThreads)
                                double *reducedVariables,
                            //      real dismxs,tmax,tmin
                                float *dismxs, float *tmax, float *tmin,
                            //      integer ndir,nvarg
                                int *ndir, int *nvarg,
                            //      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
                                float *uvxazm, float *uvyazm, float *uvzdec, float *uvhdec,
                            //      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
                                float *csatol, float *csdtol, float *bandwh, float *bandwd,
                            //      real atol(ndir)
                                float *atol,
                            //      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
                                int *ivtype, int *ivtail, int *ivhead,
                            //      real vr(maxdat,MAXVAR)
                                float *vr)
{
	float *d_x,*d_y,*d_z;
	DT *d_np,*d_dis,*d_gam,*d_hm,*d_tm;
	DT *h_np,*h_dis,*h_gam,*h_hm,*h_tm;
	float *d_uvxazm,*d_uvyazm,*d_uvzdec,*d_uvhdec,*d_csatol,*d_csdtol,*d_bandwh,*d_bandwd,*d_atol,*d_vr;
	int *d_ivtype,*d_ivtail,*d_ivhead;
    	cudaSetDevice(0);
    	dim3 threads(THREADSX,THREADSY,1);
	int frac_nd;

    	int chunk_sh_mem = 4;
	h_np = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_dis = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_gam = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_hm = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_tm = (DT*)malloc(sizeof(DT)* *mxdlv);
	int shared_mem_size;
    	int i;
    	for (i = 0; i < *mxdlv; i++){
        	h_np[i] = 0.0;
        	h_dis[i] = 0.0;
        	h_gam[i] = 0.0;
        	h_hm[i] = 0.0;
        	h_tm[i] = 0.0;
    	}
   	cudaMalloc( (void **)&d_x, sizeof(float) * (*maxdat) );
   	//Check_CUDA_Error("malloc coord");
   	cudaMalloc( (void **)&d_y, sizeof(float) * (*maxdat) );
   	//Check_CUDA_Error("malloc coord");
   	cudaMalloc( (void **)&d_z, sizeof(float) * (*maxdat) );
   	//Check_CUDA_Error("malloc coord");
   	cudaMalloc( (void **)&d_np, sizeof(DT) * (*mxdlv) );
   	//Check_CUDA_Error("malloc np, dis, gam, hm, tm");
   	cudaMalloc( (void **)&d_dis, sizeof(DT) * (*mxdlv) );
   	//Check_CUDA_Error("malloc np, dis, gam, hm, tm");
   	cudaMalloc( (void **)&d_gam, sizeof(DT) * (*mxdlv) );
   	//Check_CUDA_Error("malloc np, dis, gam, hm, tm");
   	cudaMalloc( (void **)&d_hm, sizeof(DT) * (*mxdlv) );
   	//Check_CUDA_Error("malloc np, dis, gam, hm, tm");
   	cudaMalloc( (void **)&d_tm, sizeof(DT) * (*mxdlv) );
   	//Check_CUDA_Error("malloc np, dis, gam, hm, tm");
   	cudaMalloc( (void **)&d_uvxazm, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_uvyazm, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_uvzdec, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_uvhdec, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_csatol, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_csdtol, sizeof(float) * (100) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_bandwh, sizeof(float) * (*ndir) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_bandwd, sizeof(float) * (*ndir) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_atol, sizeof(float) * (*ndir) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_vr, sizeof(float) * (*maxdat* *MAXVAR) );
   	//Check_CUDA_Error("small mallocs ");
   	cudaMalloc( (void **)&d_ivtype, sizeof(float) * (*nvarg) );
   	//Check_CUDA_Error("iv mallocs");
   	cudaMalloc( (void **)&d_ivtail, sizeof(float) * (*nvarg) );
   	//Check_CUDA_Error("iv mallocs");
   	cudaMalloc( (void **)&d_ivhead, sizeof(float) * (*nvarg) );
   	//Check_CUDA_Error("iv mallocs");
   	cudaMemcpy( d_x, x,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpy( d_y, y,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpy( d_z, z,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpy( d_np, h_np,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpy( d_dis, h_dis,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpy( d_gam, h_gam,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpy( d_hm, h_hm,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpy( d_tm, h_tm,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpy( d_uvxazm, uvxazm,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_uvyazm, uvyazm,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_uvzdec, uvzdec,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_uvhdec, uvhdec,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_csatol, csatol,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_csdtol, csdtol,sizeof(float) * (100), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_bandwh, bandwh,sizeof(float) * (*ndir), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_bandwd, bandwd,sizeof(float) * (*ndir), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_atol, atol,sizeof(float) * (*ndir), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_vr, vr,sizeof(float) * (*maxdat* *MAXVAR), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpy( d_ivtype, ivtype,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy iv var h -> d");
   	cudaMemcpy( d_ivtail, ivtail,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy iv var h -> d");
   	cudaMemcpy( d_ivhead, ivhead,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice );
   	//Check_CUDA_Error("cpy iv var h -> d");
   	cudaEvent_t start, stop;
   	float time;

    	if (MEM_OPTIMIZED){
        	printf("---------\nmaxdat %d nd %d\n---------\n",*maxdat,*nd);
        	frac_nd = *maxdat/THREADSX;
        	dim3 blocks( (frac_nd + threads.x - 1)/threads.x,(frac_nd + threads.y - 1)/threads.y,1 );
        	shared_mem_size = sizeof(DT)*(*mxdlv*5*chunk_sh_mem);
        	cudaEventCreate(&start);
        	cudaEventCreate(&stop);
        	cudaEventRecord(start, 0);
        	variogramKernelMemoryOptimized<<< blocks, threads,shared_mem_size >>>(*nd,*irepo,*maxdat,*MAXVAR,
                                            d_x,d_y,d_z,
                                            *EPSLON,
                                            *nlag,
                                            *xlag,*xltol,
                                            *mxdlv,
                                            d_np,d_dis,d_gam,d_hm,d_tm,
                                            *dismxs,*tmax,*tmin,
                                            *ndir,*nvarg,
                                            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
                                            d_csatol,d_csdtol,d_bandwh,d_bandwd,
                                            d_atol,
                                            d_ivtype,d_ivtail,d_ivhead,
                                            d_vr,chunk_sh_mem,frac_nd);
        	cudaDeviceSynchronize();
        	Check_CUDA_Error("fitness kernel");
        	cudaEventRecord(stop, 0);
        	cudaEventSynchronize(stop);
        	cudaEventElapsedTime(&time, start, stop);
                printf ("Time for the Optimized kernel: %f ms\n", time);
        	printf ("GPU time: %f\n", time/1000);
        	printf("------------------------------\n");
    	}else{
        	frac_nd = *maxdat/THREADSX;
        	dim3 blocks( (frac_nd + threads.x - 1)/threads.x,(frac_nd + threads.y - 1)/threads.y,1 );
        	shared_mem_size = sizeof(DT)*(*mxdlv*5);
        	cudaEventCreate(&start);
        	cudaEventCreate(&stop);
        	cudaEventRecord(start, 0);
        	variogramKernel<<< blocks, threads,shared_mem_size >>>(*nd,*irepo,*maxdat,*MAXVAR,
                                            d_x,d_y,d_z,
                                            *EPSLON,
                                            *nlag,
                                            *xlag,*xltol,
                                            *mxdlv,
                                            d_np,d_dis,d_gam,d_hm,d_tm,
                                            *dismxs,*tmax,*tmin,
                                            *ndir,*nvarg,
                                            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
                                            d_csatol,d_csdtol,d_bandwh,d_bandwd,
                                            d_atol,
                                            d_ivtype,d_ivtail,d_ivhead,
                                            d_vr,frac_nd);
        	cudaDeviceSynchronize();
        	Check_CUDA_Error("fitness kernel");
        	cudaEventRecord(stop, 0);
        	cudaEventSynchronize(stop);
        	cudaEventElapsedTime(&time, start, stop);
                printf ("Time for variogram kernel: %f ms\n", time);
        	printf ("GPU time: %f\n", time/1000);
        	printf("------------------------------\n");
    	}

    	cudaMemcpy( h_np, d_np,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpy( h_dis, d_dis,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpy( h_gam, d_gam,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpy( h_hm, d_hm,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpy( h_tm, d_tm,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost);
    	//Check_CUDA_Error("cpy d -> h");

   	// printf("np, dis, gam, hm, tm\n");
	//    float sum_np = 0.0;
    	for (i = 0; i < *mxdlv; i++){
        	np[i] = (double)h_np[i];
        	dis[i] = (double)h_dis[i];
        	gam[i] = (double)h_gam[i];
        	hm[i] = (double)h_hm[i];
        	tm[i] = (double)h_tm[i];
      	//  printf("%lf\t, %lf\t, %lf\t, %lf\t, %lf\n",np[i],dis[i],gam[i],hm[i],tm[i]);
    	}


    	cudaFree(d_x);
    	cudaFree(d_y);
    	cudaFree(d_z);
    	cudaFree(d_np);
    	cudaFree(d_dis);
    	cudaFree(d_gam);
    	cudaFree(d_hm);
    	cudaFree(d_tm);
    	cudaFree(d_uvxazm);
    	cudaFree(d_uvyazm);
    	cudaFree(d_uvzdec);
    	cudaFree(d_uvhdec);
    	cudaFree(d_csatol);
    	cudaFree(d_csdtol);
    	cudaFree(d_bandwh);
    	cudaFree(d_bandwd);
    	cudaFree(d_atol);
    	cudaFree(d_vr);
    	cudaFree(d_ivtype);
    	cudaFree(d_ivtail);
    	cudaFree(d_ivhead);
    	free(h_np);
    	free(h_dis);
    	free(h_gam);
    	free(h_hm);
    	free(h_tm);
	return 0;
//end routine
}







/*

TODO:

- Add documentation in each routine
- Add diagrams (UML?) of the sequence and interactions between the CPU and GPU
- Document the memory optimization proposed.

*/








