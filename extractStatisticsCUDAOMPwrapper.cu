
/*
Add description and legal texts
*/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

#define DT double
#define MAX(x,y)  ((x) >= (y) ? (x) : (y))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#define MEM_OPTIMIZED 0
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

uint64_t get_gtod_clock_time ()
{
    struct timeval tv;

    if (gettimeofday (&tv, NULL) == 0)
        return (uint64_t) (tv.tv_sec * 1000000 + tv.tv_usec);
    else
        return 0;
}

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



__host__ void computeVariogramOMP(float dxj, float dyj, float dzj, int i, int  j,const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    float *sh_np,float *sh_dis,float *sh_tm,float *sh_hm,float *sh_gam,
                                    //float *sh_loc,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr, int sh_pos, float xlaginv){

//    int half_nd = nd/2;
//    float dx,dy,dz;

    float dx,dy,dz,dxs,dys,dzs,hs,invh;
    int id,ii,jj,il,it,iv;
    int lagbeg,lagend,ilag;
    float band,dcazm,dcdec,dxy,vrh,vrhpr,vrt,vrtpr,h;
    int omni;
    int d_ivtailtmp, d_ivheadtmp;
    int index;

    //if(i<nd && j<nd){

    dx  = d_x[i] - dxj;
    dy  = d_y[i] - dyj;
    dz  = d_z[i] - dzj;
    //dx  = d_x[i] - d_x[j];
    //dy  = d_y[i] - d_y[j];
    //dz  = d_z[i] - d_z[j];
    //dx  = 0.0;//d_x[i] - d_x[j];
    //dy  = 0.0;//d_y[i] - d_y[j];
    //dz  = 0.0;//d_z[i] - d_z[j];
    dxs = dx*dx;
    dys = dy*dy;
    dzs = dz*dz;
    hs  = dxs + dys + dzs;

    if(hs <= dismxs)
    {
        if(hs < 0.0) hs = 0.0;
        h   = sqrtf(hs);
        invh=1.0f/h;


    //
    // Determine which lag this is and skip if outside the defined distance
    // tolerance:
    //
        if(h<=EPSLON){
            lagbeg = 1;
            lagend = 1;
        }
        else{
        //if(h>EPSLON){
            lagbeg = -1;
            lagend = -1;

            //for(ilag=2;ilag<=nlag+2;ilag++){
            //    if(h>=(xlag*(float)(ilag-2)-xltol) && h<=(xlag*(float)(ilag-2)+xltol)){
            //        if(lagbeg<0) lagbeg = ilag;
            //        lagend = ilag;
            //    }
            //}

            int ilag=0;
            int liminf=ceil((h-xltol)*xlaginv)+2;
            int limsup=floor((h+xltol)*xlaginv)+2;
//            for(ilag=liminf;ilag<=limsup;ilag++){
//                if(lagbeg<0)lagbeg=ilag;
//                lagend=ilag;
//            }

            lagbeg=liminf;
            lagend=limsup;


	}
        if(lagend>=0)
        {
        //			printf("dx=%f dy=%f dz=%fh=%f lagbeg=%d lagend=%d\n",dx,dy,dz,h,lagbeg,lagend);


        //
        // Definition of the direction corresponding to the current pair. All
        // directions are considered (overlapping of direction tolerance cones
        // is allowed):
        //


            dxy = sqrtf(MAX((dxs+dys),0.0));
            float invdxy=1.0f/dxy;
            for(id=0;id<ndir;id++){
            //
            // Check for an acceptable azimuth angle:
            //
                if(dxy<EPSLON){
                    dcazm = 1.0;
                }
                else{
                //if(dxy>=EPSLON){
                    //dcazm = (dx*d_uvxazm[id]+dy*d_uvyazm[id])/dxy;
                    dcazm = (dx*d_uvxazm[id]+dy*d_uvyazm[id])*invdxy;
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
                        //if(lagbeg!=1){
                            //dcdec = (dxy*d_uvhdec[id]+dz*d_uvzdec[id])/h;
                            dcdec = (dxy*d_uvhdec[id]+dz*d_uvzdec[id])*invh;

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


///////////////////
//				for(iv=0;iv<nvarg;iv++){
////
//// For this variogram, sort out which is the tail and the head value:
////
//					it = d_ivtype[iv];
//					if(dcazm>=0.0 && dcdec>=0.0){
//						ii = d_ivtail[iv]-1;
//						vrh   = d_vr[i+ii*(maxdat)];
//						ii = d_ivhead[iv]-1;
//						vrt   = d_vr[j+ii*(maxdat)];
//						if(omni || it==2){
//							ii    = d_ivhead[iv]-1;
//							vrtpr = d_vr[i+ii*(maxdat)];
//							ii    = d_ivtail[iv]-1;
//							vrhpr = d_vr[j+ii*(maxdat)];
//						}
//					}
//					else{
//						ii = d_ivtail[iv]-1;
//						vrh   = d_vr[j+ii*(maxdat)];
//						ii = d_ivhead[iv]-1;
//						vrt   = d_vr[i+ii*(maxdat)];
//						if(omni || it==2){
//							ii    = d_ivhead[iv]-1;
//							vrtpr = d_vr[j+ii*(maxdat)];
//							ii    = d_ivtail[iv]-1;
//							vrhpr = d_vr[i+ii*(maxdat)];
//						}
//					}




///////////////////

                            for(iv=0;iv<nvarg;iv++){
                    //
                    // For this variogram, sort out which is the tail and the head value:
                    //
				d_ivtailtmp=(d_ivtail[iv]-1)*maxdat;
				d_ivheadtmp=(d_ivhead[iv]-1)*maxdat;

                                it = d_ivtype[iv];
                                if(dcazm>=0.0 && dcdec>=0.0){
                                    //ii = d_ivtail[iv]-1;
                                    vrh   = d_vr[i+d_ivtailtmp];
                                    //jj = d_ivhead[iv]-1;
                                    vrt   = d_vr[j+d_ivheadtmp];
                                    if(omni || it==2){
                                        //ii    = d_ivhead[iv]-1;
                                        vrtpr = d_vr[i+d_ivheadtmp];
                                        //jj    = d_ivtail[iv]-1;
                                        vrhpr = d_vr[j+d_ivtailtmp];
                                    }
                                }
                                else{
                                    //jj = d_ivtail[iv]-1;
                                    vrh   = d_vr[j+d_ivtailtmp];
                                    //ii = d_ivhead[iv]-1;
                                    vrt   = d_vr[i+d_ivheadtmp];
                                    if(omni || it==2){
                                        //jj    = d_ivhead[iv]-1;
                                        vrtpr = d_vr[j+d_ivheadtmp];
                                        //ii    = d_ivtail[iv]-1;
                                        vrhpr = d_vr[i+d_ivtailtmp];
                                    }
                                }
                    //
                    // Reject this pair on the basis of missing values:
                    //
                                if(vrt>=tmin && vrh>=tmin && vrt<=tmax && vrh<=tmax && it!=2 || (vrtpr>=tmin && vrhpr>=tmin && vrtpr<=tmax && vrhpr<=tmax))
                                {
                                    if(it==1 || it==5 || it>=9){
					index=(id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2) - 1; 
                                        for(il=lagbeg;il<=lagend;il++){
                                            //ii = index + il;
                                            ii = (id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2) - 1 + il;
						//sh_np[ii + mxdlv*sh_pos]+=1.0;
						//sh_dis[ii + mxdlv*sh_pos]+=(h);
						//sh_tm[ii + mxdlv*sh_pos]+=(vrt);
						//sh_hm[ii + mxdlv*sh_pos]+=(vrh);
						//sh_gam[ii + mxdlv*sh_pos]+=((vrh-vrt)*(vrh-vrt));

						sh_np[ii]+=1.0;
						sh_dis[ii]+=(h);
						sh_tm[ii]+=(vrt);
						sh_hm[ii]+=(vrh);
						sh_gam[ii]+=((vrh-vrt)*(vrh-vrt));
						//sh_loc[ii]+=1.0;
						//sh_loc[ii+1]+=(h);
						//sh_loc[ii+2]+=(vrt);
						//sh_loc[ii+3]+=(vrh);
						//sh_loc[ii+4]+=((vrh-vrt)*(vrh-vrt));

                                            if(omni){
                                                if(vrtpr>=tmin && vrhpr>=tmin && vrtpr<tmax && vrhpr<tmax){

						//sh_np[ii + mxdlv*sh_pos]+=1.0;
						//sh_dis[ii + mxdlv*sh_pos]+=(h);
						//sh_tm[ii + mxdlv*sh_pos]+=(vrtpr);
						//sh_hm[ii + mxdlv*sh_pos]+=(vrhpr);
						//sh_gam[ii + mxdlv*sh_pos]+=((vrhpr-vrtpr)*(vrhpr-vrtpr));

						sh_np[ii ]+=1.0;
						sh_dis[ii]+=(h);
						sh_tm[ii ]+=(vrtpr);
						sh_hm[ii ]+=(vrhpr);
						sh_gam[ii]+=((vrhpr-vrtpr)*(vrhpr-vrtpr));
						//sh_loc[ii ]+=1.0;
						//sh_loc[ii+1]+=(h);
						//sh_loc[ii+2]+=(vrtpr);
						//sh_loc[ii+3]+=(vrhpr);
						//sh_loc[ii+4]+=((vrhpr-vrtpr)*(vrhpr-vrtpr));

                                                }
                                            }
                                        }
                                    }

                                    // The Traditional Cross Semivariogram:
                //
                                    else if(it==2){
					index=(id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2) - 1; 
                                        for(il=lagbeg;il<=lagend;il++){
                                            //ii = index + il;
                                            ii = (id)*(nvarg)*((nlag)+2)+(iv)*((nlag)+2) - 1 + il;
						sh_np[ii ]+=1.0;
						sh_dis[ii]+=(h);
						sh_tm[ii ]+=(0.5*(vrt+vrtpr));
						sh_hm[ii ]+=(0.5*(vrh+vrhpr));
						sh_gam[ii]+=((vrhpr-vrh)*(vrt-vrtpr));
						//sh_loc[ii ]+=1.0;
						//sh_loc[ii+1]+=(h);
						//sh_loc[ii+2]+=(0.5*(vrt+vrtpr));
						//sh_loc[ii+3]+=(0.5*(vrh+vrhpr));
						//sh_loc[ii+4]+=((vrhpr-vrh)*(vrt-vrtpr));

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

    //}

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

//__host__ void computePointsValuesOMP(int idx, int  idy,const int nd, const int irepo, const int maxdat, const int MAXVAR,
//                                    float *d_x, float *d_y, float *d_z,
//                                    const float EPSLON,
//                                    const int nlag,
//                                    const float xlag, const float xltol,
//                                    const int mxdlv,
//                                    float *sh_np,float *sh_dis,float *sh_tm,float *sh_hm,float *sh_gam,
//                                    const float dismxs, const float tmax, const float tmin,
//                                    const int ndir, const int nvarg,
//                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
//                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
//                                    float *d_atol,
//                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
//                                    float *d_vr, int sh_pos,int half_nd, float xlaginv){
//
//    int i,j;
//    j = idx + half_nd;
//    i = idy;
//
//    computeVariogramOMP(i,j,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,sh_pos,xlaginv);
//
//    if (idx > idy){
//        i = idy;
//        j = idx;
//
//        computeVariogramOMP(i,j,nd,irepo,maxdat,MAXVAR,
//            d_x,d_y,d_z,
//            EPSLON,nlag,xlag,xltol,
//            mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//            dismxs,tmax,tmin,ndir,nvarg,
//            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//            d_csatol, d_csdtol, d_bandwh, d_bandwd,
//            d_atol,
//            d_ivtype, d_ivtail, d_ivhead,
//            d_vr,sh_pos,xlaginv);
//
//    } else{
//        if (idx == idy){
//            i = idy;
//            j = idy;
//
//            computeVariogramOMP(i,j,nd,irepo,maxdat,MAXVAR,
//                d_x,d_y,d_z,
//                EPSLON,nlag,xlag,xltol,
//                mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//                dismxs,tmax,tmin,ndir,nvarg,
//                d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//                d_csatol, d_csdtol, d_bandwh, d_bandwd,
//                d_atol,
//                d_ivtype, d_ivtail, d_ivhead,
//                d_vr,sh_pos,xlaginv);
//        }
//        i = idx + half_nd;
//        j = idy + half_nd;
//
//        computeVariogramOMP(i,j,nd,irepo,maxdat,MAXVAR,
//            d_x,d_y,d_z,
//            EPSLON,nlag,xlag,xltol,
//            mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//            dismxs,tmax,tmin,ndir,nvarg,
//            d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//            d_csatol, d_csdtol, d_bandwh, d_bandwd,
//            d_atol,
//            d_ivtype, d_ivtail, d_ivhead,
//            d_vr,sh_pos,xlaginv);
//    }
//
//
//}
//
//




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
                                    float *d_vr,int chunks_sh_mem,int frac_nd, int thres_hybrid){

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

    //if (idx < frac_nd && idy < frac_nd){
    if (idx>=thres_hybrid && idy>=thres_hybrid && idx < frac_nd && idy < frac_nd){
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
                                    float *d_vr,int frac_nd, int thres_hybrid){

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

    //if (idx < frac_nd && idy < frac_nd){
    if (idx>=thres_hybrid && idy>=thres_hybrid && idx < frac_nd && idy < frac_nd){
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

__global__ void variogramKernelShrinked(    const int nd, const int irepo, const int maxdat, const int MAXVAR,
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
                                    float *d_vr,int frac_nd, int thres_hybrid2){

    int tidx=threadIdx.x;
    int tidy=threadIdx.y;
    int bidx=blockIdx.x;
    int bidy=blockIdx.y;
    int bdimx=blockDim.x;
    int bdimy=blockDim.y;
    //int idx = bidx*bdimx + tidx + thres_hybrid2;
    //int idy = bidy*bdimy + tidy + thres_hybrid2;
    int idx = bidx*bdimx + tidx;
    int idy = bidy*bdimy + tidy;
    int threadId = tidx + bdimx*tidy;
    int half_nd = (nd-thres_hybrid2)/2;
    //int half_nd = thres_hybrid + (nd-thres_hybrid)/2;
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
    //if (idx>=thres_hybrid && idy>=thres_hybrid && idx < frac_nd && idy < frac_nd){
        for (i = idx; i < half_nd; i += frac_nd){
            for (j = idy; j < half_nd; j += frac_nd){
                //computePointsValues(i,j,nd,irepo,maxdat,MAXVAR,
                computePointsValues(i+thres_hybrid2,j+thres_hybrid2,nd,irepo,maxdat,MAXVAR,
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



//
//__host__ void variogramKernelOMP(    const int nd, const int irepo, const int maxdat, const int MAXVAR,
//                                    float *d_x, float *d_y, float *d_z,
//                                    const float EPSLON,
//                                    const int nlag,
//                                    const float xlag, const float xltol,
//                                    const int mxdlv,
//                                    DT *h_np, DT *h_dis, DT *h_gam, DT *h_hm, DT *h_tm,
//                                    const float dismxs, const float tmax, const float tmin,
//                                    const int ndir, const int nvarg,
//                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
//                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
//                                    float *d_atol,
//                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
//                                    float *d_vr,int frac_nd, int thres_hybrid){
//    printf("Inside host kernel\n");
//    //int tidx=threadIdx.x;
//    //int tidy=threadIdx.y;
//    //int bidx=blockIdx.x;
//    //int bidy=blockIdx.y;
//    //int bdimx=blockDim.x;
//    //int bdimy=blockDim.y;
//    //int idx = bidx*bdimx + tidx;
//    //int idy = bidy*bdimy + tidy;
//
//    float xlaginv=1.0/xlag;
//
//    int idx=0;
//    int idy=0;
//    //int threadId = tidx + bdimx*tidy;
//    int threadId=0;
//    int half_nd = nd/2;
//    int i,j,ii,jj;
//    //int num_threads = bdimx*bdimy;
//    int num_threads=1;
//#pragma omp parallel
//{
//    num_threads=omp_get_num_threads();
//}
//    printf("num_threads=%d\n",num_threads);
//    //extern __shared__ float buffer[];
//    float buffer[num_threads*1*mxdlv*5] ;
//    for(i=0;i<num_threads*mxdlv*5;i++)
//        buffer[i]=0;
//    float *sh_np = &buffer[0];
//    float *sh_dis = &buffer[mxdlv*1*num_threads];
//    float *sh_gam = &buffer[2*mxdlv*1*num_threads];
//    float *sh_hm = &buffer[3*mxdlv*1*num_threads];
//    float *sh_tm = &buffer[4*mxdlv*1*num_threads];
//
//
//    int blocksx = (frac_nd + THREADSX - 1)/THREADSX;
//    int blocksy = (frac_nd + THREADSY - 1)/THREADSY;
//
//    //int init_sh_mem = threadId;
//    //while (init_sh_mem < mxdlv){
//    //    sh_np[init_sh_mem] = 0;
//    //    sh_dis[init_sh_mem] = 0.0;
//    //    sh_gam[init_sh_mem] = 0.0;
//    //    sh_hm[init_sh_mem] = 0.0;
//    //    sh_tm[init_sh_mem] = 0.0;
//    //    init_sh_mem += num_threads;
//    //}
//
//    //__syncthreads();
//
//    //if (idx < frac_nd && idy < frac_nd){
//    //if (idx<thres_hybrid || idy<thres_hybrid){
//#pragma omp parallel shared(d_x,d_y,d_z,buffer,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,d_vr)
//{
//    threadId=omp_get_thread_num();
//
//    //for(idy=0;idy<blocksy;idy++){
//    //for(idx=0;idx<thres_hybrid;idx++){
//    #pragma omp for 
//    for (jj = 0; jj < half_nd; jj += 1){
//        for (ii = 0; ii < thres_hybrid*THREADSX/2; ii += 1){
//            //for (i = ii; i < nd; i += half_nd){
//            //    for (j = jj; j < nd; j += half_nd){
//                computePointsValuesOMP(ii,jj,nd,irepo,maxdat,MAXVAR,
//                    d_x,d_y,d_z,
//                    EPSLON,nlag,xlag,xltol,
//                    mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//                    dismxs,tmax,tmin,ndir,nvarg,
//                    d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//                    d_csatol, d_csdtol, d_bandwh, d_bandwd,
//                    d_atol,
//                    d_ivtype, d_ivtail, d_ivhead,
//                    d_vr,threadId,half_nd,xlaginv);
//        //    }
//        //}
//        //    }
//        //}
//    }
//    }
//
//    //for(idx=thres_hybrid;idx<blocksx;idx++){
//    //for(idy=0;idy<thres_hybrid;idy++){
//    #pragma omp for 
//    for (ii = thres_hybrid*THREADSX/2; ii < half_nd; ii += 1){
//        for (jj = 0; jj < thres_hybrid*THREADSY/2; jj += 1){
//            //for (i = ii; i < nd; i += half_nd){
//            //    for (j = jj; j < nd; j += half_nd){
//                computePointsValuesOMP(ii,jj,nd,irepo,maxdat,MAXVAR,
//                    d_x,d_y,d_z,
//                    EPSLON,nlag,xlag,xltol,
//                    mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//                    dismxs,tmax,tmin,ndir,nvarg,
//                    d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//                    d_csatol, d_csdtol, d_bandwh, d_bandwd,
//                    d_atol,
//                    d_ivtype, d_ivtail, d_ivhead,
//                    d_vr,threadId,half_nd,xlaginv);
//        //    }
//        //}
//        //    }
//        //}
//    }
//    }
//}
//
//
//    //__syncthreads();
//
//	for(threadId=0;threadId<num_threads;threadId++){
//	for(ii=0;ii<mxdlv;ii++){
//		h_np[ii]+=sh_np[ii+threadId*mxdlv];
//		h_dis[ii]+=sh_dis[ii+threadId*mxdlv];
//		h_tm[ii]+=sh_tm[ii+threadId*mxdlv];
//		h_hm[ii]+=sh_hm[ii+threadId*mxdlv];
//		h_gam[ii]+=sh_gam[ii+threadId*mxdlv];
//	}
//	}
//
//
//
////    if (threadId < mxdlv){
////
////        atomicAdd(&d_np[threadId],sh_np[threadId]);
////        atomicAdd(&d_dis[threadId],sh_dis[threadId]);
////        atomicAdd(&d_tm[threadId],sh_tm[threadId]);
////        atomicAdd(&d_hm[threadId],sh_hm[threadId]);
////        atomicAdd(&d_gam[threadId],sh_gam[threadId]);
////    }
//}
//
//
//__host__ void variogramKernelOMPOptimized(    const int nd, const int irepo, const int maxdat, const int MAXVAR,
//                                    float *d_x, float *d_y, float *d_z,
//                                    const float EPSLON,
//                                    const int nlag,
//                                    const float xlag, const float xltol,
//                                    const int mxdlv,
//                                    DT *h_np, DT *h_dis, DT *h_gam, DT *h_hm, DT *h_tm,
//                                    const float dismxs, const float tmax, const float tmin,
//                                    const int ndir, const int nvarg,
//                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
//                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
//                                    float *d_atol,
//                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
//                                    float *d_vr,int frac_nd, int thres_hybrid){
//    printf("Inside host kernel\n");
//    //int tidx=threadIdx.x;
//    //int tidy=threadIdx.y;
//    //int bidx=blockIdx.x;
//    //int bidy=blockIdx.y;
//    //int bdimx=blockDim.x;
//    //int bdimy=blockDim.y;
//    //int idx = bidx*bdimx + tidx;
//    //int idy = bidy*bdimy + tidy;
//
//    float xlaginv=1.0/xlag;
//
//    int idx=0;
//    int idy=0;
//    //int threadId = tidx + bdimx*tidy;
//    //int threadId=0;
//    int half_nd = nd/2;
//    int i,j,ii,jj;
//    //int num_threads = bdimx*bdimy;
//    int num_threads=1;
//#pragma omp parallel
//{
//    num_threads=omp_get_num_threads();
//}
//    printf("num_threads=%d\n",num_threads);
//    //extern __shared__ float buffer[];
//    float buffer[num_threads*1*mxdlv*5] ;
//    for(i=0;i<num_threads*mxdlv*5;i++)
//        buffer[i]=0;
//    float *sh_np = &buffer[0];
//    float *sh_dis = &buffer[mxdlv*1*num_threads];
//    float *sh_gam = &buffer[2*mxdlv*1*num_threads];
//    float *sh_hm = &buffer[3*mxdlv*1*num_threads];
//    float *sh_tm = &buffer[4*mxdlv*1*num_threads];
//
//
//    int blocksx = (frac_nd + THREADSX - 1)/THREADSX;
//    int blocksy = (frac_nd + THREADSY - 1)/THREADSY;
//
//    //int init_sh_mem = threadId;
//    //while (init_sh_mem < mxdlv){
//    //    sh_np[init_sh_mem] = 0;
//    //    sh_dis[init_sh_mem] = 0.0;
//    //    sh_gam[init_sh_mem] = 0.0;
//    //    sh_hm[init_sh_mem] = 0.0;
//    //    sh_tm[init_sh_mem] = 0.0;
//    //    init_sh_mem += num_threads;
//    //}
//
//    //__syncthreads();
//
//    //if (idx < frac_nd && idy < frac_nd){
//    //if (idx<thres_hybrid || idy<thres_hybrid){
//
////    int thresTHREADSYhalf =thres_hybrid*THREADSY/2; 
////    int thresTHREADSXhalf =thres_hybrid*THREADSX/2; 
////
////    int counter=0;
////    for (idy = 0; idy < thresTHREADSYhalf ; idy++){
////        for (idx = idy; idx < nd; idx++){
////            counter=counter+1;
////        }
////    }
////    for (idx = half_nd; idx < half_nd + thresTHREADSXhalf; idx += 1){
////        for (idy = thresTHREADSYhalf; idy < idx; idy += 1){
////            counter=counter+1;
////        }
////    }
////    for (idy = half_nd; idy < half_nd + thresTHREADSYhalf; idy += 1){
////        for (idx = half_nd + thresTHREADSXhalf; idx < nd ; idx += 1){
////            counter=counter+1;
////        }
////    }
////
////    //printf("counter=%u\n",counter);
////    int indexesIDX[counter];
////    int indexesIDY[counter];
////    //printf("after counter=%u\n",counter);
////
////    counter=0;
////    for (idy = 0; idy < thresTHREADSYhalf ; idy++){
////        for (idx = idy; idx < nd; idx++){
////            //printf("idx=%d, idy=%d, counter=%u\n",idx,idy,counter);
////            indexesIDX[counter]=idx;
////            indexesIDY[counter]=idy;
////            counter=counter+1;
////        }
////    }
////    for (idx = half_nd; idx < half_nd + thresTHREADSXhalf; idx += 1){
////        for (idy = thresTHREADSYhalf; idy < idx; idy += 1){
////            indexesIDX[counter]=idx;
////            indexesIDY[counter]=idy;
////            counter=counter+1;
////        }
////    }
////    for (idy = half_nd; idy < half_nd + thresTHREADSYhalf; idy += 1){
////        for (idx = half_nd + thresTHREADSXhalf; idx < nd ; idx += 1){
////            indexesIDX[counter]=idx;
////            indexesIDY[counter]=idy;
////            counter=counter+1;
////        }
////    }    
//
//
//
////#pragma omp parallel default(none) shared(d_x,d_y,d_z,buffer,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,d_csatol,d_csdtol,d_bandwh,d_bandwd,d_atol,d_ivtype,d_ivtail,d_ivhead,d_vr,thres_hybrid,xlaginv,half_nd,indexesIDX,indexesIDY,counter) private(idy,idx,jj,ii,thresTHREADSXhalf,thresTHREADSYhalf)
//#pragma omp parallel default(none) shared(d_x,d_y,d_z,buffer,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,d_csatol,d_csdtol,d_bandwh,d_bandwd,d_atol,d_ivtype,d_ivtail,d_ivhead,d_vr,thres_hybrid,xlaginv,half_nd) private(idy,idx,jj,ii)
//{
//    int threadId=omp_get_thread_num();
//    int countPairs=0;
//    int thresTHREADSYhalf =thres_hybrid*THREADSY/2; 
//    int thresTHREADSXhalf =thres_hybrid*THREADSX/2; 
//    printf("thresTHREADSYhalf=%d\n",thresTHREADSYhalf);
//    printf("thresTHREADSXhalf=%d\n",thresTHREADSXhalf);
//    printf("nd=%d\n",nd);
//
//    float sh_np_loc[mxdlv];
//    float sh_dis_loc[mxdlv]; 
//    float sh_gam_loc[mxdlv]; 
//    float sh_hm_loc[mxdlv]; 
//    float sh_tm_loc[mxdlv]; 
//
//    for(jj=0;jj<mxdlv;jj++){
//        sh_np_loc[jj]=0.0;
//        sh_dis_loc[jj]=0.0;
//        sh_tm_loc[jj]=0.0;
//        sh_hm_loc[jj]=0.0;
//        sh_gam_loc[jj]=0.0;
//    }
//
//
//    #pragma omp for schedule(runtime) 
//    for (idx = half_nd ; idx < nd; idx++){
//        int mm = MIN(idx-half_nd,thresTHREADSYhalf);
//        for (idy = 0; idy < mm ; idy++){
//
//    countPairs++;
//    computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
//        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,0,xlaginv);
//        //d_vr,threadId,xlaginv);
//
//
//    countPairs++;
//    computeVariogramOMP(idx-half_nd,idy,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
//        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,0,xlaginv);
//        //d_vr,threadId,xlaginv);
//
//    countPairs++;
//    computeVariogramOMP(idx,idy+half_nd,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
//        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,0,xlaginv);
//        //d_vr,threadId,xlaginv);
//
//    countPairs++;
//    computeVariogramOMP(idy+half_nd,idx-half_nd,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
//        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,0,xlaginv);
//        //d_vr,threadId,xlaginv);
//
//
//
//    }
//    }
//
//
//
//
////
////
////    #pragma omp for schedule(static,THREADSX) 
////    for (idy = 0; idy < thresTHREADSYhalf ; idy++){
////        for (idx = idy; idx < nd; idx++){
////    countPairs++;
////    computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
////        d_x,d_y,d_z,
////        EPSLON,nlag,xlag,xltol,
////        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
////        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
////        dismxs,tmax,tmin,ndir,nvarg,
////        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
////        d_csatol, d_csdtol, d_bandwh, d_bandwd,
////        d_atol,
////        d_ivtype, d_ivtail, d_ivhead,
////        d_vr,0,xlaginv);
////        //d_vr,threadId,xlaginv);
////    }
////    }
////
////
////    #pragma omp for schedule(static,THREADSX) 
////    for (idx = half_nd; idx < half_nd + thresTHREADSXhalf; idx += 1){
////        for (idy = thresTHREADSYhalf; idy < idx; idy += 1){
////    //printf("Entro loop 2\n");
////    countPairs++;
////    computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
////        d_x,d_y,d_z,
////        EPSLON,nlag,xlag,xltol,
////        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
////        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
////        dismxs,tmax,tmin,ndir,nvarg,
////        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
////        d_csatol, d_csdtol, d_bandwh, d_bandwd,
////        d_atol,
////        d_ivtype, d_ivtail, d_ivhead,
////        d_vr,0,xlaginv);
////        //d_vr,threadId,xlaginv);
////    }
////    }
////
////
////    #pragma omp for schedule(static,THREADSX) 
////    for (idy = half_nd; idy < half_nd + thresTHREADSYhalf; idy += 1){
////        //for (idx = half_nd + thresTHREADSXhalf; idx < nd ; idx += 1){
////        for (idx = idy; idx < nd ; idx += 1){
////    //printf("Entro loop 3\n");
////    countPairs++;
////    computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
////        d_x,d_y,d_z,
////        EPSLON,nlag,xlag,xltol,
////        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
////        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
////        dismxs,tmax,tmin,ndir,nvarg,
////        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
////        d_csatol, d_csdtol, d_bandwh, d_bandwd,
////        d_atol,
////        d_ivtype, d_ivtail, d_ivhead,
////        d_vr,0,xlaginv);
////        //d_vr,threadId,xlaginv);
////    }
////    }
////
//
//    //for(ii=0;ii<num_threads;ii++){
//    for(jj=0;jj<mxdlv;jj++){
//        sh_np[jj+threadId*mxdlv]	=sh_np_loc[jj];;
//        sh_dis[jj+threadId*mxdlv]	=sh_dis_loc[jj];
//        sh_tm[jj+threadId*mxdlv]	=sh_tm_loc[jj];
//        sh_hm[jj+threadId*mxdlv]	=sh_hm_loc[jj];
//        sh_gam[jj+threadId*mxdlv]	=sh_gam_loc[jj];
//    }
//    //}
//
//    printf("countPairs=%d\n",countPairs);
//}
//    //__syncthreads();
//
//
//        int threadIdX=0;
//	for(threadIdX=0;threadIdX<num_threads;threadIdX++){
//	for(ii=0;ii<mxdlv;ii++){
//		h_np[ii]+=sh_np[ii+threadIdX*mxdlv];
//		h_dis[ii]+=sh_dis[ii+threadIdX*mxdlv];
//		h_tm[ii]+=sh_tm[ii+threadIdX*mxdlv];
//		h_hm[ii]+=sh_hm[ii+threadIdX*mxdlv];
//		h_gam[ii]+=sh_gam[ii+threadIdX*mxdlv];
//	}
//	}
//
////	free(indexes);
//
////    if (threadId < mxdlv){
////
////        atomicAdd(&d_np[threadId],sh_np[threadId]);
////        atomicAdd(&d_dis[threadId],sh_dis[threadId]);
////        atomicAdd(&d_tm[threadId],sh_tm[threadId]);
////        atomicAdd(&d_hm[threadId],sh_hm[threadId]);
////        atomicAdd(&d_gam[threadId],sh_gam[threadId]);
////    }
//}
//

__host__ void variogramKernelOMPOptimizedShrinked(    const int nd, const int irepo, const int maxdat, const int MAXVAR,
                                    float *d_x, float *d_y, float *d_z,
                                    const float EPSLON,
                                    const int nlag,
                                    const float xlag, const float xltol,
                                    const int mxdlv,
                                    DT *h_np, DT *h_dis, DT *h_gam, DT *h_hm, DT *h_tm,
                                    const float dismxs, const float tmax, const float tmin,
                                    const int ndir, const int nvarg,
                                    float *d_uvxazm,  float *d_uvyazm,  float *d_uvzdec,  float *d_uvhdec,
                                    float *d_csatol, float *d_csdtol, float *d_bandwh, float *d_bandwd,
                                    float *d_atol,
                                    int *d_ivtype, int *d_ivtail, int *d_ivhead,
                                    float *d_vr,int frac_nd, int thres_hybrid){
    printf("Inside host kernel\n");

    float xlaginv=1.0/xlag;

    int idx=0;
    int idy=0;
    //int threadId = tidx + bdimx*tidy;
    //int threadId=0;
    int half_nd = nd/2;
    int i,j,ii,jj;
    //int num_threads = bdimx*bdimy;
    int num_threads=1;
#pragma omp parallel
{
    num_threads=omp_get_num_threads();
}
    printf("num_threads=%d\n",num_threads);
    //extern __shared__ float buffer[];
    float buffer[num_threads*1*mxdlv*5] ;
    for(i=0;i<num_threads*mxdlv*5;i++)
        buffer[i]=0;
    float *sh_np = &buffer[0];
    float *sh_dis = &buffer[mxdlv*1*num_threads];
    float *sh_gam = &buffer[2*mxdlv*1*num_threads];
    float *sh_hm = &buffer[3*mxdlv*1*num_threads];
    float *sh_tm = &buffer[4*mxdlv*1*num_threads];


    //int blocksx = (frac_nd + THREADSX - 1)/THREADSX;
    //int blocksy = (frac_nd + THREADSY - 1)/THREADSY;

    float d_xidy, d_yidy, d_zidy;

#pragma omp parallel default(none) shared(d_x,d_y,d_z,buffer,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,d_csatol,d_csdtol,d_bandwh,d_bandwd,d_atol,d_ivtype,d_ivtail,d_ivhead,d_vr,thres_hybrid,xlaginv,half_nd) private(idy,idx,jj,ii,d_xidy, d_yidy, d_zidy)
{
    int threadId=omp_get_thread_num();
    int countPairs=0;
    int thresTHREADSYhalf =thres_hybrid*THREADSY/2; 
    int thresTHREADSXhalf =thres_hybrid*THREADSX/2; 
//    printf("thresTHREADSYhalf=%d\n",thresTHREADSYhalf);
//    printf("thresTHREADSXhalf=%d\n",thresTHREADSXhalf);
//    printf("nd=%d\n",nd);

    float sh_np_loc[mxdlv];
    float sh_dis_loc[mxdlv]; 
    float sh_gam_loc[mxdlv]; 
    float sh_hm_loc[mxdlv]; 
    float sh_tm_loc[mxdlv]; 

    //float sh_loc[5*mxdlv]; 

    for(jj=0;jj<mxdlv;jj++){
        sh_np_loc[jj]=0.0;
        sh_dis_loc[jj]=0.0;
        sh_tm_loc[jj]=0.0;
        sh_hm_loc[jj]=0.0;
        sh_gam_loc[jj]=0.0;
        //sh_loc[jj]=0.0;
        //sh_loc[jj+1]=0.0;
        //sh_loc[jj+2]=0.0;
        //sh_loc[jj+3]=0.0;
        //sh_loc[jj+4]=0.0;

    }


    #pragma omp for schedule(runtime) nowait
    for (idy = 0; idy < thres_hybrid ; idy++){
    d_xidy=d_x[idy];
    d_yidy=d_y[idy];
    d_zidy=d_z[idy];
    for (idx = idy ; idx < nd; idx++){
        //int mm = MIN(idx,thres_hybrid);
        //for (idy = 0; idy < mm ; idy++){

    countPairs++;
    //computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
    computeVariogramOMP(d_xidy, d_yidy, d_zidy, idx,idy,nd,irepo,maxdat,MAXVAR,
        d_x,d_y,d_z,
        EPSLON,nlag,xlag,xltol,
        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
        //mxdlv,sh_loc,
        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
        dismxs,tmax,tmin,ndir,nvarg,
        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
        d_csatol, d_csdtol, d_bandwh, d_bandwd,
        d_atol,
        d_ivtype, d_ivtail, d_ivhead,
        d_vr,0,xlaginv);
        //d_vr,threadId,xlaginv);


    }
//    for (; idx < nd; idx++){
//        //int mm = MIN(idx,thres_hybrid);
//        //for (idy = 0; idy < mm ; idy++){
//
//    countPairs++;
//    //computeVariogramOMP(idx,idy,nd,irepo,maxdat,MAXVAR,
//    computeVariogramOMP(d_xidy, d_yidy, d_zidy, idx,idy,nd,irepo,maxdat,MAXVAR,
//        d_x,d_y,d_z,
//        EPSLON,nlag,xlag,xltol,
//        mxdlv,sh_np_loc,sh_dis_loc,sh_tm_loc,sh_hm_loc,sh_gam_loc,
//        //mxdlv,sh_np,sh_dis,sh_tm,sh_hm,sh_gam,
//        dismxs,tmax,tmin,ndir,nvarg,
//        d_uvxazm,d_uvyazm,d_uvzdec,d_uvhdec,
//        d_csatol, d_csdtol, d_bandwh, d_bandwd,
//        d_atol,
//        d_ivtype, d_ivtail, d_ivhead,
//        d_vr,0,xlaginv);
//        //d_vr,threadId,xlaginv);
//
//    }
    }


    //for(ii=0;ii<num_threads;ii++){
    for(jj=0;jj<mxdlv;jj++){
        sh_np[jj+threadId*mxdlv]	=sh_np_loc[jj];;
        sh_dis[jj+threadId*mxdlv]	=sh_dis_loc[jj];
        sh_tm[jj+threadId*mxdlv]	=sh_tm_loc[jj];
        sh_hm[jj+threadId*mxdlv]	=sh_hm_loc[jj];
        sh_gam[jj+threadId*mxdlv]	=sh_gam_loc[jj];
    }
    //}

    //printf("countPairs=%d\n",countPairs);
}


        int threadIdX=0;
	for(threadIdX=0;threadIdX<num_threads;threadIdX++){
	for(ii=0;ii<mxdlv;ii++){
		h_np[ii]+=sh_np[ii+threadIdX*mxdlv];
		h_dis[ii]+=sh_dis[ii+threadIdX*mxdlv];
		h_tm[ii]+=sh_tm[ii+threadIdX*mxdlv];
		h_hm[ii]+=sh_hm[ii+threadIdX*mxdlv];
		h_gam[ii]+=sh_gam[ii+threadIdX*mxdlv];
	}
	}

}


extern "C" int extractstatisticscudaompwrapper_(
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
	DT *hh_np,*hh_dis,*hh_gam,*hh_hm,*hh_tm;
	float *d_uvxazm,*d_uvyazm,*d_uvzdec,*d_uvhdec,*d_csatol,*d_csdtol,*d_bandwh,*d_bandwd,*d_atol,*d_vr;
	int *d_ivtype,*d_ivtail,*d_ivhead;



    	cudaSetDevice(0);
	cudaStream_t streamid;
	cudaStreamCreate(&streamid);

	// CUDA kernel will process the first half of the data.
	float thres_factor = (float)THRES;
	//int thres_factor = 2;
	//int thres_hybrid = (int)(thres_factor*(float)(*maxdat/THREADSX));
	int thres_hybrid2 = (int)ceil(((*maxdat))*thres_factor);
	int thres_hybrid = (int)ceil(((*maxdat)/THREADSX)*(thres_factor));
	printf("thres_factor=%f\n",thres_factor);	



    	dim3 threads(THREADSX,THREADSY,1);
	int frac_nd;

    	int chunk_sh_mem = 4;
	h_np = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_dis = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_gam = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_hm = (DT*)malloc(sizeof(DT)* *mxdlv);
	h_tm = (DT*)malloc(sizeof(DT)* *mxdlv);
	hh_np = (DT*)malloc(sizeof(DT)* *mxdlv);
	hh_dis = (DT*)malloc(sizeof(DT)* *mxdlv);
	hh_gam = (DT*)malloc(sizeof(DT)* *mxdlv);
	hh_hm = (DT*)malloc(sizeof(DT)* *mxdlv);
	hh_tm = (DT*)malloc(sizeof(DT)* *mxdlv);


	int shared_mem_size;
    	int i;
    	for (i = 0; i < *mxdlv; i++){
        	h_np[i] = 0.0;
        	h_dis[i] = 0.0;
        	h_gam[i] = 0.0;
        	h_hm[i] = 0.0;
        	h_tm[i] = 0.0;
        	hh_np[i] = 0.0;
        	hh_dis[i] = 0.0;
        	hh_gam[i] = 0.0;
        	hh_hm[i] = 0.0;
        	hh_tm[i] = 0.0;
    	}
   	
	cudaEvent_t start, stop;
   	float time;

       	cudaEventCreate(&start);
       	cudaEventCreate(&stop);
       	cudaEventRecord(start, streamid);

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
   	cudaMemcpyAsync( d_x, x,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpyAsync( d_y, y,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpyAsync( d_z, z,sizeof(float) * (*maxdat), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy coords h -> d");
   	cudaMemcpyAsync( d_np, h_np,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpyAsync( d_dis, h_dis,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpyAsync( d_gam, h_gam,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpyAsync( d_hm, h_hm,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpyAsync( d_tm, h_tm,sizeof(DT) * (*mxdlv), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy np, dis, gam, hm, tm h -> d");
   	cudaMemcpyAsync( d_uvxazm, uvxazm,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_uvyazm, uvyazm,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_uvzdec, uvzdec,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_uvhdec, uvhdec,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_csatol, csatol,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_csdtol, csdtol,sizeof(float) * (100), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_bandwh, bandwh,sizeof(float) * (*ndir), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_bandwd, bandwd,sizeof(float) * (*ndir), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_atol, atol,sizeof(float) * (*ndir), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_vr, vr,sizeof(float) * (*maxdat* *MAXVAR), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy small data h -> d");
   	cudaMemcpyAsync( d_ivtype, ivtype,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy iv var h -> d");
   	cudaMemcpyAsync( d_ivtail, ivtail,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy iv var h -> d");
   	cudaMemcpyAsync( d_ivhead, ivhead,sizeof(float) * (*nvarg), cudaMemcpyHostToDevice , streamid);
   	//Check_CUDA_Error("cpy iv var h -> d");
//       	cudaStreamSynchronize(streamid);
       	cudaEventRecord(stop, streamid);
       	cudaEventSynchronize(stop);
       	cudaEventElapsedTime(&time, start, stop);
	printf ("Time for the GPU copy cpu-gpu: %f s\n", time/1000);
 



	//int num_devices;
	//cudaGetDeviceCount(&num_devices);
	//printf("num_devices=%d\n",num_devices);



	frac_nd = *maxdat/THREADSX;
       	dim3 blocks( (frac_nd + threads.x - 1)/threads.x,(frac_nd + threads.y - 1)/threads.y,1 );
       	cudaEventCreate(&start);
       	cudaEventCreate(&stop);
       	cudaEventRecord(start, streamid);
    	if (MEM_OPTIMIZED){
       		shared_mem_size = sizeof(DT)*(*mxdlv*5*chunk_sh_mem);
		//printf("Starting asynchronous CUDA kernel (mem-opt)...\n");
        	variogramKernelMemoryOptimized<<< blocks, threads,shared_mem_size,streamid >>>(*nd,*irepo,*maxdat,*MAXVAR,
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
                                            d_vr,chunk_sh_mem,frac_nd,thres_hybrid);
        	//cudaDeviceSynchronize();
    	}else{
        	shared_mem_size = sizeof(DT)*(*mxdlv*5);
		//printf("Starting asynchronous CUDA kernel...\n");
        	//variogramKernel<<< blocks, threads,shared_mem_size, streamid >>>(*nd,*irepo,*maxdat,*MAXVAR,
		frac_nd = (*maxdat-thres_hybrid2)/THREADSX;
       		dim3 blocks2( (frac_nd + threads.x - 1)/threads.x,(frac_nd + threads.y - 1)/threads.y,1 );
		if(thres_factor<1.0){
        	variogramKernelShrinked<<< blocks2, threads,shared_mem_size, streamid >>>(*nd,*irepo,*maxdat,*MAXVAR,
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
                                            //d_vr,frac_nd,thres_hybrid);
                                            d_vr,frac_nd,thres_hybrid2);
		}
     

		//cudaDeviceSynchronize();
    	}
//       	Check_CUDA_Error("fitness kernel");
//	cudaEventRecord(stop, streamid);
//       	cudaEventSynchronize(stop);
//       	cudaEventElapsedTime(&time, start, stop);
//	printf ("Time for the GPU kernel: %f s\n", time/1000);
 

      	//printf ("GPU time: %f\n", time/1000);
       	//printf("------------------------------\n");


        uint64_t prev_time_value, time_value;
        uint64_t time_diff;
        prev_time_value = get_gtod_clock_time ();



	frac_nd = (*maxdat)/THREADSX;
	variogramKernelOMPOptimizedShrinked(*nd,*irepo,*maxdat,*MAXVAR,
                                    x,y,z,
                                    *EPSLON,
                                    *nlag,
                                    *xlag,*xltol,
                                    *mxdlv,
                                    hh_np,hh_dis,hh_gam,hh_hm,hh_tm,
                                    *dismxs,*tmax,*tmin,
                                    *ndir,*nvarg,
                                    uvxazm,uvyazm,uvzdec,uvhdec,
                                    csatol,csdtol,bandwh,bandwd,
                                    atol,
                                    ivtype,ivtail,ivhead,
                                    vr,frac_nd,thres_hybrid2);
                                    //vr,frac_nd,thres_hybrid*(THREADSX));
        time_value = get_gtod_clock_time ();
        time_diff = time_value - prev_time_value;
        printf ("Time for the CPU kernel: %f s\n", (float)(time_diff)/1000000.0); 



	//cudaStreamSynchronize(0);
       	//cudaStreamSynchronize(streamid);
       	       	
	cudaEventCreate(&start);
       	cudaEventCreate(&stop);
       	cudaEventRecord(start, streamid);

       	cudaStreamSynchronize(streamid);
    	cudaMemcpyAsync( h_np, d_np,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost, streamid);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpyAsync( h_dis, d_dis,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost, streamid);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpyAsync( h_gam, d_gam,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost, streamid);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpyAsync( h_hm, d_hm,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost, streamid);
    	//Check_CUDA_Error("cpy d -> h");
    	cudaMemcpyAsync( h_tm, d_tm,sizeof(DT) * (*mxdlv),cudaMemcpyDeviceToHost, streamid);
    	//Check_CUDA_Error("cpy d -> h");
       	cudaEventRecord(stop, streamid);
       	cudaEventSynchronize(stop);
       	cudaEventElapsedTime(&time, start, stop);
	printf ("Time for the GPU copy gpu->cpu: %f s\n", time/1000);
 
   	// printf("np, dis, gam, hm, tm\n");
	//    float sum_np = 0.0;
 
	cudaStreamDestroy(streamid);


   	for (i = 0; i < *mxdlv; i++){
        	np[i] = (double)h_np[i];
        	dis[i] = (double)h_dis[i];
        	gam[i] = (double)h_gam[i];
        	hm[i] = (double)h_hm[i];
        	tm[i] = (double)h_tm[i];
      	//  printf("%lf\t, %lf\t, %lf\t, %lf\t, %lf\n",np[i],dis[i],gam[i],hm[i],tm[i]);
    	}
   	for (i = 0; i < *mxdlv; i++){
        	np[i] += (double)hh_np[i];
        	dis[i] += (double)hh_dis[i];
        	gam[i] += (double)hh_gam[i];
        	hm[i] += (double)hh_hm[i];
        	tm[i] += (double)hh_tm[i];
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
    	free(hh_np);
    	free(hh_dis);
    	free(hh_gam);
    	free(hh_hm);
    	free(hh_tm);
	return 0;
//end routine
}







/*

TODO:

- Add documentation in each routine
- Add diagrams (UML?) of the sequence and interactions between the CPU and GPU
- Document the memory optimization proposed.

*/








