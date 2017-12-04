      subroutine extractStatisticsFortran(
     + nd,irepo,maxdat,MAXVAR,x,y,z,EPSLON,nlag,xlag,xltol,
     + mxdlv,np,dis,gam,hm,tm,hv,tv,numThreads,reducedVariables,
     + dismxs,tmax,tmin,ndir,nvarg,
     + uvxazm,uvyazm,uvzdec,uvhdec,csatol,csdtol,bandwh,bandwd,
     + atol,ivtype,ivtail,ivhead,vr)

#ifdef _OPENMP
      use omp_lib
#endif 
      implicit none

      integer nd,irepo,maxdat,MAXVAR
      real x(maxdat),y(maxdat),z(maxdat)
      real EPSLON
      integer nlag
      real xlag,xltol
      integer mxdlv
      real*8 np(mxdlv),dis(mxdlv),gam(mxdlv),hm(mxdlv),
     + tm(mxdlv),hv(mxdlv),tv(mxdlv)
      integer numThreads
      real*8 reducedVariables(7,mxdlv,numThreads)
      real dismxs,tmax,tmin
      integer ndir,nvarg
      real uvxazm(100),uvyazm(100),uvzdec(100),uvhdec(100)
      real csatol(100),csdtol(100),bandwh(ndir),bandwd(ndir)
      real atol(ndir)
      integer ivtype(nvarg),ivtail(nvarg),ivhead(nvarg)
      real vr(maxdat,MAXVAR)

      integer threadId,i,j,id,ii,il,it,iv,jj
      real dx,dy,dz,dxs,dys,dzs,hs,h
      integer lagbeg,lagend,ilag
      real band,dcazm,dcdec,dxy,gamma,vrh,vrhpr,vrt,vrtpr
      logical omni
      
c------------------init extractStatistics----------------------

c$omp parallel default(firstprivate)
c$omp& shared(x,y,z,reducedVariables,vr)
#ifdef _OPENMP
      threadId = int(OMP_get_thread_num())+1 
#else
      threadId = 1
#endif

c$omp do schedule(runtime) 
      do i=1,nd
        if((int(i/irepo)*irepo).eq.i) write(*,103) i,nd
 103  format('   currently on ',i9,' of ',i9)
      do 4 j=i,nd
c
c Definition of the lag corresponding to the current pair:
c
            dx  = x(j) - x(i)
            dy  = y(j) - y(i)
            dz  = z(j) - z(i)
            dxs = dx*dx
            dys = dy*dy
            dzs = dz*dz
            hs  = dxs + dys + dzs
            if(hs.gt.dismxs) go to 4
            if(hs.lt.0.0) hs = 0.0
            h   = sqrt(hs)
c
c Determine which lag this is and skip if outside the defined distance
c tolerance:
c
            if(h.le.EPSLON) then
                  lagbeg = 1
                  lagend = 1
            else
                  lagbeg = -1
                  lagend = -1
                  do ilag=2,nlag+2
                        if(h.ge.(xlag*real(ilag-2)-xltol).and.
     +                     h.le.(xlag*real(ilag-2)+xltol)) then
                              if(lagbeg.lt.0) lagbeg = ilag 
                              lagend = ilag 
                        end if
                  end do
                  if(lagend.lt.0) go to 4
            endif

c
c Definition of the direction corresponding to the current pair. All
c directions are considered (overlapping of direction tolerance cones
c is allowed):
c
            do 5 id=1,ndir
c
c Check for an acceptable azimuth angle:
c
                  dxy = sqrt(max((dxs+dys),0.0))
                  if(dxy.lt.EPSLON) then
                        dcazm = 1.0
                  else
                        dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
                  endif
                  if(abs(dcazm).lt.csatol(id)) go to 5

c
c Check the horizontal bandwidth criteria (maximum deviation 
c perpendicular to the specified direction azimuth):
c
                  band = uvxazm(id)*dy - uvyazm(id)*dx
                  if(abs(band).gt.bandwh(id)) go to 5
c                  print *,'dxy=',dxy,' dcazm=',dcazm,' band=',band
c
c Check for an acceptable dip angle:
c
                  if(dcazm.lt.0.0) dxy = -dxy
                  if(lagbeg.eq.1) then
                        dcdec = 0.0
                  else
                        dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
                        if(abs(dcdec).lt.csdtol(id)) go to 5
                  endif
c
c Check the vertical bandwidth criteria (maximum deviation perpendicular
c to the specified dip direction):
c
                  band = uvhdec(id)*dz - uvzdec(id)*dxy
                  if(abs(band).gt.bandwd(id)) go to 5

c
c Check whether or not an omni-directional variogram is being computed:
c
                  omni = .false.
                  if(atol(id).ge.90.0) omni = .true.
c
c This direction is acceptable - go ahead and compute all variograms:
c
                  do 6 iv=1,nvarg
c
c For this variogram, sort out which is the tail and the head value:
c
                      it = ivtype(iv)
                      if(dcazm.ge.0.0.and.dcdec.ge.0.0) then
                            ii = ivtail(iv)
                            vrh   = vr(i,ii)
                            ii = ivhead(iv)
                            vrt   = vr(j,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(i,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(j,ii)
                            endif
                      else
                            ii = ivtail(iv)
                            vrh   = vr(j,ii)
                            ii = ivhead(iv)
                            vrt   = vr(i,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(j,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(i,ii)
                            endif
                      endif
c
c Reject this pair on the basis of missing values:
c
                      if(vrt.lt.tmin.or.vrh.lt.tmin.or.
     +                   vrt.gt.tmax.or.vrh.gt.tmax) go to 6
                      if(it.eq.2.and.(vrtpr.lt.tmin.or.vrhpr.lt.tmin.or.
     +                                vrtpr.gt.tmax.or.vrhpr.gt.tmax))
     +                                               go to 6
c
c             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
c
c
c The Semivariogram:
c
      if(it.eq.1.or.it.eq.5.or.it.ge.9) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble((vrhpr-vrtpr)*
     +                                              (vrhpr-vrtpr))
                     endif
               endif
            end do
c
c The Traditional Cross Semivariogram:
c
      else if(it.eq.2) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
               hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
               gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
            end do
c
c The Covariance:
c
      else if(abs(it).eq.3) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
c
c The Correlogram:
c
      else if(abs(it).eq.4) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               hv(ii)  = hv(ii)  + dble(vrh*vrh)
               tv(ii)  = tv(ii)  + dble(vrt*vrt)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
                           tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
c
c The Pairwise Relative:
c
      else if(it.eq.6) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(abs(vrt+vrh).gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(abs(vrtpr+vrhpr).gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            enddo
c
c Variogram of Logarithms:
c
      else if(it.eq.7) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = alog(vrt)-alog(vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = alog(vrt)-alog(vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            end do
c
c Madogram:
c
      else if(it.eq.8) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
                     endif
               endif
            end do
      endif
c
c Finish loops over variograms, directions, and the double data loops:
c
 6              continue
 5          continue
 4    continue

      end do
c$omp end do

#ifdef _OPENMP
      reducedVariables(1,:,threadId)=dis(:)
      reducedVariables(2,:,threadId)=gam(:)
      reducedVariables(3,:,threadId)=np(:)
      reducedVariables(4,:,threadId)=hm(:)
      reducedVariables(5,:,threadId)=tm(:)
      reducedVariables(6,:,threadId)=hv(:)
      reducedVariables(7,:,threadId)=tv(:)
#endif

c$omp end parallel

#ifdef _OPENMP
      dis(:)=0.0 
      gam(:)=0.0 
      np(:)=0.0  
      hm(:)=0.0  
      tm(:)=0.0  
      hv(:)=0.0  
      tv(:)=0.0       
      do ii=1,numThreads 
         do jj=1,mxdlv
            dis(jj) = dis(jj) + reducedVariables(1,jj,ii)
            gam(jj) = gam(jj) + reducedVariables(2,jj,ii)
            np(jj)  = np(jj)  + reducedVariables(3,jj,ii)
            hm(jj)  = hm(jj)  + reducedVariables(4,jj,ii)
            tm(jj)  = tm(jj)  + reducedVariables(5,jj,ii)
            hv(jj)  = hv(jj)  + reducedVariables(6,jj,ii)
            tv(jj)  = tv(jj)  + reducedVariables(7,jj,ii)
         end do
      end do
#endif

c------------------end extractStatistics----------------------
      return
      end
