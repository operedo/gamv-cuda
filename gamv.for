C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 2003, Statios Software and Services Incorporated.  All %
C rights reserved.                                                     %
C                                                                      %
C This program has been modified from the one distributed in 1996 (see %
C below).  This version is also distributed in the hope that it will   %
C be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
C code may be redistributed without restriction; however, this code is %
C for one developer only. Each developer or user of this source code   %
C must purchase a separate copy from Statios.                          %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Module to declare dynamic arrays in multiple subroutines:
c
c      module geostat
c
c      real,allocatable         :: x(:),y(:),z(:),vr(:,:),azm(:),atol(:),
c     +                            bandwh(:),dip(:),dtol(:),bandwd(:)
c      real*8,allocatable       :: sills(:),dis(:),gam(:),hm(:),
c     +                            tm(:),hv(:),tv(:),np(:)
c      integer,allocatable      :: ivtail(:),ivhead(:),ivtype(:)
c      character*12,allocatable :: names(:)
c
c      real      EPSLON,VERSION
c      real      xlag,xltol,tmin,tmax
c      integer   nd,nlag,ndir,nvarg,isill,test
c      character outfl*512
c
c      end module
c
c
c
      program main
c-----------------------------------------------------------------------
c
c               Variogram of Irregularly Spaced 3-D Data
c               ****************************************
c
c This is a template driver program for GSLIB's "gamv" subroutine. The
c input data must be entered with coordinates in a GEOEAS format file.
c The User's Guide can be referred to for more details.
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example gamv.par)
c
c
c
c The output file will contain each directional variogram ordered by
c direction and then variogram (the directions cycle fastest then the
c variogram number).  For each variogram there will be a one line
c description and then "nlag" lines with the following:
c
c        a) lag number (increasing from 1 to nlag)
c        b) separation distance
c        c) the "variogram" value
c        d) the number of pairs for the lag
c        e) the mean of the data contributing to the tail
c        f) the mean of the data contributing to the head
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c-----------------------------------------------------------------------
#ifdef _OPENMP
      use omp_lib
#endif
#ifdef TRACE
      use extrae_module
#endif


c geostat module
      real,allocatable         :: x(:),y(:),z(:),vr(:,:),azm(:),atol(:),
     +                            bandwh(:),dip(:),dtol(:),bandwd(:)
      real*8,allocatable       :: sills(:),dis(:),gam(:),hm(:),
     +                            tm(:),hv(:),tv(:),np(:)
      integer,allocatable      :: ivtail(:),ivhead(:),ivtype(:)
      character*12,allocatable :: names(:)

      real      EPSLON,VERSION
      real      xlag,xltol,tmin,tmax
      integer   nd,nlag,ndir,nvarg,isill,test
      character outfl*512


c readparam variables
      parameter(MV=500)

      real      var(MV),cut(MV)
      real*8    avg(MV),ssq(MV)
      integer   ivar(MV),num(MV),ivc(MV),indflag(MV)
      character datafl*512,str*512
      logical   testfl,testdat
      real,allocatable :: vrmin(:),vrmax(:)

      data      lin/1/,ncut/0/

c check variables
c      real      vrmin(*),vrmax(*)
      character title*80

c gamv variables
      parameter(PI=3.14159265)
      real      uvxazm(100),uvyazm(100),uvzdec(100),
     +          uvhdec(100),csatol(100),csdtol(100)
      logical   omni
      integer   threadId,numThreads
      real*8,allocatable ::    reducedVariables(:,:,:)
      integer extractValue
      integer :: extractstatisticscwrapper
      integer :: extractstatisticscudawrapper

c writeout variables
      character titlewr*132
      data      lout/1/
      integer*8 clock_start, clock_end, clock_rate
      real etime 
      real elapsed(2)
      real init, fin
      real total
      integer i, j

#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads()
c$omp end parallel
#else
      numThreads = 1
#endif

c      print *,numThreads


c      use geostat
      EPSLON  = 1.0e-20
      VERSION = 3.000
cc
cc Read the Parameter File:
cc
c      call readparm


c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' GAMV Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'gamv.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'gamv.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl
      write(*,*) ' columns for X,Y,Z = ',ixl,iyl,izl

      read(lin,*,err=98) nvar
      write(*,*) ' number of variables = ',nvar
      backspace lin

      read(lin,*,err=98) j,(ivar(i),i=1,nvar)
      write(*,*) ' columns = ',(ivar(i),i=1,nvar)

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nlag
      write(*,*) ' number of lags = ',nlag
      if(nlag.lt.1)       stop 'nlag is too small: check parameters'

      read(lin,*,err=98) xlag
      write(*,*) ' lag distance = ',xlag
      if(xlag.le.0.0) stop 'xlag is too small: check parameter file'

      read(lin,*,err=98) xltol
      write(*,*) ' lag tolerance = ',xltol

      read(lin,*,err=98) ndir
      write(*,*) ' number of directions = ',ndir

      allocate (azm(ndir),stat = test)
      allocate (atol(ndir),stat = test)
      allocate (bandwh(ndir),stat = test)
      allocate (dip(ndir),stat = test)
      allocate (dtol(ndir),stat = test)
      allocate (bandwd(ndir),stat = test)

      if(ndir.lt.1)       stop 'ndir is too small: check parameters'

      do i=1,ndir
            read(lin,*,err=98) azm(i),atol(i),bandwh(i),
     +                         dip(i),dtol(i),bandwd(i)
            write(*,*) ' azm, atol, bandwh = ',azm(i),atol(i),bandwh(i)
            write(*,*) ' dip, dtol, bandwd = ',dip(i),dtol(i),bandwd(i)
            if(bandwh(i).lt.0.0) then
                  write(*,*) ' Horizontal bandwidth is too small!'
                  stop
            endif
            if(bandwd(i).lt.0.0) then
                  write(*,*) ' Vertical bandwidth is too small!'
                  stop
            endif
      end do

      read(lin,*,err=98) isill
      write(*,*) ' flag to standardize sills = ',isill

      read(lin,*,err=98) nvarg
      write(*,*) ' number of variograms = ',nvarg
      if(nvarg.lt.1)      stop 'nvarg is too small: check parameters'

      mxdlv=ndir*(nlag+2)*nvarg

      allocate(reducedVariables(7,mxdlv,numThreads),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if

      allocate (dis(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (gam(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (hm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (tm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (hv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (tv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (np(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (ivtail(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (ivhead(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (ivtype(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (names(nvar+nvarg),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if

      ncut = 0
      do i=1,nvarg
            read(lin,*,err=98) ivtail(i),ivhead(i),ivtype(i)
            write(*,*) ' tail,head,type = ',
     +                   ivtail(i),ivhead(i),ivtype(i)
            if(ivtype(i).eq.9.or.ivtype(i).eq.10) then
                   ncut = ncut + 1
                   if(tmin.gt.0.0)stop'tmin interferes with indicators!'
                   if(tmax.le.1.0)stop'tmax interferes with indicators!'
                   backspace lin
                   read(lin,*,err=98) ii,jj,kk,cut(ncut)
                   if(ivtype(i).eq.9)  indflag(ncut) = 1
                   if(ivtype(i).eq.10) indflag(ncut) = 0
                   ivc(ncut) = ivtail(i)
                   ivtail(i) = nvar + ncut
                   ivhead(i) = nvar + ncut
                   write(names(nvar+ncut),140) ncut
 140               format('Indicator ',i2)
                   write(*,*) ' indicator threshold = ',cut(ncut)
            endif
      end do
      write(*,*)
      close(lin)
      MAXVAR = nvar + ncut
c
c Perform some quick error checking:
c
      if(xltol.le.0.0) then
            write(*,*) 'xltol is too small: resetting to xlag/2'
            xltol = 0.5*xlag
      endif
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
      open(lin,file=datafl,status='OLD')
c
      read(lin,*,err=99)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      maxdat = 0
 22   read(lin,*,end=44,err=99)(var(j),j=1,nvari)
      maxdat = maxdat + 1
      go to 22
 44   continue
      write(*,*)'maxdat = ',maxdat
c
      allocate (vr(maxdat,MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (vrmin(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (vrmax(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (sills(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (x(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (y(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (z(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      rewind(lin)
      read(lin,'(a40)',err=99) str
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,'(a40)',err=99) str
            do iv=1,nvar
                  j=ivar(iv)
                  if(i.eq.j) names(iv) = str(1:12)
            end do
            num(i) = 0
            avg(i) = 0.0
            ssq(i) = 0.0
      end do
c
c Read all the data until the end of the file:
c
      nd = 0
 2    continue
      read(lin,*,end=9,err=99) (var(j),j=1,nvari)
      testdat = .false.
      do iv=1,nvar
            j=ivar(iv)
            if(var(j).ge.tmin.and.var(j).lt.tmax) testdat = .true.
      end do
      if(.not.testdat) go to 2
      nd = nd + 1
c
c Acceptable data, make sure there are not too many data:
c
      do iv=1,nvar
            j=ivar(iv)
            vr(nd,iv) = var(j)
            if(var(j).ge.tmin.and.var(j).lt.tmax) then
                  num(iv) = num(iv) + 1
                  avg(iv) = avg(iv) + dble(var(j))
                  ssq(iv) = ssq(iv) + dble(var(j)*var(j))
            endif
      end do
      if(ixl.le.0) then
            x(nd) = 0.0
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = 0.0
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = 0.0
      else
            z(nd) = var(izl)
      endif
      go to 2
 9    continue
      close(lin)
c
c Compute the averages and variances as an error check for the user:
c
      do iv=1,nvar
            sills(iv) = -999.
            if(num(iv).gt.0) then
                  avg(iv)   = avg(iv)/dble(num(iv))
                  ssq(iv)   =(ssq(iv)/dble(num(iv)))-avg(iv)*avg(iv)
                  sills(iv) = ssq(iv)
                  write(*,*) 'Variable number ',iv
                  write(*,*) '  Number   = ',num(iv)
                  write(*,*) '  Average  = ',real(avg(iv))
                  write(*,*) '  Variance = ',real(ssq(iv))
            endif
      end do
c
c Construct Indicator Variables if necessary:
c
      do ic=1,ncut
            iv   = ivc(ic)
            jv   = nvar + ic
            ptot = 0.0
            p1   = 0.0
            do id=1,nd
                  if(vr(id,iv).le.tmin.or.vr(id,iv).gt.tmax) then
                        vr(id,jv) = tmin - EPSLON
                  else
                        if(indflag(ic).eq.1) then
                              if(vr(id,iv).lt.cut(ic)) then
                                    vr(id,jv) = 0.0
                              else
                                    vr(id,jv) = 1.0
                              endif
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        else
                              vr(id,jv) = 0.0
                              if(int(vr(id,iv)+0.5).eq.int(cut(ic)+0.5))
     +                        vr(id,jv) = 1.0
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        end if
                  end if
            end do
            p1        = p1 / max(ptot,1.0)
            sills(jv) = dble (p1*(1.0-p1))
      end do
c
c Establish minimums and maximums:
c
      do i=1,MAXVAR
            vrmin(i) =  1.0e21
            vrmax(i) = -1.0e21
      end do
      do id=1,nd
            do iv=1,nvar+ncut
                  if(vr(id,iv).ge.tmin.and.vr(id,iv).lt.tmax) then
                        if(vr(id,iv).lt.vrmin(iv)) vrmin(iv) = vr(id,iv)
                        if(vr(id,iv).gt.vrmax(iv)) vrmax(iv) = vr(id,iv)
                  end if
            end do
      end do
c
c Check on the variogams that were requested:
c
c      call check(vrmin,vrmax)

c
c Loop over all the variograms to be computed:
c
      write(*,*)
      do iv=1,nvarg
c
c Note the variogram type and the variables being used:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
      if(it.eq. 3) title(1:24) = 'Covariance             :'
      if(it.eq. 4) title(1:24) = 'Correlogram            :'
      if(it.eq. 5) title(1:24) = 'General Relative       :'
      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
 100  format('  tail=',a12,' head=',a12)
      write(*,101) iv,title(1:64)
 101  format(' Variogram ',i2,1x,a64)
c
c Check for possible errors or inconsistencies:
c
      if(it.eq.2) then
            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
 201        format('  WARNING: cross variogram with the same variable!')
      else if(it.eq.5) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,502)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,502)
 501        format('  WARNING: cross general relative variogram are',
     +             ' difficult to interpret!')
 502        format('  WARNING: there are both positive and negative',
     +             ' values - lag mean could be zero!')
      else if(it.eq.6) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,602)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,602)
 601        format('  WARNING: cross pairwise relative variogram are',
     +             ' difficult to interpret!')
 602        format('  WARNING: there are both positive and negative',
     +             ' values - pair means could be zero!')
      else if(it.eq.7) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0)
     +      write(*,702)
 701        format('  WARNING: cross logarithmic variograms may be',
     +             ' difficult to interpret!')
 702        format('  WARNING: there are zero or negative',
     +             ' values - logarithm undefined!')
      else if(it.eq.8) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,801)
 801        format('  WARNING: cross rodograms may be difficult to',
     +             ' interpret!')
      else if(it.eq.9) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
 901        format('  WARNING: cross madograms may be difficult to',
     +             ' interpret!')
      endif
c
c END Loop over all variograms:
c
      end do


      goto 1001
c
c Error in an Input File Somewhere:
c
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'

1001  print *,'Parameters ok.'

cc
cc Call gamv to compute the required variograms:
cc
c      call gamv

c     call system_clock(COUNT_RATE=clock_rate)

c
c Define the distance tolerance if it isn't already:
c
      if(xltol.le.0.0) xltol = 0.5 * xlag
c
c Define the angles and tolerance for each direction:
c
      do id=1,ndir
c
c The mathematical azimuth is measured counterclockwise from EW and
c not clockwise from NS as the conventional azimuth is:
c
            azmuth     = (90.0-azm(id))*PI/180.0
            uvxazm(id) = cos(azmuth)
            uvyazm(id) = sin(azmuth)
            if(atol(id).le.0.0) then
                  csatol(id) = cos(45.0*PI/180.0)
            else
                  csatol(id) = cos(atol(id)*PI/180.0)
            endif
c
c The declination is measured positive down from vertical (up) rather
c than negative down from horizontal:
c
            declin     = (90.0-dip(id))*PI/180.0
            uvzdec(id) = cos(declin)
            uvhdec(id) = sin(declin)
            if(dtol(id).le.0.0) then
                  csdtol(id) = cos(45.0*PI/180.0)
            else
                  csdtol(id) = cos(dtol(id)*PI/180.0)
            endif
      end do
c
c Initialize the arrays for each direction, variogram, and lag:
c
      nsiz = ndir*nvarg*(nlag+2)
      do i=1,nsiz
            np(i)  = 0.
            dis(i) = 0.0
            gam(i) = 0.0
            hm(i)  = 0.0
            tm(i)  = 0.0
            hv(i)  = 0.0
            tv(i)  = 0.0
      end do
      dismxs = ((real(nlag) + 0.5 - EPSLON) * xlag) ** 2
c
c MAIN LOOP OVER ALL PAIRS:
c
      irepo = max(1,min((nd/10),1000))


c     call system_clock(COUNT=clock_start)
      init = etime(elapsed)
#ifdef ANSIC
      extractValue = extractStatisticsCwrapper(
     + nd,irepo,maxdat,MAXVAR,x,y,z,EPSLON,nlag,xlag,xltol,
     + mxdlv,np,dis,gam,hm,tm,hv,tv,numThreads,reducedVariables,
     + dismxs,tmax,tmin,ndir,nvarg,
     + uvxazm,uvyazm,uvzdec,uvhdec,csatol,csdtol,bandwh,bandwd,
     + atol,ivtype,ivtail,ivhead,vr)
#else
#ifdef CUDA
      extractValue = extractStatisticsCUDAwrapper(
     + nd,irepo,maxdat,MAXVAR,x,y,z,EPSLON,nlag,xlag,xltol,
     + mxdlv,np,dis,gam,hm,tm,hv,tv,numThreads,reducedVariables,
     + dismxs,tmax,tmin,ndir,nvarg,
     + uvxazm,uvyazm,uvzdec,uvhdec,csatol,csdtol,bandwh,bandwd,
     + atol,ivtype,ivtail,ivhead,vr)
#else
#ifdef FORTRAN
      extractValue = extractStatisticsFortran(
     + nd,irepo,maxdat,MAXVAR,x,y,z,EPSLON,nlag,xlag,xltol,
     + mxdlv,np,dis,gam,hm,tm,hv,tv,numThreads,reducedVariables,
     + dismxs,tmax,tmin,ndir,nvarg,
     + uvxazm,uvyazm,uvzdec,uvhdec,csatol,csdtol,bandwh,bandwd,
     + atol,ivtype,ivtail,ivhead,vr)
#endif
#endif
#endif
      fin = etime(elapsed)

      total = fin - init
      print *, 'extract statistics total ', total
      print *, 'total ', fin
c     call system_clock(COUNT=clock_end)
c     print *, 'CLOCK=',real((real(clock_end,kind=8)- 
c     + real(clock_start,kind=8))/real(clock_rate,kind=8),kind=8)


cc------------------init extractStatistics----------------------
c
cc$omp parallel default(firstprivate)
cc$omp& shared(x,y,z,reducedVariables)
c#ifdef _OPENMP
c      threadId = int(OMP_get_thread_num())+1
c#else
c      threadId = 1
c#endif
c
cc$omp do schedule(runtime)
c      do i=1,nd
c        if((int(i/irepo)*irepo).eq.i) write(*,103) i,nd
c 103  format('   currently on ',i9,' of ',i9)
c      do 4 j=i,nd
cc
cc Definition of the lag corresponding to the current pair:
cc
c            dx  = x(j) - x(i)
c            dy  = y(j) - y(i)
c            dz  = z(j) - z(i)
c            dxs = dx*dx
c            dys = dy*dy
c            dzs = dz*dz
c            hs  = dxs + dys + dzs
c            if(hs.gt.dismxs) go to 4
c            if(hs.lt.0.0) hs = 0.0
c            h   = sqrt(hs)
cc
cc Determine which lag this is and skip if outside the defined distance
cc tolerance:
cc
c            if(h.le.EPSLON) then
c                  lagbeg = 1
c                  lagend = 1
c            else
c                  lagbeg = -1
c                  lagend = -1
c                  do ilag=2,nlag+2
c                        if(h.ge.(xlag*real(ilag-2)-xltol).and.
c     +                     h.le.(xlag*real(ilag-2)+xltol)) then
c                              if(lagbeg.lt.0) lagbeg = ilag
c                              lagend = ilag
c                        end if
c                  end do
c                  if(lagend.lt.0) go to 4
c            endif
cc
cc Definition of the direction corresponding to the current pair. All
cc directions are considered (overlapping of direction tolerance cones
cc is allowed):
cc
c            do 5 id=1,ndir
cc
cc Check for an acceptable azimuth angle:
cc
c                  dxy = sqrt(max((dxs+dys),0.0))
c                  if(dxy.lt.EPSLON) then
c                        dcazm = 1.0
c                  else
c                        dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
c                  endif
c                  if(abs(dcazm).lt.csatol(id)) go to 5
cc
cc Check the horizontal bandwidth criteria (maximum deviation
cc perpendicular to the specified direction azimuth):
cc
c                  band = uvxazm(id)*dy - uvyazm(id)*dx
c                  if(abs(band).gt.bandwh(id)) go to 5
cc
cc Check for an acceptable dip angle:
cc
c                  if(dcazm.lt.0.0) dxy = -dxy
c                  if(lagbeg.eq.1) then
c                        dcdec = 0.0
c                  else
c                        dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
c                        if(abs(dcdec).lt.csdtol(id)) go to 5
c                  endif
cc
cc Check the vertical bandwidth criteria (maximum deviation perpendicular
cc to the specified dip direction):
cc
c                  band = uvhdec(id)*dz - uvzdec(id)*dxy
c                  if(abs(band).gt.bandwd(id)) go to 5
cc
cc Check whether or not an omni-directional variogram is being computed:
cc
c                  omni = .false.
c                  if(atol(id).ge.90.0) omni = .true.
cc
cc This direction is acceptable - go ahead and compute all variograms:
cc
c                  do 6 iv=1,nvarg
cc
cc For this variogram, sort out which is the tail and the head value:
cc
c                      it = ivtype(iv)
c                      if(dcazm.ge.0.0.and.dcdec.ge.0.0) then
c                            ii = ivtail(iv)
c                            vrh   = vr(i,ii)
c                            ii = ivhead(iv)
c                            vrt   = vr(j,ii)
c                            if(omni.or.it.eq.2) then
c                                  ii    = ivhead(iv)
c                                  vrtpr = vr(i,ii)
c                                  ii    = ivtail(iv)
c                                  vrhpr = vr(j,ii)
c                            endif
c                      else
c                            ii = ivtail(iv)
c                            vrh   = vr(j,ii)
c                            ii = ivhead(iv)
c                            vrt   = vr(i,ii)
c                            if(omni.or.it.eq.2) then
c                                  ii    = ivhead(iv)
c                                  vrtpr = vr(j,ii)
c                                  ii    = ivtail(iv)
c                                  vrhpr = vr(i,ii)
c                            endif
c                      endif
c
cc                      print *,'vrh=',vrh,' vrt=',vrt
c
cc
cc Reject this pair on the basis of missing values:
cc
c                      if(vrt.lt.tmin.or.vrh.lt.tmin.or.
c     +                   vrt.gt.tmax.or.vrh.gt.tmax) go to 6
c                      if(it.eq.2.and.(vrtpr.lt.tmin.or.vrhpr.lt.tmin.or.
c     +                                vrtpr.gt.tmax.or.vrhpr.gt.tmax))
c     +                                               go to 6
cc
cc             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
cc
cc
cc The Semivariogram:
cc
c      if(it.eq.1.or.it.eq.5.or.it.ge.9) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble((vrhpr-vrtpr)*
c     +                                              (vrhpr-vrtpr))
c                     endif
c               endif
c            end do
cc
cc The Traditional Cross Semivariogram:
cc
c      else if(it.eq.2) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
c               hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
c               gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
c            end do
cc
cc The Covariance:
cc
c      else if(abs(it).eq.3) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble(vrh*vrt)
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
c                     endif
c               endif
c            end do
cc
cc The Correlogram:
cc
c      else if(abs(it).eq.4) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               hv(ii)  = hv(ii)  + dble(vrh*vrh)
c               tv(ii)  = tv(ii)  + dble(vrt*vrt)
c               gam(ii) = gam(ii) + dble(vrh*vrt)
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
c                           tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
c                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
c                     endif
c               endif
c            end do
cc
cc The Pairwise Relative:
cc
c      else if(it.eq.6) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               if(abs(vrt+vrh).gt.EPSLON) then
c                     np(ii)  = np(ii)  + 1.
c                     dis(ii) = dis(ii) + dble(h)
c                     tm(ii)  = tm(ii)  + dble(vrt)
c                     hm(ii)  = hm(ii)  + dble(vrh)
c                     gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
c                     gam(ii) = gam(ii) + dble(gamma*gamma)
c               endif
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                     if(abs(vrtpr+vrhpr).gt.EPSLON) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
c                           gam(ii) = gam(ii) + dble(gamma*gamma)
c                     endif
c                     endif
c               endif
c            enddo
cc
cc Variogram of Logarithms:
cc
c      else if(it.eq.7) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
c                     np(ii)  = np(ii)  + 1.
c                     dis(ii) = dis(ii) + dble(h)
c                     tm(ii)  = tm(ii)  + dble(vrt)
c                     hm(ii)  = hm(ii)  + dble(vrh)
c                     gamma   = alog(vrt)-alog(vrh)
c                     gam(ii) = gam(ii) + dble(gamma*gamma)
c               endif
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gamma   = alog(vrt)-alog(vrh)
c                           gam(ii) = gam(ii) + dble(gamma*gamma)
c                     endif
c                     endif
c               endif
c            end do
cc
cc Madogram:
cc
c      else if(it.eq.8) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
c                     endif
c               endif
c            end do
c      endif
cc
cc Finish loops over variograms, directions, and the double data loops:
cc
c 6              continue
c 5          continue
c 4    continue
cc 3    continue
cc      print *,'i=',i,'gam=',gam(:)
c
cc      print *,'threadId=',threadId,'/',numThreads,'i=',i,'np=',np(1:4)
c
c
c
c      end do
cc$omp end do
c
c#ifdef _OPENMP
c      reducedVariables(1,:,threadId)=dis(:)
c      reducedVariables(2,:,threadId)=gam(:)
c      reducedVariables(3,:,threadId)=np(:)
c      reducedVariables(4,:,threadId)=hm(:)
c      reducedVariables(5,:,threadId)=tm(:)
c      reducedVariables(6,:,threadId)=hv(:)
c      reducedVariables(7,:,threadId)=tv(:)
c#endif
c
cc$omp end parallel
c
c#ifdef _OPENMP
c      dis(:)=0.0
c      gam(:)=0.0
c      np(:)=0.0
c      hm(:)=0.0
c      tm(:)=0.0
c      hv(:)=0.0
c      tv(:)=0.0
c      do ii=1,numThreads
c         do jj=1,mxdlv
c            dis(jj) = dis(jj) + reducedVariables(1,jj,ii)
c            gam(jj) = gam(jj) + reducedVariables(2,jj,ii)
c            np(jj)  = np(jj)  + reducedVariables(3,jj,ii)
c            hm(jj)  = hm(jj)  + reducedVariables(4,jj,ii)
c            tm(jj)  = tm(jj)  + reducedVariables(5,jj,ii)
c            hv(jj)  = hv(jj)  + reducedVariables(6,jj,ii)
c            tv(jj)  = tv(jj)  + reducedVariables(7,jj,ii)
c         end do
c      end do
cc      print *,'np(1)=',reducedVariables(3,:,1)
cc      print *,'np(2)=',reducedVariables(3,:,2)
c#endif
c
cc------------------end extractStatistics----------------------

c
c Get average values for gam, hm, tm, hv, and tv, then compute
c the correct "variogram" measure:
c
      do 7 id=1,ndir
      do 7 iv=1,nvarg
      do 7 il=1,nlag+2
            i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
            if(np(i).le.0.) go to 7
            rnum   = np(i)
            dis(i) = dis(i) / dble(rnum)
            gam(i) = gam(i) / dble(rnum)
            hm(i)  = hm(i)  / dble(rnum)
            tm(i)  = tm(i)  / dble(rnum)
            hv(i)  = hv(i)  / dble(rnum)
            tv(i)  = tv(i)  / dble(rnum)
            it = ivtype(iv)
c
c Attempt to standardize:
c
            if(isill.eq.1) then
                  if(ivtail(iv).eq.ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it.eq.1.or.it.ge.9).and.sills(iii).gt.0.0)
     +                    gam(i) = gam(i) / sills(iii)
                  end if
            end if
c
c  1. report the semivariogram rather than variogram
c  2. report the cross-semivariogram rather than variogram
c  3. the covariance requires "centering"
c  4. the correlogram requires centering and normalizing
c  5. general relative requires division by lag mean
c  6. report the semi(pairwise relative variogram)
c  7. report the semi(log variogram)
c  8. report the semi(madogram)
c
            if(it.eq.1.or.it.eq.2) then
                  gam(i) = 0.5 * gam(i)
            else if(abs(it).eq.3) then
                  gam(i) = gam(i) - hm(i)*tm(i)
                  if(it.lt.0) then
                        if(sills(ivtail(iv)).lt.0.0.or.
     +                     sills(ivhead(iv)).lt.0.0) then
                              gam(i) = -999.0
                        else
                              variance = ( sqrt(sills(ivtail(iv)))
     +                                 *   sqrt(sills(ivhead(iv))) )
                              gam(i) = variance - gam(i)
                        end if
                  end if
            else if(abs(it).eq.4) then
                  hv(i)  = hv(i)-hm(i)*hm(i)
                  if(hv(i).lt.0.0) hv(i) = 0.0
                  hv(i)  = dsqrt(hv(i))
                  tv(i)  = tv(i)-tm(i)*tm(i)
                  if(tv(i).lt.0.0) tv(i) = 0.0
                  tv(i)  = dsqrt(tv(i))
                  if((hv(i)*tv(i)).lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                  endif
                          if(it.lt.0) gam(i) = 1.0 - gam(i)
c
c Square "hv" and "tv" so that we return the variance:
c
                  hv(i)  = hv(i)*hv(i)
                  tv(i)  = tv(i)*tv(i)
            else if(it.eq.5) then
                  htave  = 0.5*(hm(i)+tm(i))
                  htave  = htave   *   htave
                  if(htave.lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) = gam(i)/dble(htave)
                  endif
            else if(it.ge.6) then
                  gam(i) = 0.5 * gam(i)
            endif
 7    continue

cc
cc Write Results:
cc
c      call writeout

c
c Loop over all the variograms that have been computed:
c
      open(lout,file=outfl,status='UNKNOWN')
      do iv=1,nvarg
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) titlewr(1:24) = 'Semivariogram           '
      if(it.eq. 2) titlewr(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) titlewr(1:24) = 'Covariance              '
      if(it.eq. 4) titlewr(1:24) = 'Correlogram             '
      if(it.eq. 5) titlewr(1:24) = 'General Relative        '
      if(it.eq. 6) titlewr(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) titlewr(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) titlewr(1:24) = 'Semimadogram            '
      if(it.eq. 9) titlewr(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) titlewr(1:24) = 'Indicator 1/2 Variogram '
      write(titlewr(25:62),1100) names(ivtail(iv)),names(ivhead(iv))
 1100 format('tail:',a12,' head:',a12)
c
c Loop over all the directions (note the direction in the title):
c
      do id=1,ndir
            write(titlewr(62:74),1101) id
 1101       format('direction ',i2)
            write(lout,'(a74)') titlewr(1:74)
c
c Write out all the lags:
c
            do il=1,nlag+2
                  i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                  nump = int(np(i))
                  if(it.eq.4) then
                        write(lout,1102) il,dis(i),gam(i),nump,
     +                                  hm(i),tm(i),hv(i),tv(i)
                  else
                        write(lout,1102) il,dis(i),gam(i),nump,
     +                                  hm(i),tm(i)
                  endif
 1102             format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
c
c End loop over variograms
c
      end do
c
c Finished:
c
      close(lout)


cc
cc Finished:
cc

      write(*,9998) VERSION
 9998 format(/' GAMV Version: ',f5.3, ' Finished'/)
      stop
      end



c      subroutine readparm
cc-----------------------------------------------------------------------
cc
cc                  Initialization and Read Parameters
cc                  **********************************
cc
cc The input parameters and data are read in from their files. Some quick
cc error checking is performed and the statistics of all the variables
cc being considered are written to standard output.
cc
cc
cc
cc-----------------------------------------------------------------------
c      use geostat
c      parameter(MV=500)
c
c      real      var(MV),cut(MV)
c      real*8    avg(MV),ssq(MV)
c      integer   ivar(MV),num(MV),ivc(MV),indflag(MV)
c      character datafl*512,str*512
c      logical   testfl,testdat
c      real,allocatable :: vrmin(:),vrmax(:)
c
c      data      lin/1/,ncut/0/
cc
cc Note VERSION number:
cc
c      write(*,9999) VERSION
c 9999 format(/' GAMV Version: ',f5.3/)
cc
cc Get the name of the parameter file - try the default name if no input:
cc
c      do i=1,512
c            str(i:i) = ' '
c      end do
c      call getarg(1,str)
c      if(str(1:1).eq.' ')then
c            write(*,*) 'Which parameter file do you want to use?'
c            read (*,'(a)') str
c      end if
c      if(str(1:1).eq.' ') str(1:20) = 'gamv.par            '
c      inquire(file=str,exist=testfl)
c      if(.not.testfl) then
c            write(*,*) 'ERROR - the parameter file does not exist,'
c            write(*,*) '        check for the file and try again  '
c            write(*,*)
c            if(str(1:20).eq.'gamv.par            ') then
c                  write(*,*) '        creating a blank parameter file'
c                  call makepar
c                  write(*,*)
c            end if
c            stop
c      endif
c      open(lin,file=str,status='OLD')
cc
cc Find Start of Parameters:
cc
c 1    read(lin,'(a4)',end=98) str(1:4)
c      if(str(1:4).ne.'STAR') go to 1
cc
cc Read Input Parameters:
cc
c      read(lin,'(a512)',err=98) datafl
c      call chknam(datafl,512)
c      write(*,*) ' data file = ',datafl(1:40)
c
c      read(lin,*,err=98) ixl,iyl,izl
c      write(*,*) ' columns for X,Y,Z = ',ixl,iyl,izl
c
c      read(lin,*,err=98) nvar
c      write(*,*) ' number of variables = ',nvar
c      backspace lin
c
c      read(lin,*,err=98) j,(ivar(i),i=1,nvar)
c      write(*,*) ' columns = ',(ivar(i),i=1,nvar)
c
c      read(lin,*,err=98) tmin,tmax
c      write(*,*) ' trimming limits = ',tmin,tmax
c
c      read(lin,'(a512)',err=98) outfl
c      call chknam(outfl,512)
c      write(*,*) ' output file = ',outfl(1:40)
c
c      read(lin,*,err=98) nlag
c      write(*,*) ' number of lags = ',nlag
c      if(nlag.lt.1)       stop 'nlag is too small: check parameters'
c
c      read(lin,*,err=98) xlag
c      write(*,*) ' lag distance = ',xlag
c      if(xlag.le.0.0) stop 'xlag is too small: check parameter file'
c
c      read(lin,*,err=98) xltol
c      write(*,*) ' lag tolerance = ',xltol
c
c      read(lin,*,err=98) ndir
c      write(*,*) ' number of directions = ',ndir
c
c      allocate (azm(ndir),stat = test)
c      allocate (atol(ndir),stat = test)
c      allocate (bandwh(ndir),stat = test)
c      allocate (dip(ndir),stat = test)
c      allocate (dtol(ndir),stat = test)
c      allocate (bandwd(ndir),stat = test)
c
c      if(ndir.lt.1)       stop 'ndir is too small: check parameters'
c
c      do i=1,ndir
c            read(lin,*,err=98) azm(i),atol(i),bandwh(i),
c     +                         dip(i),dtol(i),bandwd(i)
c            write(*,*) ' azm, atol, bandwh = ',azm(i),atol(i),bandwh(i)
c            write(*,*) ' dip, dtol, bandwd = ',dip(i),dtol(i),bandwd(i)
c            if(bandwh(i).lt.0.0) then
c                  write(*,*) ' Horizontal bandwidth is too small!'
c                  stop
c            endif
c            if(bandwd(i).lt.0.0) then
c                  write(*,*) ' Vertical bandwidth is too small!'
c                  stop
c            endif
c      end do
c
c      read(lin,*,err=98) isill
c      write(*,*) ' flag to standardize sills = ',isill
c
c      read(lin,*,err=98) nvarg
c      write(*,*) ' number of variograms = ',nvarg
c      if(nvarg.lt.1)      stop 'nvarg is too small: check parameters'
c
c      mxdlv=ndir*(nlag+2)*nvarg
c      allocate (dis(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (gam(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (hm(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (tm(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (hv(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (tv(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (np(mxdlv),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (ivtail(nvarg),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (ivhead(nvarg),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (ivtype(nvarg),stat=test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c      allocate (names(nvar+nvarg),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
c
c      ncut = 0
c      do i=1,nvarg
c            read(lin,*,err=98) ivtail(i),ivhead(i),ivtype(i)
c            write(*,*) ' tail,head,type = ',
c     +                   ivtail(i),ivhead(i),ivtype(i)
c            if(ivtype(i).eq.9.or.ivtype(i).eq.10) then
c                   ncut = ncut + 1
c                   if(tmin.gt.0.0)stop'tmin interferes with indicators!'
c                   if(tmax.le.1.0)stop'tmax interferes with indicators!'
c                   backspace lin
c                   read(lin,*,err=98) ii,jj,kk,cut(ncut)
c                   if(ivtype(i).eq.9)  indflag(ncut) = 1
c                   if(ivtype(i).eq.10) indflag(ncut) = 0
c                   ivc(ncut) = ivtail(i)
c                   ivtail(i) = nvar + ncut
c                   ivhead(i) = nvar + ncut
c                   write(names(nvar+ncut),140) ncut
c 140               format('Indicator ',i2)
c                   write(*,*) ' indicator threshold = ',cut(ncut)
c            endif
c      end do
c      write(*,*)
c      close(lin)
c      MAXVAR = nvar + ncut
cc
cc Perform some quick error checking:
cc
c      if(xltol.le.0.0) then
c            write(*,*) 'xltol is too small: resetting to xlag/2'
c            xltol = 0.5*xlag
c      endif
cc
cc Check to make sure the data file exists, then either read in the
cc data or write an error message and stop:
cc
c      inquire(file=datafl,exist=testfl)
c      if(.not.testfl) then
c            write(*,*) 'ERROR data file ',datafl,' does not exist!'
c            stop
c      endif
cc
cc The data file exists so open the file and read in the header
cc information. Initialize the storage that will be used to summarize
cc the data found in the file:
cc
c      open(lin,file=datafl,status='OLD')
cc
c      read(lin,*,err=99)
c      read(lin,*,err=99)       nvari
c      do i=1,nvari
c            read(lin,*)
c      end do
c      maxdat = 0
c 22   read(lin,*,end=44,err=99)(var(j),j=1,nvari)
c      maxdat = maxdat + 1
c      go to 22
c 44   continue
c      write(*,*)'maxdat = ',maxdat
cc
c      allocate (vr(maxdat,MAXVAR),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (vrmin(MAXVAR),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (vrmax(MAXVAR),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (sills(MAXVAR),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (x(maxdat),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (y(maxdat),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      allocate (z(maxdat),stat = test)
c      if (test.ne.0) then
c            write(*,*) 'Error: Allocation failed due to ',
c     +                 'insufficient memory!', test
c            stop
c      end if
cc
c      rewind(lin)
c      read(lin,'(a40)',err=99) str
c      read(lin,*,err=99)       nvari
c      do i=1,nvari
c            read(lin,'(a40)',err=99) str
c            do iv=1,nvar
c                  j=ivar(iv)
c                  if(i.eq.j) names(iv) = str(1:12)
c            end do
c            num(i) = 0
c            avg(i) = 0.0
c            ssq(i) = 0.0
c      end do
cc
cc Read all the data until the end of the file:
cc
c      nd = 0
c 2    continue
c      read(lin,*,end=9,err=99) (var(j),j=1,nvari)
c      testdat = .false.
c      do iv=1,nvar
c            j=ivar(iv)
c            if(var(j).ge.tmin.and.var(j).lt.tmax) testdat = .true.
c      end do
c      if(.not.testdat) go to 2
c      nd = nd + 1
cc
cc Acceptable data, make sure there are not too many data:
cc
c      do iv=1,nvar
c            j=ivar(iv)
c            vr(nd,iv) = var(j)
c            if(var(j).ge.tmin.and.var(j).lt.tmax) then
c                  num(iv) = num(iv) + 1
c                  avg(iv) = avg(iv) + dble(var(j))
c                  ssq(iv) = ssq(iv) + dble(var(j)*var(j))
c            endif
c      end do
c      if(ixl.le.0) then
c            x(nd) = 0.0
c      else
c            x(nd) = var(ixl)
c      endif
c      if(iyl.le.0) then
c            y(nd) = 0.0
c      else
c            y(nd) = var(iyl)
c      endif
c      if(izl.le.0) then
c            z(nd) = 0.0
c      else
c            z(nd) = var(izl)
c      endif
c      go to 2
c 9    continue
c      close(lin)
cc
cc Compute the averages and variances as an error check for the user:
cc
c      do iv=1,nvar
c            sills(iv) = -999.
c            if(num(iv).gt.0) then
c                  avg(iv)   = avg(iv)/dble(num(iv))
c                  ssq(iv)   =(ssq(iv)/dble(num(iv)))-avg(iv)*avg(iv)
c                  sills(iv) = ssq(iv)
c                  write(*,*) 'Variable number ',iv
c                  write(*,*) '  Number   = ',num(iv)
c                  write(*,*) '  Average  = ',real(avg(iv))
c                  write(*,*) '  Variance = ',real(ssq(iv))
c            endif
c      end do
cc
cc Construct Indicator Variables if necessary:
cc
c      do ic=1,ncut
c            iv   = ivc(ic)
c            jv   = nvar + ic
c            ptot = 0.0
c            p1   = 0.0
c            do id=1,nd
c                  if(vr(id,iv).le.tmin.or.vr(id,iv).gt.tmax) then
c                        vr(id,jv) = tmin - EPSLON
c                  else
c                        if(indflag(ic).eq.1) then
c                              if(vr(id,iv).lt.cut(ic)) then
c                                    vr(id,jv) = 0.0
c                              else
c                                    vr(id,jv) = 1.0
c                              endif
c                              p1   = p1   + vr(id,jv)
c                              ptot = ptot + 1.0
c                        else
c                              vr(id,jv) = 0.0
c                              if(int(vr(id,iv)+0.5).eq.int(cut(ic)+0.5))
c     +                        vr(id,jv) = 1.0
c                              p1   = p1   + vr(id,jv)
c                              ptot = ptot + 1.0
c                        end if
c                  end if
c            end do
c            p1        = p1 / max(ptot,1.0)
c            sills(jv) = dble (p1*(1.0-p1))
c      end do
cc
cc Establish minimums and maximums:
cc
c      do i=1,MAXVAR
c            vrmin(i) =  1.0e21
c            vrmax(i) = -1.0e21
c      end do
c      do id=1,nd
c            do iv=1,nvar+ncut
c                  if(vr(id,iv).ge.tmin.and.vr(id,iv).lt.tmax) then
c                        if(vr(id,iv).lt.vrmin(iv)) vrmin(iv) = vr(id,iv)
c                        if(vr(id,iv).gt.vrmax(iv)) vrmax(iv) = vr(id,iv)
c                  end if
c            end do
c      end do
cc
cc Check on the variogams that were requested:
cc
c      call check(vrmin,vrmax)
c      return
cc
cc Error in an Input File Somewhere:
cc
c 98   stop 'ERROR in parameter file!'
c 99   stop 'ERROR in data file!'
c      end
c
c
c
c      subroutine gamv
cc-----------------------------------------------------------------------
cc
cc              Variogram of 3-D Irregularly Spaced Data
cc              ****************************************
cc
cc This subroutine computes a variety of spatial continuity measures of a
cc set for irregularly spaced data.  The user can specify any combination
cc of direct and cross variograms using any of eight "variogram" measures
cc
cc
cc
cc INPUT VARIABLES:
cc
cc   nd               Number of data (no missing values)
cc   x(nd)            X coordinates of the data
cc   y(nd)            Y coordinates of the data
cc   z(nd)            Z coordinates of the data
cc   nv               The number of variables
cc   vr(nd,nv)        Data values
cc   tmin,tmax        Trimming limits
cc   nlag             Number of lags to calculate
cc   xlag             Length of the unit lag
cc   xltol            Distance tolerance (if <0 then set to xlag/2)
cc   ndir             Number of directions to consider
cc   azm(ndir)        Azimuth angle of direction (measured positive
cc                      degrees clockwise from NS).
cc   atol(ndir)       Azimuth (half window) tolerances
cc   bandwh           Maximum Horizontal bandwidth (i.e., the deviation
cc                      perpendicular to the defined azimuth).
cc   dip(ndir)        Dip angle of direction (measured in negative
cc                      degrees down from horizontal).
cc   dtol(ndir)       Dip (half window) tolerances
cc   bandwd           Maximum "Vertical" bandwidth (i.e., the deviation
cc                      perpendicular to the defined dip).
cc   isill            1=attempt to standardize, 0=do not
cc   sills            the sills (variances) to standardize with
cc   nvarg            Number of variograms to compute
cc   ivtail(nvarg)    Variable for the tail of each variogram
cc   ivhead(nvarg)    Variable for the head of each variogram
cc   ivtype(nvarg)    Type of variogram to compute:
cc                      1. semivariogram
cc                      2. cross-semivariogram
cc                      3. covariance
cc                      4. correlogram
cc                      5. general relative semivariogram
cc                      6. pairwise relative semivariogram
cc                      7. semivariogram of logarithms
cc                      8. rodogram
cc                      9. indicator semivariogram (continuous)
cc                     10. indicator semivariogram (categorical)
cc
cc
cc
cc OUTPUT VARIABLES:  The following arrays are ordered by direction,
cc                    then variogram, and finally lag, i.e.,
cc                      iloc = (id-1)*nvarg*MAXLG+(iv-1)*MAXLG+il
cc
cc   np()             Number of pairs
cc   dis()            Distance of pairs falling into this lag
cc   gam()            Semivariogram, covariance, correlogram,... value
cc   hm()             Mean of the tail data
cc   tm()             Mean of the head data
cc   hv()             Variance of the tail data
cc   tv()             Variance of the head data
cc
cc
cc
cc PROGRAM NOTES:
cc
cc   1. The file "gamv.inc" contains the dimensioning parameters.
cc      These may have to be changed depending on the amount of data
cc      and the requirements to compute different variograms.
cc
cc
cc
cc Original:  A.G. Journel                                           1978
cc Revisions: K. Guertin                                             1980
cc-----------------------------------------------------------------------
c
c
c
c      use geostat
c      parameter(PI=3.14159265)
c      real      uvxazm(100),uvyazm(100),uvzdec(100),
c     +          uvhdec(100),csatol(100),csdtol(100)
c      logical   omni
c
c      integer threadId
c
cc
cc Define the distance tolerance if it isn't already:
cc
c      if(xltol.le.0.0) xltol = 0.5 * xlag
cc
cc Define the angles and tolerance for each direction:
cc
c      do id=1,ndir
cc
cc The mathematical azimuth is measured counterclockwise from EW and
cc not clockwise from NS as the conventional azimuth is:
cc
c            azmuth     = (90.0-azm(id))*PI/180.0
c            uvxazm(id) = cos(azmuth)
c            uvyazm(id) = sin(azmuth)
c            if(atol(id).le.0.0) then
c                  csatol(id) = cos(45.0*PI/180.0)
c            else
c                  csatol(id) = cos(atol(id)*PI/180.0)
c            endif
cc
cc The declination is measured positive down from vertical (up) rather
cc than negative down from horizontal:
cc
c            declin     = (90.0-dip(id))*PI/180.0
c            uvzdec(id) = cos(declin)
c            uvhdec(id) = sin(declin)
c            if(dtol(id).le.0.0) then
c                  csdtol(id) = cos(45.0*PI/180.0)
c            else
c                  csdtol(id) = cos(dtol(id)*PI/180.0)
c            endif
c      end do
cc
cc Initialize the arrays for each direction, variogram, and lag:
cc
c      nsiz = ndir*nvarg*(nlag+2)
c      do i=1,nsiz
c            np(i)  = 0.
c            dis(i) = 0.0
c            gam(i) = 0.0
c            hm(i)  = 0.0
c            tm(i)  = 0.0
c            hv(i)  = 0.0
c            tv(i)  = 0.0
c      end do
c      dismxs = ((real(nlag) + 0.5 - EPSLON) * xlag) ** 2
cc
cc MAIN LOOP OVER ALL PAIRS:
cc
c      irepo = max(1,min((nd/10),1000))
c
c
cc$omp parallel default(firstprivate)
c
c#ifdef _OPENMP
c      threadId = int(OMP_get_thread_num())
c#else
c
c#endif
c
cc$omp do
c      do 3 i=1,nd
c        if((int(i/irepo)*irepo).eq.i) write(*,103) i,nd
c 103  format('   currently on ',i9,' of ',i9)
c      do 4 j=i,nd
cc
cc Definition of the lag corresponding to the current pair:
cc
c            dx  = x(j) - x(i)
c            dy  = y(j) - y(i)
c            dz  = z(j) - z(i)
c            dxs = dx*dx
c            dys = dy*dy
c            dzs = dz*dz
c            hs  = dxs + dys + dzs
c            if(hs.gt.dismxs) go to 4
c            if(hs.lt.0.0) hs = 0.0
c            h   = sqrt(hs)
cc
cc Determine which lag this is and skip if outside the defined distance
cc tolerance:
cc
c            if(h.le.EPSLON) then
c                  lagbeg = 1
c                  lagend = 1
c            else
c                  lagbeg = -1
c                  lagend = -1
c                  do ilag=2,nlag+2
c                        if(h.ge.(xlag*real(ilag-2)-xltol).and.
c     +                     h.le.(xlag*real(ilag-2)+xltol)) then
c                              if(lagbeg.lt.0) lagbeg = ilag
c                              lagend = ilag
c                        end if
c                  end do
c                  if(lagend.lt.0) go to 4
c            endif
cc
cc Definition of the direction corresponding to the current pair. All
cc directions are considered (overlapping of direction tolerance cones
cc is allowed):
cc
c            do 5 id=1,ndir
cc
cc Check for an acceptable azimuth angle:
cc
c                  dxy = sqrt(max((dxs+dys),0.0))
c                  if(dxy.lt.EPSLON) then
c                        dcazm = 1.0
c                  else
c                        dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
c                  endif
c                  if(abs(dcazm).lt.csatol(id)) go to 5
cc
cc Check the horizontal bandwidth criteria (maximum deviation
cc perpendicular to the specified direction azimuth):
cc
c                  band = uvxazm(id)*dy - uvyazm(id)*dx
c                  if(abs(band).gt.bandwh(id)) go to 5
cc
cc Check for an acceptable dip angle:
cc
c                  if(dcazm.lt.0.0) dxy = -dxy
c                  if(lagbeg.eq.1) then
c                        dcdec = 0.0
c                  else
c                        dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
c                        if(abs(dcdec).lt.csdtol(id)) go to 5
c                  endif
cc
cc Check the vertical bandwidth criteria (maximum deviation perpendicular
cc to the specified dip direction):
cc
c                  band = uvhdec(id)*dz - uvzdec(id)*dxy
c                  if(abs(band).gt.bandwd(id)) go to 5
cc
cc Check whether or not an omni-directional variogram is being computed:
cc
c                  omni = .false.
c                  if(atol(id).ge.90.0) omni = .true.
cc
cc This direction is acceptable - go ahead and compute all variograms:
cc
c                  do 6 iv=1,nvarg
cc
cc For this variogram, sort out which is the tail and the head value:
cc
c                      it = ivtype(iv)
c                      if(dcazm.ge.0.0.and.dcdec.ge.0.0) then
c                            ii = ivtail(iv)
c                            vrh   = vr(i,ii)
c                            ii = ivhead(iv)
c                            vrt   = vr(j,ii)
c                            if(omni.or.it.eq.2) then
c                                  ii    = ivhead(iv)
c                                  vrtpr = vr(i,ii)
c                                  ii    = ivtail(iv)
c                                  vrhpr = vr(j,ii)
c                            endif
c                      else
c                            ii = ivtail(iv)
c                            vrh   = vr(j,ii)
c                            ii = ivhead(iv)
c                            vrt   = vr(i,ii)
c                            if(omni.or.it.eq.2) then
c                                  ii    = ivhead(iv)
c                                  vrtpr = vr(j,ii)
c                                  ii    = ivtail(iv)
c                                  vrhpr = vr(i,ii)
c                            endif
c                      endif
cc
cc Reject this pair on the basis of missing values:
cc
c                      if(vrt.lt.tmin.or.vrh.lt.tmin.or.
c     +                   vrt.gt.tmax.or.vrh.gt.tmax) go to 6
c                      if(it.eq.2.and.(vrtpr.lt.tmin.or.vrhpr.lt.tmin.or.
c     +                                vrtpr.gt.tmax.or.vrhpr.gt.tmax))
c     +                                               go to 6
cc
cc             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
cc
cc
cc The Semivariogram:
cc
c      if(it.eq.1.or.it.eq.5.or.it.ge.9) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble((vrhpr-vrtpr)*
c     +                                              (vrhpr-vrtpr))
c                     endif
c               endif
c            end do
cc
cc The Traditional Cross Semivariogram:
cc
c      else if(it.eq.2) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
c               hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
c               gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
c            end do
cc
cc The Covariance:
cc
c      else if(abs(it).eq.3) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble(vrh*vrt)
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
c                     endif
c               endif
c            end do
cc
cc The Correlogram:
cc
c      else if(abs(it).eq.4) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               hv(ii)  = hv(ii)  + dble(vrh*vrh)
c               tv(ii)  = tv(ii)  + dble(vrt*vrt)
c               gam(ii) = gam(ii) + dble(vrh*vrt)
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
c                           tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
c                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
c                     endif
c               endif
c            end do
cc
cc The Pairwise Relative:
cc
c      else if(it.eq.6) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               if(abs(vrt+vrh).gt.EPSLON) then
c                     np(ii)  = np(ii)  + 1.
c                     dis(ii) = dis(ii) + dble(h)
c                     tm(ii)  = tm(ii)  + dble(vrt)
c                     hm(ii)  = hm(ii)  + dble(vrh)
c                     gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
c                     gam(ii) = gam(ii) + dble(gamma*gamma)
c               endif
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                     if(abs(vrtpr+vrhpr).gt.EPSLON) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
c                           gam(ii) = gam(ii) + dble(gamma*gamma)
c                     endif
c                     endif
c               endif
c            enddo
cc
cc Variogram of Logarithms:
cc
c      else if(it.eq.7) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
c                     np(ii)  = np(ii)  + 1.
c                     dis(ii) = dis(ii) + dble(h)
c                     tm(ii)  = tm(ii)  + dble(vrt)
c                     hm(ii)  = hm(ii)  + dble(vrh)
c                     gamma   = alog(vrt)-alog(vrh)
c                     gam(ii) = gam(ii) + dble(gamma*gamma)
c               endif
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gamma   = alog(vrt)-alog(vrh)
c                           gam(ii) = gam(ii) + dble(gamma*gamma)
c                     endif
c                     endif
c               endif
c            end do
cc
cc Madogram:
cc
c      else if(it.eq.8) then
c            do il=lagbeg,lagend
c               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c               np(ii)  = np(ii)  + 1.
c               dis(ii) = dis(ii) + dble(h)
c               tm(ii)  = tm(ii)  + dble(vrt)
c               hm(ii)  = hm(ii)  + dble(vrh)
c               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
c               if(omni) then
c                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and.
c     +                  vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
c                           np(ii)  = np(ii)  + 1.
c                           dis(ii) = dis(ii) + dble(h)
c                           tm(ii)  = tm(ii)  + dble(vrtpr)
c                           hm(ii)  = hm(ii)  + dble(vrhpr)
c                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
c                     endif
c               endif
c            end do
c      endif
cc
cc Finish loops over variograms, directions, and the double data loops:
cc
c 6              continue
c 5          continue
c 4    continue
c 3    continue
cc$omp end do
cc$omp end parallel
c
cc
cc Get average values for gam, hm, tm, hv, and tv, then compute
cc the correct "variogram" measure:
cc
c      do 7 id=1,ndir
c      do 7 iv=1,nvarg
c      do 7 il=1,nlag+2
c            i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c            if(np(i).le.0.) go to 7
c            rnum   = np(i)
c            dis(i) = dis(i) / dble(rnum)
c            gam(i) = gam(i) / dble(rnum)
c            hm(i)  = hm(i)  / dble(rnum)
c            tm(i)  = tm(i)  / dble(rnum)
c            hv(i)  = hv(i)  / dble(rnum)
c            tv(i)  = tv(i)  / dble(rnum)
c            it = ivtype(iv)
cc
cc Attempt to standardize:
cc
c            if(isill.eq.1) then
c                  if(ivtail(iv).eq.ivhead(iv)) then
c                        iii = ivtail(iv)
c                        if((it.eq.1.or.it.ge.9).and.sills(iii).gt.0.0)
c     +                    gam(i) = gam(i) / sills(iii)
c                  end if
c            end if
cc
cc  1. report the semivariogram rather than variogram
cc  2. report the cross-semivariogram rather than variogram
cc  3. the covariance requires "centering"
cc  4. the correlogram requires centering and normalizing
cc  5. general relative requires division by lag mean
cc  6. report the semi(pairwise relative variogram)
cc  7. report the semi(log variogram)
cc  8. report the semi(madogram)
cc
c            if(it.eq.1.or.it.eq.2) then
c                  gam(i) = 0.5 * gam(i)
c            else if(abs(it).eq.3) then
c                  gam(i) = gam(i) - hm(i)*tm(i)
c                  if(it.lt.0) then
c                        if(sills(ivtail(iv)).lt.0.0.or.
c     +                     sills(ivhead(iv)).lt.0.0) then
c                              gam(i) = -999.0
c                        else
c                              variance = ( sqrt(sills(ivtail(iv)))
c     +                                 *   sqrt(sills(ivhead(iv))) )
c                              gam(i) = variance - gam(i)
c                        end if
c                  end if
c            else if(abs(it).eq.4) then
c                  hv(i)  = hv(i)-hm(i)*hm(i)
c                  if(hv(i).lt.0.0) hv(i) = 0.0
c                  hv(i)  = dsqrt(hv(i))
c                  tv(i)  = tv(i)-tm(i)*tm(i)
c                  if(tv(i).lt.0.0) tv(i) = 0.0
c                  tv(i)  = dsqrt(tv(i))
c                  if((hv(i)*tv(i)).lt.EPSLON) then
c                        gam(i) = 0.0
c                  else
c                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
c                  endif
c                          if(it.lt.0) gam(i) = 1.0 - gam(i)
cc
cc Square "hv" and "tv" so that we return the variance:
cc
c                  hv(i)  = hv(i)*hv(i)
c                  tv(i)  = tv(i)*tv(i)
c            else if(it.eq.5) then
c                  htave  = 0.5*(hm(i)+tm(i))
c                  htave  = htave   *   htave
c                  if(htave.lt.EPSLON) then
c                        gam(i) = 0.0
c                  else
c                        gam(i) = gam(i)/dble(htave)
c                  endif
c            else if(it.ge.6) then
c                  gam(i) = 0.5 * gam(i)
c            endif
c 7    continue
c      return
c      end
c
c
c
c      subroutine writeout
cc-----------------------------------------------------------------------
cc
cc                  Write Out the Results of GAMV3
cc                  ******************************
cc
cc An output file will be written which contains each directional
cc variogram ordered by direction and then variogram (the directions
cc cycle fastest then the variogram number).  For each variogram there
cc will be a one line description and then "nlag" lines with:
cc
cc        a) lag number (increasing from 1 to nlag)
cc        b) separation distance
cc        c) the "variogram" value
cc        d) the number of pairs for the lag
cc        e) the mean of the data contributing to the tail
cc        f) the mean of the data contributing to the head
cc        g) IF the correlogram - variance of tail values
cc        h) IF the correlogram - variance of head values
cc
cc
cc
cc
cc
cc-----------------------------------------------------------------------
c      use geostat
c      character title*132
c      data      lout/1/
cc
cc Loop over all the variograms that have been computed:
cc
c      open(lout,file=outfl,status='UNKNOWN')
c      do iv=1,nvarg
cc
cc Construct a title that reflects the variogram type and the variables
cc that were used to calculate the variogram:
cc
c      it = abs(ivtype(iv))
c      if(it.eq. 1) title(1:24) = 'Semivariogram           '
c      if(it.eq. 2) title(1:24) = 'Cross Semivariogram     '
c      if(it.eq. 3) title(1:24) = 'Covariance              '
c      if(it.eq. 4) title(1:24) = 'Correlogram             '
c      if(it.eq. 5) title(1:24) = 'General Relative        '
c      if(it.eq. 6) title(1:24) = 'Pairwise Relative       '
c      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms '
c      if(it.eq. 8) title(1:24) = 'Semimadogram            '
c      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram '
c      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram '
c      write(title(25:62),100) names(ivtail(iv)),names(ivhead(iv))
c 100  format('tail:',a12,' head:',a12)
cc
cc Loop over all the directions (note the direction in the title):
cc
c      do id=1,ndir
c            write(title(62:74),101) id
c 101        format('direction ',i2)
c            write(lout,'(a74)') title(1:74)
cc
cc Write out all the lags:
cc
c            do il=1,nlag+2
c                  i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
c                  nump = int(np(i))
c                  if(it.eq.4) then
c                        write(lout,102) il,dis(i),gam(i),nump,
c     +                                  hm(i),tm(i),hv(i),tv(i)
c                  else
c                        write(lout,102) il,dis(i),gam(i),nump,
c     +                                  hm(i),tm(i)
c                  endif
c 102              format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
c            end do
c      end do
cc
cc End loop over variograms
cc
c      end do
cc
cc Finished:
cc
c      close(lout)
c      return
c      end
c
c
c
c      subroutine check(vrmin,vrmax)
cc-----------------------------------------------------------------------
cc
cc                Error Check and Note Variogram Types
cc                ************************************
cc
cc Go through each variogram type and note the type to the screen and
cc report any possible errors.
cc
cc
cc
cc
cc
cc-----------------------------------------------------------------------
c      use geostat
c      real      vrmin(*),vrmax(*)
c      character title*80
cc
cc Loop over all the variograms to be computed:
cc
c      write(*,*)
c      do iv=1,nvarg
cc
cc Note the variogram type and the variables being used:
cc
c      it = abs(ivtype(iv))
c      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
c      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
c      if(it.eq. 3) title(1:24) = 'Covariance             :'
c      if(it.eq. 4) title(1:24) = 'Correlogram            :'
c      if(it.eq. 5) title(1:24) = 'General Relative       :'
c      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
c      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
c      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
c      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
c      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
c      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
c 100  format('  tail=',a12,' head=',a12)
c      write(*,101) iv,title(1:64)
c 101  format(' Variogram ',i2,1x,a64)
cc
cc Check for possible errors or inconsistencies:
cc
c      if(it.eq.2) then
c            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
c 201        format('  WARNING: cross variogram with the same variable!')
c      else if(it.eq.5) then
c            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
c            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
c     +            write(*,502)
c            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
c     +            write(*,502)
c 501        format('  WARNING: cross general relative variogram are',
c     +             ' difficult to interpret!')
c 502        format('  WARNING: there are both positive and negative',
c     +             ' values - lag mean could be zero!')
c      else if(it.eq.6) then
c            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
c            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
c     +            write(*,602)
c            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
c     +            write(*,602)
c 601        format('  WARNING: cross pairwise relative variogram are',
c     +             ' difficult to interpret!')
c 602        format('  WARNING: there are both positive and negative',
c     +             ' values - pair means could be zero!')
c      else if(it.eq.7) then
c            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
c            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0)
c     +      write(*,702)
c 701        format('  WARNING: cross logarithmic variograms may be',
c     +             ' difficult to interpret!')
c 702        format('  WARNING: there are zero or negative',
c     +             ' values - logarithm undefined!')
c      else if(it.eq.8) then
c            if(ivtail(iv).ne.ivhead(iv)) write(*,801)
c 801        format('  WARNING: cross rodograms may be difficult to',
c     +             ' interpret!')
c      else if(it.eq.9) then
c            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
c 901        format('  WARNING: cross madograms may be difficult to',
c     +             ' interpret!')
c      endif
cc
cc END Loop over all variograms:
cc
c      end do
cc
cc Finished:
cc
c      return
c      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='gamv.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for GAMV',/,
     +       '                  *******************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat               ',
     +       '-file with data')
      write(lun,12)
 12   format('1   2   0                         ',
     +       '-   columns for X, Y, Z coordinates')
      write(lun,13)
 13   format('2   3   4                         ',
     +       '-   number of variables,col numbers')
      write(lun,14)
 14   format('-1.0e21     1.0e21                ',
     +       '-   trimming limits')
      write(lun,15)
 15   format('gamv.out                          ',
     +       '-file for variogram output')
      write(lun,16)
 16   format('10                                ',
     +       '-number of lags')
      write(lun,17)
 17   format('5.0                               ',
     +       '-lag separation distance')
      write(lun,18)
 18   format('3.0                               ',
     +       '-lag tolerance')
      write(lun,19)
 19   format('3                                 ',
     +       '-number of directions')
      write(lun,20)
 20   format('0.0  90.0 50.0   0.0  90.0  50.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,21)
 21   format('0.0  22.5 25.0   0.0  22.5  25.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,22)
 22   format('90.  22.5 25.0   0.0  22.5  25.0  ',
     +       '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,23)
 23   format('0                                 ',
     +       '-standardize sills? (0=no, 1=yes)')
      write(lun,24)
 24   format('3                                 ',
     +       '-number of variograms')
      write(lun,25)
 25   format('1   1   1                         ',
     +       '-tail var., head var., variogram type')
      write(lun,26)
 26   format('1   2   2                         ',
     +       '-tail var., head var., variogram type')
      write(lun,27)
 27   format('2   2   1                         ',
     +       '-tail var., head var., variogram type')
      write(lun,40)
 40   format(//,'type 1 = traditional semivariogram',/,
     +          '     2 = traditional cross semivariogram',/,
     +          '     3 = covariance',/,
     +          '     4 = correlogram',/,
     +          '     5 = general relative semivariogram',/,
     +          '     6 = pairwise relative semivariogram',/,
     +          '     7 = semivariogram of logarithms',/,
     +          '     8 = semimadogram',/,
     +          '     9 = indicator semivariogram - continuous',/,
     +          '     10= indicator semivariogram - categorical')

      close(lun)
      return
      end
