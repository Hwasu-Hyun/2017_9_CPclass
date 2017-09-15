      program MCMC_SNIa_anal
c---  Cosmological parameter estimation of the FRW world model
c---  using the supernova type Ia redshift-distance data. 
c---  + CPL DE equation of state w(a)=w_0 + w_a (1-a/a_0)  2009/07/28
c---  + BAO A parameter  2009/08/06
      implicit double precision (a-h,o-z)
      !------- Metropolis Algorithm Variables ----------
      parameter (mp=5000000,nb=30)
      double precision chg(mp,nb)
      character*80     file_chains,file_array,file_clevel
      character*1      c100,c10,c1
      double precision xp(mp),yp(mp)
      !------- Likelihood Map ----------
      parameter (nx=500,nthres=5)
      real array1(0:nx),array2(0:nx)
      double precision thres(nthres),vfrac(nthres),clevel(nthres)
      double precision xleft(nthres),xright(nthres)
      parameter (nhist=50)
      double precision hist(0:nhist)
      character*80     file_hist
c     data thres/1d0,2d0,3d0,4d0,5d0/ 
c     data vfrac/0.682689492137086d0,0.954499736103642d0,
c    +           0.997300203936740d0,0.999936657516334d0,     
c    +           0.999999426696856d0/    
      ! savgol
      parameter(nmaxc=1000,ld=0,m=4)
      double precision c(nmaxc)

      pi=4d0*datan(1d0)   ! PI

      file_chains='normal_density_0127_finalv2.txt'
      file_array ='contour_array_1d_0127_final.dat'
      file_clevel='contour_array_clevel_1d_0127_final.dat'
      file_hist='histogram_0127_final.dat'

      ip=1
      open(2,file=file_chains,status='old')
10    read(2,*,end=11) (chg(ip,ib),ib=1,30)
      ip=ip+1
      goto 10
11    continue
      np=ip-1
      print *,'Number of chain elements=',np
      close(2)

      print *,'Input bin number (ib=1,30):'
      read(5,*) ib
      if(ib.lt.1.or.ib.gt.30) stop 'wrong ib.'
      xp=chg(:,ib)

      !print *,'Choose density assigning scheme: (1) Simple (2) CIC '
      !read(5,*) iassign
      iassign=2
      if(iassign.lt.1.or.iassign.gt.2) stop 'wrong assigning scheme.'

      xmin=+1d30
      xmax=-1d30
      xmean=0d0
      xstd=0d0
      ic=0
      do i=1,np
         if(xp(i).ge.0.d0) then
            if(xp(i).le.xmin) xmin=xp(i)
            if(xp(i).gt.xmax) xmax=xp(i)
            xmean=xmean+xp(i)
            xstd=xstd+xp(i)**2
            ic=ic+1
         endif
      enddo
      xmean=xmean/dble(ic)
      xstd=dsqrt(xstd/dble(ic)-xmean**2)
      print *,'X(chain): min & max:',xmin,xmax
      print *,'X(chain): mean & std:',xmean,xstd

      f10=(xmax-xmin)/10d0
      xmin=xmin-f10
      xmax=xmax+f10

!--------------------------------------------------------
!     histogram
!----------------------------------------------------------
      nh=30
      do ih=0,nh
         hist(ih)=0.0
      enddo
      dx_hist=(xmax-xmin)/dble(nh) 
      do ip=1,np
         if(xp(ip).ge.0.d0) then
            ih=nint((xp(ip)-xmin)/dx_hist)
            hist(ih)=hist(ih)+1.0
          endif
      enddo
      open(10,file=file_hist)
      do ih=0,nh
         xcrd=dble(ih)*dx_hist+xmin
         write(10,*) ih,xcrd,hist(ih)
      enddo
      close(10)

c---  Save Likelihood map and determine contour levels
      dx=(xmax-xmin)/dble(nx) ! delta in x dimension
      print *,'X(grid): min & max = ',xmin,xmax
      print *,'dx = ',dx

      ! save pixelized data into array 
      do ix=0,nx
         array1(ix)=0.0
      enddo

      if(iassign.eq.1) then 
      ! particle assignment into gridded array
      ! using simple method
         do ip=1,np
            if(xp(ip).ge.0.d0) then
               ix=nint((xp(ip)-xmin)/dx)
               array1(ix)=array1(ix)+1.0
            endif
         enddo
      else if(iassign.eq.2) then
      ! particle assignment into gridded array
      ! using Cloud-In-Cell (CIC) scheme
         do ip=1,np
            if(xp(ip).ge.0.d0) then
               xcrd=(xp(ip)-xmin)/dx ! x coordinate in the array
               ix=nint(xcrd)
               xsep=xcrd-dble(ix) ! separation wrt ix (major cell)
               ixx=ix+int(dsign(1d0,xsep))
               x_wt=abs(xsep) ! weight
               array1(ix)=array1(ix)+real(1d0-x_wt)
               array1(ixx)=array1(ixx)+real(x_wt)
            endif
         enddo 
      endif

      !goto 999
      ! applying savgol filtering algorithm
      iwx=20
      array2=0.0
      do ix=0,nx
         nl=iwx
         nr=iwx
         if(ix-nl.lt.0) nl=ix
         if(ix+nr.gt.nx) nr=nx-ix
         nps=nl+nr+1
         call savgol(c,nps,nl,nr,ld,m)
         xsmoo1=0d0
         do j=nl+1,1,-1
            xsmoo1=xsmoo1+c(j)*dble(array1(ix-j+1))
         enddo
         xsmoo2=0d0
         do j=1,nr
            xsmoo2=xsmoo2+c(nps-j+1)*dble(array1(ix+j))
         enddo
         array2(ix)=real(xsmoo1+xsmoo2)
      enddo
      array1=array2
      !-----------------------------------------------------
999   continue 

      amin=+1.0d40
      amax=-1.0d40
      amax_pos=0d0
      do ix=0,nx
         if(array1(ix).gt.amax) then
            amax=array1(ix)
            amax_xcrd=dble(ix)*dx+xmin
            ipos_max=ix
         endif
         if(array1(ix).lt.amin) amin=array1(ix)
         !if(array1(ix).ne.0.0) print *,ix,array1(ix)
      enddo
      !print *,'Min Max in Array : ',amin,amax
      write(*,fmt='(4f10.5)') amax_xcrd,xstd

      open(7,file=file_array)
      do ix=0,nx
         xcrd=dble(ix)*dx+xmin
         write(7,*) ix,xcrd,array1(ix)
      enddo
      close(7)

      !== Setting contour levels  ===
      ! set confidence limits using the error function defined as
      ! erf(x)=(2/sqrt{pi}) \int_0^x exp(-t^2)dt
      !    =2 (1/sqrt{2pi}) \int_0^{sqrt{2} x} exp(-t^2/2) dt
      !    =2 \int_0^{sqrt{2} x} [Gaussian func(t), m=0,s=1]dt
      ! where x=threshold level / sqrt{2}.
      do ilevel=1,nthres
         thres(ilevel)=dble(ilevel) ! simply 1 sigma, 2 sigma, ...
         vfrac(ilevel)=erf(thres(ilevel)/dsqrt(2d0))
         !write(*,*) ilevel,thres(ilevel),erf(thres(ilevel)/dsqrt(2d0))
      enddo

      ntot=0
      vtot=0d0
      do ix=0,nx
         ntot=ntot+1
         vtot=vtot+dble(array1(ix))
      enddo
      print *,'Total useful pixels :',ntot
      print *,'Total volume of prob. mountain :',vtot

      open(7,file=file_clevel)
      write(7,fmt='(6(1pe15.5))') xmin,xmax,dx
      write(7,fmt='(3(1pe15.5))') amax_xcrd,amax

      !-- determine area fraction threshold levels (clevel)
      do ilevel=1,nthres
         vol_thr=vtot*vfrac(ilevel)
         iup=np
         idn=0
50       continue
         iav=(iup+idn)/2
         iavp1=iav+1
         c_up=0d0
         c_dn=0d0
         c_av=0d0
         c_avp1=0d0
         do ix=0,nx
            if(array1(ix).ge.iup) c_up=c_up+dble(array1(ix))
            if(array1(ix).ge.idn) c_dn=c_dn+dble(array1(ix))
            if(array1(ix).ge.iav) c_av=c_av+dble(array1(ix))
            if(array1(ix).ge.iavp1) c_avp1=c_avp1+dble(array1(ix))
         enddo
         if(abs(iup-iav).eq.0.or.abs(idn-iav).eq.0) then
            wt=(c_av-vol_thr)/(c_av-c_avp1)
            clevel(ilevel)=(1d0-wt)*dble(iav)+wt*dble(iavp1)
            goto 51
         endif
         if(c_av.lt.vol_thr) iup=iav
         if(c_av.ge.vol_thr) idn=iav
         goto 50
51    continue
      write(7,fmt="(i6,f10.2,f15.9,f15.4)")
     +          ilevel,thres(ilevel),vfrac(ilevel),clevel(ilevel)
      write(*,fmt="(i6,f10.2,f15.9,f15.4)")
     +          ilevel,thres(ilevel),vfrac(ilevel),clevel(ilevel)
      enddo

c     print *,'Input min & max of X coord. in supermongo:'
c     read(5,*) smin,smax
c     print *,nint((smin-xmin)/dx),nint((smax-xmin)/dx)
c     print *,'Input min & max of Y coord. in supermongo:'
c     read(5,*) smin,smax
c     print *,nint((smin-ymin)/dy),nint((smax-ymin)/dy)

      do ilevel=1,3
         do ix=ipos_max,0,-1
            if(array1(ix).lt.clevel(ilevel)) then
               ileft=ix
               iright=ix+1
               cleft=array1(ix)
               cright=array1(ix+1)
               xleft(ilevel)=xmin+(dble(ix)+(clevel(ilevel)-cleft)
     &               /(cright-cleft))*dx
               goto 60
            endif
         enddo
60       continue
         do ix=ipos_max,nx,+1
            if(array1(ix).lt.clevel(ilevel)) then
               iright=ix
               ileft=ix-1
               cright=array1(ix)
               cleft=array1(ix-1)
               xright(ilevel)=xmin+(dble(ix)-(clevel(ilevel)-cright)
     &                /(cleft-cright))*dx
               goto 61
            endif
         enddo
61       continue
      enddo
      write(*,*) '-3s,-2s,-1s,mean,+1s,+2s,+3s'
      write(*,fmt='(30(1pe15.5))')
     &   xleft(3),xleft(2),xleft(1),xmean,
     &   xright(1),xright(2),xright(3)
      write(*,*) '-3s,-2s,-1s,mode,+1s,+2s,+3s'
      write(*,fmt='(30(1pe15.5))')
     &   xleft(3),xleft(2),xleft(1),amax_xcrd,
     &   xright(1),xright(2),xright(3)

      write(7,fmt='(30(1pe15.5))')
     &   xleft(3),xleft(2),xleft(1),xmean,
     &   xright(1),xright(2),xright(3)
      close(7)

c     print *,'Total elapsed time (sec): ',etime(cpu)-tzero

      end
c-----------------------------------------------------------------------------
      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      DOUBLE PRECISION c(np)
      PARAMETER (MMAX=6)
CU    USES lubksb,ludcmp
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      DOUBLE PRECISION d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+
     *1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m)
     *pause 'bad args in savgol'
      do 14 ipj=0,2*m
        sum=0.d0
        if(ipj.eq.0)sum=1.d0
        do 11 k=1,nr
          sum=sum+ dble(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+ dble(-k)**ipj
12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.d0
15    continue
      b(ld+1)=1.d0
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.d0
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.d0
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0
