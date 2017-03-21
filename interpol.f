!--------------------------------------------------------
! interpol.f: This subroutine is for interpolations of the given data
! using trigonometric interplation and Lagrange interpolation
!--------------------------------------------------------

!****************************************************
!     Trigonometric interpolation of y1 using Cosine functions
!     See Eq. (8.1) in Berrut and Trefethen, SIAM 46 (2004) 501
      subroutine IntTriCos(nin,x,y1,nout,
     1     xx,yy1,eps)

!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1  - values by interpolation of y1 in terms of x at xx

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,nin,nout,isame
 
      real *8  x(*), y1(*), xx(*), yy1(*), w(nin)
      real *8 ltemp, ltemp2

!     find the weights for Barycentric weight 
!     using nin-1 points (except x(nin) by assuming x(1)=x(nin))
!      call barycentric_wtriC(nin-1,x,w)
      call barycentric_wtriC(nin,x,w)
c      write(*,*) 'X',X(1:nin)
c      write(*,*) 'XX',XX(1:nout)
c      write(*,*) 'W',W(1:nin)
      do i=1,nout
         isame = 0
         j=0
         
         ! determine if there is a value in x to be equal to xx 
         do while ((j.le.nin) .and. (isame.ne.1))
            j=j+1
            if (dabs(x(j)-xx(i)) < eps) then
               yy1(i) = y1(j) ! If xx is same as x, yy=y at x
               isame = 1
               exit
            end if
         end do
  
         ! If xx is different from all points at x, interpolate them
         if (isame.ne.1) then
          
            ltemp2 = 0.0d0
            yy1(i) = 0.0d0
            do j=1, nin !j=1, nin-1
               ltemp=w(j)/(dcos(xx(i))-dcos(x(j)))
               ltemp2=ltemp2+ltemp
               yy1(i)=yy1(i)+y1(j)*ltemp 
            end do
            yy1(i)=yy1(i)/ltemp2      
         end if
      end do
c      write(*,*) 'YY',YY1(1:nout)
      
      return
      end 

!****************************************************
!     Trigonometric interpolation of y1 using Sine functions
      subroutine IntTriSin(nin,x,y1,nout,
     1     xx,yy1,eps)

!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1  - values by interpolation of y1 in terms of x at xx


      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,nin,nout,isame
 
      real *8  x(*), y1(*), xx(*), yy1(*), w(nin-1)
      real *8 ltemp, ltemp2

!     find the weights for Barycentric weight 
!     using nin-1 points (except x(nin) by assuming x(1)=x(nin))
      call barycentric_wtriS(nin-1,x,w)

      do i=1,nout
         isame = 0
         j=0
         do while ((j.lt.nin) .and. (isame.ne.1))
            j=j+1
            if (dabs(x(j)-xx(i)) < eps) then
               yy1(i) = y1(j)
               isame = 1
            end if
         end do         
         if (isame.ne.1) then
            ltemp2 = 0.0d0
            yy1(i) = 0.0d0
            do j=1, nin-1
               ltemp=w(j)/(dsin(xx(i))-dsin(x(j)))
               ltemp2=ltemp2+ltemp
               yy1(i)=yy1(i)+y1(j)*ltemp 
            end do
            yy1(i)=yy1(i)/ltemp2      
         end if
      end do

      return
      end 

!********************************************************
!     First derivative of Trigoonometric interpolated function
!     in terms of x at the all given x. Use cosine function
      subroutine IntTriDerCos(nin,x,y1,dyy1)

!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   - known values of the dependent variable
!     
!     output:
!     dyy1 - dy/dx evalauted by 
!            the trigonometric interpolation of y1

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,nin
      real *8  x(*), y1(*), w(nin-1)
      real *8  dyy1(*)
      real *8  ltemp, ltemp2
    

      dyy1(1:nin) = 0.0d0

!     find the weights for Barycentric weight 
!     using nin-1 points (except x(nin) by assuming x(1)=x(nin))
!      call barycentric_wtriC(nin-1,x,w)
      call barycentric_wtriC(nin,x,w)
      do i=1,nin
         ltemp2 = 0.0d0
         dyy1(i)= 0.0d0
         do j=1, nin
            if (j .ne. i) then 
               ltemp=w(j)/(dcos(x(i))-dcos(x(j)))
               ltemp2=ltemp2+ltemp
               dyy1(i)=dyy1(i)+y1(j)*ltemp 
            end if
         end do
         dyy1(i)=dyy1(i)-y1(i)*ltemp2       
         dyy1(i)=dyy1(i)*(-dsin(x(i)))/w(i)
      end do
      dyy1(nin)=dyy1(1)
      return
      end 

********************************************************
!     First derivative of Trigoonometric interpolated function
!     in terms of x at the all given x. Use sine function
      subroutine IntTriDerSin(nin,x,y1,dyy1)

!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   - known values of the dependent variable
!     
!     output:
!     dyy1 - dy/dx evalauted by 
!            the trigonometric interpolation of y1

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,nin
      real *8  x(*), y1(*), w(nin-1)
      real *8  dyy1(*)
      real *8  ltemp, ltemp2
    

      dyy1(1:nin) = 0.0d0

!     find the weights for Barycentric weight 
!     using nin-1 points (except x(nin) by assuming x(1)=x(nin))
      call barycentric_wtriS(nin-1,x,w)

      do i=1,nin-1
         ltemp2 = 0.0d0
         dyy1(i)= 0.0d0
         do j=1, nin-1
            if (j .ne. i) then
               ltemp=w(j)/(dsin(x(i))-dsin(x(j)))
               ltemp2=ltemp2+ltemp
               dyy1(i)=dyy1(i)+y1(j)*ltemp 
            end if
         end do
         dyy1(i)=dyy1(i)-y1(i)*ltemp2       
         dyy1(i)=dyy1(i)*(dcos(x(i)))/w(i)
      end do
      dyy1(nin)=dyy1(1)
      return
      end 


!****************************************************
!     Langrange interpolation of two periodic data y1 and y2
      subroutine IntLag(nin,x,y1,y2,norder,nout,
     1     xx,yy1,yy2,dyy1,dyy2,eps)

!     x and xx need to be increasing order
!     y1 and y2 are periodic (i.e. y1(1)=y1(nin) and y2(1)=y2(nin))
!     
!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable y1
!     y2   -  known values of the dependent variable y2
!     norder - order of Lagrange interpolation.
!              number of points used for the interpolation of a point
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1 - interpolation of y1 in terms of x at xx
!     yy2 - interpolation of y2 in terms of x at xx
!     dyy1 - interpolation of dy1/dx in terms of x at xx
!     dyy2 - interpolation of dy2/dx in terms of x at xx
      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin,ix,ij,nlorder, nrorder
      integer *4 nsame, ndiff, error
!      integer :: norder
      real *8  x(*), y1(*), y2(*), xx(*)

      parameter (nout2=1000000)
      real *8  xxs(nout2),yy1s(nout2),yy2s(nout2),dyy1s(nout2)
      real *8  dyy2s(nout2) ,dyy2d(nout2) 
      real *8  xxd(nout2),yy1d(nout2),yy2d(nout2),dyy1d(nout2)
      real *8  yy1(*), yy2(*), dyy1(*), dyy2(*)
      real *8 ltemp, ltemp2, ltemp3
      integer * 8 isame(nout2),jsame(nin),idiff(nout2),it


!     find if there is x which is same as xx
      call sort_samevalue(nout,xx,nin,x,nsame,isame
     1    ,jsame ,ndiff,idiff,eps)
c      write(*,*) 'nout2',nout
      write(*,*) 'nsame,ndiff', nsame, ndiff
c      write(*,*) 'isame,idiff', isame(1:nsame), idiff(1:ndiff)
c      write(*,*) 'xsame,xdiff', xx(isame(1:nsame)), xx(idiff(1:ndiff))
      if(nsame .ne. 0) then

      ! for xx which values are the same as a value in x

         xxs(1:nsame)=xx(isame(1:nsame))
         ! calculate dy/dx for the xx
         call IntLag_samepts2(nin,x,y1,y2, norder, nsame,
     1        xxs, yy1s, yy2s, dyy1s,dyy2s,eps)
         yy1(isame(1:nsame)) = yy1s(1:nsame)
         yy2(isame(1:nsame)) = yy2s(1:nsame)
         dyy1(isame(1:nsame)) = dyy1s(1:nsame)
         dyy2(isame(1:nsame)) = dyy2s(1:nsame)


      end if

      if(ndiff .ne. 0) then  

      ! for xx which values are different from the values in x
         xxd(1:ndiff)=xx(idiff(1:ndiff))
         call IntLag_diffpts2(nin,x,y1,y2, norder,ndiff,
     1        xxd, yy1d, yy2d, dyy1d,dyy2d)
         yy1(idiff(1:ndiff)) = yy1d(1:ndiff)
         yy2(idiff(1:ndiff)) = yy2d(1:ndiff)
         dyy1(idiff(1:ndiff)) = dyy1d(1:ndiff)
         dyy2(idiff(1:ndiff)) = dyy2d(1:ndiff)
      end if
      return

      end 

!****************************************************
!     Langrange interpolation of the periodic data for xx 
!     , which values are different from the values in x
      subroutine IntLag_diffpts2(nin,x,y1,y2,norder,nout,
     1     xx,yy1,yy2,dyy1,dyy2)

!     x and xx need to be increasing order
!     y1 and y2 are periodic (i.e. y1(1)=y1(nin) and y2(1)=y2(nin))
!     
!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable y1
!     y2   -  known values of the dependent variable y2
!     norder - order of Lagrange interpolation.
!              number of points used for the interpolation of a point
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1 - interpolation of y1 in terms of x at xx
!     yy2 - interpolation of y2 in terms of x at xx
!     dyy1 - interpolation of dy1/dx in terms of x at xx
!     dyy2 - interpolation of dy2/dx in terms of x at xx

      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder, nout
      real *8 x(*), y1(*), y2(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*), yy2(*), dyy1(*), dyy2(*)
      real *8  xtdiff(norder+nin) 
      real *8  ytdiff1(norder+nin), ytdiff2(norder+nin) 
      real *8 ltemp, ltemp2, ltemp3
      
!      nin = size(x)
!      nout = size(xx)

      ! use norder values for interpolation 
      ! (nlorder values in x, which are less than xx
      !   + nrorder values in x, which are greater than xx)
      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder
 
      ! y is assumed to be periodic y(1)=y(nin)
      ! expand y1 and y2 for first few points and last few points
      ! to constuct periodic data

      ytdiff1(1:nlorder)=y1(nin-nlorder:nin-1)
      ytdiff2(1:nlorder)=y2(nin-nlorder:nin-1)
      xtdiff(1:nlorder)=x(nin-nlorder:nin-1)-x(nin)
    
      ytdiff1(nlorder+1:nin+nlorder)=y1(1:nin)
      ytdiff2(nlorder+1:nin+nlorder)=y2(1:nin)
      xtdiff(nlorder+1:nin+nlorder)=x(1:nin)
   
      ytdiff1(nin+nlorder+1:)=y1(2:nrorder+1)
      ytdiff2(nin+nlorder+1:)=y2(2:nrorder+1)
      xtdiff(nin+nlorder+1:)=x(nin)+x(2:nrorder+1)

      yy1(1:nout) = 0.0
      dyy1(1:nout) = 0.0
      yy2(1:nout) = 0.0
      dyy2(1:nout) = 0.0

!
      ix = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while ((x(ix+1)<xx(i)) .and. ix<nin)  
            ix=ix+1
         end do

!     find the weights for Barycentric weight 
!     using norder points around xx 
         call barycentric_w(norder,xtdiff(ix:ix+norder-1),w)
         ltemp=0.0
         ltemp2=0.0
 
         do j=1, norder
            ij=ix+j-1
            yy1(i)=yy1(i)+ytdiff1(ij)*w(j)/(xx(i)-xtdiff(ij))
            dyy1(i)=dyy1(i)+ytdiff1(ij)*w(j)/(xx(i)-xtdiff(ij))**2
            yy2(i)=yy2(i)+ytdiff2(ij)*w(j)/(xx(i)-xtdiff(ij))
            dyy2(i)=dyy2(i)+ytdiff2(ij)*w(j)/(xx(i)-xtdiff(ij))**2

            ltemp=ltemp+w(j)/(xx(i)-xtdiff(ij))
            ltemp2=ltemp2+w(j)/(xx(i)-xtdiff(ij))**2
         end do
         dyy1(i)=(ltemp2/ltemp*yy1(i)-dyy1(i))/ltemp
         yy1(i)=yy1(i)/ltemp
         dyy2(i)=(ltemp2/ltemp*yy2(i)-dyy2(i))/ltemp
         yy2(i)=yy2(i)/ltemp
      end do

      end 


!****************************************************
!     derivative of the periodic data for xx
!     , which values are the same as the values in x.
!     See Lagrange formula given in Berrut and Trefethen, SIAM Review 2004

      subroutine IntLag_samepts2(nin, x,y1,y2,norder,nout,
     1     xx,yy1, yy2, dyy1,dyy2,eps)

!     x and xx need to be increasing order
!     y1 and y2 are periodic (i.e. y1(1)=y1(nin) and y2(1)=y2(nin))
!     
!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable y1
!     y2   -  known values of the dependent variable y2
!     norder - order of Lagrange interpolation.
!              number of points used for the interpolation of a point
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1 - interpolation of y1 in terms of x at xx
!     yy2 - interpolation of y2 in terms of x at xx
!     dyy1 - interpolation of dy1/dx in terms of x at xx
!     dyy2 - interpolation of dy2/dx in terms of x at xx


      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder
      real *8 x(*), y1(*), y2(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*), yy2(*), dyy1(*), dyy2(*)
      real *8  xtsame(norder+nin) 
      real *8  ytsame1(norder+nin), ytsame2(norder+nin) 
      real *8  ltemp, ltemp2, ltemp3
      

      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder

      ytsame1(1:nlorder)=y1(nin-nlorder:nin-1)
      ytsame2(1:nlorder)=y2(nin-nlorder:nin-1)
      xtsame(1:nlorder)=x(nin-nlorder:nin-1)-x(nin)
    
      ytsame1(nlorder+1:nin+nlorder)=y1(1:nin)
      ytsame2(nlorder+1:nin+nlorder)=y2(1:nin)
      xtsame(nlorder+1:nin+nlorder)=x(1:nin)
   
      ytsame1(nin+nlorder+1:)=y1(2:nrorder+1)
      ytsame2(nin+nlorder+1:)=y2(2:nrorder+1)
      xtsame(nin+nlorder+1:)=x(nin)+x(2:nrorder+1)

      dyy1(1:nout) = 0.0
      dyy2(1:nout) = 0.0

      ix = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while (abs(x(ix)-xx(i)) > eps)  
            ix=ix+1
         end do
         yy1(i) = y1(ix)
         yy2(i) = y2(ix)
         call barycentric_w(norder,xtsame(ix:ix+norder-1),w)
         ltemp2 = 0.0
         do j=1, norder
            if (j .ne. nlorder+1) then 
               ij=ix+j-1
               ltemp=w(j)/(xx(i)-xtsame(ij))
               ltemp2=ltemp2+ltemp
               dyy1(i)=dyy1(i)+ytsame1(ij)*ltemp 
               dyy2(i)=dyy2(i)+ytsame2(ij)*ltemp
            end if
         end do
         dyy1(i)=dyy1(i)-ytsame1(ix+nlorder)*ltemp2      
         dyy1(i)=dyy1(i)/w(nlorder+1)
         dyy2(i)=dyy2(i)-ytsame2(ix+nlorder)*ltemp2      
         dyy2(i)=dyy2(i)/w(nlorder+1)
      end do

      end 

!****************************************************
!     Langrange interpolation using all given data
      subroutine IntLag_sameptsAll(nin,x,y,
     1     xx,yy,dyy,eps)

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin, isame
      real *8  x(*), y(*), xx, yy, dyy, eps, w(nin)
      real *8 ltemp2,ltemp

      isame=1
      ix=1
      do while (dabs(x(ix)-xx) > eps)  
         ix=ix+1
         if(ix.gt.nin) then
            isame=0
            exit
         end if
      end do 

      if (isame.ne.0) then

         yy = y(ix)
         call barycentric_w(nin,x,w)
         ltemp2 = 0.0d0
         do j=1, nin
            if (j .ne. ix) then 
              
               ltemp=w(j)/(xx-x(j))
               ltemp2=ltemp2+ltemp
               dyy =dyy + y(j)*ltemp 
            end if
         end do
         dyy=dyy-yy*ltemp2      
         dyy=dyy/w(ix)
      end if 
      write(*,*) 'yy,dyy',yy,dyy
      return

      end 

!****************************************************
!     check if there is a integer value in x, which is the same as xx
      subroutine CheckSameInt(nin,x,xx,isame)

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin, isame
      integer *4  x(*), xx

      isame=0
      do j=1, nin
         if (x(j).eq.xx) then
            isame=j
            return
         end if
      end do
      return

      end 
!****************************************************
!     check if there is a real value in x, which is the same as xx
!     If abs(xx-x) < eps, they are considered to be the same.
      subroutine CheckSame(nin,x,xx,eps,isame)

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin, isame
      real *8  x(*), xx

      isame=0
      do j=1, nin
         if (dabs(x(j)-xx) < eps) then
            isame=j
            return
         end if
      end do
      return

      end 

!****************************************************
!     Langrange interpolation using all given data
      subroutine IntLagAll(nin,x,y,
     1     xx,yy,w,eps,isame)

      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin, isame
      real *8  x(*), y(*), xx, yy, eps, w(*)

      do j=1, nin
         if (dabs(x(j)-xx) < eps) then
            yy=y(j)
            isame=j
            return
         end if
      end do

      call barycentric_w(nin,x,w)
      temp=0.0d0
      yy=0.0d0
      do j=1, nin
         yy=yy + y(j)*w(j)/(xx-x(j))
         temp=temp+w(j)/(xx-x(j))          
      end do
      yy = yy/temp
      isame=0
      return

      end 

!****************************************************
!     Langrange interpolation using all given data. 
!     Only last data points in x and y are new, 
!     so modify the existing weight w for the new data

      subroutine IntLagAdd(nin,x,y,
     1     xx,yy,w,eps,isame)

!     x(1:nin-1) and y(1:nin-1) are used for the weight w(1:nin-1)
!     x(nin) and y(nin) are added.

!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!           
!     y    -  known values of the dependent variable
!     w    -  w(1:nin-1) is precomputed weight using x(1:nin-1) and y(1:nin-1)
!     xx   - value of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     isame - index of x if there is a point in x which is the same as xx.
!             Otherwise, isame= 0
!     yy    - value by interpolation of y1 in terms of x at xx
      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin, isame
      real *8  x(*), y(*), w(*),xx, yy, eps

      do j=1, nin
         if (dabs(x(j)-xx) < eps) then
            yy=y(j)
            isame=j
            return
         end if
      end do

      call barycentric_wAdd(nin,x,w)
      temp=0.0d0
      yy=0.0d0
      do j=1, nin
         yy=yy + y(j)*w(j)/(xx-x(j))
         temp=temp+w(j)/(xx-x(j))          
      end do
      yy = yy/temp
      isame=0
      return

      end 

!****************************************************
!     Langrange interpolation of a non-periodic data y1
      subroutine IntLag1np(nin,x,y1,norder,nout,
     1     xx,yy1,dyy1,eps)

!     x and xx need to be increasing order
!     
!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable y1
!     norder - order of Lagrange interpolation.
!              number of points used for the interpolation of a point
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1 - interpolation of y1 in terms of x at xx
!     dyy1 - interpolation of dy1/dx in terms of x at xx


      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin,ix,ij,nlorder, nrorder
      integer *4 nsame, ndiff, error

      real *8  x(*), y1(*), xx(*)

      parameter (nout2=1000000)
      real *8  xxs(nout2),yy1s(nout2),dyy1s(nout2)
      real *8  xxd(nout2),yy1d(nout2),dyy1d(nout2)
      real *8  yy1(*), dyy1(*)
      real *8 ltemp, ltemp2, ltemp3
      integer * 8 isame(nout2),jsame(nin),idiff(nout2),it


      call sort_samevalue(nout,xx,nin,x,nsame,isame
     1     ,jsame,ndiff,idiff,eps)
c      write(*,*) 'nsame,ndiff', nsame, ndiff
      if(nsame .ne. 0) then

         xxs(1:nsame)=xx(isame(1:nsame))
         call IntLag_samepts1np(nin,x,y1, norder, nsame,
     1        xxs, yy1s,  dyy1s,eps)
         yy1(isame(1:nsame)) = yy1s(1:nsame)
         dyy1(isame(1:nsame)) = dyy1s(1:nsame)

      end if
      if(ndiff .ne. 0) then  

         xxd(1:ndiff)=xx(idiff(1:ndiff))
!
         call IntLag_diffpts1np(nin,x,y1, norder,ndiff,
     1        xxd, yy1d,  dyy1d)
         yy1(idiff(1:ndiff)) = yy1d(1:ndiff)
         dyy1(idiff(1:ndiff)) = dyy1d(1:ndiff)
      end if
      return

      end 
!****************************************************
!     Langrange interpolation of the non-periodic data
      subroutine IntLag_diffpts1np(nin,x,y1,norder,nout,
     1     xx,yy1,dyy1)

      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder, nout,il,ir
      real *8 x(*), y1(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*), dyy1(*)
      real *8 ltemp, ltemp2, ltemp3
      
      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder
 
 
      yy1(1:nout) = 0.0
      ix = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while ((x(ix)<xx(i)) .and. ix<nin)  
            ix=ix+1
         end do

         il=ix-nlorder
         ir=ix+nrorder
         if (il<1) then
            il=1
            ir=norder
         else if (ir>nin) then
            ir=nin
            il=nin-norder+1
         end if

         call barycentric_w(norder,x(il:ir),w)
         ltemp=0.0
         ltemp2=0.0
 
         do j=1, norder
            ij=il+j-1
          
            yy1(i)=yy1(i)+y1(ij)*w(j)/(xx(i)-x(ij))
            dyy1(i)=dyy1(i)+y1(ij)*w(j)/(xx(i)-x(ij))**2
 
            ltemp=ltemp+w(j)/(xx(i)-x(ij))
            ltemp2=ltemp2+w(j)/(xx(i)-x(ij))**2
           
         end do
         dyy1(i)=(ltemp2/ltemp*yy1(i)-dyy1(i))/ltemp
         yy1(i)=yy1(i)/ltemp
      end do

      end 
!****************************************************
!     derivative of the non periodic data 
!     using Lagrange formula given in Berrut and Trefethen, SIAM Review 2004
      subroutine IntLag_samepts1np(nin, x,y1,norder,nout,
     1     xx,yy1,  dyy1,eps)
      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder
      real *8 x(*), y1(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*),  dyy1(*)
      real *8  ltemp, ltemp2, ltemp3
      
      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder

      dyy1(1:nout) = 0.0

      ix = 1
      isame = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while (abs(x(ix)-xx(i)) > eps)  
            ix=ix+1
            if(ix.gt.nin) then
               isame=0
               exit
            end if
         end do 

         if (isame.eq.1) then
  
            il=ix-nlorder
            ir=ix+nrorder
       
            if (il<1) then
               il=1
               ir=norder
               ixx=ix
            else if (ir>nin) then
               ir=nin
               il=nin-norder+1
               ixx=ix-il+1
            else 
               ixx=nlorder+1
            end if
c            write(*,*) 'il,ir,ix,ixx',il,ir,ix,ixx
            yy1(i) = y1(ix)

            call barycentric_w(norder,x(il:ir),w)
            ltemp2 = 0.0
            do j=1, norder
               ij=il+j-1
               if (ij .ne. ix) then 
                  
                  ltemp=w(j)/(xx(i)-x(ij))
                  ltemp2=ltemp2+ltemp
                  dyy1(i)=dyy1(i)+y1(ij)*ltemp 
               end if
            end do
            dyy1(i)=dyy1(i)-y1(ix)*ltemp2      
            dyy1(i)=dyy1(i)/w(ixx)
         end if
      end do

      return
      end 

!****************************************************
!     Langrange interpolation of a periodic data y1
      subroutine IntLag1(nin,x,y1,norder,nout,
     1     xx,yy1,dyy1,eps)

!     x and xx need to be increasing order
!     y1 and y2 are periodic (i.e. y1(1)=y1(nin) and y2(1)=y2(nin))
!     
!     input:
!     nin  - number of known values used for interpolation
!     x    - known values of the indepedent variable
!     y1   -  known values of the dependent variable y1
!     norder - order of Lagrange interpolation.
!              number of points used for the interpolation of a point
!     nout - number of points to be interpolated
!     xx   - values of the independent variable to be interpolated 
!            (min(x) <= xx <= max(x))
!     eps  - small value used to determine if xx is equal to be x

!     output:
!     yy1 - interpolation of y1 in terms of x at xx
!     dyy1 - interpolation of dy1/dx in terms of x at xx


      implicit real *8 (a-h,o-z)
      save
      integer *4 i,j,k,nin,ix,ij,nlorder, nrorder
      integer *4 nsame, ndiff, error
      real *8  x(*), y1(*), xx(*)

      parameter (nout2=1000000)
      real *8  xxs(nout2),yy1s(nout2),dyy1s(nout2)
      real *8  xxd(nout2),yy1d(nout2),dyy1d(nout2)
      real *8  yy1(*), dyy1(*)
      real *8 ltemp, ltemp2, ltemp3
      integer * 8 isame(nout2),jsame(nin),idiff(nout2),it


      call sort_samevalue(nout,xx,nin,x,nsame,isame
     1     ,jsame,ndiff,idiff,eps)

      if(nsame .ne. 0) then

         xxs(1:nsame)=xx(isame(1:nsame))
         call IntLag_samepts1(nin,x,y1, norder, nsame,
     1        xxs, yy1s,  dyy1s,eps)
         yy1(isame(1:nsame)) = yy1s(1:nsame)
         dyy1(isame(1:nsame)) = dyy1s(1:nsame)

      end if
      if(ndiff .ne. 0) then  

         xxd(1:ndiff)=xx(idiff(1:ndiff))
!
         call IntLag_diffpts1(nin,x,y1, norder,ndiff,
     1        xxd, yy1d, dyy1d)
         yy1(idiff(1:ndiff)) = yy1d(1:ndiff)
         dyy1(idiff(1:ndiff)) = dyy1d(1:ndiff)

      end if
      return

      end 

!****************************************************
!     Langrange interpolation of the periodic data
      subroutine IntLag_diffpts1(nin,x,y1,norder,nout,
     1     xx,yy1,dyy1)

      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder, nout
      real *8 x(*), y1(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*),dyy1(*)
      real *8  xtdiff(norder+nin) 
      real *8  ytdiff1(norder+nin)
      real *8 ltemp, ltemp2, ltemp3
      
      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder
 
      ! y is assumed to be periodic y(1)=y(nin)

      ytdiff1(1:nlorder)=y1(nin-nlorder:nin-1)
      xtdiff(1:nlorder)=x(nin-nlorder:nin-1)-x(nin)
    
      ytdiff1(nlorder+1:nin+nlorder)=y1(1:nin)
      xtdiff(nlorder+1:nin+nlorder)=x(1:nin)
   
      ytdiff1(nin+nlorder+1:)=y1(2:nrorder+1)
      xtdiff(nin+nlorder+1:)=x(nin)+x(2:nrorder+1)

      yy1(1:nout) = 0.0
      dyy1(1:nout) = 0.0
!
      ix = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while ((x(ix+1)<xx(i)) .and. ix<nin)  
            ix=ix+1
         end do

         call barycentric_w(norder,xtdiff(ix:ix+norder-1),w)
         ltemp=0.0
         ltemp2=0.0
 
         do j=1, norder
            ij=ix+j-1
          
            yy1(i)=yy1(i)+ytdiff1(ij)*w(j)/(xx(i)-xtdiff(ij))
            dyy1(i)=dyy1(i)+ytdiff1(ij)*w(j)/(xx(i)-xtdiff(ij))**2
 
            ltemp=ltemp+w(j)/(xx(i)-xtdiff(ij))
            ltemp2=ltemp2+w(j)/(xx(i)-xtdiff(ij))**2
           
         end do
         dyy1(i)=(ltemp2/ltemp*yy1(i)-dyy1(i))/ltemp
         yy1(i)=yy1(i)/ltemp
      end do

      end 


!****************************************************
!     derivative of the periodic data 
!     using Lagrange formula given in Berrut and Trefethen, SIAM Review 2004
      subroutine IntLag_samepts1(nin, x,y1,norder,nout,
     1     xx,yy1,  dyy1,eps)
      implicit real *8 (a-h,o-z)
      save
      integer i,j,k,ix,ij,nlorder, nrorder
      real *8 x(*), y1(*), xx(*)
      real *8  w(norder)
      real *8  yy1(*),  dyy1(*)
      real *8  xtsame(norder+nin) 
      real *8  ytsame1(norder+nin)
      real *8  ltemp, ltemp2, ltemp3
      
      nlorder = floor(norder/2.0)
      nrorder = norder-nlorder

      ytsame1(1:nlorder)=y1(nin-nlorder:nin-1)
      xtsame(1:nlorder)=x(nin-nlorder:nin-1)-x(nin)
    
      ytsame1(nlorder+1:nin+nlorder)=y1(1:nin)
      xtsame(nlorder+1:nin+nlorder)=x(1:nin)
   
      ytsame1(nin+nlorder+1:)=y1(2:nrorder+1)
      xtsame(nin+nlorder+1:)=x(nin)+x(2:nrorder+1)
      dyy1(1:nout) = 0.0

      ix = 1
      isame = 1
      do i=1,nout
         !find the index of the given data which is close to the target 
         do while (abs(x(ix)-xx(i)) > eps)  
            ix=ix+1
            if(ix.gt.nin) then
               isame=0
               exit
            end if
         end do 

         if (isame.eq.1) then

         yy1(i) = y1(ix)
         call barycentric_w(norder,xtsame(ix:ix+norder-1),w)
         ltemp2 = 0.0
         do j=1, norder
            if (j .ne. nlorder+1) then 
               ij=ix+j-1
               ltemp=w(j)/(xx(i)-xtsame(ij))
               ltemp2=ltemp2+ltemp
               dyy1(i)=dyy1(i)+ytsame1(ij)*ltemp 
            end if
         end do
         dyy1(i)=dyy1(i)-ytsame1(ix+nlorder)*ltemp2      
         dyy1(i)=dyy1(i)/w(nlorder+1)
         end if
      end do

      return
      end 

!****************************************************
!     find the weights for Barycentric weight of Langrange interpolation
      subroutine barycentric_w(n,x,w)

      implicit real *8 (a-h,o-z)

      integer j,k,n
      real *8 x(*), w(*)
c      write(*,*) 'x', x(1:n)
      w(1:n) = 1.0d0
      do j=2,n 
         do k=1,j-1
            w(k)=(x(k)-x(j))*w(k)
            w(j)=(x(j)-x(k))*w(j)
         end do
      end do
      do j=1,n
         w(j)=1.0d0/w(j)
      end do
c      write(*,*) 'bary_w', w(1:n)
      return
      end 

!****************************************************
!     Modify the weights for Barycentric weight of Langrange interpolation
!     using a new point x(n). 
!     w(1:n-1) was precomputed using x(1:n-1)

      subroutine barycentric_wAdd(n,x,w)

      implicit real *8 (a-h,o-z)

      integer j,k,n
      real *8 x(*), w(*)
      
!      n=size(x)-1
      w(n)=1.0d0
      do j=1,n-1 
         dx=x(j)-x(n)
         w(j)=w(j)/dx
         w(n)=-w(n)/dx
      end do
      
      return
      end 

!****************************************************
!     find the weights for Barycentric weight of trigonometric interpolation
!     Cosine function is used
      subroutine barycentric_wtriC(n,x,w)

      implicit real *8 (a-h,o-z)

      integer j,k,n
      real *8 x(*), w(*)
      w(1:n) = 1.0d0
      do j=2,n 
         do k=1,j-1
            w(k)=(dcos(x(k))-dcos(x(j)))*w(k)
            w(j)=(dcos(x(j))-dcos(x(k)))*w(j)
         end do
      end do
      do j=1,n
         w(j)=1.0d0/w(j)
      end do


      return
      end 
!****************************************************
!     find the weights for Barycentric weight of trigonometric interpolation
!     Sine function is used
      subroutine barycentric_wtriS(n,x,w)

      implicit real *8 (a-h,o-z)

      integer j,k,n
      real *8 x(*), w(*)
      w(1:n) = 1.0d0
      do j=2,n 
         do k=1,j-1
            w(k)=(dsin(x(k))-dsin(x(j)))*w(k)
            w(j)=(dsin(x(j))-dsin(x(k)))*w(j)
         end do
      end do
      do j=1,n
         w(j)=1.0d0/w(j)
      end do

      return
      end 

!****************************************************
!     find the indexes of x1 which values are the same as the values of x2
!     within the error of eps
      subroutine sort_samevalue(n1,x1,n2,x2,
     1     nsame,isame,jsame, ndiff, idiff, eps)
!           
!     x1 is to be sorted
!     x2 is given data to be compared  
      implicit real *8 (a-h,o-z)
      save
      integer i, j,k,n1,n2
!      real *8 eps
      real *8 x1(*), x2(*)
      integer *8 isame(*),jsame(*), idiff(*)


      nsame = 0
      ndiff = 0
      do i=1,n1
         do j=1, n2
            if (dabs(x1(i)-x2(j)) < eps) then
               nsame = nsame +1
               isame(nsame)=i 
               jsame(nsame)=j 
               exit
            else if (j.eq.n2) then
               ndiff = ndiff +1
               idiff(ndiff)=i
            end if
         end do
      end do

      return
      end 

!****************************************************
!     test of the Lagrange interpolation
      subroutine testintlag(nin,nout,norder,eps)

      implicit real*8 (a-h,o-z)
      save
      real *8 xin(nin),y1(nin), y2(nin)
      real *8 xtest(nout),yo1(nout), yex1(nout),yo2(nout),yex2(nout)
      real *8 dyo1(nout), dyex1(nout),dyo2(nout),dyex2(nout)
      integer i

      pi = 4.0d0*datan(1.0d0)
      do i=1,nin
         xin(i)=(i-1)*2*pi/(nin-1)
         y1(i) = sin(xin(i))
         y2(i) = sin(xin(i))**3
      end do
c
      do i=1,nout
         xtest(i)=(i-1)*2*pi/(nout-1)
         yex1(i) = sin(xtest(i))
         yex2(i) = sin(xtest(i))**3
         dyex1(i) = cos(xtest(i))
         dyex2(i) = 3*sin(xtest(i))**2*cos(xtest(i))
      end do

      call IntLag(nin,xin,y1,y2,norder,nout,
     1     xtest,yo1,yo2,dyo1,dyo2,eps)


c
c     check error
c
      err1 = 0.0d0
      errmax1 = 0.0d0
      errdy1 = 0.0d0
      sol1 = 0.0d0
      soldy1 = 0.0d0
      err2 = 0.0d0
      errmax2 = 0.0d0
      errdy2 = 0.0d0
      sol2 = 0.0d0
      soldy2 = 0.0d0



      do i = 1,nout

         errtmp = abs(yo1(i)-yex1(i))
         errdytmp = abs(dyo1(i)-dyex1(i))
          err1 = err1 + errtmp**2
         if (errmax1.lt.errtmp) then
            errmax1= errtmp
         end if
         errdy1 = errdy1 +errdytmp**2 
 
         sol1 = sol1 + (yex1(i))**2
         soldy1 = soldy1 + (dyex1(i))**2

         errtmp = abs(yo2(i)-yex2(i))
         errdytmp = abs(dyo2(i)-dyex2(i))
          err2 = err2 + errtmp**2
         if (errmax2.lt.errtmp) then
            errmax2= errtmp
         end if
         errdy2 = errdy2 +errdytmp**2 
 
         sol2 = sol2 + (yex2(i))**2
         soldy2 = soldy2 + (dyex2(i))**2
   
      enddo
     
      relerr1 = sqrt(err1/sol1)
      relerrdy1 = sqrt(errdy1/soldy1)
      relerr2 = sqrt(err2/sol2)
      relerrdy2 = sqrt(errdy2/soldy2)


c      write (*,*) 'error in Lagrange interpolation order of ', norder
c      write (*,*) 'err1, errmax1', relerr1, errmax1
c      write (*,*) 'err2, errmax2', relerr2, errmax2
c      write (*,*) 'err1dy, err2dy', relerrdy1, relerrdy2

      return
      end

      recursive subroutine QSortdec(nA, A, iA)
 
      real *8 A(nA)
      integer *4 iA(nA)

      integer *4 left, right, marker, itemp
      real*8 random, pivot, temp

      save
c      write(*,*) 'nA,A',nA!,A
      if (nA .gt. 1) then
 
        call random_number(random)
     
        pivot = A(int(random*real(nA-1))+1)   ! random pivor (not best performance, but avoids worst-case)
c        write(*,*) 'pivot',pivot,random
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right) < pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left) > pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
                itemp = iA(left)
                iA(left) = iA(right)
                iA(right) = itemp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QSortdec(marker-1,A(1:marker-1),iA(1:marker-1))
        call QSortdec(nA-marker+1,A(marker:nA),iA(marker:nA))
 
      end if
 
      end subroutine QSortdec

      recursive subroutine QSortinc(nA, A, iA)
 
      real *8 A(nA)
      integer *4 iA(nA)

      integer *4 left, right, marker, itemp, indx
      real*8 random, pivot, temp

      save
c      write(*,*) 'nA,A',nA!,A
      if (nA .gt. 1) then
 
        call random_number(random)
        indx=int(random*real(nA-1))+1
        pivot = A(indx)   ! random pivor (not best performance, but avoids worst-case)
c        write(*,*) 'pivot',pivot,indx
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left) < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
                itemp = iA(left)
                iA(left) = iA(right)
                iA(right) = itemp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QSortinc(marker-1,A(1:marker-1),iA(1:marker-1))
        call QSortinc(nA-marker+1,A(marker:nA),iA(marker:nA))
 
      end if
 
      end subroutine QSortinc


      subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
      implicit none
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer i, j, gap
      double precision h

      gap = n-1
!     check input
      if ( n < 2 ) return
      if ( n < 3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1)) ! linear interpolation
         c(1) = 0.
         d(1) = 0.
         b(2) = b(1)
         c(2) = 0.
         d(2) = 0.
         return
      end if
!
! step 1: preparation
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
         d(i) = x(i+1) - x(i)
         b(i) = 2.0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
!
! step 2: end conditions 
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if(n /= 3) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
!
! step 3: forward elimination 
!
      do i = 2, n
         h = d(i-1)/b(i-1)
         b(i) = b(i) - h*d(i-1)
         c(i) = c(i) - h*c(i-1)
      end do
!
! step 4: back substitution
!
      c(n) = c(n)/b(n)
      do j = 1, gap
         i = n-j
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!
! step 5: compute spline coefficients
!
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
      do i = 1, gap
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      end subroutine spline

      subroutine ispline(u, x, y, b, c, d, n, uy)
!======================================================================
! subroutine ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! uy = interpolated value at point u
!=======================================================================
      implicit none
      integer n
      double precision  u, uy, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      double precision dx

      save

! if u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
         dx = u - x(1)
         uy = y(1) + dx*(b(1) + dx*(c(1) + dx*d(1)))
!         uy = y(1)
         return
      end if
      if(u >= x(n)) then
         dx = u - x(n)
         uy = y(n) + dx*(b(n) + dx*(c(n) + dx*d(n)))
!         uy = y(n)
         return
      end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
      i = 1
      j = n+1
      do while (j > i+1)
         k = (i+j)/2
         if(u < x(k)) then
            j=k
         else
            i=k
         end if
      end do
!*
!  evaluate spline interpolation
!*
      dx = u - x(i)
      uy = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      end subroutine ispline
