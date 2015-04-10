!--------------------------------------------------------
! timefile.f: This subroutine is to measure the elapsed run time and 
!--------------------------------------------------------

      subroutine runtime(sloc)

      use arrays, only: ici, icr, ic_max, ftime
      
      implicit double precision(a-h, o-z)
      character (20) sloc

      call system_clock(ict, icr, ic_max)
      if (ict.gt.ici) then 
         tsec=real(ict-ici)/real(icr)
      else
         tsec=real(ict-ici+ic_max)/real(icr) !Assumed one trip of clock
      end if

      write(ftime,*) 'CPU time until '//sloc//' : t(sec)= ',tsec

      return
      end
     
      subroutine initfiles

      use arrays, only: nt2, kcheb, nsub, iecom, itrinity
      use arrays, only: fiter, ftime, fout
      use arrays, only: fiter1, fiter2, fiter3,fiter4,fiter5
      use arrays, only: fiter6,fiter7,fiter22
  
      character(3) ko, ns
      character(4) nt

      ftime =18
      fiter =19
      fout  =20
      fiter1=11;fiter2=12;fiter3=13;fiter4=14
      fiter5=15;fiter6=16;fiter22=22
      write (nt,'(i4.4)') nt2
      write (ko,'(i3.3)') kcheb
      write (ns,'(i3.3)') nsub
      open(ftime, file='time_nt'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out')
      open(fiter, file='iter_nt'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out')
      open(fout, file='ecom.out')

      if (iecom.ne.0) then
         open(fiter1,file='fpol1.out')
         open(fiter2,file='fpol2.out')
         open(fiter3,file='fpol3.out' )
         open(fiter4,file='fpol4.out' )
         open(fiter5,file='fpol5.out')
         open(fiter6,file='fpol6.out')
      end if

      if (itrinity.ne.0) then
         open(fiter22,file='trinity.out')
      end if

      return
      end

      subroutine closefiles

      use arrays, only: fiter, ftime, fout
      use arrays, only: fiter1, fiter2, fiter3,fiter4,fiter5
      use arrays, only: fiter6,fiter7,fiter22
     
      close(ftime)
      close(fiter)
      close(fout)
      if (iecom.ne.0) then
         close(fiter1)
         close(fiter2)
         close(fiter3)
         close(fiter4)
         close(fiter5)
         close(fiter6)
      end if

      if (itrinity.ne.0) then
         close(fiter22)
      end if
      return
      end
