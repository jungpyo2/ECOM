!--------------------------------------------------------
! quaplot.f: This is a plotting routine for conformal mapping 
!--------------------------------------------------------
c
c
c
c
c
        subroutine quaplot(iw,x,y,n,itype1,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=1
c
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
        subroutine quaplot2(iw,x,y,n,itype1,x2,y2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=2
c
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
        subroutine quaplot3(iw,x,y,n,itype1,x2,y2,n2,itype2,
     1      x3,y3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=3
c
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
        subroutine quagraph(iw,x,y,n,itype1,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=1
c
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
        subroutine quagraph2(iw,x,y,n,itype1,x2,y2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=2
c
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
c
        subroutine quagraph3(iw,x,y,n,itype1,x2,y2,n2,itype2,
     1      x3,y3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=3
c
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c
        return
        end
c
c
c
c
c
        subroutine quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c
c        convert the user-specified Fortran unit number to 
c        character format
c
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c
 2000 format(1a8)
        read(dummy,2000) anum8
c
c        construct the file name on which the Gnuplot instructions
c        are to be written
c
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c
        write(iun,2250)
c
c        generate the title for the plot
c
        line(1)=blank
        line(2)=blank
        line(3)=quo
c
        call quamesslen(title,nchar,line(4))
c
        line(nchar+4)=quo
 2300 format('  set title ',80a1)
        write(iun,2300) (line(i),i=1,nchar+4)
c
 2350 format('   show title')
        write(iun,2350)
c
c        write the instructions 
c
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c
 2800 format('plot ','"',1a8,'"     ','notitle  with dots')
 2803 format('plot ','"',1a8,'"     ','notitle  with points')
 2805 format('plot ','"',1a8,'"     ','notitle  with lines')
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) file8
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) file8
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) file8

 2830 format('plot ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2831 format('plot ','"',1a8,'"     ','notitle  with points, ',1a1)
 2832 format('plot ','"',1a8,'"     ','notitle  with lines, ',1a1)
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) file8,backslash
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1       write(iun,2831) file8,backslash
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2832) file8,backslash
c
c        write the first data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x(i),y(i),i=1,n)        
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the second data file
c
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return

        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c
 2850 format('     ','"',1a8,'"     ','notitle  with dots')
 2851 format('     ','"',1a8,'"     ','notitle  with points')
 2852 format('     ','"',1a8,'"     ','notitle  with lines')
c
        if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1      write(iun,2850) file8
c
        if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1      write(iun,2851) file8
c
        if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1      write(iun,2852) file8
c
c
        if( (inumgr .eq. 3) .and. (itype2 .eq. 1) )
     1      write(iun,2855) file8,backslash
c
        if( (inumgr .eq. 3) .and. (itype2 .eq. 2) )
     1      write(iun,2856) file8,backslash
c
        if( (inumgr .eq. 3) .and. (itype2 .eq. 3) )
     1      write(iun,2857) file8,backslash
c
 2855 format('     ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2856 format('     ','"',1a8,'"     ','notitle  with points, ',1a1)
 2857 format('     ','"',1a8,'"     ','notitle  with lines, ',1a1)
c
c        write the second data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the third data file
c
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return

        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c
        if(itype3 .eq. 1) write(iun,2870) file8
        if(itype3 .eq. 2) write(iun,2871) file8
        if(itype3 .eq. 3) write(iun,2872) file8
c
 2870 format('     ','"',1a8,'"     ','notitle  with dots')
 2871 format('     ','"',1a8,'"     ','notitle  with points')
 2872 format('     ','"',1a8,'"     ','notitle  with lines')
c
c        write the third data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
 3000 format(2x,e11.5,2x,e11.5)
c
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c
        close(iun22)
        close(iun)
c
        return
        end
c
c
c
c
c
        subroutine quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c
c        convert the user-specified Fortran unit number to 
c        character format
c
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c
 2000 format(1a8)
        read(dummy,2000) anum8
c
c        construct the file name on which the Gnuplot instructions
c        are to be written
c
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c
        write(iun,2250)
c
c        generate the title for the plot
c
        line(1)=blank
        line(2)=blank
        line(3)=quo
c
        call quamesslen(title,nchar,line(4))
c
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c
        write(iun,2270) (line(i),i=1,nchar+4)
c
 2280 format('   show title')
        write(iun,2280)
c
c        find the limits for both x and y
c
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c
        do 2300 i=1,n
        if(x(i) .lt. xmin) xmin=x(i)
        if(y(i) .lt. ymin) ymin=y(i)
        if(x(i) .gt. xmax) xmax=x(i)
        if(y(i) .gt. ymax) ymax=y(i)
 2300 continue
c
        if(inumgr .eq. 1) goto 2340
c
        do 2310 i=1,n2
        if(x2(i) .lt. xmin) xmin=x2(i)
        if(y2(i) .lt. ymin) ymin=y2(i)
        if(x2(i) .gt. xmax) xmax=x2(i)
        if(y2(i) .gt. ymax) ymax=y2(i)
 2310 continue
        if(inumgr .eq. 2) goto 2340
c
        do 2320 i=1,n3
        if(x3(i) .lt. xmin) xmin=x3(i)
        if(y3(i) .lt. ymin) ymin=y3(i)
        if(x3(i) .gt. xmax) xmax=x3(i)
        if(y3(i) .gt. ymax) ymax=y3(i)
 2320 continue
c
 2340 continue
c
        xcenter=(xmin+xmax)/2 
        ycenter=(ymax+ymin)/2
c
        xsize=(xmax-xmin)
        ysize=(ymax-ymin)
        size=xsize
        if(ysize .gt. size) size=ysize
        size=size*1.1
c
        xmin=xcenter-size/2
        xmax=xcenter+size/2
        ymin=ycenter-size/2
        ymax=ycenter+size/2
c
c        set the size of the stupid thing
c
cccc 2350 format(2x,' set size 0.75,1.0')
 2350 format(2x,' set size 1.0, 1.0')
        write(iun,2350) 
c
c        write the instructions 
c
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c
 2800 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots')

 2803 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points')

 2805 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines')
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) xmin,xmax,ymin,ymax,file8
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) xmin,xmax,ymin,ymax,file8
c
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) xmin,xmax,ymin,ymax,file8
c
 2830 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots, ',1a1)
c
 2833 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points, ',1a1)
c
 2835 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines, ',1a1)
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) xmin,xmax,ymin,ymax,file8,backslash
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2833) xmin,xmax,ymin,ymax,file8,backslash
c
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2835) xmin,xmax,ymin,ymax,file8,backslash
c
c        write the first data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x(i),y(i),i=1,n)        
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the second data file
c
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return
c
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c
 2850 format('     ', '"',1a8,'" ','notitle with dots')
 2851 format('     ', '"',1a8,'" ','notitle with points')
 2852 format('     ', '"',1a8,'" ','notitle with lines')
c
         if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2850) file8
c
         if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2851) file8
c
         if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2852) file8
c
 2855 format('     ', '"',1a8,'" ','notitle with dots, ',1a1)
 2856 format('     ', '"',1a8,'" ','notitle with points, ',1a1)
 2857 format('     ', '"',1a8,'" ','notitle with lines, ',1a1)
c
         if( (inumgr .ne. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2855) file8,backslash
c
         if( (inumgr .ne. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2856) file8,backslash
c
         if( (inumgr .ne. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2857) file8,backslash
c
c        write the second data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the third data file
c
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return

        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c
         if(itype3 .eq. 1) write(iun,2870) file8
         if(itype3 .eq. 2) write(iun,2871) file8
         if(itype3 .eq. 3) write(iun,2872) file8
 2870 format('     ','"',1a8,'"','notitle  with dots')
 2871 format('     ','"',1a8,'"','notitle  with points')
 2872 format('     ','"',1a8,'"','notitle  with lines')
c
c        write the third data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
 3000 format(2x,e11.5,2x,e11.5)
c
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c
        close(iun22)
        close(iun)
c
        return
        end
c
c
c
c
c
        SUBROUTINE quamesslen(MES,nchar,line)
        CHARACTER *1 MES(1),AST,line(1)
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c
        nchar=i1
        do 1800 i=1,nchar
        line(i)=mes(i)
 1800 continue
         RETURN
         END



