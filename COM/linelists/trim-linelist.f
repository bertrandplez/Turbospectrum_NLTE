      program trim_linelist

! B. Plez 2022-12-06

      implicit none

      character*128 filein,fileout
      character*256 oneline,bla
      integer i,nline,nwrite,ion
      real*8 wave,wmin,wmax

      print*,' linelist to read and trim'
      read(5,10) filein
      open(10,file=filein,status='old')

      print*,' output linelist'
      read(5,10) fileout
      open(20,file=fileout,status='new')
10    format(a)

      print*,'min and max wavelength for trimming'
      read(5,*) wmin,wmax
      if (wmin.gt.wmax) then
        wave=wmax
        wmax=wmin
        wmin=wave
      endif

      open(15,status='scratch',form='unformatted')

      do while (.true.)
        rewind(15)
        nwrite=0
        read(10,10,end=99) oneline
        read(oneline,*) bla,ion,nline
        write(15) oneline
        read(10,10) oneline
        write(15) oneline
        do i=1,nline
          read(10,10) oneline
          read(oneline,*) wave
          if (wave.le.wmax.and.wave.ge.wmin) then
            nwrite=nwrite+1
            write(15) oneline
          endif
        enddo
        if (nwrite.gt.0) then
          rewind(15)
          read(15) oneline
          read(oneline,*) bla,ion
          write(20,20) trim(bla),ion,nwrite
20        format('''',a,'''',1x,i2,1x,i7)
          read(15) oneline
          write(20,10) trim(oneline)
          do i=1,nwrite
            read(15) oneline
            write(20,10) trim(oneline)
          enddo
        endif
      enddo
99    continue
      end
