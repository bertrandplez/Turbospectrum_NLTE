!
      subroutine read_modelatom(iunit,modelatomfile,maxn,modnlevel,
     &                          modenergy,modg,modion,modid,
     &                          species)
!
! read NLTE model atom
! created by BPlez 15/04-2020
!
! unit number to open, file name, max number of levels that can be accomodated,
! number of level read, energy, statistical weight, ionisation stage, and 
! identification of the levels
!
      implicit none
      integer maxn
      character modelatomfile*256,oneline*256,species*2
      character modid(maxn)*40
      integer iunit,modnlevel,modion(maxn),i
      real modenergy(maxn),abundance,mass,modg(maxn)
      logical header

      open(iunit,file=modelatomfile,status='old')
      header=.true.
      do while (header)
        read(iunit,10,end=99) oneline
        if (oneline(1:1).ne.'*') then
          read(oneline,10) species
          read(iunit,*) abundance,mass
          read(iunit,*) modnlevel
          header=.false.
        endif
      enddo
      i=0
      do while (i.lt.modnlevel)
        read(iunit,10,end=99) oneline
        if (oneline(1:1).ne.'*') then
          i=i+1
          read(oneline,*) modenergy(i),modg(i),modid(i),modion(i)
        endif
      enddo

      print*,'read_modelatom ',species,modnlevel,' levels'
      do i=1,modnlevel
        print*,i,modion(i),modid(i),modenergy(i),modg(i)
      enddo
10    format(a)

      close(iunit)

      return

99    print*,'unexpected end of file in read_modelatom'
      stop 'stop in read_modelatom'

      end
