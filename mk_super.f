CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Code to make a supercell POSCAR dfile from a a CONTCAR unit cell
C       written by Blair Tuttle
C       2005 updated 2017
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      character test*1 , title*26, atom(10)*2 
      real x(10,1000,3),na(10)
      real a1vec(3),a2vec(3),a3vec(3),d(3,990),b(3,990)
      real b1vec(3),b2vec(3),b3vec(3),c(3,990) 
      real cc,a_latt,b_latt 
      real xl(3,0:26),junk1,junk2
      real xx,yy,zz,dist(3,3)
      integer jj,nx,ny,nz,nn,ii,i,j,k,i1,i2,i3, lbeta 
        integer ntot,ntp,ntype(1000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C               Read CONTCAR file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(unit=1, file='CONTCAR')
      read(1,*) title  
        write(*,*) title
      read(1,*) a_latt 
      read(1,*) a1vec(1),a1vec(2),a1vec(3)
      read(1,*) a2vec(1),a2vec(2),a2vec(3)
      read(1,*) a3vec(1),a3vec(2),a3vec(3)
        write(*,*) a3vec(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(*,*) 'How many species in your CONTCAR file?'
        read(*,*)        ntp   
        read(1,*) (atom(i), i = 1,ntp)
        write(*,*) (atom(i), i = 1,ntp)
      read(1,*) (na(i),i=1,ntp)     
C  only needed if selective dynamics chosen      read(1,*) test
      read(1,*) test  
        if (test.eq.'D') then
        else
          write(*,*) 'Must use Direct coordinates'
          Stop
        endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	do jj = 1,ntp  
          nn = na(jj)
	  do i = 1,nn  
	    read(1,*) x(jj,i,1),x(jj,i,2),x(jj,i,3)
	  enddo
        enddo
        close(unit = 1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                  now make POSCAR supercell
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(*,*) 'periodic integer multiples for x , y, z:  '
        read(*,*)  nx, ny, nz
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(unit=1, file='POSCAR')
      write(1,*) 'Supercell of ', title 
      write(1,*) a_latt
      write(1,*) nx*a1vec(1),nx*a1vec(2),nx*a1vec(3)
      write(1,*) ny*a2vec(1),ny*a2vec(2),ny*a2vec(3)
      write(1,*) nz*a3vec(1),nz*a3vec(2),nz*a3vec(3)
      write(1,*) (atom(i), i = 1,ntp)
      write(1,*) (nx*ny*nz*na(i),i=1,ntp)
C  only needed if selective dynamics chosen      read(1,*) test
      write(1,*) test
        do jj = 1,ntp
	do i = 1,nx 
	do j = 1,ny 
	do k = 1,nz 
        nn = na(jj)
	do ii = 1,nn 	
	   xx = (x(jj,ii,1) + i-1)/nx
           yy = (x(jj,ii,2) + j-1)/ny  
           zz = (x(jj,ii,3) + k-1)/nz   
	   write(1,*) xx,yy,zz
	enddo
	enddo
	enddo
	enddo
	enddo
	end 
