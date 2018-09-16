    program anlyze
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C program analyzes nearest neighbor distances from a POSCAR VASP file
C  written by  Blair Tuttle 
C	2017
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      character test*26, atom*4 
      real a1vec(3),a2vec(3),a3vec(3),d(3,990),b(3,990)
      real b1vec(3),b2vec(3),b3vec(3),c(3,990) 
      real cc,a_latt,b_latt 
      real xl(3,0:26),junk1,junk2
      real dist(3,3)
      integer nn,ii,i,j,k,i1,i2,i3, lbeta 
	integer ntot,ntp,na(2),ntype(1000)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=1, file='CONTCAR')
      read(1,*) test 
      read(1,*) a_latt 
      read(1,*) a1vec(1),a1vec(2),a1vec(3)
      read(1,*) a2vec(1),a2vec(2),a2vec(3)
      read(1,*) a3vec(1),a3vec(2),a3vec(3)
      read(1,*) test
	write(*,*) 'How many species in your CONTCAR file?'
	read(*,*) 	 ntp   
      read(1,*) (na(i),i=1,ntp)     
C  only needed if selective dynamics chosen      read(1,*) test
      read(1,*) test
      i=0
	do ii = 1,ntp
      do j=1,na(ii)  
        i=i+1
        read(1,*) b(1,i),b(2,i),b(3,i)
        d(1,i) = b(1,i)*a1vec(1)+b(2,i)*a2vec(1)+b(3,i)*a3vec(1)
        d(2,i) = b(1,i)*a1vec(2)+b(2,i)*a2vec(2)+b(3,i)*a3vec(2)
        d(3,i) =b(1,i)*a1vec(3)+b(2,i)*a2vec(3)+ b(3,i)*a3vec(3)
        d(1,i) = d(1,i)*a_latt
        d(2,i) = d(2,i)*a_latt
        d(3,i) = d(3,i)*a_latt
        b(1,i) = d(1,i)
        b(2,i) = d(2,i)
        b(3,i) = d(3,i)
	ntype(i) = ii  
      enddo
      enddo
	ntot = i
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i =1,3
          a1vec(i) =a1vec(i)*a_latt
          a2vec(i) = a2vec(i)*a_latt
          a3vec(i) = a3vec(i)*a_latt
        enddo 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	write(*,*) 'choose three nearest neighbor distances in Angstroms:   '
	do i = 1, ntp
		do j=i,ntp	
	write(*,*) 'The distance between atom ',i,' and atom ',j,' is : '
	read(*,*) dist(i,j) 
	dist(j,i) = dist(i,j)
		enddo
	enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 51 i1=-1,1
      do 52 i2=-1,1
      do 53 i3=-1,1
          lbeta=9*(i3+1)+3*(i2+1)+(i1+1)
      do 55 k=1,3
          xl(k,lbeta)=i1*a1vec(k)+i2*a2vec(k)+i3*a3vec(k)
55    continue
53    continue
52    continue
51    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	open(unit = 2, file = 'NN_DIST.DAT')
	open(unit=3,file = 'NN_LIST.DAT')
	write(*,*) 'Neighbor distances writeen in file NN_DIST.DAT'
        write(*,*) 'Number of neighbors writeen in file NN_LIST.DAT'
        do 9 i=1,ntot     
        nn = 0
c loop over all possible neighbors
        write(2,111) i,ntype(i),b(1,i),b(2,i),b(3,i)
        do 11 j=1,ntot 
         if (j .eq. i) go to 12
        do 13 lbeta=0,26 
          c(i,j)=(b(1,i)-(xl(1,lbeta)+b(1,j)))**2 +
     1    (b(2,i)-(xl(2,lbeta)+b(2,j)))**2 +
     2    (b(3,i)-(b(3,j) +xl(3,lbeta)      )    )**2
          if (c(i,j).le.dist(ntype(i),ntype(j))**2  ) then 
            cc=sqrt(c(i,j) )
            write(2,121) i,ntype(i),j,ntype(j), cc 
	    nn = nn + 1
            goto 12
          endif
13    continue
12    continue
11    continue
10    continue
	write(3,*) i, nn
9     continue 
111   format('##',2i3,3f8.4)
121   format(4i4,1f8.3)
      stop
      end

