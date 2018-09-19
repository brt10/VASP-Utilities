c    program anlyze
c program to analyze nearest neighbor structures Blair Tuttle
      implicit none
      character froot*6,fname*100,froot1*7,test*6, atom(3)*14
      real a1vec(3),a2vec(3),a3vec(3),d(3,990),b(3,990)
      real b1vec(3),b2vec(3),b3vec(3),c(3,990)
      real x2(3),x1(3),x,y,z,cc,a_latt,b_latt,dist
      real dh2,xl(3,0:26),junk1,junk2
      integer nn,natom(900),ii,jj,i,j,k,i1,i2,i3, lbeta
      integer ntp,kk,nstep, na(10)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        dist = 6.9
        write(*,*) 'Input : two POSCAR files: POSCAR00 and POSCAR17'
        write(*,*) 'Output : 16 POSCAR files with positions scaled for
     c   VASP NEB SIMULATION'
        write(*,*) 'Works for Orho-rhombic supercells only'
        nstep  = 16
        write(*,*) '========================================='
        write(*,*) 'How many species in your POSCAR files ? '
        read(*,*) ntp

        open(unit=1,file = 'POSCAR00' )
      read(1,*) test
      read(1,*) a_latt
      read(1,*) a1vec(1),a1vec(2),a1vec(3)
      read(1,*) a2vec(1),a2vec(2),a2vec(3)
      read(1,*) a3vec(1),a3vec(2),a3vec(3)
      read(1,*) (atom(i),i=1,ntp)  
      read(1,*) (na(i), i = 1,ntp)  
      read(1,*) test
        open(unit=2,file = 'POSCAR17' )
      read(2,*) test
      read(2,*) a_latt
      read(2,*) a1vec(1),a1vec(2),a1vec(3)
      read(2,*) a2vec(1),a2vec(2),a2vec(3)
      read(2,*) a3vec(1),a3vec(2),a3vec(3)
      read(2,*) (atom(i),i=1,ntp)
      read(2,*) (na(i), i = 1,ntp)
      read(2,*) test
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        i = 0
        do ii = 1, ntp
      do jj =1,na(ii) 
        i = i + 1
        read(1,*) b(1,i),b(2,i),b(3,i)
        read(2,*) c(1,i),c(2,i),c(3,i)
        c(1,i) = c(1,i)*a_latt
        c(2,i) = c(2,i)*a_latt
        c(3,i) = c(3,i)*a_latt
        b(1,i) = b(1,i)*a_latt
        b(2,i) = b(2,i)*a_latt
        b(3,i) = b(3,i)*a_latt
      enddo
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i =1,3
          a1vec(i) = a1vec(i)*a_latt
          a2vec(i) = a2vec(i)*a_latt
          a3vec(i) = a3vec(i)*a_latt
        enddo
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do kk = 1,nstep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        froot = 'POSCAR'
        write(fname,'(a,i2.2)')froot,kk
      open(unit=3,file = fname)
      write(3,*) 'neb'
      write(3,*) a_latt
      write(3,*) 1.0,0.0,0.0
      write(3,*) 0.0,1.0,0.0
      write(3,*) 0.0,0.0,1.0
      write(1,*) (atom(i),i=1,ntp)
      write(1,*) (na(i), i = 1,ntp)
      write(3,*) test
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 8  ii = 1,ntp
        do 9 jj=1,na(ii) 
          i = i + 1
          j = i
        do 13 lbeta=0,26
          d(1,j) = (xl(1,lbeta)+b(1,j))
          d(2,j) = (xl(2,lbeta)+b(2,j))
          d(3,j) = (xl(3,lbeta)+b(3,j))
          cc=(c(1,i)-(xl(1,lbeta)+b(1,j)))**2 +
     1    (c(2,i)-(xl(2,lbeta)+b(2,j)))**2 +
     2    (c(3,i)-(b(3,j) +xl(3,lbeta)  )    )**2
        if (cc.le.dist*dist  ) then
        x = d(1,i) + kk*(c(1,j) - d(1,j))/(nstep+1)
        y = d(2,i) + kk*(c(2,j) - d(2,j))/(nstep+1)
        z = d(3,i) + kk*(c(3,j) - d(3,j))/(nstep+1)
        write(3,*) x/a_latt,y/a_latt,z/a_latt
        endif
13    continue
9     continue
8       continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      stop
      end

