	implicit none
        real chg_ave,ra(1000,3),chg(500,500,500),rho(10)
        real x1,x2,x3,y1,y2,y3,z1,z2,z3,rlatt, xitems
        real chg_zave,xlines,chg_max
	real x,y,z,xdiff,ydiff,zdiff
        integer nlines,nt,ii,i,j,k,nitems, nx,ny,nz
        integer ixmax,iymax,izmax,ntot_lrg, ilrg(1000)
        integer nextra,ix,iy,iz
        integer n,kk,kn,in,jn
	integer ia(1000,3),na(10),nline 
	integer jj,ia_max
	character tmp*50
	write(*,*) 'INPUT: fort.1 : CHGCAR / LOCPOT / PARCHG'
	write(*,*) 'OUTPUT: XY_AVE.DAT, XZ_AVE.DAT, YZ_AVE.DAT' 
	write(*,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
C	write(*,*) 'How many species in your file? '
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	read(*,*) nt
	nt = 3
C	write(*,*) 'How many items per line in your file (fort.1) 5 or 10?'
CCCCCCCCCCCCCCCCCCCCCCC	read(*,*) nitems
	nitems =10
	read(1,*) tmp
	read(1,*) rlatt
	read(1,*) x1,y1,z1
	read(1,*) x2,y2,z2
	read(1,*) x3,y3,z3
	read(1,*) tmp
	read(1,*) (na(i), i = 1,nt)
	read(1,*) tmp
	ii = 0
	do i= 1,nt
	   do j = 1,na(i)
	   	ii = ii +1
	   	read(1,*) ra(ii,1),ra(ii,2),ra(ii,3) 
	  enddo
	enddo
        read(1,*) nx,ny,nz
	write(*,*) nx,ny,nz
	nlines = nx*ny*nz/nitems    
	xlines =  nx*ny*nz
	xitems = nitems
	xlines = xlines/xitems 
	write(*,*) nlines
        do ii = 1,nlines 
           read(1,*)  ( rho(k), k = 1,nitems   )
           do jj = 1, nitems   
            write(2,*) rho(jj)
           enddo
        enddo
	nextra = Nint( (xlines - nlines)*10) 
	write(*,*) xlines, nlines, nextra
	if (nextra.gt.0) then
	read(1,*)  ( rho(k), k = 1,nextra    )
	do ii = 1,nextra
		write(2,*) rho(ii)
	enddo
	endif
        close(unit =2 )
	open(unit=2)
        do k = 1, nz 
           do j = 1, ny
              do i = 1, nx
                read(2,*) chg(i,j,k)
              enddo
           enddo
        enddo
	close(unit=2)
	close(unit=1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        open(unit =3, file = 'XY_AVE.DAT')
        do k = 1, nz 
           chg_ave = 0
           do i = 1, nx
              do j = 1, ny
                chg_ave = chg_ave + chg(i,j,k)
              enddo
           enddo
          write(3,*) k*rlatt*z3/nz, chg_ave/(nx*ny)
        enddo
	close(unit =3 )
        open(unit =3, file = 'YZ_AVE.DAT')
        do i = 1, nx
           chg_ave = 0
           do j = 1, ny
              do k = 1, nz
                chg_ave = chg_ave + chg(i,j,k)
              enddo
           enddo
        write(3,*) i*rlatt*x1/nx, chg_ave/(nz*ny)
        enddo
        close(unit =3 )
        open(unit =3, file = 'XZ_AVE.DAT')
        do j = 1, ny
           chg_ave = 0
           do i = 1, nx
              do k = 1, nz
                chg_ave = chg_ave + chg(i,j,k)
              enddo
           enddo
        write(3,*) j*rlatt*y2/ny, chg_ave/(nz*nx)
        enddo
	close(unit=3)



        end
