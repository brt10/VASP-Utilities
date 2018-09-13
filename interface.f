         program www 
c==============================================================
c       program performs WWW bond switching algorithm     
c=================================================================
	implicit none 
	include 'parameter.h' 
c================================================================
c  values below from input file
c==============================================================
	integer ibond					!ibond = 1 read bonds , =2 find bonds
	real kt          				!temperature of bond switches 
	real l_vec(3,3)					!lattice vectors 
	integer nsw     				!number of switches 
	integer nsp 					!number of species 
	integer na(nsp_max)				!number of atoms per spec, total atoms
	real a_vec(nsp_max,na_max,3)  			! a_vec = atom vectors
c==============================================================
c   values above from input file
c=====================================================
	integer isw,count,nout 
	integer nbtot(nsp_max,na_max)  			!number of bonds per atom 
        real b_vec(nsp_max,na_max,nbond_max,3)  	!nn atom vectors
	integer bi(nsp_max,na_max,nbond_max,3)   	! index of nn
        real force(nsp_max,na_max,3)			! force on atoms
	real Ei,Ef,PE
	real a_vec_o(nsp_max,na_max,3)    		! old a_vec = atom vectors
	real b_vec_o(nsp_max,na_max,nbond_max,3)  	!old nn atom vectors
	integer bi_o(nsp_max,na_max,nbond_max,3)   	!old  index of nn
	real xfrac,prob
c=============================================================
	nout =10
	call input(kt,l_vec,nsw,nsp,na,a_vec,ibond,b_vec,bi,nbtot,PE)
	Ei = PE
	count = 0
	write(*,*) count,count,Ei
c===========================================================================
c   begin main loop
c===========================================================================
        do isw = 1,nsw
	call save(nsp,na,nbtot,a_vec,b_vec,bi,a_vec_o,b_vec_o,bi_o,1) 
	call bond_switch(nsp,na,l_vec,a_vec,b_vec,bi,nbtot)
	call relax(nsp,na,nbtot,a_vec,b_vec,bi,l_vec,PE)
        Ef = PE 
	prob = exp(-(Ef-Ei)/kt)
	xfrac =  rand()
	if (prob.lt.xfrac) then 	!go back to old coords 
           PE = Ei
           call save(nsp,na,nbtot,a_vec,b_vec,bi,a_vec_o,b_vec_o,bi_o,2)
	else   				!keep new coords
	   count = count + 1
	   write(*,*) isw,count,Ef 
	   Ei = Ef
	   PE = Ei
           call re_align(nsp,na,a_vec,nbtot,bi,b_vec,l_vec)
	   if (count.eq.nout) then
             call  output(nsp,na,nbtot,bi,b_vec,a_vec,l_vec)
	     nout = nout + 10
	   endif 
	endif 
        enddo	
c==========================================================================
c   end main loop
c==========================================================================
	call relax(nsp,na,nbtot,a_vec,b_vec,bi,l_vec,PE)
	call re_align(nsp,na,a_vec,nbtot,bi,b_vec,l_vec)
	call output(nsp,na,nbtot,bi,b_vec,a_vec,l_vec)
	 end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine relax(nsp,na,nbtot,a_vec,b_vec,bi,l_vec,PE)
c==============================================================================
        include 'parameter.h'
        implicit none
c=========================passed variables
        integer nsp                     	!number of species
        integer na(nsp_max)            	 	!number of atoms per spec
        integer nbtot(nsp_max,na_max)   	!number of bonds for each atom
	real l_vec(3,3)
        real a_vec(nsp_max,na_max,3)    	! a_vec = atom vectors
	real a_vec_p(nsp_max,na_max,3)            !previous step  a_vec = atom vectors
        real b_vec(nsp_max,na_max,nbond_max,3)  !nn atom vectors
        integer bi(nsp_max,na_max,nbond_max,3)  ! index of nn
c===================================== eventually to be put in parameter.h 
c============================== but now are found below
	real mass(nsp_max)                      !mass from parameter.h in units of amu
	real force_tol				!from parameter.h
	real  time_step
	integer max_time_step			!from parameter.h
c=========================new variables
	real s_vec(nsp_max,na_max,nbond_max,3)  ! shell atom vectors
        integer si(nsp_max,na_max,nbond_max,2)  ! index of shell atoms 
        integer nstot(nsp_max,na_max)           !number of shell atoms for each atom
        real force(nsp_max,na_max,3)		!forces on atoms
	real d_a_vec(nsp_max,na_max,3)		!delta x for each atom
	real f_max				!max to determine if tolerance reached
	real f_tol
	real dtf	!conversion factor:  dt^2/ mass time F(Ev/ A) so result in A 
	real deltaX				!temp variable for position
	integer it,i,j,k,ii,jj,kk 
	real test
	real Eo,Ei,Ef,PE
	real DE,DEtol,xfrac,prob 
	integer irelax 
c===================================================================
	f_tol = 0.2
	time_step = 0.5 
	dtf = time_step**2   
	mass(1) = 28.08*931.5/9 	!convert amu -> eV fs fs / ang ang 
        max_time_step =200 
	DEtol = 28    !eV
	irelax = 1 	! for now 1 is the only option 
c===================================================================
c============ initializing positions ===============================  
c===================================================================
        call find_bonds(nsp,na,l_vec,a_vec,b_vec,bi,nbtot,
     1                                        s_vec,si,nstot,3)	
c=================================================================
       Eo = PE
       call Energy(nsp,na,nbtot,a_vec,b_vec,bi,
     1                                        s_vec,si,nstot,Ei)
	DE = Ei - Eo
	xfrac = rand()
	prob = 1.5**(-(DE-DEtol))
	if (prob.lt.xfrac) then
	   PE = Ei
	   write(8,*) 'DE too large'
	   return 
	endif
c==========================================================================
	do it = 1,max_time_step   	!loop it 
	   call Calc_Force(nsp,na,nbtot,a_vec,b_vec,bi,
     1                    s_vec,si,nstot,force) 
	   f_max = 0
	   do i = 1, nsp  		!loop i 
	      do j = 1,na(i)		!loop j 
	         do k = 1,3		!loop k 
	          if (irelax.eq.1) then   ! simple steepest descent 
		    d_a_vec(i,j,k) = force(i,j,k)*dtf/mass(i)  
                    a_vec(i,j,k) = a_vec(i,j,k) + d_a_vec(i,j,k) 
                    if (force(i,j,k).gt.f_max) then
                       f_max = force(i,j,k)
                    endif
	          else
		   STOP 'irelax must be 1 for now'
	          endif
		 enddo			!loop k 
	      enddo			!loop j  
	   enddo			!loop i
c===============================================================
c========move b_vec by amount from a_vec 
c===============================================================
	  do i = 1, nsp
	    do j = 1 , na(i)
	      do k = 1,nbtot(i,j)
		ii = bi(i,j,k,1)
		jj = bi(i,j,k,2)
	        b_vec(i,j,k,1) = b_vec(i,j,k,1) + d_a_vec(ii,jj,1)
		b_vec(i,j,k,2) = b_vec(i,j,k,2) + d_a_vec(ii,jj,2)
                b_vec(i,j,k,3) = b_vec(i,j,k,3) + d_a_vec(ii,jj,3)   	
             enddo
	    enddo
	  enddo	
          do i = 1, nsp
            do j = 1 , na(i)
              do k = 1,nstot(i,j)
                ii = si(i,j,k,1)
                jj = si(i,j,k,2)
                s_vec(i,j,k,1) = s_vec(i,j,k,1) + d_a_vec(ii,jj,1)
                s_vec(i,j,k,2) = s_vec(i,j,k,2) + d_a_vec(ii,jj,2)
                s_vec(i,j,k,3) = s_vec(i,j,k,3) + d_a_vec(ii,jj,3)
              enddo
            enddo
          enddo  
	   if (f_max.lt.f_tol) then
	      goto 400			!exit loop it 
	   endif
           if ( f_max.lt.10) then
	      dtf = 25 
           endif  
	   if (f_max.lt.5) then 
	      dtf = 100 
	   endif
	enddo				!loop it 
c=============================== begin keep atoms in unit cell
400 	continue
c============================== now find energy of relaxed system 
        call Energy(nsp,na,nbtot,a_vec,b_vec,bi,
     1                                        s_vec,si,nstot,Ef)
c================================================
c	write(9,*) it, dtf,Ef
	PE = Ef
	return          
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine re_align(nsp,na,a_vec,nbtot,bi,b_vec,l_vec) 
c==================================================================
        implicit none
        include 'parameter.h'
c================================================================
        integer nsp                     !number of species
        integer na(nsp_max)             !number of atoms per spec
        integer nbtot(nsp_max,na_max)           !number of bonds per atom
        integer bi(nsp_max,na_max,nbond_max,3)   ! index of nn
        integer i1, i2, i3, lbeta 
	integer i,j,k,ii,jj,kk
        real a_vec(nsp_max,na_max,3)    ! a_vec = atom vectors
        real b_vec(nsp_max,na_max,nbond_max,3)  !nn atom vectors
        real xl(3,30)
	real l_vec(3,3)  
	real ds2
	real rs
	real x,y,z
c============================================================
c               bring atoms back into unit cell
c============================================================
       do i = 1, nsp
           do j = 1 , na(i)
              do kk = 1,3
                if (a_vec(i,j,kk).lt.0.0) then
                   a_vec(i,j,kk) = a_vec(i,j,kk) + l_vec(kk,kk) 
                elseif (a_vec(i,j,kk).gt.l_vec(kk,kk) ) then
                   a_vec(i,j,kk) = a_vec(i,j,kk) - l_vec(kk,kk)  
                endif
              enddo
          enddo
        enddo
c============================================================
c           re align neighbor atoms 
c=============================================================
      do 51 i1=-1,1
      do 52 i2=-1,1
      do 53 i3=-1,1
          lbeta=9*(i3+1)+3*(i2+1)+(i1+1)
      do 55 k=1,3
          xl(k,lbeta)=i1*l_vec(1,k)+i2*l_vec(2,k)+i3*l_vec(3,k)
55    continue
53    continue
52    continue
51    continue
c=========================================================
      do i = 1,nsp
           do j = 1,na(i)
              do k = 1,nbtot(i,j)
		 ds2 = 2500 
                 ii = bi(i,j,k,1)
                 jj = bi(i,j,k,2)
                 do lbeta=0,26
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.ds2) then
                      b_vec(i,j,k,1) = x
                      b_vec(i,j,k,2) = y
                      b_vec(i,j,k,3) = z
                      ds2 = rs 
                   endif
                 enddo
              enddo 
           enddo 
        enddo   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine save(nsp,na,nbtot,a_vec,
     1                          b_vec,bi,a_vec_o,b_vec_o,bi_o,flag)   
c====================================================================================
c  subroutine saves old coods if flag = 1 or returns old coords to new if flag = 2 
c====================================================================================
        implicit none
        include 'parameter.h'
c================================================================ 
        integer nsp                     !number of species
        integer na(nsp_max)             !number of atoms per spec
        real a_vec(nsp_max,na_max,3)    ! a_vec = atom vectors
        integer nbtot(nsp_max,na_max)           !number of bonds per atom
        real b_vec(nsp_max,na_max,nbond_max,3)  !nn atom vectors
        integer bi(nsp_max,na_max,nbond_max,3)   ! index of nn
        real a_vec_o(nsp_max,na_max,3)    ! old a_vec = atom vectors
        real b_vec_o(nsp_max,na_max,nbond_max,3)  !old nn atom vectors
        integer bi_o(nsp_max,na_max,nbond_max,3)   !old  index of nnu 
	integer flag
	integer i,j,k,jj
c=========================================================================
	if (flag.eq.1) then
c==================================== now save bonding and position vectors
        do i = 1, nsp
           do j = 1,na(i)
              do k = 1,3
                 a_vec_o(i,j,k) = a_vec(i,j,k)
              enddo
           enddo
        enddo
        do i = 1,nsp
           do j = 1,na(i)
              do jj = 1,nbtot(i,j)
                do k =1 ,3
                   b_vec_o(i,j,jj,k) = b_vec(i,j,jj,k)
                enddo
              enddo
           enddo
        enddo
        do i = 1,nsp
           do j = 1,na(i)
              do jj = 1,nbtot(i,j)
                do k =1 ,3
                   bi_o(i,j,jj,k) = bi(i,j,jj,k)
                enddo
              enddo
           enddo   
        enddo
	return
	elseif (flag.eq.2) then
        do i = 1, nsp
           do j = 1,na(i)
              do k = 1,3
                 a_vec(i,j,k) = a_vec_o(i,j,k)
              enddo
           enddo
        enddo
        do i = 1,nsp
           do j = 1,na(i)
              do jj = 1,nbtot(i,j)
                do k =1 ,3
                   b_vec(i,j,jj,k) = b_vec_o(i,j,jj,k)
                enddo
              enddo
           enddo
        enddo
        do i = 1,nsp
           do j = 1,na(i)
              do jj = 1,nbtot(i,j)
                do k =1 ,2
                   bi(i,j,jj,k) = bi_o(i,j,jj,k)
                enddo
              enddo
           enddo  
        enddo
        return
	else
	  STOP 'iflag in save must be 1 or 2'
	endif
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine Energy(nsp,na,nbtot,a_vec,b_vec,bi,
     1                                        s_vec,si,nstot,PE)
c==============================================================================
c        calculates the potential energy for a keating-style potential see Tersoff et al.
c==============================================================================
        include 'parameter.h'
        implicit none
c======================================================passed variables
        integer nsp                     !number of species
        integer na(nsp_max)             !number of atoms per spec
        integer nbtot(nsp_max,na_max)   !number of bonds for each atom
        real a_vec(nsp_max,na_max,3)    ! a_vec = atom vectors
        real b_vec(nsp_max,na_max,nbond_max,3)  !nn atom vectors
        integer bi(nsp_max,na_max,nbond_max,3)   ! index of nn
        integer nstot(nsp_max,na_max)   !number of non bonds for each atom
        real s_vec(nsp_max,na_max,nbond_max,3)  !non bond atom vectors
        integer si(nsp_max,na_max,nbond_max,2)   ! index of non bonds
        real PE,PE1,PE2 
c===================================================new variables
        integer i,j,k,li,lj,mi,mj,m
        real dx1,dy1,dz1,dr1,dx2,dy2,dz2
        real dr1_dr2,dr2,DR,r1dot2
        real kqtmp,cos_q,cos_qo
        real d_cos_q(3)                         ! derivative of cos(q)
c================================== force constants
        real kr(nsp_max,nsp_max)
        real ro(nsp_max,nsp_max)
        real kq(nsp_max,nsp_max,nsp_max)
        real qo(nsp_max,nsp_max,nsp_max)
	real ds
	real gamma
c==================================================================
        kr(1,1) = 9.08          ! units = eV/Ang^2
        ro(1,1) = 2.3513        ! units = Ang
        kq(1,1,1) = 3.58        ! units = eV
        qo(1,1,1) = 1.9106      ! units = radians
        ds = 3.8 
        gamma = 0.5
c=================================================================
c       values for long and short bonds to simulate si-sio2
c=================================================================
C	kr(1,2) = 1.89		! units = eV / ang^2
C	ro(1,2) = 3.04 		! units = Ang
CCCCCCCCCCCCCCC  Above from Vanderbilt paper PRB, 1996
CCCCCCCCCCCCCCC  Below has new guesses
	kr(1,2)  = 3.00 
	ro(1,2) = 3.00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  just a guess
	kq(1,2,2) = 4.03	! units = eV
	kq(1,2,1) = 3.81	! units = eV
	kq(1,1,2) = 3.81	! units = eV
	qo(1,2,1) = 1.9106	! units = radians
	qo(1,1,2) = 1.9106	! units = radians 
	qo(1,2,2) = 1.9106	! units = radians
c=================================================================
c========================     calc energy for each atom
c=================================================================
        PE = 0.0
	PE1 = 0.0
	PE2 = 0.00
        do i = 1,nsp
           do j = 1,na(i)
              do k = 1,nbtot(i,j)
                 dx1 = a_vec(i,j,1) - b_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - b_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - b_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 DR = dr1-ro(i,bi(i,j,k,3))
                 PE = PE + 0.25*kr(i,bi(i,j,k,3))*DR**2 !0.25 since double counting terms
              enddo
              do k = 1,nbtot(i,j) - 1
                 dx1 = a_vec(i,j,1) - b_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - b_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - b_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 do m = k+1,nbtot(i,j)
                    dx2 = a_vec(i,j,1) - b_vec(i,j,m,1)
                    dy2 = a_vec(i,j,2) - b_vec(i,j,m,2)
                    dz2 = a_vec(i,j,3) - b_vec(i,j,m,3)
                    dr2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                    cos_q = (dx1*dx2+dy1*dy2+dz1*dz2)/(dr1*dr2)
                    cos_qo = cos(qo(i,bi(i,j,k,3),bi(i,j,m,3)))
                    kqtmp = kq(i,bi(i,j,k,3),bi(i,j,m,3)  )
                    PE = PE + 0.5*kqtmp*(cos_q - cos_qo)**2
                 enddo
             enddo
          enddo  ! na loop
        enddo    ! nsp loop
        do i = 1, nsp
           do j = 1,na(i)
              do k = 1,nstot(i,j)
                 dx1 = a_vec(i,j,1) - s_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - s_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - s_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 if (dr1.lt.ds) then
                    DR = (ds - dr1)
                 else
                    DR = 0
                 endif 
                 PE = PE +gamma*DR*DR*DR 
              enddo
           enddo
        enddo
c==================================================================================
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine Calc_Force(nsp,na,nbtot,a_vec,b_vec,bi,
     1                    s_vec,si,nstot,force)
c==============================================================================
c        calculates the force for a keating-style potential see Tersoff et al.
c==============================================================================
	include 'parameter.h'
	implicit none
c======================================================passed variables
        integer nsp                     !number of species
        integer na(nsp_max)             !number of atoms per spec
        integer nbtot(nsp_max,na_max)   !number of bonds for each atom
        real a_vec(nsp_max,na_max,3)    ! a_vec = atom vectors
        real b_vec(nsp_max,na_max,nbond_max,3)  !nn atom vectors
        integer bi(nsp_max,na_max,nbond_max,3)   ! index of nn
        real s_vec(nsp_max,na_max,nbond_max,3)  ! shell atom vectors
        integer si(nsp_max,na_max,nbond_max,2)  ! index of shell atoms
        integer nstot(nsp_max,na_max)           !number of shell atoms for each atom
c===================================================new variables 
	integer i,j,k,li,lj,mi,mj,m
	real dx1,dy1,dz1,dr1,dx2,dy2,dz2
	real dr1_dr2,dr2,DR,r1dot2 
	real kqtmp,cos_q,cos_qo
	real force(nsp_max,na_max,3), f_tmp	! force on atoms
      	real d_cos_q(3)                         ! derivative of cos(q)      
c===============   variables from parameter.h 
        real kr(nsp_max,nsp_max)
        real ro(nsp_max,nsp_max)
        real kq(nsp_max,nsp_max,nsp_max)
        real qo(nsp_max,nsp_max,nsp_max)
	real ds
	real gamma
c==================================================================
        kr(1,1) = 9.08          ! units = eV/Ang^2
        ro(1,1) = 2.3513        ! units = Ang
        kq(1,1,1) = 3.58        ! units = eV
        qo(1,1,1) = 1.9106      ! units = radians
	ds = 4.5	        !  ds = 3.8 for si, ds = 4.5 for sio2 
	gamma = 0.5
c=================================================================
c       values for long and short bonds to simulate si-sio2
c=================================================================
        kr(1,2) = 1.89          ! units = eV / ang^2
        ro(1,2) = 3.04          ! units = Ang
        kq(1,2,2) = 4.03        ! units = eV
        kq(1,2,1) = 3.81        ! units = eV
        kq(1,1,2) = 3.81        ! units = eV
        qo(1,2,1) = 1.9106      ! units = radians
        qo(1,1,2) = 1.9106      ! units = radians 
        qo(1,2,2) = 1.9106      ! units = radians
c=================================================================
c============================================ radial loop
c==================================================================
	do i = 1, nsp
	   do j = 1,na(i)
	      force(i,j,1) = 0
	      force(i,j,2) = 0
	      force(i,j,3) = 0
              do k = 1,nbtot(i,j)
                 dx1 = a_vec(i,j,1) - b_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - b_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - b_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 DR = (dr1-ro(i,bi(i,j,k,3)))
		 force(i,j,1) = force(i,j,1) - kr(i,bi(i,j,k,3))*DR*dx1/dr1   
		 force(i,j,2) = force(i,j,2) - kr(i,bi(i,j,k,3))*DR*dy1/dr1
		 force(i,j,3) = force(i,j,3) - kr(i,bi(i,j,k,3))*DR*dz1/dr1
	      enddo
	   enddo
	enddo
c===================================================== angular loop 
c========================================== atom at center of angle
        do i = 1,nsp          ! i loop
           do j = 1,na(i)     ! j loop
              do k = 1,nbtot(i,j) - 1      !k loop
                 dx1 = - a_vec(i,j,1) + b_vec(i,j,k,1)
                 dy1 = - a_vec(i,j,2) + b_vec(i,j,k,2)
                 dz1 = - a_vec(i,j,3) + b_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
		 do m = k+1,nbtot(i,j)     !m loop
                    dx2 = - a_vec(i,j,1) + b_vec(i,j,m,1)
                    dy2 = - a_vec(i,j,2) + b_vec(i,j,m,2)
                    dz2 = - a_vec(i,j,3) + b_vec(i,j,m,3)
                    dr2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                    dr1_dr2 = dr1*dr2
                    r1dot2 = dx1*dx2 + dy1*dy2 + dz1*dz2
                    cos_q = (dx1*dx2+dy1*dy2+dz1*dz2)/(dr1_dr2)
                    cos_qo = cos(qo(i,bi(i,j,k,3),bi(i,j,m,3)))
                    kqtmp = kq(i,bi(i,j,k,3),bi(i,j,m,3)  )
         d_cos_q(1) =  (- dx1 - dx2 +dx2*r1dot2/dr2**2
     1                             + dx1*r1dot2/dr1**2)/(dr1_dr2)
         d_cos_q(2) =  (- dy1 + dy2*r1dot2/dr2**2
     1                             - dy2 + dy1*r1dot2/dr1**2)/(dr1_dr2)
         d_cos_q(3) =  (- dz1 + dz2*r1dot2/dr2**2
     1                             - dz2 + dz1*r1dot2/dr1**2)/(dr1_dr2)
c======================================================================
        force(i,j,1) = force(i,j,1) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(1)
        force(i,j,2) = force(i,j,2) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(2)
        force(i,j,3) = force(i,j,3) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(3)
                  enddo  ! m loop
                enddo    ! k loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    enddo	 ! j loop
 	enddo		 ! i loop
c================================================== atom not at center 
        do i = 1,nsp          ! i loop
           do j = 1,na(i)     ! j loop
              do k = 1,nbtot(i,j)       !k loop
                 dx1 =   a_vec(i,j,1) - b_vec(i,j,k,1)
                 dy1 =   a_vec(i,j,2) - b_vec(i,j,k,2)
                 dz1 =   a_vec(i,j,3) - b_vec(i,j,k,3)
		 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
	         mi = bi(i,j,k,1)
	 	 mj = bi(i,j,k,2)
                 do m = 1,nbtot(mi,mj)     !m loop
	if ((i.eq.bi(mi,mj,m,1)).and.(j.eq.bi(mi,mj,m,2))) then
		    goto 110
	else
                    dx2 = - a_vec(mi,mj,1) + b_vec(mi,mj,m,1)
                    dy2 = - a_vec(mi,mj,2) + b_vec(mi,mj,m,2)
                    dz2 = - a_vec(mi,mj,3) + b_vec(mi,mj,m,3)
                    dr2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
                    dr1_dr2 = dr1*dr2
                    r1dot2 = dx1*dx2 + dy1*dy2 + dz1*dz2
                    cos_q = (dx1*dx2+dy1*dy2+dz1*dz2)/(dr1_dr2)
                    cos_qo = cos(qo(bi(i,j,k,3),i,bi(mi,mj,m,3)))
                    kqtmp = kq( bi(i,j,k,3),i,bi(mi,mj,m,3)  )
        d_cos_q(1) =  (+dx2 - dx1*r1dot2/dr1**2)/(dr1_dr2)
        d_cos_q(2) =  (+dy2 - dy1*r1dot2/dr1**2)/(dr1_dr2)
        d_cos_q(3) =  (+dz2 - dz1*r1dot2/dr1**2)/(dr1_dr2)
        force(i,j,1) = force(i,j,1) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(1)
        force(i,j,2) = force(i,j,2) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(2)
        force(i,j,3) = force(i,j,3) -
     1                    kqtmp*(cos_q - cos_qo)*d_cos_q(3)
	endif
110              enddo  ! m loop
		enddo   ! k loop
           enddo        ! j loop
        enddo           ! i loop
c=============================================================================
c        repulsive loop  
c=============================================================================
        do i = 1, nsp
           do j = 1,na(i)
              do k = 1,nstot(i,j)
                 dx1 = a_vec(i,j,1) - s_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - s_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - s_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 if (dr1.lt.ds) then
		    DR = (ds - dr1)
		 else
		    DR = 0
		 endif	
                 force(i,j,1) = force(i,j,1) +3*gamma*DR*DR*dx1/dr1
                 force(i,j,2) = force(i,j,2) +3*gamma*DR*DR*dy1/dr1
                 force(i,j,3) = force(i,j,3) +3*gamma*DR*DR*dz1/dr1
              enddo
           enddo
        enddo
c==============================================================================
100	return	
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine bond_switch(nsp, na, l_vec, a_vec, b_vec, bi,nbtot)
c=========================================================================
c  Subprogram takes a random atom and considers a random bonding atom.  Then, from the remaining
c  bonding atoms it chooses one at random. Then, it chooses another bonding atom of the second atom.
c  This last atom is then is made a bond with the first bond atom.  
c  Then, to preserve number of bonds, a bond atom from the first bond atom is switched to become a bond
c  of the third atom. See below and  PRL vol.85 pg. 1984.
c=====================================
	include 'parameter.h'
        integer nsp
	integer nox       !  atom number associated with the oxide
        integer na(nsp_max)
	integer nbtot(nsp_max,na_max) 
	real l_vec(3,3)
        real a_vec(nsp_max,na_max,3) 
	real b_vec(nsp_max,na_max,nbond_max,3)
        integer bi(nsp_max, na_max, nbond_max,3)
	integer bitmp
        real rs ,x, y, z
        real xl(3,30)
        integer i1, i2, i3, lbeta , k
        integer i, j , ii, jj
	integer n1,n2,n3,n4,n5
	integer isp,iat,n1sp,n1at,n2sp,n2at,n3sp,n3at
c================================= do bond switch =================
c    First, define atoms to be switched:  
c    atom i is bonded to n1 and n2 . atom n1 also is bonded to n3
c    the bonds between n1-n3 & i-n2 switch so in the end
c    the bonds are i-n3 and n1-n2    
c===================================================
600	continue    !allows a start over if three fold rings occur
	isp = 1 + nsp*rand()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	nox = 41 
	iat = nox + (na(isp) - nox - 4 )*rand()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Above is special code to only switch oxide silicons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  	n1 = 1 + nbtot(isp,iat)*rand()
	n2 = 1 + nbtot(isp,iat)*rand()
        if (n2.eq.n1) then
          n2 = n1 + 1 + (nbtot(isp,iat)-2)*rand()
          if (n2.gt.nbtot(isp,iat)) then
              n2 = n2 - nbtot(isp,iat)
          endif
        endif 
        n1sp = bi(isp,iat,n1,1)
	n1at = bi(isp,iat,n1,2)
        n2sp = bi(isp,iat,n2,1)
	n2at = bi(isp,iat,n2,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c==  make sure n2 and iat do not share neighbor 
c==  dissallowing this prevents three fold rings
c====================================================
	do i = 1,nbtot(isp,iat)
	   do j = 1,nbtot(n2sp,n2at)
	        if ( (bi(isp,iat,i,1).eq.bi(n3sp,n3at,j,1)).and.
     1             (bi(isp,iat,i,2).eq.bi(n3sp,n3at,j,2))) then
	        goto 600
	      endif
	   enddo
	enddo 	
c========================== no three fold rings allowed
	n3 = 1 + nbtot(n1sp,n1at)*rand()
	n3sp = bi(n1sp,n1at,n3,1)
	n3at = bi(n1sp,n1at,n3,2)
	if ((n3sp.eq.isp).and.(n3at.eq.iat)) then
	  n3 = n3 + 1 + ( nbtot(n1sp,n1at) -2 )*rand() 
	  if (n3.gt.nbtot(n1sp,n1at)) then
	      n3 = n3 - nbtot(n1sp,n1at) 
	  endif
          n3sp = bi(n1sp,n1at,n3,1)
          n3at = bi(n1sp,n1at,n3,2)
        endif
c================== again disallow 3-fold rings
        do i = 1,nbtot(isp,iat) 
           do j = 1,nbtot(n2sp,n2at)    
              if ((bi(n3sp,n3at,i,1).eq.bi(n2sp,n2at,j,1)).and.
     1         (bi(n3sp,n3at,i,2).eq.bi(n2sp,n2at,j,2))) then   
                goto 600        
              endif     
           enddo        
        enddo   
c===================== 3 fold rings disallowed
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Below is special code to only switch oxide silicons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        if ( (n1at.lt.nox).or.(n1at.gt.na(isp)-4) ) then
          goto 600
        endif
        if ( (n2at.lt.nox).or.(n2at.gt.na(isp) -4) ) then
          goto 600
        endif
        if ( (n3at.lt.nox).or.(n3at.gt.na(isp)-4) ) then
          goto 600
        endif
C	write(*,*) iat,n1at,n2at,n3at
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c================================================
c	 now switch n2 and n3  ' 
c===============================================
        do i = 1,3
           bitmp = bi(isp,iat,n2,i)
           bi(isp,iat,n2,i) = bi(n1sp,n1at,n3,i)
           bi(n1sp,n1at,n3,i) = bitmp 
        enddo
c=================================================
c   update nn b_vec s
c=================================================
      do 51 i1=-1,1
      do 52 i2=-1,1
      do 53 i3=-1,1
          lbeta=9*(i3+1)+3*(i2+1)+(i1+1)
      do 55 k=1,3
          xl(k,lbeta)=i1*l_vec(1,k)+i2*l_vec(2,k)+i3*l_vec(3,k)
55    continue
53    continue
52    continue
51    continue
c==================================================================
      i=isp 
      j = iat  
      ii = n3sp  
      jj = n3at   
      d2 = 1000   ! choose any initial value greater than nn dist
                 do lbeta=0,26         
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
                      b_vec(i,j,n2,1) = x
                      b_vec(i,j,n2,2) = y
                      b_vec(i,j,n2,3) = z
                      d2 = rs 
                   endif
                 enddo                   
c===========================================================
      i=n1sp   
      j =n1at   
      ii = n2sp  
      jj = n2at    
      d2 = 1000
                   do lbeta=0,26
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
                      b_vec(i,j,n3,1) = x
                      b_vec(i,j,n3,2) = y
                      b_vec(i,j,n3,3) = z
                      d2 = rs
                   endif
                 enddo                    ! 5 loop
c================================================================
c   now re-set bond of n2 ' 
c=========================================================
	do i = 1, nbtot(n2sp,n2at)
	  if((bi(n2sp,n2at,i,1).eq.isp).and.
     1         (bi(n2sp,n2at,i,2).eq.iat)) then 
         n4 = i
	      goto 20
	  endif
	enddo
20	continue
        bi(n2sp,n2at,n4,1) = n1sp
	bi(n2sp,n2at,n4,2) = n1at 
c===========================================================
      i=n2sp
      j =n2at
      ii = n1sp
      jj = n1at
      d2 = 1000
		 do lbeta=0,26
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
                      b_vec(i,j,n4,1) = x
                      b_vec(i,j,n4,2) = y
                      b_vec(i,j,n4,3) = z
                      d2 = rs
                   endif
                 enddo                    ! 5 loop
c============================================================
c   now resest bond for n3  ' 
c=============================================
        do i = 1, nbtot(n3sp,n3at)
          if( (bi(n3sp,n3at,i,1).eq.n1sp).
     1             and.(bi(n3sp,n3at,i,2).eq.n1at)) then 
              n5 = i
              goto 30
          endif
        enddo
30      continue
        bi(n3sp,n3at,n5,1) = isp
        bi(n3sp,n3at,n5,2) = iat 
c========================================================
      i=n3sp
      j =n3at
      ii = isp
      jj = iat
      d2 = 1000
              do lbeta=0,26
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
                      b_vec(i,j,n5,1) = x
                      b_vec(i,j,n5,2) = y
                      b_vec(i,j,n5,3) = z
                      d2 = rs
                   endif
                 enddo                    ! 5 loop
c====================================================================
c=====================done with bond switch  ========================
c===================================================================
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine output(nsp,na,nbtot,bi,b_vec,a_vec,l_vec)
	include 'parameter.h'
	implicit none
	integer nsp
	integer na(nsp_max)
	integer nbtot(nsp_max,na_max)
	integer bi(nsp_max,na_max,nbond_max,3)
	integer i,j,k
	integer ii,jj
	integer m,n,i1,i2,i3 
	integer lbeta
	integer count 
	real d2
	real dist(nsp_max,na_max,na_max)
	real xl(3,30)
        real l_vec(3,3)                                 !lattice vectors
        real a_vec(nsp_max,na_max,3)                    ! a_vec = atom vectors
	real b_vec(nsp_max,na_max,nbond_max,3)
	real x,y,z
	real rs,r,rave,rsq 
	real angle,ang(nsp_max,na_max,nbond_max,nbond_max)
	real atmp,a_ave,a_sdev
	real dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2 
c===========================================================================
c         write atom vectors and neighbor lists
c===========================================================================
	open(unit = 1, file = 'out.dat')
	write(1,*) l_vec(1,1), l_vec(2,2),l_vec(3,3)
	write(1,*) nsp
	write(1,*)  (na(i),i=1,nsp)
	write(1,*) 'A'
	do i = 1,nsp
	   do j = 1,na(i)
	write(1,111) i,j,a_vec(i,j,1),a_vec(i,j,2),a_vec(i,j,3) 
	   enddo
	enddo
111     format(2i4,3f10.5)
        do i = 1,nsp
           do j = 1,na(i)
        write(1,*) i,j, nbtot(i,j),
     1              (bi(i,j,k,1),bi(i,j,k,2),
     2                 bi(i,j,k,3), k = 1,nbtot(i,j) )
           enddo
        enddo
	close(unit = 1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       nearest neighbor file
ccccccccccccccccccccccccccccccccccccccccc
	open(unit = 2,file = 'out.nn')
        do i = 1,nsp
	   count = 0
           do j = 1,na(i)
              do k = 1,nbtot(i,j)
                 rs = (a_vec(i,j,1) - b_vec(i,j,k,1) )**2 +
     1             (a_vec(i,j,2) - b_vec(i,j,k,2) )**2 +
     2             (a_vec(i,j,3) - b_vec(i,j,k,3) )**2
	         rave = rave + sqrt(rs)
	         count = count + 1 
                 write(2,*) i,j,k,sqrt(rs) 
              enddo
           enddo
	   rave = rave/count
         enddo
	do i = 1,nsp
           count = 0
           do j = 1,na(i)
              do k = 1,nbtot(i,j)
                 rs = (a_vec(i,j,1) - b_vec(i,j,k,1) )**2 +
     1             (a_vec(i,j,2) - b_vec(i,j,k,2) )**2 +
     2             (a_vec(i,j,3) - b_vec(i,j,k,3) )**2
                 rsq = rsq + sqrt((rave - sqrt(rs))**2)
                 count = count + 1
              enddo
           enddo
           rsq = rsq/count
         enddo
	write(2,*) 'The average bond length and the rms deviation'
	write(2,*) rave,rsq
	close(unit = 2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	open(unit = 3, file = 'out.xyz')
        do i = 1,nsp
           do j = 1,na(i)
                    x =  a_vec(i,j,1) 
                    y =  a_vec(i,j,2) 
                    z =  a_vec(i,j,3) 
                    write(3,330) x,y,z
           enddo
        enddo
330     format('Si',3f8.4)
        do i = 1,nsp
           do j = 1,na(i)
              do k = 1,nbtot(i,j)
                 if (bi(i,j,k,3).eq.2) then
                    x = ( a_vec(i,j,1) + b_vec(i,j,k,1) ) /2
                    y =  ( a_vec(i,j,2) + b_vec(i,j,k,2) ) /2
                    z =   ( a_vec(i,j,3) + b_vec(i,j,k,3) ) /2
                    write(3,333) x,y,z
                  endif
                enddo
           enddo
        enddo
333     format('O',3f8.4)
	close(unit = 3 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        average bond angle 
cccccccccccccccccccccccccccccccccccc
       count = 0
       do i = 1,nsp
           do j = 1,na(i)
              do k = 1,nbtot(i,j) - 1
                 dx1 = a_vec(i,j,1) - b_vec(i,j,k,1)
                 dy1 = a_vec(i,j,2) - b_vec(i,j,k,2)
                 dz1 = a_vec(i,j,3) - b_vec(i,j,k,3)
                 dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                 do m = k+1,nbtot(i,j)
                    dx2 = a_vec(i,j,1) - b_vec(i,j,m,1)
                    dy2 = a_vec(i,j,2) - b_vec(i,j,m,2)
                    dz2 = a_vec(i,j,3) - b_vec(i,j,m,3)
                    dr2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
          atmp =acos( (dx1*dx2+dy1*dy2+dz1*dz2)/(dr1*dr2) )
	            ang(i,j,k,m) = atmp * 180 / 3.14159
	            a_ave = a_ave + ang(i,j,k,m)   
		    count = count + 1
              enddo
             enddo
          enddo  ! na loop
        enddo    ! nsp loop
	a_ave = a_ave/count
        count = 0
        do i = 1,nsp
           do j = 1,na(i)
              do k = 1,nbtot(i,j) - 1
                 do m = k+1,nbtot(i,j)
		   a_sdev =  a_sdev + sqrt((a_ave - ang(i,j,k,m))**2)                    
                   count = count + 1
                 enddo
             enddo
          enddo  ! na loop
        enddo    ! nsp loop
	a_sdev = a_sdev/count
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           calc angle distribution
ccccccccccccccccccccccccccccccccccccccccccccc
        open(unit = 4, file = 'out.adf')
cccccccccccccccccccccccccccccccccccccccccccc
        angle = 0
        do ii = 1,360   !outer loop
           angle = angle + 5 
           n = 0
           do i = 1,nsp
             do j = 1,na(i)
               do k = 1,nbtot(i,j) - 1
                 do m = k+1,nbtot(i,j)
                   if ( ( ang(i,j,k,m).gt.angle-5).and.
     1                    (ang(i,j,k,m).lt.angle) )   then
                     n = n+1
                   endif
                 enddo
	        enddo
              enddo
            enddo
          write(4,128) angle,n
        enddo      !outer loop
        write(4,*) 'The average bond angle and the rms deviation'
        write(4,*) a_ave,a_sdev
128     format(f10.3,i5)
	close(unit = 4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          calculate radial distribution function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        open(unit = 3, file = 'out.rdf')
cccccccccccccccccccccccccccccccccccccccccc
      do 51 i1=-1,1
      do 52 i2=-1,1
      do 53 i3=-1,1
          lbeta=9*(i3+1)+3*(i2+1)+(i1+1)
      do 55 k=1,3
          xl(k,lbeta)=i1*l_vec(1,k)+i2*l_vec(2,k)+i3*l_vec(3,k)
55    continue
53    continue
52    continue
51    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        d2 = 25
        do i = 1,nsp
           do j = 1,na(i)-1 
              ii = i  
              do jj = j+1, na(i) 
                 do lbeta=0,26
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
		      dist(i,j,jj) = sqrt(rs)
		      if (sqrt(rs).lt.(1.5)) then
	                write(2,*) j,jj,sqrt(rs) 
			write(2,*) j, (bi(i,j,k,2),k=1,nbtot(i,j))
		        write(2,*) jj,(bi(i,jj,k,2),k=1,nbtot(i,jj))
                      endif
		      goto 10
                   endif
                 enddo
10               continue
              enddo
           enddo
        enddo
	r = 0
	do ii = 1,100   !outer loop  
	   r = r + 0.05
	   n = 0
	   do i = 1,nsp
             do j = 1,na(i)-1 
               do jj = j+1, na(i)
                   if ( ( dist(i,j,jj).gt.r-0.05).and.
     1                    (dist(i,j,jj).lt.r))   then
	             n = n+1
		   endif 
                 enddo
              enddo
            enddo
	  write(3,125) r,n
	enddo      !outer loop  
125	format(f6.3,i5)  
	close(unit = 3)
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine 
     1   input(kt,l_vec,nsw,nsp,na,a_vec,ibond,b_vec,bi,nbtot,PE) 
c=========================================================
c  reads input from in.dat
c
c        take care to get in.dat in correct format 
c        there currently no checks for errors
c
c=========================================================
	include 'parameter.h'
	implicit none
	integer ibond
        real kt
        real xlatt,l_vec(3,3)
	integer irand
        integer i,j, nsw, iii, jjj, k 
        integer isp, nsp
        integer na(nsp_max)
        character latt_var*1
        real a_vec(nsp_max,na_max,3)  ! a_vec = atom vectors
        integer nbtot(nsp_max,na_max)                   !number of bonds per atom
        real b_vec(nsp_max,na_max,nbond_max,3)          !nn atom vectors
        integer bi(nsp_max,na_max,nbond_max,3)     
        integer nstot(nsp_max,na_max)   !number of non bonds for each atom
        real s_vec(nsp_max,na_max,nbond_max,3)  !non bond atom vectors
        integer si(nsp_max,na_max,nbond_max,2)   ! index of non bonds
	real x,y,z
	real PE
c==========================================================
	do i = 1,3
	   do j = 1,3
	      l_vec(i,j) = 0
	   enddo
	enddo
c==========================================================
	open(unit = 1,file = 'in.dat')
	read(1,*) irand
	call srand(irand)
	read(1,*) ibond
	if ( (ibond.ne.1).and.(ibond.ne.2)) then
	 write(*,*) 'ibond initially must be 1 or 2 & ibond = ',ibond
	 STOP 
	endif
	read(1,*) kt         
 	read(1,*) nsw
	read(1,*) l_vec(1,1),l_vec(2,2),l_vec(3,3)
	read(1,*) nsp
	read(1,*) (na(isp), isp =1,nsp)
	read(1,*) latt_var 
	if (latt_var.eq.'L') then 
	   do i =1,nsp
	      do j = 1,na(i)
		read(1,*) iii,jjj,x,y,z
                if ( (iii.ne.i).or.(jjj.ne.j)) then
                   STOP  ' Error in read file species index not correct'
                endif
		a_vec(i,j,1) = x * l_vec(1,1) 
		a_vec(i,j,2) = y * l_vec(2,2) 
		a_vec(i,j,3) = z * l_vec(3,3) 
	      enddo
	   enddo
	elseif (latt_var.eq.'A') then
           do i =1,nsp
              do j = 1,na(i)
                read(1,*) iii,jjj,x,y,z
                if ( (iii.ne.i).or.(jjj.ne.j)) then
                   STOP  ' Error in read file species index not correct' 
                endif
                a_vec(i,j,1) = x
                a_vec(i,j,2) = y
                a_vec(i,j,3) = z
              enddo
           enddo
	else
	   STOP 'error in atom flag: it must be L or A'
	endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccc  find nearest & close non second nearest neighbors
        call find_bonds(nsp, na,l_vec,a_vec,b_vec,bi,nbtot,
     1      s_vec,si,nstot,ibond)
ccccccccccccccccccc  calculate initial energy
	call Energy(nsp,na,nbtot,a_vec,b_vec,bi,
     1                                        s_vec,si,nstot,PE)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	close(unit = 1)
	return
	end  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine  find_bonds(nsp,na,l_vec,a_vec,b_vec,bi,nbtot,
     1                                        s_vec,si,nstot,ibond)
c==================================================================
c=======  finds neighbor bonds and their closest image vectors
c==================================================================
        implicit none
	include 'parameter.h'
        real l_vec(3,3)
        real a_vec(nsp_max, na_max, 3) 
        real b_vec(nsp_max, na_max,nbond_max,3)   !closest image positions of bond atom 
        real s_vec(nsp_max, na_max,nbond_max,3)   !closest image positions of shell atom
	real d2, dist(nsp_max,nsp_max,2)
	real rs ,x, y, z 
        real xl(3,30)
	integer nsp
        integer bi(nsp_max, na_max, nbond_max,3)    !integer of bond atom
        integer si(nsp_max, na_max, nbond_max,2)    !integer of shell atom
	integer na(nsp_max) 
        integer i1, i2, i3, lbeta , k
        integer i, j , ii, jj , iii , jjj , kkk, ik, jk , l, il ,jl 
 	integer nbtemp
	integer nbtot(nsp_max,na_max)   ! number of bond atoms
	integer nstot(nsp_max,na_max)   ! number of shell atoms 
	integer ibond
c===================================================================
	real ds2 
c====================================================================
c  important nn dist parameters to be changed as needed  !!!!!!!!!!!!!
c=================================================================
	dist(1,1,1) = 2.8    ! unit = Ang & near. neigh. dist between specy 1 and 1
	dist(1,1,2) = 3.5   
c========================================================================
c    brt -  move above into a  parameterss file
c=======================================================================
      do 51 i1=-1,1
      do 52 i2=-1,1
      do 53 i3=-1,1
          lbeta=9*(i3+1)+3*(i2+1)+(i1+1)
      do 55 k=1,3
          xl(k,lbeta)=i1*l_vec(1,k)+i2*l_vec(2,k)+i3*l_vec(3,k)
55    continue
53    continue
52    continue
51    continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ibond.eq.1) then 
ccccccc flag to find close non second nearest neighbors
      ibond = 3
ccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i = 1,nsp
	   do j = 1,na(i)
	      read(1,*) iii,jjj,kkk,
     1            (bi(i,j,k,1),bi(i,j,k,2),bi(i,j,k,3),k=1,kkk) 
	      if ((i.ne.iii).or.(j.ne.jjj) ) then
	        STOP ' error reading out.nnlist '
	      endif
	      nbtot(i,j) = kkk
	      do k = 1,nbtot(i,j)
                 ii = bi(i,j,k,1) 
		 jj = bi(i,j,k,2)
	         d2 = 1000
	         do lbeta=0,26          
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then
                      b_vec(i,j,k,1) = x
                      b_vec(i,j,k,2) = y
                      b_vec(i,j,k,3) = z
	              d2 = rs
                   endif
                 enddo                  
              enddo                       
           enddo  
        enddo    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     routine to test that bond types match
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     need to add this
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	close (unit =1 )         
      elseif(ibond.eq.2) then
ccccccc flag to find close non second nearest neighbors
      ibond = 3   
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      d2 = 500 
      do i=1,nsp                   ! 1 loop 
         do j = 1,na(i)			! 2 loop 
	    nbtemp =0
	    do ii = 1,nsp			!3 loop 
		do jj = 1,na(ii)                   ! 4 loop 
                 if ( (i.eq.ii).and.(j.eq.jj)) then
	            goto 11
      		 endif
                 do lbeta=0,26        	!	5 loop 
	           x = xl(1,lbeta)+ a_vec(ii,jj,1) 
	           y = xl(2,lbeta)+ a_vec(ii,jj,2)
	           z = xl(3,lbeta)+ a_vec(ii,jj,3) 
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.d2) then 
                      nbtemp = nbtemp + 1
                      b_vec(i,j,nbtemp,1) = x
 	              b_vec(i,j,nbtemp,2) = y
	              b_vec(i,j,nbtemp,3) = z
	              bi(i,j,nbtemp,1) = ii
	              bi(i,j,nbtemp,2) = jj
	              d2 = rs 
	           endif
	         enddo                    ! 5 loop 
11		continue
	        enddo                    ! 4 loop
	    enddo                            ! 3 loop 
         nbtot(i,j) = nbtemp
	 enddo   ! 2 loop  
      enddo     ! 1 loop 
      endif
c=================================== find third nearest neighbors
      if (ibond.eq.3) then
	ds2 = 4.2*4.2 
      do i=1,nsp                   ! 1 loop
         do j = 1,na(i)                 ! 2 loop
            nbtemp =0
            do ii = 1,nsp                       !3 loop
                do jj = 1,na(ii)                   ! 4 loop
                 if ( (i.eq.ii).and.(j.eq.jj)) then
                    goto 12
                 endif
c=======================================================
c               check that ii,jj not nearest neighbor
c=======================================================
	         do k = 1,nbtot(i,j)
	           ik = bi(i,j,k,1)
		   jk = bi(i,j,k,2)
	           if (  (ii.eq.ik).and.(jj.eq.jk)  ) then
	            goto 12 
		   endif
		   do l = 1,nbtot(ik,jk) 
	              il = bi(ik,jk,l,1) 
		      jl = bi(ik,jk,l,2)
		      if (   (ii.eq.il).and.(jj.eq.jl) ) then
			goto 12
		      endif
		   enddo
		 enddo		 
c===========================================================
                 do lbeta=0,26          !       5 loop
                   x = xl(1,lbeta)+ a_vec(ii,jj,1)
                   y = xl(2,lbeta)+ a_vec(ii,jj,2)
                   z = xl(3,lbeta)+ a_vec(ii,jj,3)
                   rs=(a_vec(i,j,1) -x)**2 +
     1                      (a_vec(i,j,2) -y)**2 +
     2                      (a_vec(i,j,3) -z)**2
                   if (rs.le.ds2) then
                      nbtemp = nbtemp + 1
                      s_vec(i,j,nbtemp,1) = x
                      s_vec(i,j,nbtemp,2) = y
                      s_vec(i,j,nbtemp,3) = z
                      si(i,j,nbtemp,1) = ii
                      si(i,j,nbtemp,2) = jj
                      goto 12
                   endif
                 enddo                    ! 5 loop
12               continue
                enddo                    ! 4 loop
            enddo                            ! 3 loop
         nstot(i,j) = nbtemp
        enddo   ! 2 loop
      enddo     ! 1 loop
      else
	 STOP 'ibond can only be 1,2 or 3 '
      endif
	return          	
        end  
