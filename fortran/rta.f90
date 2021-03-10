! ## File: rta.f90
! ## - main program: running the RTA quasistationary method to simulate SIS dynamics,
! ## See README.md for more information and use
!-----------------------------------------------------------------------------
! SIS epidemic model algorithm based on the article
!           "Simple quasistationary method for simulations of epidemic
!           processes with localized states"
! Copyright (C) 2021 Guilherme S. Costa, Silvio C. Ferreira
!
! Please cite the above cited paper ()
! as reference to our code.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
! Author    : Guilherme S. Costa
! Email     : guilherme.h.costa@ufv.br
! Date      : 06 Mar 2021
! Version   : 1.0
!-----------------------------------------------------------------------------
! See README.md for more details

program sis
implicit none

    !INPUTS
    real*8 :: lb				!infection rate
    integer :: t_ave				!averaging time
    integer :: t_rlx				!relaxation time
    integer :: iseed				!initial seed
    character*199 :: file_net			!network file

    !NETWORK VARIABLES
    integer :: n_stubs,n_vert,n_con,deg_min,deg_max
    integer :: a,b,iErr,ax,vt_in,vt_out
    integer, allocatable :: adj(:),ini(:),deg(:),deg2(:),edg(:,:)

    !RAT METHOD VARIABLES
    integer :: box,sinf(100),t_max(100)
    integer, parameter :: n_boxes = 100
    integer, allocatable :: t_i(:),tbox(:),boxes(:,:),lab_boxes(:)
    real*8, allocatable :: tn_i(:)
    real*8 :: sinft,total_t

    !QS QUANTITIES
    real*8 :: rho, chi, tau,norm,rho2
    real*8, allocatable :: his(:)

    !SIS SIMULATION VARIABLES
    integer*1, allocatable :: sigma(:)
    integer, allocatable :: list(:)
    integer :: tser,cont,edg_inf
    real*8 :: dt,t,p

    !AUXILIAR VARIABLES AND COUNTERS
    integer :: i,j,rnd
    real*8 ::z1,z2,z

    read*,lb,t_ave,t_rlx,file_net,iseed

    call read_net()

    print*,"Read the network"

	allocate(sigma(n_vert),list(n_vert))
	allocate(his(0:n_vert),t_i(n_vert),tn_i(n_vert))
	allocate(boxes(n_boxes,n_vert/n_boxes),tbox(n_boxes),lab_boxes(n_vert))

	call build_boxes()
	call evolve_sis()
	call calc_qs()

    print*,"Finished simulation"

 contains

subroutine evolve_sis()

    call initial_condition()

    tser = 1
    !####################RELAXATION###########################
    do while ( t .le. t_rlx)

        p = 1d0*cont/(1d0*cont+1d0*lb*edg_inf)
        z = int(ran2(iseed)*cont + 1)
        z1 = ran2(iseed)

        if(z1 .le. p) then
            sigma(list(z)) = 0
			edg_inf = edg_inf - deg(list(z))

			t_i(list(z)) = t_i(list(z)) + 1
			box = lab_boxes(list(z))
			t_max(box) = max(t_max(box),t_i(list(z)))
			sinft = sinft+1
			sinf(box) = sinf(box) + 1

			list(z) = list(cont)
			cont = cont -1

			if(cont == 0) call qs_method_rta()
        else
            do
                z = int(ran2(iseed)*cont +1d0)
                z2 = ran2(iseed)
                if(z2 .le. 1d0*deg(list(z))/deg_max) exit
            end do
            rnd = int(ran2(iseed)*deg(list(z)))
            if(sigma(adj(ini(list(z))+rnd)) == 0) then
                sigma(adj(ini(list(z))+ rnd)) = 1
                cont = cont + 1
                list(cont) = adj(ini(list(z))+rnd)
                edg_inf = edg_inf + deg(adj(ini(list(z))+rnd))
            endif
		endif

        if(t .ge. tser) then
            if(mod(int(tser),t_rlx/10) == 0)  print*,"Time elapsed (Relaxation)::", t
			tser = tser+1
        endif
        dt = 1d0*p/cont
        t = t + dt
        total_t = t
	end do

    t = 0 ; tser = 1
    !####################RELAXATION###########################
    do while ( t .le. t_ave)

	    p = 1d0*cont/(cont + lb*edg_inf)
        z = int(ran2(iseed)*cont + 1)
        z1 = ran2(iseed)


        if(z1 .le. p) then
			sigma(list(z)) = 0
			edg_inf = edg_inf - deg(list(z))

			t_i(list(z)) = t_i(list(z)) + 1
			box = lab_boxes(list(z))
			t_max(box) = max(t_max(box),t_i(list(z)))
			sinft = sinft+1
			sinf(box) = sinf(box) + 1

			list(z) = list(cont)
			cont = cont -1

            if(cont == 0) call qs_method_rta()
        else

			do
				z = int(ran2(iseed)*cont +1d0)
                z2 = ran2(iseed)
                if(z2 .le. 1d0*deg(list(z))/deg_max) exit
            end do

            rnd = int(ran2(iseed)*deg(list(z)))
            if(sigma(adj(ini(list(z))+rnd)) == 0) then
                sigma(adj(ini(list(z))+ rnd)) = 1
                cont = cont + 1
                list(cont) = adj(ini(list(z))+rnd)
                edg_inf = edg_inf + deg(adj(ini(list(z))+rnd))
            endiF

        endif
        if(t .ge. tser) then

			tser = tser+1
			if(mod(int(tser),t_ave/10) == 0)  print*,"Time elapsed (Average):", tser
        endif

        dt = 1d0*p/cont
        t = t + dt
        total_t = total_t + dt
        his(cont) = his(cont) + dt

    end do

end subroutine

!SUBROUTINE TO CALCULATE THE QUASISTATIONARY QUANTITIES

subroutine calc_qs()
    norm = sum(his)
	his = 1d0*his/norm

    open(121,file='pn.dat')

	do i = 1,n_vert
		rho =  rho + 1d0*i*his(i)
		rho2 = rho2 + 1d0*i*i*his(i)
		write(121,*)i,his(i)
	end do
    close(121)

	tau = 1d0/his(1)
	rho = 1d0*rho/n_vert
	rho2 = 1d0*rho2/(1d0*n_vert)/(1d0*n_vert)
	chi = 1d0*n_vert*(rho2-rho*rho)/rho

	if(tau .ge. 1e10) tau = t_ave

	print*,"Density of infected nodes:", rho
	print*,"Dynamical susceptibility:", chi
	print*,"Lifespam:", tau

end subroutine
!###########################################################

!SET INITIAL CONDITION
subroutine initial_condition()
    rho = 0 ; rho2 = 0 ; his = 0 ;
    cont = 0 ; edg_inf = 0 ; sigma = 0 ;
    sinft =0 ; sinf = 0;

    do i = 1,n_vert
        sigma(i)=1
        cont=cont+1
        list(cont) = i
        edg_inf = edg_inf + deg(i)
    end do

    t = 0 ; t_i = 0 ; t_max = 0 ;

 end subroutine
!##############################################################


subroutine qs_method_rta()
	real*8 :: propo,pdeci,time1,time2
	integer ::nativ,n_reativ,tgt,maxtbox

	maxtbox = maxval(sinf)
	propo = 1d0*sinft/total_t
	pdeci = propo - int(propo)
	n_reativ = propo
	z1 = ran2(iseed)

	if ( z1 .le. pdeci) n_reativ = n_reativ + 1

	cont = 0

	do nativ = 1,n_reativ

	      do
                z = int(ran2(iseed)*n_boxes +1)
                z1 = ran2(iseed)
                if(z1 .le. 1d0*sinf(z)/maxtbox) exit
	      end do

	      do
                z1 = int(ran2(iseed)*n_vert/n_boxes + 1)
                z2 = ran2(iseed)
                tgt = boxes(z,z1)
                if(z2 .le. 1d0*t_i(tgt)/t_max(z) .and. sigma(tgt) == 0) exit
	      end do

	      sigma(tgt) = 1
	      cont = cont+1
	      list(cont) = tgt
	      edg_inf = edg_inf + deg(tgt)

	end do

end subroutine

subroutine build_boxes()

	do j = 1,n_boxes
		do i = 1,n_vert/n_boxes

			boxes(j,i) = (n_vert/n_boxes)*(j-1) + i
			lab_boxes((n_vert/n_boxes)*(j-1) +i) = j

		end do
	end do

end subroutine

!READ THE NETWORK FILE AND BUILD THE ADJACENCY AND AUXILIARY LISTS
subroutine read_net()

    open(100,file = file_net)

    n_stubs = 0 ; n_con = 0
    ax=0

    do
        read(100,*,IOSTAT = iErr)a,b

        n_con = n_con + 1
        if(a .gt. ax)ax = a
        if(b .gt. ax)ax = b

        if(iErr .ne. 0)exit
    end do
    close(100)
    n_con = n_con -1
    n_vert = ax

    allocate(ini(n_vert),deg(n_vert),edg(n_con,2),deg2(n_vert))
    allocate(adj(n_con))

    ini = 0 ; deg = 0 ; edg = 0 ; adj = 0

    open(100,file = file_net)

    !BUILDING THE DEGREES MATRIX

    do i = 1,n_con
        read(100,*)edg(i,1),edg(i,2)
        deg(edg(i,1)) = deg(edg(i,1))+1
    end do
    deg2 = deg

    !BUILDING THE INI MATRIX
    ini(1) = 1
    do i = 2,n_vert
        ini(i) = ini(i-1)+deg(i-1)
    end do

    !BUILDING THE ADJACENCY LIST

    do i = 1,n_con
        vt_in = edg(i,1)
        vt_out= edg(i,2)
        if (vt_out .gt. vt_in)then
            adj(ini(vt_in)+deg2(vt_in)-1) = vt_out
            adj(ini(vt_out)+deg2(vt_out)-1) = vt_in
            deg2(vt_in) = deg2(vt_in)-1
            deg2(vt_out)= deg2(vt_out)-1
        endif
    end do
    deg_min = minval(deg)
    deg_max = maxval(deg)
    deallocate(deg2)
end subroutine
!###########################################################


!PSEUDO NUMBER GENERATOR

FUNCTION ran2(idum)
    INTEGER :: idum
    REAL*8 ::ran2
    Integer, PARAMETER ::IM1=2147483563,IM2=2147483399,IMM1=IM1-1,&
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
         NTAB=32,NDIV=1+IMM1/NTAB
    real*8,parameter   ::EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1
    INTEGER            ::idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
END Function ran2
!############################################################################

end program
