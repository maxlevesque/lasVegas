module histogram
  implicit none
  private
  double precision, public :: h ! distance between 2 bins
  integer, public :: bin(180) ! 100 bins of length h
  public :: init
contains
  subroutine init(len)
    implicit none
    double precision, intent(in) :: len
    double precision :: maxrange
    maxrange = len*sqrt(3.)/2.
    if(size(bin)<=0) stop "size(bin) must not be <=0"
    h = maxrange/size(bin)
  end subroutine init
end module histogram



program lasvegas
  use histogram, only: h, bin, histogram_init=>init
  implicit none
  integer, parameter :: N=128 ! number of molecules
  double precision, parameter :: rO(3)=[0.,0.,0.]
  double precision, parameter :: rH1(3)=[ 0.816495,0.,0.5773525]
  double precision, parameter :: rH2(3)=[-0.816495,0.,0.5773525]
  double precision, parameter :: len=15.6664 ! Supercell length in angstroms
  double precision, parameter :: targetacceptanceratio = 0.4
  integer, parameter :: mcstepmax=10**5 ! number of MonteCarlo steps
  integer :: i, j, mcstep, k, l
  double precision :: rx(N), ry(N), rz(N) ! cartesian coordinates of the center of mass of all N particles
  double precision :: rxH1(N), ryH1(N), rzH1(N), rxH2(N), ryH2(N), rzH2(N)
  double precision, parameter :: eps=0.65, sig=3.166 ! lennard jones
  double precision, parameter :: sig6=sig**6, sig12=sig**12
  double precision :: r2, r6, r12, U, dx, dy, dz, Uij, dr, drmax, rand, du, dub, ratio, dl, r, intensity
  double precision :: rxinew, ryinew, rzinew, rsqrt
  integer :: ntrial, naccpt, ibin
  integer, parameter :: nadjst=10**3
  double precision, parameter :: boltzmanncst = 8.3144621e-3 ! (kJ/mol)/K  (joules per kelvin)
  double precision, parameter :: temperature = 300 ! K
  double precision, parameter :: beta=1./(boltzmanncst*temperature)
  double precision, parameter :: halflen=len/2.
  double precision, parameter :: pi=acos(-1.d0)
  double precision :: q0, q1, q2, q3 ! quaternion used to easily define the rotation matrix
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33

  call print_header
  call random_seed()
  drmax = sig/2.
  ntrial = 0
  naccpt = 0
  call histogram_init(len)

  ! init positions the farther possible in cubic lattice
  l=0
  dl = len/real(nint(N**(1./3.))+1)
  do i=1,nint(N**(1./3.))+1
    do j=1,nint(N**(1./3.))+1
      do k=1,nint(N**(1./3.))+1
        l=l+1
        if (l>N) exit
        rx(l)=(i-1)*dl
        ry(l)=(j-1)*dl
        rz(l)=(k-1)*dl
        if( any([rx(l),ry(l),rz(l)] <0) .or. any([rx(l),ry(l),rz(l)] >=len)) stop "rxinew etc must be >=0 and <len"
      end do
    end do
  end do

  ! init position of hydrogen atoms
  do i=1,N
    rxH1(i)=rx(i)+rH1(1)
    ryH1(i)=ry(i)+rH1(2)
    rzH1(i)=rz(i)+rH1(3)
    rxH2(i)=rx(i)+rH2(1)
    ryH2(i)=ry(i)+rH2(2)
    rzH2(i)=rz(i)+rH2(3)
  end do

  ! Total potential energy of the system
  ! ====================================
  ! 1/ Lennard-Jones
  u = 0
  do i=1,N-1
    do j=i+1,N
      dx = abs(rx(i)-rx(j))
      if( dx>halflen ) dx=len-dx ! pbc
      dy = abs(ry(i)-ry(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rz(i)-rz(j))
      if( dz>halflen ) dz=len-dz
      r2 = dx**2 + dy**2 + dz**2
      r6 = r2**3
      r12 = r6**2
      Uij = 4*eps*( sig12/r12 - sig6/r6 )
      U = U + Uij
    end do
  end do
  ! 2/ Electrostatics


  open(43,file="internal-energy.dat")
  write(43,*) u

  ! Monte Carlo Metropolis steps
  ! ============================
  do mcstep=1,mcstepmax

    ! randomnly select a particle i between 1 and N
    call random_number(rand) ! u real in [0,1[
    i = 1 + FLOOR(N*rand)    ! i integer in {1,...,N}

    ! displace x coordinate of particle i by a random umount, dx, which is given by dx = dR*u, where u is a uniform random variate in [-0.5:0.5]
    ! x
    call random_number(rand)
    dr = drmax*(2*rand-1)
    rxinew   = modulo(rx(i)  +dr,len)
    rxinewh1 = modulo(rxh1(i)+dr,len)
    rxinewh2 = modulo(rxh2(i)+dr,len)
    ! y
    call random_number(rand)
    dr = drmax*(2*rand-1)
    ryinew   = modulo(ry(i)  +dr,len)
    ryinewh1 = modulo(ryh1(i)+dr,len)
    ryinewh2 = modulo(ryh2(i)+dr,len)
    ! z
    call random_number(rand)
    dr = drmax*(2*rand-1)
    rzinew   = rz(i)  +dr
    rzinewh1 = rzh1(i)+dr
    rzinewh2 = rzh2(i)+dr

    ! rotate molecule by a random quantity
    ! ====================================
    call generate_random_quaternions (q0,q1,q2,q3)
    ! rotation matrix
    a11=q0**2+q1**2-q2**2-q3**2; a12=2*(q1*q2+q0*q3)        ; a13=2*(-q0*q2+q1*q3)
    a21=2*(q1*q2-q0*q3)        ; a22=q0**2-q1**2+q2**2-q3**2; a23=2*(q0*q1+q2*q3)
    a31=2*(q0*q2+q1*q3)        ; a32=2*(-q0*q1+q2*q3)       ; a33=q0**2-q1**2-q2**2+q3**2
    ! apply rotation matrix
    rxinewrot   = a11*rxinew + a12*ryinew + a13*rzinew
    ryinewrot   = a21*rxinew + a22*ryinew + a23*rzinew
    rzinewrot   = a31*rxinew + a32*ryinew + a33*rzinew
    rxinew = modulo(rxinewrot,len)
    ryinew = modulo(ryinewrot,len)
    rzinew = modulo(rzinewrot,len)
    rxinewh1rot = a11*rxinewh1 + a12*ryinewh1 + a13*rzinewh1
    ryinewh1rot = a21*rxinewh1 + a22*ryinewh1 + a23*rzinewh1
    rzinewh1rot = a31*rxinewh1 + a32*ryinewh1 + a33*rzinewh1
    rxinewh1 = modulo(rxinewh1rot,len)
    ryinewh1 = modulo(ryinewh1rot,len)
    rzinewh1 = modulo(rzinewh1rot,len)
    rxinewh2rot = a11*rxinewh2 + a12*ryinewh2 + a13*rzinewh2
    ryinewh2rot = a21*rxinewh2 + a22*ryinewh2 + a23*rzinewh2
    rzinewh2rot = a31*rxinewh2 + a32*ryinewh2 + a33*rzinewh2
    rxinewh2 = modulo(rxinewh2rot,len)
    ryinewh2 = modulo(ryinewh2rot,len)
    rzinewh2 = modulo(rzinewh2rot,len)

    ! compute variation in energy
    du = 0
    do j=1,N
      ! add new contribution to LJ
      if( j==i ) cycle
      dx = abs(rxinew-rx(j))
      if( dx>halflen ) dx=len-dx
      dy = abs(ryinew-ry(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rzinew-rz(j))
      if( dz>halflen ) dz=len-dz
      r2 = dx**2 + dy**2 + dz**2
      r6 = r2**3
      r12 = r6**2
      Uij = 4*eps*( sig12/r12 - sig6/r6 )
      du = du + Uij
      ! remove old contribution
      dx = abs(rx(i)-rx(j))
      if( dx>halflen ) dx=len-dx
      dy = abs(ry(i)-ry(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rz(i)-rz(j))
      if( dz>halflen ) dz=len-dz
      r2 = dx**2 + dy**2 + dz**2
      r6 = r2**3
      r12 = r6**2
      Uij = 4*eps*( sig12/r12 - sig6/r6 )
      du = du -Uij
    end do

    dub = du * beta
    if( dub <= 0 ) then ! always accept
      ! print*,"dub<=0  =>accpt"
      u = u + du
      rx(i) = rxinew
      ry(i) = ryinew
      rz(i) = rzinew
      naccpt = naccpt +1
    else
      call random_number(rand)
      if( exp(-dub) > rand ) then
        u = u + du
        rx(i) = rxinew
        ry(i) = ryinew
        rz(i) = rzinew
        naccpt = naccpt +1
      end if
    end if

    if(modulo(ntrial,100)==0) write(43,*) u

    ntrial = ntrial + 1
    if( modulo(ntrial,nadjst)==0 ) then
      ratio = real(naccpt)/real(nadjst)
      if( ratio > targetacceptanceratio ) then
        drmax = drmax * 1.05
      else
        drmax = drmax * 0.95
      end if
      naccpt = 0
      print*,"progress(%)",floor(real(ntrial)/real(mcstepmax)*100)
      if(abs((ratio-targetacceptanceratio)/targetacceptanceratio)>0.2) print*,"WARNING: accpt ratio=",&
        real(ratio,4),"only. Target=>0.3"
    end if

    ! histogram for the center of mass
    if(ntrial > mcstepmax/10) then ! trick to not accumulate stats before "melting". Much smarter tricks could be used.
      do i=1,N
        do j=1,N
          if (j/=i) then
            dx = abs(rx(i)-rx(j))
            if( dx>halflen ) dx=len-dx
            dy = abs(ry(i)-ry(j))
            if( dy>halflen ) dy=len-dy
            dz = abs(rz(i)-rz(j))
            if( dz>halflen ) dz=len-dz
            r2 = dx**2 + dy**2 + dz**2
            rsqrt = sqrt(r2)
            ibin = int(rsqrt/h) +1
            bin(ibin) = bin(ibin)+1
          end if
        end do
      end do
    end if

  end do ! mcstep

  close(43) ! internal energy


  ! print radial distribution function
  ! ==================================
  block
    double precision :: const, rupper, rlower, rho, nideal, gr(size(bin)), grsm(size(bin))
    integer :: imax
    rho = real(N)/len**3
    const = 4.*pi*rho/3.
    open(77,file="rdf.dat")
    imax=0
    do i=1,size(bin)
      r = (i-0.5)*h
      if( r > len/2.) cycle
      if(imax<i) imax=i
      rlower = real(i-1)*h
      rupper = rlower+h
      nideal = const *(rupper**3-rlower**3)
      gr(i) = real(bin(i)) / real(ntrial) / real(N) / nideal
      write(77,*) r, gr(i)
    end do
    close(77)
    ! write a smoothed g(r) see Tildesley p.204
    ! third degree, five point smoothing
    do i=1,imax
      if(i==1) then
        grsm(i)=(69*gr(i)+4* gr(i+1)-6*gr(i+2)+4*gr(i+3)-gr(i+4))/70.
      else if(i==2) then
        grsm(i)=(2* gr(i-1)+27*gr(i)+12*gr(i+1)-8*gr(i+2)+2*gr(i+3))/35.
      else if(i>2 .and. i<imax-1) then
        grsm(i)=(-3*gr(i-2)+12*gr(i-1)+17*gr(i)+12*gr(i+1)-3*gr(i+2))/35.
      else if(i==imax-1) then
        grsm(i)=(2* gr(i+1)+27*gr(i)+12*gr(i-1)-8*gr(i-2)+2*gr(i-3))/35.
      else if(i==imax) then
        grsm(i)=(69*gr(i)+4*gr(i-1)-6*gr(i-2)+4*gr(i-3)-gr(i-4))/70.
      end if
    end do
    open(78,file="smoothed-rdf.dat")
    do i=1,imax
      r=(i-0.5)*h
      write(78,*)r,grsm(i)
    end do
    close(78)
  end block

  ! print positions
  open(39,file="pos.dat",form="unformatted")
  write(39)rx,ry,rz
  close(39)




contains

  ! Generate random quaternion
  ! Method by G. Marsaglia, Choosing a point from the surface of a sphere, Ann. Math. Stat. 43, 645–646 (1972).
  subroutine generate_random_quaternions (q0, q1, q2, q3)
    implicit none
    double precision, intent(out) :: q0, q1, q2, q3
    double precision :: x1, x2, y1, y2, s1, s2, t
    s1=10
    do while(s1>1)
      call random_number(x1)
      call random_number(y1)
      s1=x1**2+y1**2
    end do
    s2=10
    do while(s2>1)
      call random_number(x2)
      call random_number(y2)
      s2=x2**2+y2**2
    end do
    q0=x1
    q1=y1
    t=sqrt((1-s1)/s2)
    q2=x2*t
    q3=y2*t
  end subroutine generate_random_quaternions

end program lasvegas
