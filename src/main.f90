module histogram
  implicit none
  private
  double precision, public :: h ! distance between 2 bins
  integer, public :: bin(200) ! 100 bins of length h
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
  integer, parameter :: N=512
  double precision, parameter :: len=24.8689 ! Angstrom , box length
  integer, parameter :: mcstepmax=10**5
  integer :: i, j, mcstep, k, l
  double precision :: rx(N), ry(N), rz(N) ! positions of the N particles
  double precision, parameter :: eps=0.65, sig=3.166 ! lennard jones
  double precision, parameter :: sig6=sig**6, sig12=sig**12
  double precision :: r2, r6, r12, U, dx, dy, dz, Uij, drmax, rand, du, dub, ratio, dl, r, intensity
  double precision :: rxinew, ryinew, rzinew, rsqrt
  integer :: ntrial, naccpt, ibin
  integer, parameter :: nadjst=10**3
  double precision, parameter :: boltzmanncst = 8.3144621e-3 ! (kJ/mol)/K  (joules per kelvin)
  double precision, parameter :: temperature = 300 ! K
  double precision, parameter :: beta=1./(boltzmanncst*temperature)
  double precision, parameter :: halflen=len/2.

  call print_header
  call random_seed()
  drmax = sig/2.
  ntrial = 0
  naccpt = 0
  call histogram_init(len)

  ! init positions the farther possible
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
        IF(ANY([rx(l),ry(l),rz(l)]<0).or.ANY([rx(l),ry(l),rz(l)]>=len)) stop"rxinew etc must be >=0 and <len"
      end do
    end do
  end do

  ! total potential energy of the system
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

  open(43,file="internal-energy.dat")
  write(43,*) u

  do mcstep=1,mcstepmax

    ! randomnly select a particle i between 1 and N
    call random_number(rand) ! u real in [0,1[
    i = 1 + FLOOR(N*rand)    ! i integer in {1,...,N}

    ! displace x coordinate of particle i by a random umount, dx, which is given by dx = dR*u, where u is a uniform random variate in [-0.5:0.5]
    call random_number(rand)
    rxinew = modulo( rx(i) + drmax*(2*rand-1) ,len)
    call random_number(rand)
    ryinew = modulo( ry(i) + drmax*(2*rand-1) ,len)
    call random_number(rand)
    rzinew = modulo( rz(i) + drmax*(2*rand-1) ,len)

    IF(ANY([rxinew,ryinew,rzinew]<0).or.ANY([rxinew,ryinew,rzinew]>=len)) stop"rxinew etc must be >=0 and <len"

    ! compute variation in energy
    du = 0
    do j=1,N
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
      du = du +(-Uij)
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
      if( ratio > 0.3 ) then
        drmax = drmax * 1.05
      else
        drmax = drmax * 0.95
      end if
      naccpt = 0
      print*,"progress(%)",floor(real(ntrial)/real(mcstepmax)*100)
      if(abs((ratio-0.3)/0.3)>0.2) print*,"WARNING: accpt ratio=",real(ratio,4),"only. Target=>0.3"
    end if

    ! histogram
    if(ntrial > mcstepmax/10) then
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
  open(77,file="rdf.dat")
  do i=1,size(bin)
    r = (i-0.5)*h
    if( r <= len/2.) then
      intensity = real(bin(i))/real(ntrial)
      r2 = r**2
      write(77,*) r, intensity/(4*acos(-1.)*r2)
    end if
  end do
  close(77)

end program lasvegas
