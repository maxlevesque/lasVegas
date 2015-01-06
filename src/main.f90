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
  use ewald, only: ewald_init, ewald_energy
  implicit none
  integer, parameter :: N=128 ! number of molecules
  integer, parameter :: mcstepmax=10**5 ! number of MonteCarlo steps
  double precision, parameter :: len=15.6664 ! Supercell length in angstroms
  double precision, parameter :: targetacceptanceratio = 0.4
  integer :: i, j, mcstep, k, l
  double precision :: rxO(N), ryO(N), rzO(N), rxH1(N), ryH1(N), rzH1(N), rxH2(N), ryH2(N), rzH2(N) ! cartesian coordinates of all N molecules
  ! WATER MODEL: SPC/E by Berendsen et al. JCP 1987
  double precision, parameter :: rO_molecularframe(3)=[0.,0.,0.]
  double precision, parameter :: rH1_molecularframe(3)=[ 0.816495,0.,0.5773525]
  double precision, parameter :: rH2_molecularframe(3)=[-0.816495,0.,0.5773525]
  double precision, parameter :: qH=0.42380, qO=-2.*qH
  double precision, parameter :: eps=0.65 ! kJ/mol
  double precision, parameter :: sig=3.166 ! Angstroms

  double precision, parameter :: sig6=sig**6, sig12=sig**12
  double precision :: r6, r12, dx, dy, dz, Uij, dr, drmax, rand, rand3(3), du, dub, ratio, r
  double precision :: rO(3), rH1(3), rH2(3), rO_tr(3), rH1_tr(3), rH2_tr(3), rO_tr_rot(3), rH1_tr_rot(3), rH2_tr_rot(3)
  integer :: ntrial, naccpt, ibin
  integer, parameter :: nadjst=10**3
  double precision, parameter :: temperature = 300 ! K
  double precision, parameter :: boltzmanncst = 8.3144621e-3 ! kJ/(mol K)
  double precision, parameter :: beta=1./(boltzmanncst*temperature)
  double precision, parameter :: halflen=len/2.
  ! CODATA 2010 physical constants
  double precision, parameter :: pi=acos(-1.d0)
  double precision :: q0, q1, q2, q3 ! quaternion used to easily define the rotation matrix
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33
  double precision :: q(3*N), rx(3*N), ry(3*N), rz(3*N)
  double precision :: Uewald, Ulj, Utot, UewaldBEFORE, UewaldAFTER, dUewald


  call print_header
  call random_seed()
  drmax = sig/2.
  ntrial = 0
  naccpt = 0
  call histogram_init(len)
  call positions_init


  ! Total potential energy of the system
  ! ====================================
  ! 1/ Lennard-Jones
  Ulj = 0
  do i=1,N-1
    do j=i+1,N
      dx = abs(rxO(i)-rxO(j))
      if( dx>halflen ) dx=len-dx ! pbc
      dy = abs(ryO(i)-ryO(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rzO(i)-rzO(j))
      if( dz>halflen ) dz=len-dz
      r6 = (dx**2 + dy**2 + dz**2)**3
      r12 = r6**2
      Ulj = Ulj + 4*eps*( sig12/r12 - sig6/r6 )
    end do
  end do
  Ulj = Ulj ! /kb ?
print*,"Ulj initial=",Ulj

  ! 2/ Electrostatics
  call fill_q
  call fill_positions
  call ewald_init (len)
  call ewald_energy (N,len,q,rx,ry,rz,Uewald)
print*,"Uewald initial =",Uewald

  open(43,file="internal-energy.dat")
  write(43,*) Ulj+Uewald, Ulj, Uewald

  ! Monte Carlo Metropolis steps
  ! ============================
  do mcstep=1,mcstepmax

    ! randomnly select a particle i between 1 and N
    call random_number(rand) ! u real in [0,1[
    i = 1 + floor(N*rand)    ! i integer in {1,...,N}
    rO  = [rxO(i),ryO(i),rzO(i)]
    rH1 = [rxH1(i),ryH1(i),rzH1(i)]
    rH2 = [rxH2(i),ryH2(i),rzH2(i)]

    ! displace x coordinate of particle i by a random umount, dx
    ! dx = drmax*u, where u is a uniform random variate in [-0.5:0.5]
    call random_number(rand3)
    rO_tr  = rO  + drmax*(2*rand3-1)! _tr means translated

    ! rotate molecule by a random quantity
    ! ====================================
    call generate_random_quaternion (q0,q1,q2,q3)
    ! rotation matrix
    a11=q0**2+q1**2-q2**2-q3**2; a12=2*(q1*q2+q0*q3)        ; a13=2*(-q0*q2+q1*q3)
    a21=2*(q1*q2-q0*q3)        ; a22=q0**2-q1**2+q2**2-q3**2; a23=2*(q0*q1+q2*q3)
    a31=2*(q0*q2+q1*q3)        ; a32=2*(-q0*q1+q2*q3)       ; a33=q0**2-q1**2-q2**2+q3**2
    ! apply rotation matrix (1/ translate center of mass to origin of lab frame, 2/ apply rotation matrix, 3/ translate back to original position)
    rO_tr_rot     = modulo( rO_tr, len)
    rH1_tr_rot(1) = a11*rH1(1) + a12*rH1(2) + a13*rH1(3)
    rH1_tr_rot(2) = a21*rH1(1) + a22*rH1(2) + a23*rH1(3)
    rH1_tr_rot(3) = a31*rH1(1) + a32*rH1(2) + a33*rH1(3)
    rH2_tr_rot(1) = a11*rH2(1) + a12*rH2(2) + a13*rH2(3)
    rH2_tr_rot(2) = a21*rH2(1) + a22*rH2(2) + a23*rH2(3)
    rH2_tr_rot(3) = a31*rH2(1) + a32*rH2(2) + a33*rH2(3)

    ! compute variation in energy
    du = 0
    do j=1,N
      ! add new contribution to LJ
      if( j==i ) cycle
      dx = abs(rO_tr_rot(1)-rxO(j))
      if( dx>halflen ) dx=len-dx
      dy = abs(rO_tr_rot(2)-ryO(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rO_tr_rot(3)-rzO(j))
      if( dz>halflen ) dz=len-dz
      r6 = (dx**2 + dy**2 + dz**2)**3
      r12 = r6**2
      Uij = 4*eps*( sig12/r12 - sig6/r6 )
      du = du + Uij
      ! remove old contribution
      dx = abs(rxO(i)-rxO(j))
      if( dx>halflen ) dx=len-dx
      dy = abs(ryO(i)-ryO(j))
      if( dy>halflen ) dy=len-dy
      dz = abs(rzO(i)-rzO(j))
      if( dz>halflen ) dz=len-dz
      r6 = (dx**2 + dy**2 + dz**2)**3
      r12 = r6**2
      Uij = 4*eps*( sig12/r12 - sig6/r6 )
      du = du -Uij
    end do

    j=0
    do l=1,N
      j=j+1
      rx(j)=rxO(l)
      ry(j)=ryO(l)
      rz(j)=rzO(l)
      j=j+1
      rx(j)=rxO(l)+rxH1(l)
      ry(j)=ryO(l)+ryH1(l)
      rz(j)=rzO(l)+rzH1(l)
      j=j+1
      rx(j)=rxO(l)+rxH2(l)
      ry(j)=ryO(l)+ryH2(l)
      rz(j)=rzO(l)+rzH2(l)
    end do
    call ewald_energy (N,len,q,rx,ry,rz,UewaldBEFORE)

    j=0
    do l=1,N
      if (l/=i) then
        j=j+1
        rx(j)=rxO(l)
        ry(j)=ryO(l)
        rz(j)=rzO(l)
        j=j+1
        rx(j)=rxO(l)+rxH1(l)
        ry(j)=ryO(l)+ryH1(l)
        rz(j)=rzO(l)+rzH1(l)
        j=j+1
        rx(j)=rxO(l)+rxH2(l)
        ry(j)=ryO(l)+ryH2(l)
        rz(j)=rzO(l)+rzH2(l)
      else if(l==i) then
        j=j+1
        rx(j)=rO_tr_rot(1)
        ry(j)=rO_tr_rot(2)
        rz(j)=rO_tr_rot(3)
        j=j+1
        rx(j)=rO_tr_rot(1)+rH1_tr_rot(1)
        ry(j)=rO_tr_rot(2)+rH1_tr_rot(2)
        rz(j)=rO_tr_rot(3)+rH1_tr_rot(3)
        j=j+1
        rx(j)=rO_tr_rot(1)+rH2_tr_rot(1)
        ry(j)=rO_tr_rot(2)+rH2_tr_rot(2)
        rz(j)=rO_tr_rot(3)+rH2_tr_rot(3)
      end if
    end do
    call ewald_energy (N,len,q,rx,ry,rz,UewaldAFTER)

    dub = (du+UewaldAFTER-UewaldBEFORE) * beta
    ! Metropolis algorithm
    if( dub <= 0 ) then ! always accept
      ! print*,"dub<=0  =>accpt"
      Utot = Utot + du
      rxO(i) = rO_tr_rot(1)
      ryO(i) = rO_tr_rot(2)
      rzO(i) = rO_tr_rot(3)
      rxH1(i) = rH1_tr_rot(1)
      ryH1(i) = rH1_tr_rot(2)
      rzH1(i) = rH1_tr_rot(3)
      rxH2(i) = rH2_tr_rot(1)
      ryH2(i) = rH2_tr_rot(2)
      rzH2(i) = rH2_tr_rot(3)
      naccpt = naccpt +1
    else
      call random_number(rand)
      if( exp(-dub) > rand ) then
        Utot = Utot + du
        rxO(i) = rO_tr_rot(1)
        ryO(i) = rO_tr_rot(2)
        rzO(i) = rO_tr_rot(3)
        rxH1(i) = rH1_tr_rot(1)
        ryH1(i) = rH1_tr_rot(2)
        rzH1(i) = rH1_tr_rot(3)
        rxH2(i) = rH2_tr_rot(1)
        ryH2(i) = rH2_tr_rot(2)
        rzH2(i) = rH2_tr_rot(3)
        naccpt = naccpt +1
      end if
    end if

    if(modulo(ntrial,100)==0) write(43,*) Utot, Ulj, Uewald

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

    ! print positions once system has been equilibrated mcstepmax/10 times
    if(ntrial == mcstepmax/10) then
      open(38,file="positions_equilibrated.dat")
      write(38,*) 3*N ! N, which is an integer parameter that should not be read
      write(38,*)"configuration after equilibration. Box length is",len
      do i=1,N ! 3 sites per molecule
        write(38,*)"O",rxO(i),ryO(i),rzO(i)
        write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH1(i),ryH1(i),rzH1(i)]
        write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH2(i),ryH2(i),rzH2(i)]
      end do
      close(38)
      print*,"positions_equilibrated.dat written"
    end if

    ! update histogram for the center of mass
    if(ntrial > mcstepmax/10) then ! trick to not accumulate stats before "melting". Much smarter tricks could be used.
      do i=1,N
        do j=1,N
          if (j/=i) then
            dx = abs(rxO(i)-rxO(j))
            if( dx>halflen ) dx=len-dx
            dy = abs(ryO(i)-ryO(j))
            if( dy>halflen ) dy=len-dy
            dz = abs(rzO(i)-rzO(j))
            if( dz>halflen ) dz=len-dz
            ibin = int(sqrt(dx**2 + dy**2 + dz**2)/h) +1
            bin(ibin) = bin(ibin)+1
          end if
        end do
      end do
    end if

  end do ! mcstep

  close(43) ! internal energy
  call print_final_positions
  call rdf_print



contains
  subroutine print_final_positions
    implicit none
    open(38,file="positions_final.dat")
    write(38,*) 3*N ! N, which is an integer parameter that should not be read
    write(38,*)"final configuration. Box length is",len
    do i=1,N ! 3 sites per molecule
      write(38,*)"O",rxO(i),ryO(i),rzO(i)
      write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH1(i),ryH1(i),rzH1(i)]
      write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH2(i),ryH2(i),rzH2(i)]
    end do
    close(38)
    print*,"positions_final.dat written"
  end subroutine print_final_positions

  ! Generate random quaternion
  ! Method by G. Marsaglia, Choosing a point from the surface of a sphere, Ann. Math. Stat. 43, 645–646 (1972).
  subroutine generate_random_quaternion (q0, q1, q2, q3)
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
  end subroutine generate_random_quaternion

  ! if file positions_init.dat exists, then read positions from this file.
  ! else, start with cubic lattice
  subroutine positions_init
    implicit none
    double precision :: dl
    logical :: file_exists
    integer :: i,j,k,l
    character(len=1) :: string
    inquire(file="positions_init.dat", EXIST=file_exists)
    if( file_exists ) then
      open(38,file="positions_init.dat")
      read(38,*)i ! N, which is an integer parameter that should not be read
      read(38,*)string
      do i=1,N ! 3 sites per molecule
        read(38,*)string,rxO(i),ryO(i),rzO(i)
        read(38,*)string,rxH1(i),ryH1(i),rzH1(i)
        rxH1(i) = rxH1(i)-rxO(i)
        ryH1(i) = ryH1(i)-ryO(i)
        rzH1(i) = rzH1(i)-rzO(i)
        read(38,*)string,rxH2(i),ryH2(i),rzH2(i)
        rxH2(i) = rxH2(i)-rxO(i)
        ryH2(i) = ryH2(i)-ryO(i)
        rzH2(i) = rzH2(i)-rzO(i)
      end do
      close(38)
    else
      ! init positions the farther possible in cubic lattice
      l=0
      dl = len/real(nint(N**(1./3.))+1)
      do i=1,nint(N**(1./3.))+1
        do j=1,nint(N**(1./3.))+1
          do k=1,nint(N**(1./3.))+1
            l=l+1
            if (l>N) exit
            rxO(l)=(i-1)*dl
            ryO(l)=(j-1)*dl
            rzO(l)=(k-1)*dl
            if( any([rxO(l),ryO(l),rzO(l)] <0) .or. any([rxO(l),ryO(l),rzO(l)] >=len)) stop "rxinew etc must be >=0 and <len"
          end do
        end do
      end do
      ! init position of hydrogen atoms
      do i=1,N
        rxH1(i) = rH1_molecularframe(1)
        ryH1(i) = rH1_molecularframe(2)
        rzH1(i) = rH1_molecularframe(3)
        rxH2(i) = rH2_molecularframe(1)
        ryH2(i) = rH2_molecularframe(2)
        rzH2(i) = rH2_molecularframe(3)
      end do
    end if

    ! print initial positions
    open(38,file="positions_init.dat")
    write(38,*)3*N ! number of molecules
    write(38,*)len
    do i=1,N ! 3 sites per molecule
      write(38,*)"O",rxO(i),ryO(i),rzO(i)
      write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH1(i),ryH1(i),rzH1(i)]
      write(38,*)"H",[rxO(i),ryO(i),rzO(i)]+[rxH2(i),ryH2(i),rzH2(i)]
    end do
    close(38)
  end subroutine positions_init


  ! print radial distribution function
  ! ==================================
  subroutine rdf_print
    implicit none
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
  end subroutine rdf_print

  subroutine fill_q
    implicit none
    integer :: i, j
    j=0
    do i=1,N
      j=j+1
      q(j)=qO
      j=j+1
      q(j)=qH
      j=j+1
      q(j)=qH
    end do
  end subroutine fill_q

  subroutine fill_positions
    implicit none
    integer :: i,j
    j=0
    do i=1,N
      j=j+1
      rx(j)=rxO(i)
      ry(j)=ryO(i)
      rz(j)=rzO(i)
      j=j+1
      rx(j)=rxO(i)+rxH1(i)
      ry(j)=ryO(i)+ryH1(i)
      rz(j)=rzO(i)+rzH1(i)
      j=j+1
      rx(j)=rxO(i)+rxH2(i)
      ry(j)=ryO(i)+ryH2(i)
      rz(j)=rzO(i)+rzH2(i)
    end do
  end subroutine fill_positions

end program lasvegas
