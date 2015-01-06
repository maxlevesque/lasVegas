module ewald
  implicit none
  double precision, parameter :: qe=1.602176565d-19, eps0=8.854187817D-12, pi=acos(-1.D0), sqrtpi=sqrt(pi), pi2=pi**2,&
      twopi=2.d0*pi, qfact=1389.3545783908448 !kJ/mol qfact=(qe**2)/(4.D0*pi*eps0*1D-10)!    ! in Joules
  complex, parameter :: twopii=cmplx(0.d0,twopi)
  double precision, parameter :: rcut = 10.d0
  integer, parameter :: kmax=5, nkmax=27
  double precision :: alpha, alpha2, pi2Ialpha2

  public

contains
  subroutine ewald_init(len)
    implicit none
    double precision, intent(in) :: len
    alpha=5.6d0/len
    alpha2=alpha**2
    pi2Ialpha2=pi2/alpha2
  end subroutine

  subroutine ewald_energy (N,len,q,x,y,z,Uewald)
    implicit none
    integer, intent(in) :: N ! number of sites
    double precision, intent(in) :: len
    double precision, intent(out) :: Uewald
    double precision, intent(in) :: q(3*N), x(3*N), y(3*N), z(3*N)
    double precision :: Ereal, Eself, r, k2lmn, Efourier, Eintra, symfac
    complex :: Efourier_part2
    integer :: i,j,l,m,o
    ! open(46,file="/home/levesque/Recherche/src/libewald/spce_sample_config_periodic1.txt")
    ! read(46,*)
    ! read(46,*)
    ! do i=1,3*N
    !   read(46,*)j,x(i),y(i),z(i)
    ! end do
    ! close(46)
    Ereal=0.d0
    Eintra=0.d0
    do i=1,3*N-1
      do j=i+1,3*N
        r=rij(i,j)
        if(moleculeofsite(i)==moleculeofsite(j)) then
          Eintra = Eintra -q(j)*q(i)*erf(alpha*r)/r
        else
          if(r>rcut) cycle
          Ereal = Ereal +q(j)*q(i)*erfc(alpha*r)/r
        end if
      end do
    end do
    Ereal= Ereal*qfact
    Eintra= Eintra*qfact
    Eself = -alpha/sqrtpi *sum(q**2) *qfact

    Efourier = 0.d0
    do l=0,kmax
      select case(l); case(0); symfac = 1.d0; case default; symfac = 2.d0; end select
      do m=-kmax,kmax
        do o=-kmax,kmax
          if(l==0 .and. m==0 .and. o==0) cycle ! E(k=0)=0
          if(l**2+m**2+o**2 >= nkmax) cycle
          k2lmn = k2(l,m,o)
          Efourier = Efourier + symfac*Exp(-k2lmn*pi2Ialpha2)/k2lmn *abs(sum(  q*Exp(twopii/len*(l*x+m*y+o*z))  ))**2
        end do
      end do
    end do
    Efourier = Efourier *1.d0/(twopi*len**3) *qfact

    Uewald = Ereal + Eintra + Eself + Efourier

  contains

    pure function rij(i,j)
      implicit none
      double precision :: rij
      integer, intent(in) :: i,j
      double precision :: dx, dy, dz
      dx = abs( x(i)-x(j) )
      if(dx>len/2.) dx=len-dx
      dy = abs( y(i)-y(j) )
      if(dy>len/2.) dy=len-dy
      dz = abs( z(i)-z(j) )
      if(dz>len/2.) dz=len-dz
      rij=sqrt(dx**2+dy**2+dz**2)
    end function rij

    pure function moleculeofsite(i)
      integer :: moleculeofsite
      integer, parameter :: sitepermolecule=3 ! spc like water
      integer, intent(in) :: i
      moleculeofsite = (i-1)/sitepermolecule +1 ! division between integers: that's what we want here.
    end function moleculeofsite

    pure function kvec(l,m,n)
      double precision :: kvec(3)
      integer, intent(in) :: l,m,n
      kvec = real([l,m,n])/len
    end function kvec

    pure function k2(l,m,n)
      double precision :: k2
      integer, intent(in) :: l,m,n
      k2 = (l**2+m**2+n**2)/len**2
    end function k2

  end subroutine


end module ewald
