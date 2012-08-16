!********************************************************************************************
! This program solves Henry's problem using the Henry Method
! Psi = sum{sum{A_{m,n}*sin(m*pi*y)*cos(n*pi*x/xi), n=0..infinity}, m=1..infinity}
! Psi = sum{sum{B_{r,s}*sin(s*pi*x/xi)*cos(r*pi*y), s=1..infinity}, r=0..infinity}
! A_{g,h} = frac{sum{B_{r,h}*h*N(g,r),r=0..infinity} + 4/pi*W(g,h)}
!      {epsilon_2*a*pi**2*[g**2+h**2/xi**2]*xi}
! B_{g,h} = frac{pi/4*sum{sum{sum{sum{A_{m,n}*B_{r,s}*(m*s*L*R - n*r*F*G),
!      s=1..infinity},r=0..infinity},n=0..infinity},m=0..infinity}
!      + sum{A_{g,n}*g*N(h,n), n=0..infinity}
!      + epsilon_1*sum{B_{g,s}*s*N(h,s),s=1..infinity} + 4/pi*W(h,g)}
!      {epsilon_1*b*pi**2*[g**2+h**2/xi**2]*xi}
! Where A and B are both two-dimensional matrices.
! A                             : Fourier coefficients for Psi
! B                             : Fourier coefficients for C
! a=Q/(k_1*d)                   : dimensionless
! b=D/Q                         : dimensionless
! d                             : thickness of aquifer, m
! D                             : dispersion coefficient, m**2/sec
! E                             : dimensionless coefficient relating the concentration to
!                               : the density of the solution
! k                             : permeability of sand, m**2
! kbar=k*rho_0*g/mu             : transmission coefficient of aquifer, m/sec
! k_1=kbar*(rho_s-rho_o)/rho_o  : transmission coefficient times the density-difference
!                               : ratio
! l                             : length of aquifer, m
! Q                             : net freshwater discharge per unit length of beach,
!                               : m**2/sec
! rho                           : density of solution, kg/m**3
! rho_o                         : density of freshwater, kg/m**3
! rho_s                         : density of saltwater, kg/m**3
! mu                            : viscosity of water, N-sec/m**2
! xi=l/d                        : aspect ratio
!********************************************************************************************
PROGRAM Henrys_Problem
  
  USE Henry_funcs
  USE Henry_sums
  USE LU
  IMPLICIT NONE

  REAL, PARAMETER :: epsilon = 5E-4

  REAL, PARAMETER :: a = 1/pi2, b=1/pi2 !0.263, b = 0.1 
  REAL, PARAMETER :: d = 1.0, l = 2.0, dx = 0.05, dy= 0.05 
  REAL :: xi, xi2, api2, bpi2

  REAL, DIMENSION(:, :), ALLOCATABLE :: A_Matrix, B_Matrix, LHS 
  REAL, DIMENSION(:), ALLOCATABLE :: RHS, RHS_old 
  REAL, DIMENSION(:, :), ALLOCATABLE :: Psi, C 
  REAL, DIMENSION(:), ALLOCATABLE :: x, y 

  !LU decomp stuff
  INTEGER, DIMENSION(:), ALLOCATABLE :: indx
  REAL :: d0
  
  REAL :: error

  INTEGER :: h_a(3), h_b(0:2) !Size of h depends on both i and which coefficient A or B
  INTEGER :: i_a, i_b, j_a, j_b, i_y, j_x, linearsize, i, j, n 
  
  INTEGER :: AllocateStatus !Status variable for ALLOCATE
  INTEGER :: g,h

100     FORMAT("error=",ES8.2)
  OPEN(30, FILE='Psi0_1.txt')
  OPEN(31,FILE='C0_1.txt')

  h_a = 2!(/ 15, 10, 5, 3, 2 /)
  h_b = 3!(/ 20, 10, 5, 3, 2 /)

  i_a = SIZE(h_a,1)
  i_b = SIZE(h_b,1) - 1
  j_a = MAXVAL(h_a)
  j_b = MAXVAL(h_b)
  i_y = NINT(d/dy)
  j_x = NINT(l/dx)
  xi = 1!l/d

  xi2 = xi**2
  api2 = a*pi2
  bpi2 = b*pi2

  linearsize = i_a*(j_a+1)+(i_b+1)*j_b

  ALLOCATE(A_Matrix(i_a, 0:j_a), B_Matrix(0:i_b, j_b), Psi(0:i_y, 0:j_x), &
    C(0:i_y, 0:j_x), x(0:j_x), y(0:i_y), LHS(linearsize, linearsize), &
    RHS(linearsize), RHS_old(linearsize), indx(linearsize), &
    STAT = AllocateStatus) 
  IF(AllocateStatus /= 0) STOP "*** NOT ENOUGH MEMORY ***"

  A_Matrix = 0. 
  B_Matrix = 0. 
  RHS_old = 0.

  CALL initialize()
  
  error = 1.
  WRITE (*,*) 'b=', b


  DO WHILE(error > epsilon)
    CALL build_system(LHS,RHS) 
!    DO i=1,linearsize
!      IF(i .LE. i_a*(j_a+1)) THEN
!        CALL AgANDh(i,g,h)
!      ELSE
!        CALL BgANDh(i,g,h)
!      END IF
!      WRITE(*,*) g,h,(LHS(i,j), j=1,linearsize), "|", pid4*RHS(i)
!    END DO
!    STOP
 
    !A(0,h) is an independent system
!    CALL ludcmp(LHS(:j_a+1,:j_a+1),j_a+1,j_a+1,indx(:j_a+1),d0)
!    CALL lubksb(LHS(:j_a+1,:j_a+1),j_a+1,j_a+1,indx(:j_a+1),RHS(:j_a+1))

    !B(g,0) is an independent system
!    CALL ludcmp(LHS(i_a*(j_a+1)+1::j_b,i_a*(j_a+1)+1::j_b),j_b,j_b,indx(i_a*(j_a+1)+1::j_b),d0)
!    CALL lubksb(LHS(i_a*(j_a+1)+1::j_b,i_a*(j_a+1)+1::j_b),j_b,j_b,indx(i_a*(j_a+1)+1::j_b),RHS(i_a*(j_a+1)+1::j_b))
!    STOP

    !Replaces LHS with its LU decomposition
    CALL ludcmp(LHS, linearsize, linearsize, indx, d0) 
    !Solves A*x=B using ludcmp
    CALL lubksb(LHS, linearsize, linearsize, indx, RHS) 

    !Translates RHS into A and B matrices
    CALL AssignAandB(RHS) 

    error = L2Norm(RHS - RHS_old)/L2Norm(RHS_old)
    WRITE(*,100) error
    RHS_old = RHS 
  END DO

  !Assigns Psi and C based on A and B
  CALL PsiAndC() 
  CALL WriteToFile()

  CLOSE (30)
  CLOSE (31)

  DEALLOCATE(A_Matrix, B_Matrix, Psi, C, x, y, LHS, RHS)

  CONTAINS

  !This will create the Left Hand Side and Right Hand Side for the Henry Method
  SUBROUTINE build_system(LHS, RHS)
    REAL, INTENT(INOUT) :: LHS(:,:), RHS(:)
    INTEGER :: i, g, h,start, finish, stride
    
    LHS = 0.
    RHS = 0.

    DO i=1, linearsize
      !Entries associated with A(g,h)
      IF(i .LE. i_a*(j_a+1)) THEN
        CALL AgANDh(i,g,h)
!        WRITE(*,*) "A(",g, ",", h,")"
        LHS(i,i) = eps(h)*api2*(g**2 + h**2/xi2)*xi

        !This part is for the sum of B(r,h) terms
        IF(h/=0) THEN
          start = i_a*(j_a+1)+1
          stride = j_b
          LHS(i,start::stride) = -Br(g,h)
        END IF
        !Nonlinear and linear terms which are considered constants and so are on the RHS
        RHS(i) = fourdpi*W_FUNC(g,h)
      !Entries associated with B(g,h)
      ELSE
        CALL BgANDh(i,g,h)
!        WRITE(*,*) "B(",g, ",", h,")"
        LHS(i,i) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi

        !This is the sum of A(g,n) terms
        IF(g/=0) THEN
          start = g*(j_a+1)
          finish = g*(j_a+1)+j_a
          LHS(i,start:finish) = -An(g,h)
        END IF
        !This is the sum of B(g,s) terms
!        WRITE(*,*) eps(g)*Bs(h)
        start = i_a*(j_a+1) + g*j_b+1
        finish = i_a*(j_a+1) + (g+1)*j_b
!        WRITE(*,*) g, h, eps(g), Bs(h)
        LHS(i,start:finish) = LHS(i,start:finish) - eps(g)*Bs(h)

        !Nonlinear and linear terms which are considered constants and so are on the RHS
        RHS(i) = pid4*QuadVal(g,h) + fourdpi*W_FUNC(h,g) 
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE AgANDh(i, g, h) !Determines g and h given i inside the A section of LHS or RHS
    INTEGER*4, INTENT(IN) :: i
    INTEGER*4, INTENT(OUT) :: g, h
    INTEGER*4 :: q

    h = MOD(i, j_a + 1)
    g = (i - h)/(j_a+1) + 1
    h = h - 1

    IF (h < 0) THEN
      h = j_a
      g = g - 1
    END IF

  END SUBROUTINE AgANDh

  SUBROUTINE BgANDh(i, g, h) !Determines g and h given i inside the B section of LHS or RHS
    INTEGER*4, INTENT(IN) :: i
    INTEGER*4, INTENT(OUT) :: g, h
    INTEGER*4 :: q, r

    r = i - i_a*(j_a + 1)
    h = MOD(r, j_b)

    IF (h == 0) THEN
      h = j_b
    END IF
    g = (r - h)/j_b 
  END SUBROUTINE BgANDh

  !Create the Nonlinear term for a particular g,h 
  FUNCTION QuadVal(g,h)
    INTEGER, INTENT(IN) :: g, h
    INTEGER :: m, n, r, s
    REAL :: QuadVal

    QuadVal = 0.
    DO m=1,i_a
      DO n=0,j_a
        DO r=0,i_b
          DO s=1,j_b
            QuadVal = QuadVal + A_Matrix(m,n)*B_Matrix(n,s) &
              *(m*s*L_FUNC(m,r,g)*R_FUNC(h,n,s) &
              - n*r*F_FUNC(m,r,g)*G_FUNC(h,n,s)) 
          END DO
        END DO
      END DO
    END DO
  END FUNCTION QuadVal

  !Create the vector associated with the sum of the B(r,h) terms
  FUNCTION Br(g,h)
    REAL :: Br(0:i_b)
    INTEGER, INTENT(IN) :: g, h
    INTEGER :: r

    DO r=0,i_b
      Br(r) = h*N_FUNC(g,r)
    END DO
  END FUNCTION

  !Create the vector associated with the sum of the A(g,n) terms
  FUNCTION An(g,h)
    REAL :: An(0:j_a)
    INTEGER, INTENT(IN) :: g, h
    INTEGER :: n

    DO n=0,j_a
      An(n) = g*N_FUNC(h,n)
    END DO
  END FUNCTION

  !Create the vector associated with the sum of the B(g,s) terms
  FUNCTION Bs(h)
    REAL :: Bs(j_b)
    INTEGER, INTENT(IN) :: h
    INTEGER :: s

    DO s=1,j_b
      Bs(s) = s*N_FUNC(h,s)
    END DO
  END FUNCTION

  SUBROUTINE AssignAandB(Sol)
    REAL, DIMENSION(linearsize), INTENT(IN) :: Sol
    INTEGER*4 :: g, h, i

    DO i = 1, linearsize
      IF (i <= i_a*(j_a + 1)) THEN !A portion of Sol
        CALL AGANDH(i, g, h) !Finds g and h corresponding to counter i
        A_Matrix(g, h) = Sol(i)
      ELSE !B portion of Sol
        CALL BGANDH(i, g, h) !Finds g and h corresponding to counter i
        B_Matrix(g, h) = Sol(i)
      END IF
    END DO
  END SUBROUTINE AssignAandB

  FUNCTION L2Norm(F)
    REAL, INTENT(IN) :: F(:)
    REAL :: L2Norm
    INTEGER :: i,n

    n = SIZE(F)
    L2Norm = 0.

    DO i = 1, n
      L2Norm = L2Norm + (F(i))**2
    END DO

    L2Norm = L2Norm**0.5
  END FUNCTION L2NORM
  
  SUBROUTINE initialize()
    INTEGER :: i

    x_loop0: DO i = 0, j_x     !Initializing the values for x
      x(i) = dx*i
    END DO x_loop0

    y_loop0: DO i = 0, i_y     !Initializing the values for y
      y(i) = dy*i
    END DO y_loop0
  END SUBROUTINE initialize

  SUBROUTINE PsiAndC()
    INTEGER :: g, h

    DO h = 0, j_x     !Initializing the values for x
      x(h) = dx*h
    END DO

    DO g = 0, i_y     !Initializing the values for y
            y(g) = dy*g
    END DO

    DO g = 0, i_y !Find Psi(g, h) and C(g, h)
      DO h = 0, j_x
        IF (g == 0 .AND. h == 0) THEN
          Psi(g, h) = 0.
          C(g, h) = 0.
        ELSE
          Psi(g, h) = Sum_A(A_Matrix, x, y, xi, g, h, i_a, j_a, i_y, j_x) + y(g)/d
          C(g, h) = Sum_B(B_Matrix, x, y, xi, g, h, i_b, j_b, i_y, j_x) + x(h)/(d*xi)
        END IF
      END DO
    END DO
  END SUBROUTINE PSIANDC

  SUBROUTINE WriteToFile()
    INTEGER :: g, h

101     FORMAT(1x, 3(F7.3,3x))

    DO h = 0, j_x !Write Psi and C out to file
      DO g = 0, i_y
        WRITE (30, 101) x(h), y(g), Psi(g, h)
        WRITE (31, 101) x(h), y(g), C(g, h)
      END DO
      WRITE (30,*)
      WRITE (31,*)
    END DO
  END SUBROUTINE

END PROGRAM Henrys_Problem
