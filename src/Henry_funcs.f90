MODULE Henry_funcs
!This provides the various functions associated with the Henry problem
  
  IMPLICIT none
  REAL*8, PARAMETER :: pi = 4*ATAN(1.)
  REAL*8, PARAMETER :: pi2 = pi**2
  REAL*8, PARAMETER :: pid4 = pi/4., fourdpi = 4./pi
  
  CONTAINS

  FUNCTION F_Func(m, r, g)
    INTEGER, INTENT(IN) :: m, r, g
    INTEGER :: F_Func

    F_Func = delta(m - r, g) + delta(r - m, g) - delta(m + r, g)
  END FUNCTION F_Func
        
  FUNCTION L_Func(m, r, g)
    INTEGER, INTENT(IN) :: m, r, g
    INTEGER :: L_Func
              
    L_Func = delta(m - r, g) + delta(r - m, g) + delta(m + r, g)
  END FUNCTION L_Func

  FUNCTION G_Func(h, n, s)
    INTEGER, INTENT(IN) :: h, n, s
    REAL :: G_Func
    
    G_Func=T_Func(h+n-s)+T_Func(h-n+s)-T_Func(h+n+s)-T_Func(h-n-s)
  END FUNCTION G_Func

  FUNCTION R_Func(h, n, s)
    INTEGER, INTENT(IN) :: h, n, s
    REAL :: R_Func

    R_Func=T_Func(h+n-s)+T_Func(h-n+s)+T_Func(h+n+s)+T_Func(h-n-s)
  END FUNCTION R_Func
  
  FUNCTION T_Func(i)
    INTEGER, INTENT(IN) :: i
    REAL :: T_Func

    IF(i == 0) THEN
      T_Func = 0.
    ELSE
      T_Func = REAL((-1)**i-1)/REAL(i)
    END IF
  END FUNCTION T_Func

  FUNCTION N_Func(h, n)
    INTEGER, INTENT(IN) :: h, n
    REAL :: N_Func

    N_Func = T_Func(h + n) + T_Func(h - n)
  END FUNCTION N_Func

  FUNCTION W_Func(h, g) !((-1)**h - 1)/h if g=0 and 0 otherwise
    INTEGER, INTENT(IN) :: h, g
    REAL :: W_Func

    W_Func = T_Func(h)*delta(g, 0)
  END FUNCTION W_Func

  FUNCTION eps(e)
    INTEGER, INTENT(IN) :: e
    INTEGER :: eps

    eps = delta(e, 0) + 1  ! 2 if e=0 and 1 otherwise
  END FUNCTION eps

  FUNCTION delta(i, j) !Kronecker delta function 1 if i=j and 0 otherwise
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: delta
                
    IF (i == j) THEN
      delta = 1
    ELSE
      delta = 0
    END IF
  END FUNCTION delta

END MODULE Henry_funcs
