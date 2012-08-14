MODULE Henry_funcs
  
  IMPLICIT none
  REAL*8, PARAMETER :: pi = 4*ATAN(1.)
  REAL*8, PARAMETER :: pi2 = pi**2
  
  CONTAINS

  FUNCTION F_FUNC(m, r, g)
              
    INTEGER, INTENT(IN) :: m, r, g
    INTEGER :: F_FUNC

    F_FUNC = delta(m - r, g) + delta(r - m, g) - delta(m + r, g)
              
  END FUNCTION F_FUNC
        
  FUNCTION L_FUNC(m, r, g)
              
    INTEGER, INTENT(IN) :: m, r, g
    INTEGER :: L_FUNC
              
    L_FUNC = delta(m - r, g) + delta(r - m, g) + delta(m + r, g)

  END FUNCTION L_FUNC

  FUNCTION G_FUNC(h, n, s)

    INTEGER, INTENT(IN) :: h, n, s
    REAL :: G_FUNC
    
    G_FUNC=T_FUNC(h+n-s)+T_FUNC(h-n+s)-T_FUNC(h+n+s)-T_FUNC(h-n-s)

  END FUNCTION G_FUNC

  FUNCTION R_FUNC(h, n, s)

    INTEGER, INTENT(IN) :: h, n, s
    REAL :: R_FUNC

    R_FUNC=T_FUNC(h+n-s)+T_FUNC(h-n+s)+T_FUNC(h+n+s)+T_FUNC(h-n-s)

  END FUNCTION R_FUNC
  
  FUNCTION T_FUNC(i)

    INTEGER, INTENT(IN) :: i
    REAL :: T_FUNC

    IF(i == 0) THEN
      T_FUNC = 0.
    ELSE
      T_FUNC = REAL((-1)**i-1)/REAL(i)
    END IF

  END FUNCTION T_FUNC

  FUNCTION N_FUNC(h, n)

    INTEGER, INTENT(IN) :: h, n
    REAL :: N_FUNC

    N_FUNC = T_FUNC(h + n) + T_FUNC(h - n)

  END FUNCTION N_FUNC

  FUNCTION W_FUNC(h, g)                !((-1)**h - 1)/h if g=0 and 0 otherwise

    INTEGER, INTENT(IN) :: h, g
    REAL :: W_FUNC

    W_FUNC = T_FUNC(h)*delta(g, 0)

  END FUNCTION W_FUNC

  FUNCTION eps(e)

    INTEGER, INTENT(IN) :: e
    INTEGER :: eps

    eps = delta(e, 0) + 1                ! 2 if e=0 and 1 otherwise

  END FUNCTION eps

  FUNCTION delta(i, j)                !Kronecker delta function 1 if i=j and 0 otherwise

    INTEGER, INTENT(IN) :: i, j
    INTEGER :: delta
                
    IF (i == j) THEN
      delta = 1
    ELSE
      delta = 0
    END IF

  END FUNCTION delta

END MODULE Henry_funcs
