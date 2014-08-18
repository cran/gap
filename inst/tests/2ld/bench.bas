CLS
REM program to generate tests for 2LD
REM JH Zhao 14/2/01
REM Reference Zapata et al. (2001) table 1
m = 4
n = 4
p = .25

REM alternative coordinates
FOR i = 1 TO m
FOR j = 1 TO n
   IF INT((i + j) / 2) * 2 = i + j THEN
      PRINT USING "  #/# "; i; j;
   ELSE
      PRINT USING " [#/#]"; i; j;
   END IF
NEXT: PRINT
NEXT

REM set D'
DIM d(m, n), dp(m, n), dmax(m, n)
PRINT
FOR i = 1 TO m
FOR j = 1 TO n
   IF INT((i + j) / 2) * 2 = i + j THEN v = .2 ELSE v = -.6
   dp(i, j) = v
   IF dp(i, j) < 0 THEN
      dmax(i, j) = p ^ 2
   ELSE
      dmax(i, j) = p * (1 - p)
   END IF
   d(i, j) = dp(i, j) * dmax(i, j)
   PRINT USING " #.#"; v;
NEXT: PRINT
NEXT

REM print d
PRINT
FOR i = 1 TO m
FOR j = 1 TO n
   PRINT USING " #.#####"; d(i, j);
NEXT: PRINT
NEXT

REM print h
PRINT
FOR i = 1 TO m
FOR j = 1 TO n
   PRINT USING " #.#####"; d(i, j) + p ^ 2;
NEXT: PRINT
NEXT

END

