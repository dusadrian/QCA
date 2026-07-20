# tests for superSubset() have the same output

    Code
      Kro.ss
    Output
      
                   inclN   RoN   covN  
      -------------------------------- 
      1  ES + LP   1.000  0.636  0.733 
      2  WS + LP   0.909  0.917  0.909 
      3  ~WM + LP  0.909  0.583  0.667 
      -------------------------------- 
      

---

    Code
      head(Kro.coms)
    Output
         ES+LP WS+LP ~WM+LP
      SE     1     1      1
      FI     1     1      1
      NO     1     1      1
      DK     1     1      1
      NL     1     1      1
      ES     1     0      0

---

    Code
      HK.ss
    Output
      
                                    inclN   RoN   covN  
      ------------------------------------------------- 
      1  C[1] + F[2] + T[2] + V[0]  1.000  0.333  0.913 
      ------------------------------------------------- 
      

---

    Code
      HK.ss$coms[1:10, , drop = FALSE]
    Output
         C[1]+F[2]+T[2]+V[0]
      AO                   1
      BJ                   1
      BW                   1
      BF                   1
      BI                   1
      CF                   1
      CM                   1
      CV                   1
      TD                   1
      KM                   1

---

    Code
      superSubset(d.jobsecurity, outcome = JSR, relation = "suf", incl.cut = 0.9,
        cov.cut = 0.4)
    Output
      
                 inclS   PRI   covS  
      ------------------------------ 
      1  L*R     0.972  0.943  0.582 
      2  L*P     0.979  0.955  0.506 
      3  R*V     0.950  0.904  0.418 
      4  P*~V    0.950  0.883  0.417 
      5  S*~L*R  0.990  0.982  0.401 
      6  S*R*~V  0.990  0.983  0.402 
      7  C*R*P   0.932  0.836  0.410 
      ------------------------------ 
      

---

    Code
      superSubset(d.jobsecurity, outcome = ~JSR, relation = "suf", incl.cut = 0.9,
        cov.cut = 0.4)
    Output
      
                   inclS   PRI   covS  
      -------------------------------- 
      1  ~S*~C     0.936  0.906  0.519 
      2  ~C*~R     0.943  0.923  0.586 
      3  ~S*~L*~P  1.000  1.000  0.450 
      4  ~L*~R*~P  0.987  0.981  0.491 
      -------------------------------- 
      

---

    Code
      ssLC
    Output
      
                      inclN   RoN   covN  
      ----------------------------------- 
      1  DEV          1.000  0.800  0.800 
      2  LIT          1.000  0.500  0.615 
      3  STB          1.000  0.700  0.727 
      4  DEV*LIT      1.000  0.800  0.800 
      5  DEV*STB      1.000  1.000  1.000 
      6  LIT*STB      1.000  0.900  0.889 
      7  DEV*LIT*STB  1.000  1.000  1.000 
      8  ~URB + IND   1.000  0.000  0.444 
      ----------------------------------- 
      

---

    Code
      superSubset(LM, SURV)
    Output
      
                                   inclN   RoN   covN  
      ------------------------------------------------ 
      1  LIT[1]                    1.000  0.500  0.615 
      2  STB[1]                    1.000  0.700  0.727 
      3  LIT[1]*STB[1]             1.000  0.900  0.889 
      4  DEV[1] + IND[1]           1.000  0.800  0.800 
      5  URB[0] + IND[1]           1.000  0.000  0.444 
      6  DEV[2] + URB[1] + IND[0]  1.000  0.100  0.471 
      ------------------------------------------------ 
      

---

    Code
      ssCVF
    Output
      
                                                   inclN   RoN   covN  
      ---------------------------------------------------------------- 
       1  GEOCON                                   0.904  0.492  0.624 
       2  DEMOC + ETHFRACT + ~GEOCON               0.930  0.470  0.626 
       3  DEMOC + ~ETHFRACT + POLDIS               0.918  0.506  0.637 
       4  DEMOC + ETHFRACT + POLDIS                0.906  0.502  0.630 
       5  DEMOC + ~ETHFRACT + ~NATPRIDE            0.905  0.527  0.641 
       6  DEMOC + ETHFRACT + ~NATPRIDE             0.935  0.530  0.656 
       7  DEMOC + ~GEOCON + POLDIS                 0.920  0.539  0.654 
       8  DEMOC + ~GEOCON + ~NATPRIDE              0.908  0.584  0.671 
       9  DEMOC + POLDIS + ~NATPRIDE               0.916  0.596  0.682 
      10  ~ETHFRACT + POLDIS + ~NATPRIDE           0.911  0.554  0.657 
      11  ~DEMOC + ETHFRACT + POLDIS + ~NATPRIDE   0.913  0.532  0.647 
      12  ETHFRACT + ~GEOCON + POLDIS + ~NATPRIDE  0.911  0.613  0.688 
      ---------------------------------------------------------------- 
      

---

    Code
      ssCVF$coms$GEOCON
    Output
       [1] 0.95 0.35 0.35 0.78 0.35 0.78 0.78 0.78 0.78 0.05 0.78 0.35 0.95 0.95 0.35
      [16] 0.95 0.78 0.35 0.95 0.35 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95

---

    Code
      superSubset(CVF, outcome = ~PROTEST, incl.cut = 0.9, cov.cut = 0.6)
    Output
      
                              inclN   RoN   covN  
      ------------------------------------------- 
      1  NATPRIDE             0.932  0.622  0.693 
      2  ~DEMOC + ~ETHFRACT   0.951  0.548  0.663 
      3  ~ETHFRACT + ~POLDIS  0.927  0.443  0.603 
      ------------------------------------------- 
      

---

    Code
      ssLC4
    Output
      
                                            inclN   RoN   covN  
      --------------------------------------------------------- 
      1  Developed                          1.000  0.800  0.800 
      2  High literacy                      1.000  0.500  0.615 
      3  Stable                             1.000  0.700  0.727 
      4  Developed*High literacy            1.000  0.800  0.800 
      5  Developed*Stable                   1.000  1.000  1.000 
      6  High literacy*Stable               1.000  0.900  0.889 
      7  Developed*High literacy*Stable     1.000  1.000  1.000 
      8  Low urbanization + Industrialized  1.000  0.000  0.444 
      --------------------------------------------------------- 
      

