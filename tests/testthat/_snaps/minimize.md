# tests for minimize() have the same output

    Code
      minimize(d.represent, outcome = "WNP")
    Output
      
      M1: ~ES*~WS*WM*LP + ES*QU*~WS*WM + ES*QU*~WS*LP + ES*WS*~WM*~LP + ES*WS*WM*LP
          <-> WNP
      

---

    Code
      minimize(d.represent, outcome = "~WNP")
    Output
      
      M1: ~ES*~WS*WM*~LP + QU*~WS*~WM*~LP + ES*~QU*~WS*~WM*LP <-> ~WNP
      

---

    Code
      Kro.sp
    Output
      
      M1: WS + ES*WM + QU*LP + (~ES*LP) <-> WNP 
      M2: WS + ES*WM + QU*LP + (WM*LP) <-> WNP 
      
                                      ------------------- 
                 inclS   PRI   covS   covU   (M1)   (M2)   cases 
      ----------------------------------------------------------------------------- 
      1      WS  1.000  1.000  0.455  0.182  0.182  0.182  FI; DK; SE; NO,IS 
      2   ES*WM  1.000  1.000  0.545  0.091  0.091  0.091  DK; ES; NL,BE; NO,IS 
      3   QU*LP  1.000  1.000  0.545  0.091  0.091  0.091  DE; AT; NL,BE; NO,IS 
      ----------------------------------------------------------------------------- 
      4  ~ES*LP  1.000  1.000  0.182  0.000  0.091         NZ; DE 
      5   WM*LP  1.000  1.000  0.636  0.000         0.091  NZ; DE; DK; NL,BE; NO,IS 
      ----------------------------------------------------------------------------- 
             M1  1.000  1.000  1.000 
             M2  1.000  1.000  1.000 
      

---

    Code
      Kro.sp$PIchart
    Output
      
                    4  12 21 24 26 27 28 29 32
      WS            -  -  x  x  -  -  -  x  x 
      ~ES*LP        x  x  -  -  -  -  -  -  - 
      ES*WM         -  -  -  x  -  x  x  -  x 
      QU*LP         -  x  -  -  x  -  x  -  x 
      WM*LP         x  x  -  x  -  -  x  -  x 
      ES*~QU*~LP    -  -  x  -  -  -  -  -  - 
      ~QU*~WM*~LP   -  -  x  -  -  -  -  -  - 
      

---

    Code
      Kro.sp$SA
    Output
      $M1
         ES QU WS WM LP
      2   0  0  0  0  1
      5   0  0  1  0  0
      6   0  0  1  0  1
      7   0  0  1  1  0
      8   0  0  1  1  1
      10  0  1  0  0  1
      13  0  1  1  0  0
      14  0  1  1  0  1
      15  0  1  1  1  0
      16  0  1  1  1  1
      19  1  0  0  1  0
      20  1  0  0  1  1
      22  1  0  1  0  1
      23  1  0  1  1  0
      30  1  1  1  0  1
      31  1  1  1  1  0
      
      $M2
         ES QU WS WM LP
      5   0  0  1  0  0
      6   0  0  1  0  1
      7   0  0  1  1  0
      8   0  0  1  1  1
      10  0  1  0  0  1
      13  0  1  1  0  0
      14  0  1  1  0  1
      15  0  1  1  1  0
      16  0  1  1  1  1
      19  1  0  0  1  0
      20  1  0  0  1  1
      22  1  0  1  0  1
      23  1  0  1  1  0
      30  1  1  1  0  1
      31  1  1  1  1  0
      

---

    Code
      print(minimize(cbind(Kro.sp$SA[[1]], FO = 1), outcome = "FO"))
    Output
      
      M1: ~ES*WS + ~ES*~WM*LP + WS*~WM*LP + WS*WM*~LP + ES*~QU*~WS*WM <-> FO
      

---

    Code
      print(minimize(cbind(Kro.sp$SA[[2]], FO = 1), outcome = "FO"))
    Output
      
      M1: ~ES*WS + WS*~WM*LP + WS*WM*~LP + ~ES*QU*~WM*LP + ES*~QU*~WS*WM <-> FO
      

---

    Code
      Kro.sc
    Output
      
      M1: ~ES*~WS*WM*LP + ES*QU*~WS*WM + ES*QU*~WS*LP + ES*WS*~WM*~LP + ES*WS*WM*LP
          <-> WNP
      

---

    Code
      minimize(d.jobsecurity, outcome = "JSR", incl.cut = 0.9, include = "?",
        details = TRUE)
    Output
      
      M1: S*R + (L*P) -> JSR 
      M2: S*R + (P*~V) -> JSR 
      
                                    ------------------- 
               inclS   PRI   covS   covU   (M1)   (M2)  
      ------------------------------------------------- 
      1   S*R  0.871  0.821  0.610  0.231  0.256  0.335 
      ------------------------------------------------- 
      2   L*P  0.979  0.955  0.506  0.014  0.152        
      3  P*~V  0.950  0.883  0.417  0.004         0.142 
      ------------------------------------------------- 
           M1  0.883  0.821  0.762 
           M2  0.874  0.813  0.752 
      

---

    Code
      Emm.si
    Output
      
      From C1P1, C1P2: 
      
      M1:    S*R*~V + S*C*R*P + S*L*R*P + C*L*P*~V -> JSR 
      
                   inclS   PRI   covS   covU  
      --------------------------------------- 
      1    S*R*~V  0.990  0.983  0.402  0.152 
      2   S*C*R*P  0.965  0.921  0.277  0.041 
      3   S*L*R*P  1.000  1.000  0.354  0.027 
      4  C*L*P*~V  0.964  0.872  0.297  0.138 
      --------------------------------------- 
               M1  0.965  0.941  0.685 
      

---

    Code
      pof(Emm.si$i.sol$C1P1$pims, 1 - d.jobsecurity$JSR, relation = "suf")
    Output
      
                   inclS   PRI   covS   covU  
      --------------------------------------- 
      1    S*R*~V  0.438  0.017  0.198  0.031 
      2   S*C*R*P  0.557  0.000  0.178  0.000 
      3   S*L*R*P  0.492  0.000  0.193  0.000 
      4  C*L*P*~V  0.756  0.128  0.259  0.108 
      --------------------------------------- 
      

---

    Code
      Emm.si$i.sol$C1P1$PIchart
    Output
      
                 27 37 47 48 56 64
      S*R*~V     -  x  x  -  -  - 
      S*C*R*P    -  -  -  -  x  x 
      S*L*R*P    -  -  x  x  -  x 
      C*L*P*~V   x  -  -  -  -  - 
      

---

    Code
      EC1
    Output
         S C L R P V
      31 0 1 1 1 1 0
      39 1 0 0 1 1 0
      45 1 0 1 1 0 0
      53 1 1 0 1 0 0
      55 1 1 0 1 1 0
      59 1 1 1 0 1 0
      61 1 1 1 1 0 0
      63 1 1 1 1 1 0

---

    Code
      EC2
    Output
         S C L R P V
      31 0 1 1 1 1 0
      39 1 0 0 1 1 0
      45 1 0 1 1 0 0
      53 1 1 0 1 0 0
      55 1 1 0 1 1 0
      59 1 1 1 0 1 0
      61 1 1 1 1 0 0
      63 1 1 1 1 1 0

---

    Code
      minimize(cbind(Emm.si$i.sol$C1P1$EC, FO = 1), outcome = "FO")
    Output
      
      M1: S*C*R*~V + S*C*L*P*~V + S*~L*R*P*~V + S*L*R*~P*~V + C*L*R*P*~V <-> FO
      

---

    Code
      Emm.si$i.sol$C1P1$pims
    Output
         S*R*~V S*C*R*P S*L*R*P C*L*P*~V
      AU   0.00    0.00    0.00     0.00
      AT   0.33    0.67    0.57     0.33
      BE   0.33    0.67    0.43     0.33
      CA   0.00    0.00    0.00     0.00
      DK   0.00    0.00    0.00     0.40
      FI   0.40    0.40    0.40     0.40
      FR   0.67    0.20    0.20     0.20
      DE   0.00    0.60    0.43     0.00
      IE   0.33    0.00    0.14     0.00
      IT   0.67    0.33    0.57     0.33
      NL   0.00    0.00    0.00     0.29
      NZ   0.00    0.00    0.00     0.00
      NO   0.00    0.00    0.00     0.60
      PT   1.00    0.00    0.20     0.00
      ES   0.33    0.00    0.60     0.00
      SE   0.00    0.00    0.00     0.20
      CH   0.00    0.00    0.00     0.00
      GB   0.00    0.00    0.00     0.00
      US   0.00    0.00    0.00     0.00

---

    Code
      minimize(d.jobsecurity, outcome = "~JSR", incl.cut1 = 0.9, incl.cut0 = 0.4,
        explain = "C", include = "?", details = TRUE)
    Output
      
      M1: P + (~S*C) -> ~JSR 
      M2: P + (~S*L*~V) -> ~JSR 
      M3: P + (~S*~R*~V) -> ~JSR 
      
                                        -------------------------- 
                   inclS   PRI   covS   covU   (M1)   (M2)   (M3)  
      ------------------------------------------------------------ 
      1         P  0.438  0.119  0.409  0.142  0.158  0.217  0.264 
      ------------------------------------------------------------ 
      2      ~S*C  0.542  0.262  0.281  0.000  0.030               
      3   ~S*L*~V  0.596  0.265  0.259  0.004         0.067        
      4  ~S*~R*~V  0.656  0.342  0.203  0.004                0.059 
      ------------------------------------------------------------ 
               M1  0.400  0.133  0.439 
               M2  0.427  0.161  0.476 
               M3  0.451  0.167  0.468 
      

---

    Code
      HK.sp
    Output
      
      M1: C[1] + T[2] + T[1]*V[0] + (C[0]) <-> PB[1] 
      M2: C[1] + T[2] + T[1]*V[0] + (F[2]) <-> PB[1] 
      
                                         ------------------- 
                    inclS   PRI   covS   covU   (M1)   (M2)  
      ------------------------------------------------------ 
      1       C[1]  1.000  1.000  0.405  0.048  0.071  0.048 
      2       T[2]  0.941  0.941  0.762  0.095  0.214  0.167 
      3  T[1]*V[0]  1.000  1.000  0.143  0.048  0.048  0.048 
      ------------------------------------------------------ 
      4       C[0]  1.000  1.000  0.333  0.000  0.024        
      5       F[2]  1.000  1.000  0.619  0.000         0.024 
      ------------------------------------------------------ 
                M1  0.955  0.955  1.000 
                M2  0.955  0.955  1.000 
      

---

    Code
      HK.sp$pims
    Output
         C[0] C[1] F[2] T[2] T[1]*V[0]
      AO    1    0    1    0         0
      BJ    0    1    1    0         1
      BW    0    0    0    0         0
      BF    0    1    1    1         0
      BI    1    0    1    1         0
      CF    0    1    0    1         0
      CM    0    1    0    1         0
      CV    1    0    1    1         0
      TD    0    1    1    1         0
      KM    0    1    1    1         0
      CD    1    0    1    1         0
      CG    0    1    1    0         0
      CI    0    1    0    1         0
      DJ    0    1    1    1         0
      GQ    1    0    1    1         0
      ER    1    0    1    1         0
      ET    1    0    1    1         0
      GA    0    1    1    1         0
      GM    0    0    0    1         0
      GH    0    0    1    1         0
      GW    1    0    1    1         0
      GN    0    1    1    1         0
      KE    0    0    0    1         0
      LS    0    0    1    1         0
      LR    1    0    1    1         0
      MG    0    1    0    0         1
      MW    0    0    0    0         1
      ML    0    1    0    0         0
      MR    0    1    1    1         0
      MU    0    0    0    0         0
      MZ    1    0    1    0         1
      NA    0    0    0    0         1
      NG    0    0    0    1         0
      NE    0    1    1    0         1
      RW    1    0    0    1         0
      ST    1    0    0    1         0
      SN    0    1    0    0         0
      SC    1    0    0    1         0
      SL    0    0    0    1         0
      SO    1    0    1    1         0
      ZA    0    0    0    0         0
      SD    0    0    1    1         0
      SZ    0    0    1    1         0
      TZ    0    0    0    1         0
      TG    0    1    0    1         0
      UG    0    0    1    1         0
      ZM    0    0    0    1         0
      ZW    0    0    0    0         0

---

    Code
      rownames(d.partybans[d.partybans$T == 2 & d.partybans$PB != 1, ])
    Output
      [1] "KE" "ZM"

---

    Code
      minimize(d.partybans, outcome = "PB[1]", conditions = "C, F, T, V", incl.cut0 = 0.4,
        explain = "C")
    Output
      
      M1: C[2]*F[1]*T[2]*V[1]
      

---

    Code
      HK.si
    Output
      
      From C1P1: 
      
      M1:    C[1] + F[2] + T[1]*V[0] + T[2]*V[0] + C[0]*F[1]*T[2] -> PB[1] 
      
                         inclS   PRI   covS   covU  
      --------------------------------------------- 
      1            C[1]  1.000  1.000  0.405  0.048 
      2            F[2]  1.000  1.000  0.619  0.143 
      3       T[1]*V[0]  1.000  1.000  0.143  0.048 
      4       T[2]*V[0]  1.000  1.000  0.571  0.048 
      5  C[0]*F[1]*T[2]  1.000  1.000  0.071  0.024 
      --------------------------------------------- 
                     M1  1.000  1.000  0.952 
      

---

    Code
      SA.si
    Output
      
      From C1P1: 
      
      M1:    FED[0]*FIN[1]*HIS[1] + FED[2]*URB[0]*GER[1] + FIN[0]*GER[1]*HIS[0] +
             FIN[1]*GER[1]*HIS[1] -> ACC[1] 
      
                               inclS   PRI   covS   covU  
      --------------------------------------------------- 
      1  FED[0]*FIN[1]*HIS[1]  1.000  1.000  0.250  0.083 
      2  FED[2]*URB[0]*GER[1]  1.000  1.000  0.167  0.083 
      3  FIN[0]*GER[1]*HIS[0]  1.000  1.000  0.167  0.083 
      4  FIN[1]*GER[1]*HIS[1]  1.000  1.000  0.417  0.250 
      --------------------------------------------------- 
                           M1  1.000  1.000  0.750 
      

---

    Code
      minimize(d.graduate, outcome = "REC", details = TRUE, show.cases = TRUE)
    Output
      
      M1: P*E*S + P*E*A*EBA + E*A*S*~EBA <-> REC
      
                     inclS   PRI   covS   covU   cases 
      ---------------------------------------------------- 
      1       P*E*S  1.000  1.000  0.571  0.143  f; c; a,b 
      2   P*E*A*EBA  1.000  1.000  0.571  0.286  d,e; a,b 
      3  E*A*S*~EBA  1.000  1.000  0.286  0.143  m; c 
      ---------------------------------------------------- 
                 M1  1.000  1.000  1.000 
      

---

    Code
      Bau.cna
    Output
      
      M1: ~D*L -> U
      
      M1: ~U*L -> D
      
      M1: U + D <-> L
      
      M1: ~L*E -> G 
      M2: ~U*~D*E -> G 
      
      M1: L + G <-> E
      

---

    Code
      print(Bau.cna$E, details = TRUE, show.cases = TRUE)
    Output
      
      M1: L + G <-> E
      
             inclS   PRI   covS   covU   cases 
      --------------------------------------------------- 
      1   L  1.000  1.000  0.857  0.429  f; e; d; c; b; a 
      2   G  1.000  1.000  0.571  0.143  g; e; c; a 
      --------------------------------------------------- 
         M1  1.000  1.000  1.000 
      

---

    Code
      head(d.mmv)
    Output
        A B C D E
      a 2 2 0 2 3
      b 0 2 1 1 2
      c 0 2 0 2 3
      d 1 2 0 2 3
      e 1 1 0 3 0
      f 1 1 2 1 1

---

    Code
      mmv.s
    Output
      
      M1: B[2]*C[0]*E[3] <-> D[2]
      
      M1: B[2]*C[0]*D[2] + A[2]*B[0]*C[1]*D[3] <-> E[3]
      

---

    Code
      print(mmv.s$"E[3]", details = TRUE, show.cases = TRUE)
    Output
      
      M1: B[2]*C[0]*D[2] + A[2]*B[0]*C[1]*D[3] <-> E[3]
      
                              inclS   PRI   covS   covU   cases 
      ----------------------------------------------------------- 
      1       B[2]*C[0]*D[2]  1.000  1.000  0.750  0.750  c; d; a 
      2  A[2]*B[0]*C[1]*D[3]  1.000  1.000  0.250  0.250  g 
      ----------------------------------------------------------- 
                          M1  1.000  1.000  1.000 
      

---

    Code
      mmv.e3
    Output
      
      M01: B[1] + D[0] + D[1] <-> ~E[3] 
      M02: B[1] + D[0] + A[0]*C[1] <-> ~E[3] 
      M03: B[1] + D[0] + B[2]*C[1] <-> ~E[3] 
      M04: B[1] + D[1] + B[0]*C[0] <-> ~E[3] 
      M05: B[1] + A[0]*C[1] + B[0]*C[0] <-> ~E[3] 
      M06: B[1] + B[0]*C[0] + B[2]*C[1] <-> ~E[3] 
      M07: D[0] + D[1] + A[1]*D[3] <-> ~E[3] 
      M08: D[0] + D[1] + C[0]*D[3] <-> ~E[3] 
      M09: D[1] + A[1]*D[3] + B[0]*C[0] <-> ~E[3] 
      M10: D[1] + B[0]*C[0] + C[0]*D[3] <-> ~E[3] 
      

---

    Code
      ttLC
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
      32    1   1   1   1   1     1     4  1.000 1.000
      22    1   0   1   0   1     1     2  1.000 1.000
      24    1   0   1   1   1     1     2  1.000 1.000
       1    0   0   0   0   0     0     3  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       5    0   0   1   0   0     0     2  0.000 0.000
       6    0   0   1   0   1     0     1  0.000 0.000
      23    1   0   1   1   0     0     1  0.000 0.000
      31    1   1   1   1   0     0     1  0.000 0.000
      

---

    Code
      cLC
    Output
      
      M1: DEV*~URB*LIT*STB + DEV*LIT*IND*STB <-> SURV
      

---

    Code
      minimize(ttLC, details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV*~URB*LIT*STB + DEV*LIT*IND*STB <-> SURV
      
                           inclS   PRI   covS   covU   cases 
      ------------------------------------------------------------------- 
      1  DEV*~URB*LIT*STB  1.000  1.000  0.500  0.250  FI,IE; FR,SE 
      2   DEV*LIT*IND*STB  1.000  1.000  0.750  0.500  FR,SE; BE,CZ,NL,UK 
      ------------------------------------------------------------------- 
                       M1  1.000  1.000  1.000 
      

---

    Code
      minimize(ttLCn)
    Output
      
      M1: ~DEV*~URB*~IND + DEV*LIT*IND*~STB <-> ~SURV
      

---

    Code
      minimize(ttLCn, include = "?", details = TRUE, show.cases = TRUE)
    Output
      
      M1: ~DEV + ~STB <-> ~SURV
      
               inclS   PRI   covS   covU   cases 
      --------------------------------------------------------------- 
      1  ~DEV  1.000  1.000  0.800  0.300  GR,PT,ES; IT,RO; HU,PL; EE 
      2  ~STB  1.000  1.000  0.700  0.200  GR,PT,ES; HU,PL; AU; DE 
      --------------------------------------------------------------- 
           M1  1.000  1.000  1.000 
      

---

    Code
      pLC
    Output
      
      M1: DEV*STB <-> SURV
      
                  inclS   PRI   covS   covU   cases 
      ----------------------------------------------------------------- 
      1  DEV*STB  1.000  1.000  1.000    -    FI,IE; FR,SE; BE,CZ,NL,UK 
      ----------------------------------------------------------------- 
              M1  1.000  1.000  1.000 
      

---

    Code
      pLC$SA
    Output
      $M1
         DEV URB LIT IND STB
      18   1   0   0   0   1
      20   1   0   0   1   1
      26   1   1   0   0   1
      28   1   1   0   1   1
      30   1   1   1   0   1
      

---

    Code
      minimize(ttLM, details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV[2]*LIT[1]*IND[1] + DEV[1]*URB[0]*LIT[1]*IND[0] -> SURV[1]
      
                                      inclS   PRI   covS   covU   cases 
      --------------------------------------------------------------------------- 
      1         DEV[2]*LIT[1]*IND[1]  1.000  1.000  0.625  0.625  FR,SE; BE,NL,UK 
      2  DEV[1]*URB[0]*LIT[1]*IND[0]  1.000  1.000  0.250  0.250  FI,IE 
      --------------------------------------------------------------------------- 
                                  M1  1.000  1.000  0.875 
      

---

    Code
      minimize(ttLM, include = "?", details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV[2] + DEV[1]*IND[0] -> SURV[1]
      
                        inclS   PRI   covS   covU   cases 
      ------------------------------------------------------------- 
      1         DEV[2]  1.000  1.000  0.625  0.625  FR,SE; BE,NL,UK 
      2  DEV[1]*IND[0]  1.000  1.000  0.250  0.250  FI,IE 
      ------------------------------------------------------------- 
                    M1  1.000  1.000  0.875 
      

---

    Code
      minimize(ttLMn, details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV[0]*URB[0]*IND[0] + DEV[1]*URB[0]*LIT[1]*IND[1] -> ~SURV[1]
      
                                      inclS   PRI   covS   covU  
      ---------------------------------------------------------- 
      1         DEV[0]*URB[0]*IND[0]  1.000  1.000  0.800  0.800 
      2  DEV[1]*URB[0]*LIT[1]*IND[1]  1.000  1.000  0.100  0.100 
      ---------------------------------------------------------- 
                                  M1  1.000  1.000  0.900 
      
                                      cases 
      ------------------------------------- 
      1         DEV[0]*URB[0]*IND[0]  GR,IT,PT,RO,ES; EE,HU,PL
      2  DEV[1]*URB[0]*LIT[1]*IND[1]  AU
      ------------------------------------- 
      

---

    Code
      minimize(ttLMn, include = "?", details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV[0] + DEV[1]*URB[0]*IND[1] -> ~SURV[1]
      
                               inclS   PRI   covS   covU   cases 
      ----------------------------------------------------------------------------- 
      1                DEV[0]  1.000  1.000  0.800  0.800  GR,IT,PT,RO,ES; EE,HU,PL 
      2  DEV[1]*URB[0]*IND[1]  1.000  1.000  0.100  0.100  AU 
      ----------------------------------------------------------------------------- 
                           M1  1.000  1.000  0.900 
      

---

    Code
      minimize(ttLF, details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV*~URB*LIT*STB + DEV*LIT*IND*STB -> SURV
      
                           inclS   PRI   covS   covU   cases 
      ------------------------------------------------------------------- 
      1  DEV*~URB*LIT*STB  0.809  0.761  0.433  0.196  FI,IE; FR,SE 
      2   DEV*LIT*IND*STB  0.843  0.821  0.622  0.385  FR,SE; BE,CZ,NL,UK 
      ------------------------------------------------------------------- 
                       M1  0.871  0.851  0.818 
      

---

    Code
      minimize(ttLF, include = "?", details = TRUE, show.cases = TRUE)
    Output
      
      M1: DEV*STB -> SURV
      
                  inclS   PRI   covS   covU   cases 
      ----------------------------------------------------------------- 
      1  DEV*STB  0.869  0.848  0.824    -    FI,IE; FR,SE; BE,CZ,NL,UK 
      ----------------------------------------------------------------- 
              M1  0.869  0.848  0.824 
      

---

    Code
      minimize(ttLF, include = "?", details = TRUE, show.cases = TRUE, dir.exp = "1,1,1,1,1")
    Output
      
      From C1P1: 
      
      M1:    DEV*LIT*STB -> SURV 
      
                      inclS   PRI   covS   covU   cases 
      --------------------------------------------------------------------- 
      1  DEV*LIT*STB  0.869  0.848  0.824    -    FI,IE; FR,SE; BE,CZ,NL,UK 
      --------------------------------------------------------------------- 
                  M1  0.869  0.848  0.824 
      

---

    Code
      pCVF
    Output
      
      M1: ~NATPRIDE + DEMOC*GEOCON*POLDIS + (~DEMOC*ETHFRACT*POLDIS +
          DEMOC*ETHFRACT*GEOCON) -> PROTEST 
      M2: ~NATPRIDE + DEMOC*GEOCON*POLDIS + (~DEMOC*ETHFRACT*POLDIS +
          DEMOC*ETHFRACT*~POLDIS) -> PROTEST 
      M3: ~NATPRIDE + DEMOC*GEOCON*POLDIS + (DEMOC*ETHFRACT*GEOCON +
          ETHFRACT*GEOCON*POLDIS) -> PROTEST 
      M4: ~NATPRIDE + DEMOC*GEOCON*POLDIS + (DEMOC*ETHFRACT*~POLDIS +
          ETHFRACT*GEOCON*POLDIS) -> PROTEST 
      
      --------------------------------------------------------------------------------- 
                                 inclS   PRI   covS   covU   (M1)   (M2)   (M3)   (M4)  
      --------------------------------------------------------------------------------- 
      1               ~NATPRIDE  0.899  0.807  0.597  0.121  0.132  0.122  0.136  0.126 
      2     DEMOC*GEOCON*POLDIS  0.906  0.805  0.342  0.065  0.065  0.070  0.065  0.065 
      --------------------------------------------------------------------------------- 
      3  ~DEMOC*ETHFRACT*POLDIS  0.842  0.718  0.299  0.000  0.040  0.040               
      4   DEMOC*ETHFRACT*GEOCON  0.935  0.826  0.480  0.000  0.085         0.085        
      5  DEMOC*ETHFRACT*~POLDIS  0.932  0.773  0.417  0.000         0.085         0.085 
      6  ETHFRACT*GEOCON*POLDIS  0.869  0.786  0.365  0.005                0.045  0.045 
      --------------------------------------------------------------------------------- 
      
                                 cases 
      -------------------------------- 
      1               ~NATPRIDE  CrimRussiansUkr,RussiansUkraine; HungariansYugo,KosovoAlbanians;
                              RussiansLatvia; BasquesSpain; AlbaniansFYROM
      2     DEMOC*GEOCON*POLDIS  HungariansRom,CatholicsNIreland; AlbaniansFYROM; RussiansEstonia
      -------------------------------- 
      3  ~DEMOC*ETHFRACT*POLDIS  HungariansYugo,KosovoAlbanians; GagauzMoldova
      4   DEMOC*ETHFRACT*GEOCON  BasquesSpain; SerbsFYROM,CatalansSpain; AlbaniansFYROM; RussiansEstonia
      5  DEMOC*ETHFRACT*~POLDIS  BasquesSpain; SerbsFYROM,CatalansSpain
      6  ETHFRACT*GEOCON*POLDIS  HungariansYugo,KosovoAlbanians; GagauzMoldova; AlbaniansFYROM;
                              RussiansEstonia
      -------------------------------- 
      

---

    Code
      pCVF$PIchart
    Output
      
                               5  15 16 24 27 29 30 31 32
      ~NATPRIDE                x  x  -  -  x  x  -  x  - 
      ~DEMOC*ETHFRACT*POLDIS   -  x  x  -  -  -  -  -  - 
      DEMOC*ETHFRACT*GEOCON    -  -  -  -  -  x  x  x  x 
      DEMOC*ETHFRACT*~POLDIS   -  -  -  -  -  x  x  -  - 
      DEMOC*GEOCON*POLDIS      -  -  -  x  -  -  -  x  x 
      ETHFRACT*GEOCON*POLDIS   -  x  x  -  -  -  -  x  x 
      

---

    Code
      minimize(ttCVF, include = "?", row.dom = TRUE, details = TRUE, show.cases = TRUE)
    Output
      
      M1: ~NATPRIDE + DEMOC*ETHFRACT*GEOCON + DEMOC*GEOCON*POLDIS +
          ETHFRACT*GEOCON*POLDIS -> PROTEST
      
                                 inclS   PRI   covS   covU  
      ----------------------------------------------------- 
      1               ~NATPRIDE  0.899  0.807  0.597  0.136 
      2   DEMOC*ETHFRACT*GEOCON  0.935  0.826  0.480  0.085 
      3     DEMOC*GEOCON*POLDIS  0.906  0.805  0.342  0.065 
      4  ETHFRACT*GEOCON*POLDIS  0.869  0.786  0.365  0.045 
      ----------------------------------------------------- 
                             M1  0.879  0.782  0.810 
      
                                 cases 
      -------------------------------- 
      1               ~NATPRIDE  CrimRussiansUkr,RussiansUkraine; HungariansYugo,KosovoAlbanians;
                                 RussiansLatvia; BasquesSpain; AlbaniansFYROM
      2   DEMOC*ETHFRACT*GEOCON  BasquesSpain; SerbsFYROM,CatalansSpain; AlbaniansFYROM;
                                 RussiansEstonia
      3     DEMOC*GEOCON*POLDIS  HungariansRom,CatholicsNIreland; AlbaniansFYROM; RussiansEstonia
      4  ETHFRACT*GEOCON*POLDIS  HungariansYugo,KosovoAlbanians; GagauzMoldova; AlbaniansFYROM;
                                 RussiansEstonia
      -------------------------------- 
      

---

    Code
      minimize(ttCVF, include = "?", use.letters = TRUE)
    Output
      
      M1: ~E + A*C*D + (~A*B*D + A*B*C) -> PROTEST 
      M2: ~E + A*C*D + (~A*B*D + A*B*~D) -> PROTEST 
      M3: ~E + A*C*D + (A*B*C + B*C*D) -> PROTEST 
      M4: ~E + A*C*D + (A*B*~D + B*C*D) -> PROTEST 
      

---

    Code
      minimize(RS, outcome = "REC", details = TRUE, show.cases = TRUE)
    Output
      
      M1: P*E*S + P*E*A*EBA + E*A*~EBA*S <-> REC
      
                     inclS   PRI   covS   covU   cases 
      ---------------------------------------------------- 
      1       P*E*S  1.000  1.000  0.571  0.143  6; 3; 1,2 
      2   P*E*A*EBA  1.000  1.000  0.571  0.286  4,5; 1,2 
      3  E*A*~EBA*S  1.000  1.000  0.286  0.143  13; 3 
      ---------------------------------------------------- 
                 M1  1.000  1.000  1.000 
      

---

    Code
      minimize(ttLM2, details = TRUE, include = "?")
    Output
      
      M1: LIT[1]*STB[1] <-> SURV[1]
      
                        inclS   PRI   covS   covU   cases 
      -------------------------------------------------------------------- 
      1  LIT[1]*STB[1]  0.889  0.889  1.000    -    FI,IE; FR,SE; BE,NL,UK 
      -------------------------------------------------------------------- 
                    M1  0.889  0.889  1.000 
      

