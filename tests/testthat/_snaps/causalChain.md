# tests for causalChain() have the same output

    Code
      cc$E$IC
    Output
                                  ------------------- 
             inclS   PRI   covS   covU   (M1)   (M2)  
      ----------------------------------------------- 
      1   G  1.000  1.000  0.571  0.143  0.143  0.143 
      ----------------------------------------------- 
      2   U  1.000  1.000  0.571  0.000         0.143 
      3   D  1.000  1.000  0.571  0.000         0.143 
      4   L  1.000  1.000  0.857  0.000  0.429        
      ----------------------------------------------- 
         M1  1.000  1.000  1.000 
         M2  1.000  1.000  1.000 
      

---

    Code
      causalChain(d.women)
    Output
      
      M1: WS + ~ES*LP + ES*WM + QU*LP <-> WNP
      M2: WS + ES*WM + QU*LP + WM*LP <-> WNP
      

---

    Code
      causalChain(d.pban, ordering = "C, F, T, V < PB", sol.cov = 0.95)
    Output
      
      M01: C[1] + F[2] + C[0]*F[1] + C[2]*V[0] <-> PB[1]
      M02: C[1] + F[2] + C[0]*T[2] + C[2]*V[0] <-> PB[1]
      M03: C[1] + F[2] + C[0]*F[1] + C[2]*F[0] + F[1]*V[0] <-> PB[1]
      M04: C[1] + F[2] + C[0]*F[1] + C[2]*T[1] + T[2]*V[0] <-> PB[1]
      M05: C[1] + F[2] + C[0]*F[1] + T[1]*V[0] + T[2]*V[0] <-> PB[1]
      M06: C[1] + F[2] + C[0]*T[2] + C[2]*F[0] + F[1]*V[0] <-> PB[1]
      M07: C[1] + F[2] + C[0]*T[2] + C[2]*T[1] + T[2]*V[0] <-> PB[1]
      M08: C[1] + F[2] + C[0]*T[2] + T[1]*V[0] + T[2]*V[0] <-> PB[1]
      M09: C[1] + F[2] + C[0]*F[1] + C[2]*F[0] + F[1]*T[1] + T[2]*V[0] <-> PB[1]
      M10: C[1] + F[2] + C[0]*F[1] + C[2]*T[1] + F[0]*T[2] + F[1]*V[0] <-> PB[1]
      M11: C[1] + F[2] + C[0]*F[1] + F[0]*T[2] + F[1]*V[0] + T[1]*V[0] <-> PB[1]
      M12: C[1] + F[2] + C[0]*T[2] + C[2]*F[0] + F[1]*T[1] + T[2]*V[0] <-> PB[1]
      M13: C[1] + F[2] + C[0]*T[2] + C[2]*T[1] + F[0]*T[2] + F[1]*V[0] <-> PB[1]
      M14: C[1] + F[2] + C[0]*T[2] + F[0]*T[2] + F[1]*V[0] + T[1]*V[0] <-> PB[1]
      

---

    Code
      causalChain(d.pban, ordering = "C, F, T, V < PB", pi.cons = 0.93, sol.cons = 0.95)
    Output
      
      M1: C[1] + F[2] + T[2] + C[2]*T[1] <-> PB[1]
      M2: C[1] + F[2] + T[2] + C[2]*F[0] + F[1]*T[1] <-> PB[1]
      

---

    Code
      causalChain(dat2, ordering = "AU", sol.cons = 0.9, pi.cons = 0.85, sol.cov = 0.85)
    Output
      
      M1: ~RE*CN + RE*~CN <-> AU
      M2: ~RE*DE + ~CN*DE <-> AU
      

