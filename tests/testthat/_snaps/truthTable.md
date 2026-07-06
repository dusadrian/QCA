# tests for truthTable() have the same output

    Code
      truthTable(d.represent, outcome = WNP)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           ES QU WS WM LP   OUT    n  incl  PRI  
       3   0  0  0  1  0     0     2  0.000 0.000
       4   0  0  0  1  1     1     1  1.000 1.000
       9   0  1  0  0  0     0     1  0.000 0.000
      11   0  1  0  1  0     0     4  0.000 0.000
      12   0  1  0  1  1     1     1  1.000 1.000
      18   1  0  0  0  1     0     1  0.000 0.000
      21   1  0  1  0  0     1     1  1.000 1.000
      24   1  0  1  1  1     1     1  1.000 1.000
      25   1  1  0  0  0     0     3  0.000 0.000
      26   1  1  0  0  1     1     1  1.000 1.000
      27   1  1  0  1  0     1     1  1.000 1.000
      28   1  1  0  1  1     1     2  1.000 1.000
      29   1  1  1  0  0     1     1  1.000 1.000
      32   1  1  1  1  1     1     2  1.000 1.000
      

---

    Code
      truthTable(d.represent, outcome = WNP, complete = TRUE, show.cases = TRUE,
        sort.by = "incl = TRUE, n = FALSE")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           ES QU WS WM LP   OUT    n  incl  PRI   cases      
       4   0  0  0  1  1     1     1  1.000 1.000 NZ         
      12   0  1  0  1  1     1     1  1.000 1.000 DE         
      21   1  0  1  0  0     1     1  1.000 1.000 FI         
      24   1  0  1  1  1     1     1  1.000 1.000 DK         
      26   1  1  0  0  1     1     1  1.000 1.000 AT         
      27   1  1  0  1  0     1     1  1.000 1.000 ES         
      29   1  1  1  0  0     1     1  1.000 1.000 SE         
      28   1  1  0  1  1     1     2  1.000 1.000 NL,BE      
      32   1  1  1  1  1     1     2  1.000 1.000 NO,IS      
       9   0  1  0  0  0     0     1  0.000 0.000 IT         
      18   1  0  0  0  1     0     1  0.000 0.000 LU         
       3   0  0  0  1  0     0     2  0.000 0.000 CA,US      
      25   1  1  0  0  0     0     3  0.000 0.000 CH,PT,GR   
      11   0  1  0  1  0     0     4  0.000 0.000 AU,GB,FR,IE
       1   0  0  0  0  0     ?     0    -     -              
       2   0  0  0  0  1     ?     0    -     -              
       5   0  0  1  0  0     ?     0    -     -              
       6   0  0  1  0  1     ?     0    -     -              
       7   0  0  1  1  0     ?     0    -     -              
       8   0  0  1  1  1     ?     0    -     -              
      10   0  1  0  0  1     ?     0    -     -              
      13   0  1  1  0  0     ?     0    -     -              
      14   0  1  1  0  1     ?     0    -     -              
      15   0  1  1  1  0     ?     0    -     -              
      16   0  1  1  1  1     ?     0    -     -              
      17   1  0  0  0  0     ?     0    -     -              
      19   1  0  0  1  0     ?     0    -     -              
      20   1  0  0  1  1     ?     0    -     -              
      22   1  0  1  0  1     ?     0    -     -              
      23   1  0  1  1  0     ?     0    -     -              
      30   1  1  1  0  1     ?     0    -     -              
      31   1  1  1  1  0     ?     0    -     -              
      

---

    Code
      truthTable(d.represent, outcome = WNP, complete = TRUE, show.cases = TRUE,
        sort.by = "incl, n", decreasing = "TRUE, FALSE")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           ES QU WS WM LP   OUT    n  incl  PRI   cases      
      28   1  1  0  1  1     1     2  1.000 1.000 NL,BE      
      32   1  1  1  1  1     1     2  1.000 1.000 NO,IS      
       4   0  0  0  1  1     1     1  1.000 1.000 NZ         
      12   0  1  0  1  1     1     1  1.000 1.000 DE         
      21   1  0  1  0  0     1     1  1.000 1.000 FI         
      24   1  0  1  1  1     1     1  1.000 1.000 DK         
      26   1  1  0  0  1     1     1  1.000 1.000 AT         
      27   1  1  0  1  0     1     1  1.000 1.000 ES         
      29   1  1  1  0  0     1     1  1.000 1.000 SE         
      11   0  1  0  1  0     0     4  0.000 0.000 AU,GB,FR,IE
      25   1  1  0  0  0     0     3  0.000 0.000 CH,PT,GR   
       3   0  0  0  1  0     0     2  0.000 0.000 CA,US      
       9   0  1  0  0  0     0     1  0.000 0.000 IT         
      18   1  0  0  0  1     0     1  0.000 0.000 LU         
       1   0  0  0  0  0     ?     0    -     -              
       2   0  0  0  0  1     ?     0    -     -              
       5   0  0  1  0  0     ?     0    -     -              
       6   0  0  1  0  1     ?     0    -     -              
       7   0  0  1  1  0     ?     0    -     -              
       8   0  0  1  1  1     ?     0    -     -              
      10   0  1  0  0  1     ?     0    -     -              
      13   0  1  1  0  0     ?     0    -     -              
      14   0  1  1  0  1     ?     0    -     -              
      15   0  1  1  1  0     ?     0    -     -              
      16   0  1  1  1  1     ?     0    -     -              
      17   1  0  0  0  0     ?     0    -     -              
      19   1  0  0  1  0     ?     0    -     -              
      20   1  0  0  1  1     ?     0    -     -              
      22   1  0  1  0  1     ?     0    -     -              
      23   1  0  1  1  0     ?     0    -     -              
      30   1  1  1  0  1     ?     0    -     -              
      31   1  1  1  1  0     ?     0    -     -              
      

---

    Code
      Kro.tt
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           ES QU WS WM LP   OUT    n  incl  PRI   cases      
       3   0  0  0  1  0     0     2  0.000 0.000 CA,US      
      11   0  1  0  1  0     0     4  0.000 0.000 AU,GB,FR,IE
      25   1  1  0  0  0     0     3  0.000 0.000 CH,PT,GR   
      28   1  1  0  1  1     1     2  1.000 1.000 NL,BE      
      32   1  1  1  1  1     1     2  1.000 1.000 NO,IS      
      

---

    Code
      Kro.tt$removed
    Output
           ES QU WS WM LP   OUT    n  incl  PRI   cases
       4   0  0  0  1  1     1     1  1.000 1.000  NZ  
       9   0  1  0  0  0     0     1  0.000 0.000  IT  
      12   0  1  0  1  1     1     1  1.000 1.000  DE  
      18   1  0  0  0  1     0     1  0.000 0.000  LU  
      21   1  0  1  0  0     1     1  1.000 1.000  FI  
      24   1  0  1  1  1     1     1  1.000 1.000  DK  
      26   1  1  0  0  1     1     1  1.000 1.000  AT  
      27   1  1  0  1  0     1     1  1.000 1.000  ES  
      29   1  1  1  0  0     1     1  1.000 1.000  SE  
      

---

    Code
      truthTable(d.jobsecurity, outcome = JSR, incl.cut = c(0.9, 0.5))
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           S  C  L  R  P  V    OUT    n  incl  PRI  
       2   0  0  0  0  0  1     0     2  0.198 0.000
       5   0  0  0  1  0  0     C     1  0.581 0.000
      10   0  0  1  0  0  1     0     1  0.494 0.000
      20   0  1  0  0  1  1     C     2  0.716 0.500
      25   0  1  1  0  0  0     C     2  0.839 0.699
      27   0  1  1  0  1  0     1     1  0.940 0.825
      33   1  0  0  0  0  0     0     2  0.203 0.000
      37   1  0  0  1  0  0     1     2  0.977 0.966
      47   1  0  1  1  1  0     1     1  1.000 1.000
      48   1  0  1  1  1  1     1     1  1.000 1.000
      56   1  1  0  1  1  1     1     2  1.000 1.000
      57   1  1  1  0  0  0     C     1  0.717 0.000
      64   1  1  1  1  1  1     1     1  1.000 1.000
      

---

    Code
      truthTable(d.jobsecurity, outcome = ~JSR, incl.cut = c(0.9, 0.5))
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           S  C  L  R  P  V    OUT    n  incl  PRI  
       2   0  0  0  0  0  1     1     2  1.000 1.000
       5   0  0  0  1  0  0     1     1  1.000 1.000
      10   0  0  1  0  0  1     1     1  1.000 1.000
      20   0  1  0  0  1  1     C     2  0.716 0.500
      25   0  1  1  0  0  0     C     2  0.627 0.301
      27   0  1  1  0  1  0     C     1  0.714 0.175
      33   1  0  0  0  0  0     1     2  1.000 1.000
      37   1  0  0  1  0  0     0     2  0.352 0.034
      47   1  0  1  1  1  0     0     1  0.458 0.000
      48   1  0  1  1  1  1     0     1  0.459 0.000
      56   1  1  0  1  1  1     C     2  0.571 0.000
      57   1  1  1  0  0  0     1     1  0.950 0.824
      64   1  1  1  1  1  1     C     1  0.612 0.000
      

---

    Code
      HK.tt
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           C  F  T  V    OUT    n  incl  PRI  
      11   0  1  2  0     1     2  1.000 1.000
      12   0  1  2  1     1     1  1.000 1.000
      15   0  2  1  0     1     1  1.000 1.000
      16   0  2  1  1     1     1  1.000 1.000
      17   0  2  2  0     1     6  1.000 1.000
      18   0  2  2  1     1     3  1.000 1.000
      19   1  0  0  0     1     1  1.000 1.000
      27   1  1  1  0     1     1  1.000 1.000
      28   1  1  1  1     1     1  1.000 1.000
      29   1  1  2  0     1     4  1.000 1.000
      33   1  2  1  0     1     2  1.000 1.000
      34   1  2  1  1     1     1  1.000 1.000
      35   1  2  2  0     1     7  1.000 1.000
      37   2  0  0  0     0     2  0.000 0.000
      38   2  0  0  1     0     1  0.000 0.000
      39   2  0  1  0     1     1  1.000 1.000
      40   2  0  1  1     0     1  0.000 0.000
      41   2  0  2  0     1     1  1.000 1.000
      45   2  1  1  0     1     1  1.000 1.000
      47   2  1  2  0     1     1  1.000 1.000
      48   2  1  2  1     C     4  0.500 0.500
      53   2  2  2  0     1     3  1.000 1.000
      54   2  2  2  1     1     2  1.000 1.000
      

---

    Code
      HK.tt$noflevels
    Output
      [1] 3 3 3 2

---

    Code
      HK.tt$tt[which(HK.tt$tt$n > 2), ]
    Output
         C F T V OUT n incl PRI                cases
      17 0 2 2 0   1 6    1   1    CV,GQ,ER,GW,LR,SO
      18 0 2 2 1   1 3    1   1             BI,CD,ET
      29 1 1 2 0   1 4    1   1          CF,CM,CI,TG
      35 1 2 2 0   1 7    1   1 BF,TD,KM,DJ,GA,GN,MR
      48 2 1 2 1   C 4  0.5 0.5          KE,NG,TZ,ZM
      53 2 2 2 0   1 3    1   1             GH,LS,SZ

---

    Code
      truthTable(d.partybans, outcome = PB, conditions = c(C, F, T), incl.cut = c(0.9,
        0.4), show.cases = TRUE, inf.test = "binom, 0.1")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      pval1: p-value for alternative hypothesis inclusion > 0.9
      pval0: p-value for alternative hypothesis inclusion > 0.4
      
           C  F  T    OUT    n  incl  PRI   pval1 pval0 cases                     
       6   0  1  2     C     3  1.000 1.000 0.729 0.064 RW,ST,SC                  
       8   0  2  1     0     2  1.000 1.000 0.810 0.160 AO,MZ                     
       9   0  2  2     C     9  1.000 1.000 0.387 0.000 BI,CV,CD,GQ,ER,ET,GW,LR,SO
      10   1  0  0     0     1  1.000 1.000 0.900 0.400 SN                        
      14   1  1  1     0     2  1.000 1.000 0.810 0.160 MG,ML                     
      15   1  1  2     C     4  1.000 1.000 0.656 0.026 CF,CM,CI,TG               
      17   1  2  1     C     3  1.000 1.000 0.729 0.064 BJ,CG,NE                  
      18   1  2  2     C     7  1.000 1.000 0.478 0.002 BF,TD,KM,DJ,GA,GN,MR      
      19   2  0  0     0     3  0.000 0.000 1.000 1.000 BW,MU,ZW                  
      20   2  0  1     0     2  0.500 0.500 0.990 0.640 NA,ZA                     
      21   2  0  2     0     1  1.000 1.000 0.900 0.400 GM                        
      23   2  1  1     0     1  1.000 1.000 0.900 0.400 MW                        
      24   2  1  2     0     5  0.600 0.600 0.991 0.317 KE,NG,SL,TZ,ZM            
      27   2  2  2     C     5  1.000 1.000 0.590 0.010 GH,LS,SD,SZ,UG            
      

---

    Code
      truthTable(d.graduate, outcome = REC)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           P  E  A  S  EBA   OUT    n  incl  PRI  
       3   0  0  0  0   -     0     3  0.000 0.000
      15   0  1  0  0   -     0     1  0.000 0.000
      22   0  1  1  1   0     1     1  1.000 1.000
      27   1  0  0  0   -     0     1  0.000 0.000
      30   1  0  0  1   -     0     3  0.000 0.000
      36   1  0  1  1   -     0     2  0.000 0.000
      42   1  1  0  1   -     1     1  1.000 1.000
      44   1  1  1  0   1     1     2  1.000 1.000
      46   1  1  1  1   0     1     1  1.000 1.000
      47   1  1  1  1   1     1     2  1.000 1.000
      

---

    Code
      ttLC
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
       1    0   0   0   0   0     0     3  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       5    0   0   1   0   0     0     2  0.000 0.000
       6    0   0   1   0   1     0     1  0.000 0.000
      22    1   0   1   0   1     1     2  1.000 1.000
      23    1   0   1   1   0     0     1  0.000 0.000
      24    1   0   1   1   1     1     2  1.000 1.000
      31    1   1   1   1   0     0     1  0.000 0.000
      32    1   1   1   1   1     1     4  1.000 1.000
      

---

    Code
      print(ttLC, show.cases = TRUE)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI   cases      
       1    0   0   0   0   0     0     3  0.000 0.000 GR,PT,ES   
       2    0   0   0   0   1     0     2  0.000 0.000 IT,RO      
       5    0   0   1   0   0     0     2  0.000 0.000 HU,PL      
       6    0   0   1   0   1     0     1  0.000 0.000 EE         
      22    1   0   1   0   1     1     2  1.000 1.000 FI,IE      
      23    1   0   1   1   0     0     1  0.000 0.000 AU         
      24    1   0   1   1   1     1     2  1.000 1.000 FR,SE      
      31    1   1   1   1   0     0     1  0.000 0.000 DE         
      32    1   1   1   1   1     1     4  1.000 1.000 BE,CZ,NL,UK
      

---

    Code
      print(ttLC, show.cases = TRUE, complete = TRUE)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI   cases      
       1    0   0   0   0   0     0     3  0.000 0.000 GR,PT,ES   
       2    0   0   0   0   1     0     2  0.000 0.000 IT,RO      
       3    0   0   0   1   0     ?     0    -     -              
       4    0   0   0   1   1     ?     0    -     -              
       5    0   0   1   0   0     0     2  0.000 0.000 HU,PL      
       6    0   0   1   0   1     0     1  0.000 0.000 EE         
       7    0   0   1   1   0     ?     0    -     -              
       8    0   0   1   1   1     ?     0    -     -              
       9    0   1   0   0   0     ?     0    -     -              
      10    0   1   0   0   1     ?     0    -     -              
      11    0   1   0   1   0     ?     0    -     -              
      12    0   1   0   1   1     ?     0    -     -              
      13    0   1   1   0   0     ?     0    -     -              
      14    0   1   1   0   1     ?     0    -     -              
      15    0   1   1   1   0     ?     0    -     -              
      16    0   1   1   1   1     ?     0    -     -              
      17    1   0   0   0   0     ?     0    -     -              
      18    1   0   0   0   1     ?     0    -     -              
      19    1   0   0   1   0     ?     0    -     -              
      20    1   0   0   1   1     ?     0    -     -              
      21    1   0   1   0   0     ?     0    -     -              
      22    1   0   1   0   1     1     2  1.000 1.000 FI,IE      
      23    1   0   1   1   0     0     1  0.000 0.000 AU         
      24    1   0   1   1   1     1     2  1.000 1.000 FR,SE      
      25    1   1   0   0   0     ?     0    -     -              
      26    1   1   0   0   1     ?     0    -     -              
      27    1   1   0   1   0     ?     0    -     -              
      28    1   1   0   1   1     ?     0    -     -              
      29    1   1   1   0   0     ?     0    -     -              
      30    1   1   1   0   1     ?     0    -     -              
      31    1   1   1   1   0     0     1  0.000 0.000 DE         
      32    1   1   1   1   1     1     4  1.000 1.000 BE,CZ,NL,UK
      

---

    Code
      truthTable(LC, SURV, complete = TRUE)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
       1    0   0   0   0   0     0     3  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       3    0   0   0   1   0     ?     0    -     -  
       4    0   0   0   1   1     ?     0    -     -  
       5    0   0   1   0   0     0     2  0.000 0.000
       6    0   0   1   0   1     0     1  0.000 0.000
       7    0   0   1   1   0     ?     0    -     -  
       8    0   0   1   1   1     ?     0    -     -  
       9    0   1   0   0   0     ?     0    -     -  
      10    0   1   0   0   1     ?     0    -     -  
      11    0   1   0   1   0     ?     0    -     -  
      12    0   1   0   1   1     ?     0    -     -  
      13    0   1   1   0   0     ?     0    -     -  
      14    0   1   1   0   1     ?     0    -     -  
      15    0   1   1   1   0     ?     0    -     -  
      16    0   1   1   1   1     ?     0    -     -  
      17    1   0   0   0   0     ?     0    -     -  
      18    1   0   0   0   1     ?     0    -     -  
      19    1   0   0   1   0     ?     0    -     -  
      20    1   0   0   1   1     ?     0    -     -  
      21    1   0   1   0   0     ?     0    -     -  
      22    1   0   1   0   1     1     2  1.000 1.000
      23    1   0   1   1   0     0     1  0.000 0.000
      24    1   0   1   1   1     1     2  1.000 1.000
      25    1   1   0   0   0     ?     0    -     -  
      26    1   1   0   0   1     ?     0    -     -  
      27    1   1   0   1   0     ?     0    -     -  
      28    1   1   0   1   1     ?     0    -     -  
      29    1   1   1   0   0     ?     0    -     -  
      30    1   1   1   0   1     ?     0    -     -  
      31    1   1   1   1   0     0     1  0.000 0.000
      32    1   1   1   1   1     1     4  1.000 1.000
      

---

    Code
      truthTable(LC, SURV, complete = TRUE, sort.by = "incl, n")
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
       3    0   0   0   1   0     ?     0    -     -  
       4    0   0   0   1   1     ?     0    -     -  
       7    0   0   1   1   0     ?     0    -     -  
       8    0   0   1   1   1     ?     0    -     -  
       9    0   1   0   0   0     ?     0    -     -  
      10    0   1   0   0   1     ?     0    -     -  
      11    0   1   0   1   0     ?     0    -     -  
      12    0   1   0   1   1     ?     0    -     -  
      13    0   1   1   0   0     ?     0    -     -  
      14    0   1   1   0   1     ?     0    -     -  
      15    0   1   1   1   0     ?     0    -     -  
      16    0   1   1   1   1     ?     0    -     -  
      17    1   0   0   0   0     ?     0    -     -  
      18    1   0   0   0   1     ?     0    -     -  
      19    1   0   0   1   0     ?     0    -     -  
      20    1   0   0   1   1     ?     0    -     -  
      21    1   0   1   0   0     ?     0    -     -  
      25    1   1   0   0   0     ?     0    -     -  
      26    1   1   0   0   1     ?     0    -     -  
      27    1   1   0   1   0     ?     0    -     -  
      28    1   1   0   1   1     ?     0    -     -  
      29    1   1   1   0   0     ?     0    -     -  
      30    1   1   1   0   1     ?     0    -     -  
      

---

    Code
      truthTable(LC, SURV, complete = TRUE, sort.by = "incl, n = FALSE")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
       6    0   0   1   0   1     0     1  0.000 0.000
      23    1   0   1   1   0     0     1  0.000 0.000
      31    1   1   1   1   0     0     1  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       5    0   0   1   0   0     0     2  0.000 0.000
       1    0   0   0   0   0     0     3  0.000 0.000
      22    1   0   1   0   1     1     2  1.000 1.000
      24    1   0   1   1   1     1     2  1.000 1.000
      32    1   1   1   1   1     1     4  1.000 1.000
       3    0   0   0   1   0     ?     0    -     -  
       4    0   0   0   1   1     ?     0    -     -  
       7    0   0   1   1   0     ?     0    -     -  
       8    0   0   1   1   1     ?     0    -     -  
       9    0   1   0   0   0     ?     0    -     -  
      10    0   1   0   0   1     ?     0    -     -  
      11    0   1   0   1   0     ?     0    -     -  
      12    0   1   0   1   1     ?     0    -     -  
      13    0   1   1   0   0     ?     0    -     -  
      14    0   1   1   0   1     ?     0    -     -  
      15    0   1   1   1   0     ?     0    -     -  
      16    0   1   1   1   1     ?     0    -     -  
      17    1   0   0   0   0     ?     0    -     -  
      18    1   0   0   0   1     ?     0    -     -  
      19    1   0   0   1   0     ?     0    -     -  
      20    1   0   0   1   1     ?     0    -     -  
      21    1   0   1   0   0     ?     0    -     -  
      25    1   1   0   0   0     ?     0    -     -  
      26    1   1   0   0   1     ?     0    -     -  
      27    1   1   0   1   0     ?     0    -     -  
      28    1   1   0   1   1     ?     0    -     -  
      29    1   1   1   0   0     ?     0    -     -  
      30    1   1   1   0   1     ?     0    -     -  
      

---

    Code
      truthTable(LM, SURV, sort.by = "incl")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
      22    1   0   1   0   1     1     2  1.000 1.000
      32    1   1   1   1   1     1     1  1.000 1.000
      40    2   0   1   1   1     1     2  1.000 1.000
      48    2   1   1   1   1     1     3  1.000 1.000
       1    0   0   0   0   0     0     3  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       5    0   0   1   0   0     0     2  0.000 0.000
       6    0   0   1   0   1     0     1  0.000 0.000
      23    1   0   1   1   0     0     1  0.000 0.000
      31    1   1   1   1   0     0     1  0.000 0.000
      

---

    Code
      ttLM
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV URB LIT IND STB   OUT    n  incl  PRI  
      22    1   0   1   0   1     1     2  1.000 1.000
      40    2   0   1   1   1     1     2  1.000 1.000
      48    2   1   1   1   1     1     3  1.000 1.000
       1    0   0   0   0   0     0     3  0.000 0.000
       2    0   0   0   0   1     0     2  0.000 0.000
       5    0   0   1   0   0     0     2  0.000 0.000
      

---

    Code
      ttLM$removed
    Output
           DEV URB LIT IND STB   OUT    n  incl  PRI  
       6    0   0   1   0   1     0     1  0.000 0.000
      23    1   0   1   1   0     0     1  0.000 0.000
      31    1   1   1   1   0     0     1  0.000 0.000
      32    1   1   1   1   1     1     1  1.000 1.000
      

---

    Code
      truthTable(CVF, PROTEST, incl.cut = 0.8, sort.by = "incl")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEMOC ETHFRACT GEOCON POLDIS NATPRIDE   OUT    n  incl  PRI  
      31     1      1       1      1       0        1     1  0.959 0.891
      29     1      1       1      0       0        1     1  0.957 0.858
      27     1      1       0      1       0        1     1  0.938 0.783
      15     0      1       1      1       0        1     2  0.929 0.839
      30     1      1       1      0       1        1     2  0.922 0.667
      24     1      0       1      1       1        1     2  0.904 0.754
      32     1      1       1      1       1        1     1  0.892 0.717
       5     0      0       1      0       0        1     2  0.874 0.584
      16     0      1       1      1       1        1     1  0.821 0.591
      28     1      1       0      1       1        0     1  0.768 0.225
       8     0      0       1      1       1        0     1  0.749 0.248
      10     0      1       0      0       1        0     1  0.741 0.086
       6     0      0       1      0       1        0     1  0.732 0.148
      22     1      0       1      0       1        0     4  0.725 0.281
      20     1      0       0      1       1        0     1  0.711 0.162
      14     0      1       1      0       1        0     3  0.653 0.098
       2     0      0       0      0       1        0     4  0.648 0.182
      

---

    Code
      truthTable(CVF, ~PROTEST, incl.cut = 0.8, sort.by = "incl")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEMOC ETHFRACT GEOCON POLDIS NATPRIDE   OUT    n  incl  PRI  
      10     0      1       0      0       1        1     1  0.976 0.914
       6     0      0       1      0       1        1     1  0.953 0.852
      28     1      1       0      1       1        1     1  0.933 0.775
      14     0      1       1      0       1        1     3  0.929 0.816
       2     0      0       0      0       1        1     4  0.922 0.818
      20     1      0       0      1       1        1     1  0.902 0.714
       8     0      0       1      1       1        1     1  0.902 0.705
      22     1      0       1      0       1        1     4  0.846 0.597
      30     1      1       1      0       1        1     2  0.843 0.333
       5     0      0       1      0       0        1     2  0.823 0.416
      27     1      1       0      1       0        0     1  0.776 0.217
      16     0      1       1      1       1        0     1  0.742 0.409
      29     1      1       1      0       0        0     1  0.738 0.142
      32     1      1       1      1       1        0     1  0.726 0.283
      24     1      0       1      1       1        0     2  0.699 0.229
      31     1      1       1      1       0        0     1  0.662 0.109
      15     0      1       1      1       0        0     2  0.629 0.161
      

---

    Code
      truthTable(CVF, PROTEST, incl.cut = c(0.8, 0.75), sort.by = "incl")
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEMOC ETHFRACT GEOCON POLDIS NATPRIDE   OUT    n  incl  PRI  
      31     1      1       1      1       0        1     1  0.959 0.891
      29     1      1       1      0       0        1     1  0.957 0.858
      27     1      1       0      1       0        1     1  0.938 0.783
      15     0      1       1      1       0        1     2  0.929 0.839
      30     1      1       1      0       1        1     2  0.922 0.667
      24     1      0       1      1       1        1     2  0.904 0.754
      32     1      1       1      1       1        1     1  0.892 0.717
       5     0      0       1      0       0        1     2  0.874 0.584
      16     0      1       1      1       1        1     1  0.821 0.591
      28     1      1       0      1       1        C     1  0.768 0.225
       8     0      0       1      1       1        0     1  0.749 0.248
      10     0      1       0      0       1        0     1  0.741 0.086
       6     0      0       1      0       1        0     1  0.732 0.148
      22     1      0       1      0       1        0     4  0.725 0.281
      20     1      0       0      1       1        0     1  0.711 0.162
      14     0      1       1      0       1        0     3  0.653 0.098
       2     0      0       0      0       1        0     4  0.648 0.182
      

---

    Code
      truthTable(RS, REC)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           P  E  A  EBA S    OUT    n  incl  PRI  
       5   0  0  0   -  0     0     3  0.000 0.000
      17   0  1  0   -  0     0     1  0.000 0.000
      20   0  1  1   0  1     1     1  1.000 1.000
      29   1  0  0   -  0     0     1  0.000 0.000
      30   1  0  0   -  1     0     3  0.000 0.000
      36   1  0  1   -  1     0     2  0.000 0.000
      42   1  1  0   -  1     1     1  1.000 1.000
      44   1  1  1   0  1     1     1  1.000 1.000
      45   1  1  1   1  0     1     2  1.000 1.000
      46   1  1  1   1  1     1     2  1.000 1.000
      

---

    Code
      truthTable(LC, SURV, incl.cut = 0.8, use.labels = TRUE)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV         URB               LIT           IND                STB     
       1   Undeveloped Low urbanization  Low literacy  Not industrialized Unstable
       2   Undeveloped Low urbanization  Low literacy  Not industrialized Stable  
       5   Undeveloped Low urbanization  High literacy Not industrialized Unstable
       6   Undeveloped Low urbanization  High literacy Not industrialized Stable  
      22   Developed   Low urbanization  High literacy Not industrialized Stable  
      23   Developed   Low urbanization  High literacy Industrialized     Unstable
      24   Developed   Low urbanization  High literacy Industrialized     Stable  
      31   Developed   High urbanization High literacy Industrialized     Unstable
      32   Developed   High urbanization High literacy Industrialized     Stable  
             OUT         n  incl  PRI  
       1     Collapse    3  0.000 0.000
       2     Collapse    2  0.000 0.000
       5     Collapse    2  0.000 0.000
       6     Collapse    1  0.000 0.000
      22     Survival    2  1.000 1.000
      23     Collapse    1  0.000 0.000
      24     Survival    2  1.000 1.000
      31     Collapse    1  0.000 0.000
      32     Survival    4  1.000 1.000
      

---

    Code
      truthTable(LC3, SURV, incl.cut = 0.8, use.labels = TRUE)
    Output
      
        OUT: output value
          n: number of cases in configuration
       incl: sufficiency inclusion score
        PRI: proportional reduction in inconsistency
      
           DEV         URB               LIT           IND                STB     
       1   Undeveloped Low urbanization  Low literacy  Not industrialized Unstable
       2   Undeveloped Low urbanization  Low literacy  Not industrialized Stable  
       5   Undeveloped Low urbanization  High literacy Not industrialized Unstable
       6   Undeveloped Low urbanization  High literacy Not industrialized Stable  
      22   Developed   Low urbanization  High literacy Not industrialized Stable  
      23   Developed   Low urbanization  High literacy Industrialized     Unstable
      24   Developed   Low urbanization  High literacy Industrialized     Stable  
      31   Developed   High urbanization High literacy Industrialized     Unstable
      32   Developed   High urbanization High literacy Industrialized     Stable  
             OUT         n  incl  PRI  
       1     Collapse    3  0.000 0.000
       2     Collapse    2  0.000 0.000
       5     Collapse    2  0.000 0.000
       6     Collapse    1  0.000 0.000
      22     Survival    2  1.000 1.000
      23     Collapse    1  0.000 0.000
      24     Survival    2  1.000 1.000
      31     Collapse    1  0.000 0.000
      32     Survival    4  1.000 1.000
      

