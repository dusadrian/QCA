# tests for solveChart() have the same output

    Code
      chart1
    Output
      
           A~BCD   A~BC~D  ~AB~C~D ~ABCD  
      A       x       x       -       -   
      B       -       -       x       x   
      ~C      -       -       x       -   
      D       x       -       -       x   
      

---

    Code
      solveChart(chart1)
    Output
           [,1]
      [1,]    1
      [2,]    2

---

    Code
      solveChart(chart1, type = "hybrid")
    Output
           [,1]
      [1,]    1
      [2,]    2

---

    Code
      solveChart(chart1, all.sol = TRUE)
    Output
           [,1] [,2]
      [1,]    1    1
      [2,]    2    3
      [3,]    0    4

---

    Code
      chart2
    Output
      
           ABC   A~B~C A~BC  ~AB~C
      A      x     x     x     -  
      B      x     -     -     x  
      ~C     -     x     -     x  
      

---

    Code
      solveChart(chart2)
    Output
           [,1] [,2]
      [1,]    1    1
      [2,]    2    3

---

    Code
      solveChart(chart3)
    Output
           [,1] [,2] [,3] [,4]
      [1,]    1    1    2    2
      [2,]    3    3    3    3
      [3,]    4    4    4    4
      [4,]    5    6    5    6

---

    Code
      solveChart(chart4)
    Output
           [,1]
      [1,]    1
      [2,]    3

