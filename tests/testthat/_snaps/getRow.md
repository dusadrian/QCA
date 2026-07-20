# tests for getRow() have the same output

    Code
      mat
    Output
         A B C
      2  0 0 1
      4  0 1 0
      5  0 1 1
      7  0 2 0
      8  0 2 1
      10 1 0 0
      11 1 0 1
      13 1 1 0
      14 1 1 1
      16 1 2 0
      17 1 2 1

---

    Code
      getRow(row.no = 2, noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]   -1   -1    0

---

    Code
      getRow(row.no = rows2, noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]   -1    0    0
      [2,]   -1    1    0
      [3,]    0   -1    0
      [4,]    0    0    0
      [5,]    0    1    0
      [6,]    1   -1    0
      [7,]    1    0    0
      [8,]    1    1    0

---

    Code
      getRow(row.no = rows2_20, noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]   -1    0    0
      [2,]   -1    1    0
      [3,]    0   -1    0
      [4,]    0    0    0
      [5,]    0    1    0
      [6,]    1   -1    0

---

    Code
      getRow(row.no = rows879, noflevels = rep(3, 4)) - 1
    Output
            [,1] [,2] [,3] [,4]
       [1,]   -1    0    1    0
       [2,]   -1    1    1    0
       [3,]    0   -1    1    0
       [4,]    0    0    1    0
       [5,]    0    1    1    0
       [6,]    1   -1    1    0
       [7,]    1    0    1    0
       [8,]    1    1    1    0
       [9,]    1    1    1    1

---

    Code
      getRow(row.no = 14, noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]    0    0    0

---

    Code
      getRow(row.no = rows14, noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]   -1   -1    0
      [2,]   -1    0   -1
      [3,]   -1    0    0
      [4,]    0   -1   -1
      [5,]    0   -1    0
      [6,]    0    0   -1
      [7,]    0    0    0

---

    Code
      getRow(row.no = c(14, 17), noflevels = noflevels + 1) - 1
    Output
           [,1] [,2] [,3]
      [1,]    0    0    0
      [2,]    0    1    0

---

    Code
      getRow(row.no = rows1417, noflevels = noflevels + 1) - 1
    Output
            [,1] [,2] [,3]
       [1,]   -1   -1    0
       [2,]   -1    0   -1
       [3,]   -1    0    0
       [4,]   -1    1   -1
       [5,]   -1    1    0
       [6,]    0   -1   -1
       [7,]    0   -1    0
       [8,]    0    0   -1
       [9,]    0    0    0
      [10,]    0    1   -1
      [11,]    0    1    0

