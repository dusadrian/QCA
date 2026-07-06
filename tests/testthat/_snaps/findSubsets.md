# tests for findSubsets() have the same output

    Code
      findSubsets(input = 2, noflevels = noflevels + 1)
    Output
      [1]  5  8 11 14 17 20 23 26

---

    Code
      findSubsets(input = 2, noflevels = noflevels + 1, stop = 20)
    Output
      [1]  5  8 11 14 17 20

---

    Code
      findSubsets(input = c(8, 79), noflevels = rep(3, 4))
    Output
      [1] 17 26 35 44 53 62 71 80 81

