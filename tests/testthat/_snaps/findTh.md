# tests for findTh() have the same output

    Code
      findTh(x)
    Output
      [1] 66.5

---

    Code
      findTh(x, groups = 3)
    Output
      [1] 19.0 66.5

---

    Code
      findTh(gdp)
    Output
      [1] 8700

---

    Code
      findTh(gdp, n = 2)
    Output
      [1]  8700 18000

---

    Code
      findTh(gdp, n = 2, hclustm = "complete", distm = "euclidean")
    Output
      [1]  8700 18000

