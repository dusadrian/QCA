#include "commonlib.h"
#include "lp_lib.h"
#include "lp_pricePSE.h"

STATIC void freePricer(lprec *lp)
{
  if(lp != NULL)
    FREE(lp->edgeVector);
}

STATIC void simplexPricer(lprec *lp, MYBOOL isdual)
{
  (void) lp;
  (void) isdual;
}

STATIC MYBOOL initPricer(lprec *lp)
{
  if(lp != NULL)
    FREE(lp->edgeVector);
  return(TRUE);
}

STATIC MYBOOL resizePricer(lprec *lp)
{
  (void) lp;
  return(TRUE);
}

STATIC REAL getPricer(lprec *lp, int item, MYBOOL isdual)
{
  (void) lp;
  (void) item;
  (void) isdual;
  return(1.0);
}

STATIC MYBOOL restartPricer(lprec *lp, MYBOOL isdual)
{
  (void) lp;
  (void) isdual;
  return(TRUE);
}

STATIC MYBOOL updatePricer(lprec *lp, int rownr, int colnr, REAL *pcol, REAL *prow, int *nzprow)
{
  (void) lp;
  (void) rownr;
  (void) colnr;
  (void) pcol;
  (void) prow;
  (void) nzprow;
  return(TRUE);
}

STATIC MYBOOL verifyPricer(lprec *lp)
{
  (void) lp;
  return(TRUE);
}
