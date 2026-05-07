#include "lp_crash.h"
#include "lp_lib.h"

MYBOOL crash_basis(lprec *lp)
{
  if(lp == NULL)
    return(FALSE);

  default_basis(lp);
  return(TRUE);
}

MYBOOL guess_basis(lprec *lp, REAL *guessvector, int *basisvector)
{
  (void) lp;
  (void) guessvector;
  (void) basisvector;

  return(FALSE);
}
