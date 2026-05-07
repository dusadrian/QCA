#include <string.h>

#include "commonlib.h"
#include "lp_lib.h"
#include "lp_presolve.h"
#include "lp_crash.h"
#include "lp_report.h"

STATIC MYBOOL presolve_createUndo(lprec *lp)
{
  if(lp->presolve_undo != NULL)
    presolve_freeUndo(lp);

  lp->presolve_undo = (presolveundorec *) calloc(1, sizeof(presolveundorec));
  if(lp->presolve_undo == NULL)
    return(FALSE);

  lp->presolve_undo->lp = lp;
  return(TRUE);
}

STATIC MYBOOL inc_presolve_space(lprec *lp, int delta, MYBOOL isrows)
{
  int i, ii, rowcolsum, oldrowalloc, oldcolalloc, oldrowcolalloc;
  presolveundorec *psundo = lp->presolve_undo;

  if(psundo == NULL) {
    if(!presolve_createUndo(lp))
      return(FALSE);
    psundo = lp->presolve_undo;
  }

  oldrowalloc = lp->rows_alloc - delta;
  oldcolalloc = lp->columns_alloc - delta;
  oldrowcolalloc = lp->sum_alloc - delta;
  rowcolsum = lp->sum_alloc + 1;

  if(isrows)
    allocREAL(lp, &psundo->fixed_rhs, lp->rows_alloc + 1, AUTOMATIC);
  else
    allocREAL(lp, &psundo->fixed_obj, lp->columns_alloc + 1, AUTOMATIC);
  allocINT(lp, &psundo->var_to_orig, rowcolsum, AUTOMATIC);
  allocINT(lp, &psundo->orig_to_var, rowcolsum, AUTOMATIC);

  if((psundo->var_to_orig == NULL) || (psundo->orig_to_var == NULL) ||
     (isrows && (psundo->fixed_rhs == NULL)) ||
     (!isrows && (psundo->fixed_obj == NULL)))
    return(FALSE);

  ii = isrows ? oldrowalloc + 1 : oldcolalloc + 1;
  for(i = oldrowcolalloc + 1; i < rowcolsum; i++, ii++) {
    psundo->var_to_orig[i] = 0;
    psundo->orig_to_var[i] = 0;
    if(isrows)
      psundo->fixed_rhs[ii] = 0;
    else
      psundo->fixed_obj[ii] = 0;
  }

  return(TRUE);
}

STATIC MYBOOL presolve_setOrig(lprec *lp, int orig_rows, int orig_cols)
{
  presolveundorec *psundo = lp->presolve_undo;

  if(psundo == NULL)
    return(FALSE);

  psundo->orig_rows = orig_rows;
  psundo->orig_columns = orig_cols;
  psundo->orig_sum = orig_rows + orig_cols;

  if(lp->wasPresolved)
    presolve_fillUndo(lp, orig_rows, orig_cols, FALSE);

  return(TRUE);
}

STATIC MYBOOL presolve_fillUndo(lprec *lp, int orig_rows, int orig_cols, MYBOOL setOrig)
{
  int i;
  presolveundorec *psundo = lp->presolve_undo;

  if(psundo == NULL)
    return(FALSE);

  if((psundo->var_to_orig == NULL) || (psundo->orig_to_var == NULL) ||
     (psundo->fixed_rhs == NULL) || (psundo->fixed_obj == NULL))
    return(FALSE);

  for(i = 0; i <= orig_rows; i++) {
    psundo->var_to_orig[i] = i;
    psundo->orig_to_var[i] = i;
    psundo->fixed_rhs[i] = 0;
  }
  for(i = 1; i <= orig_cols; i++) {
    psundo->var_to_orig[orig_rows + i] = i;
    psundo->orig_to_var[orig_rows + i] = i;
    psundo->fixed_obj[i] = 0;
  }

  if(setOrig)
    presolve_setOrig(lp, orig_rows, orig_cols);

  return(TRUE);
}

STATIC MYBOOL presolve_rebuildUndo(lprec *lp, MYBOOL isprimal)
{
  (void) lp;
  (void) isprimal;
  return(TRUE);
}

STATIC MYBOOL presolve_freeUndo(lprec *lp)
{
  presolveundorec *psundo = lp->presolve_undo;

  if(psundo == NULL)
    return(FALSE);

  FREE(psundo->orig_to_var);
  FREE(psundo->var_to_orig);
  FREE(psundo->fixed_rhs);
  FREE(psundo->fixed_obj);
  if(psundo->deletedA != NULL)
    freeUndoLadder(&(psundo->deletedA));
  if(psundo->primalundo != NULL)
    freeUndoLadder(&(psundo->primalundo));
  if(psundo->dualundo != NULL)
    freeUndoLadder(&(psundo->dualundo));
  FREE(lp->presolve_undo);

  return(TRUE);
}

STATIC int presolve(lprec *lp)
{
  if(!lp->basis_valid) {
    crash_basis(lp);
    if(!lp->basis_valid)
      report(lp, DETAILED, "presolve: Had to repair broken basis.\n");
  }
  if(!lp->basis_valid)
    lp->var_basic[0] = AUTOMATIC;

  lp->wasPresolved = FALSE;
  lp->timepresolved = timeNow();
  return(RUNNING);
}

STATIC MYBOOL postsolve(lprec *lp, int status)
{
  if(lp->lag_status != RUNNING) {
    int itemp;

    if(status == PRESOLVED)
      status = OPTIMAL;

    if((status == OPTIMAL) || (status == SUBOPTIMAL)) {
      itemp = check_solution(lp, lp->columns, lp->best_solution,
                             lp->orig_upbo, lp->orig_lowbo, lp->epssolution);
      if((itemp != OPTIMAL) && (lp->spx_status == OPTIMAL))
        lp->spx_status = itemp;
      else if((itemp == OPTIMAL) &&
              ((status == SUBOPTIMAL) || (lp->spx_status == PRESOLVED)))
        lp->spx_status = status;
    }
    else if(status != PRESOLVED) {
      report(lp, NORMAL, "lp_solve unsuccessful after %.0f iter and a last best value of %g\n",
             (double) get_total_iter(lp), lp->best_solution[0]);
      if(lp->bb_totalnodes > 0)
        report(lp, NORMAL, "lp_solve explored %.0f nodes before termination\n",
               (double) get_total_nodes(lp));
    }
    else
      lp->spx_status = OPTIMAL;

    presolve_rebuildUndo(lp, TRUE);
  }

  if(varmap_canunlock(lp))
    lp->varmap_locked = FALSE;

  return(TRUE);
}
