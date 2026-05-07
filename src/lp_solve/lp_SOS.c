#include <stdlib.h>
#include "lp_SOS.h"

STATIC SOSgroup *create_SOSgroup(lprec *lp)
{
  return NULL;
}

STATIC void resize_SOSgroup(SOSgroup *group)
{
}

STATIC int append_SOSgroup(SOSgroup *group, SOSrec *SOS)
{
  return 0;
}

STATIC int clean_SOSgroup(SOSgroup *group, MYBOOL forceupdatemap)
{
  return 0;
}

STATIC void free_SOSgroup(SOSgroup **group)
{
  if(group != NULL)
    *group = NULL;
}

STATIC SOSrec *create_SOSrec(SOSgroup *group, char *name, int type, int priority, int size, int *variables, REAL *weights)
{
  return NULL;
}

STATIC MYBOOL delete_SOSrec(SOSgroup *group, int sosindex)
{
  return FALSE;
}

STATIC int append_SOSrec(SOSrec *SOS, int size, int *variables, REAL *weights)
{
  return 0;
}

STATIC void free_SOSrec(SOSrec *SOS)
{
}

STATIC int make_SOSchain(lprec *lp, MYBOOL forceresort)
{
  return 0;
}

STATIC int SOS_member_updatemap(SOSgroup *group)
{
  return 0;
}

STATIC MYBOOL SOS_member_sortlist(SOSgroup *group, int sosindex)
{
  return FALSE;
}

STATIC MYBOOL SOS_shift_col(SOSgroup *group, int sosindex, int column, int delta, LLrec *usedmap, MYBOOL forceresort)
{
  return FALSE;
}

int SOS_member_delete(SOSgroup *group, int sosindex, int member)
{
  return 0;
}

int SOS_get_type(SOSgroup *group, int sosindex)
{
  return 0;
}

int SOS_infeasible(SOSgroup *group, int sosindex)
{
  return 0;
}

int SOS_member_index(SOSgroup *group, int sosindex, int member)
{
  return 0;
}

int SOS_member_count(SOSgroup *group, int sosindex)
{
  return 0;
}

int SOS_memberships(SOSgroup *group, int column)
{
  return 0;
}

int *SOS_get_candidates(SOSgroup *group, int sosindex, int column, MYBOOL excludetarget, REAL *upbound, REAL *lobound)
{
  return NULL;
}

int SOS_is_member(SOSgroup *group, int sosindex, int column)
{
  return 0;
}

MYBOOL SOS_is_member_of_type(SOSgroup *group, int column, int sostype)
{
  return FALSE;
}

MYBOOL SOS_set_GUB(SOSgroup *group, int sosindex, MYBOOL state)
{
  return FALSE;
}

MYBOOL SOS_is_GUB(SOSgroup *group, int sosindex)
{
  return FALSE;
}

MYBOOL SOS_is_marked(SOSgroup *group, int sosindex, int column)
{
  return FALSE;
}

MYBOOL SOS_is_active(SOSgroup *group, int sosindex, int column)
{
  return FALSE;
}

MYBOOL SOS_is_full(SOSgroup *group, int sosindex, int column, MYBOOL activeonly)
{
  return FALSE;
}

MYBOOL SOS_can_activate(SOSgroup *group, int sosindex, int column)
{
  return FALSE;
}

MYBOOL SOS_set_marked(SOSgroup *group, int sosindex, int column, MYBOOL asactive)
{
  return FALSE;
}

MYBOOL SOS_unmark(SOSgroup *group, int sosindex, int column)
{
  return FALSE;
}

int SOS_fix_unmarked(SOSgroup *group, int sosindex, int variable, REAL *bound, REAL value,
                     MYBOOL isupper, int *diffcount, DeltaVrec *changelog)
{
  return 0;
}

int SOS_fix_list(SOSgroup *group, int sosindex, int variable, REAL *bound,
                 int *varlist, MYBOOL isleft, DeltaVrec *changelog)
{
  return 0;
}

int SOS_is_satisfied(SOSgroup *group, int sosindex, REAL *solution)
{
  return SOS_COMPLETE;
}

MYBOOL SOS_is_feasible(SOSgroup *group, int sosindex, REAL *solution)
{
  return TRUE;
}
