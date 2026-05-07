#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "lp_lib.h"
#include "lp_report.h"

#include <R.h>
#include <R_ext/Print.h>

#if defined _MSC_VER
# define vsnprintf _vsnprintf
#endif

char * __VACALL explain(lprec *lp, char *format, ...)
{
  char buff[DEF_STRBUFSIZE + 1];
  va_list ap;

  va_start(ap, format);
  vsnprintf(buff, DEF_STRBUFSIZE, format, ap);
  va_end(ap);

  allocCHAR(lp, &(lp->ex_status), (int) strlen(buff), AUTOMATIC);
  strcpy(lp->ex_status, buff);
  return(lp->ex_status);
}

void __VACALL report(lprec *lp, int level, char *format, ...)
{
  (void) lp;
  (void) level;
  (void) format;
}

void debug_print(lprec *lp, char *format, ...)
{
  va_list ap;

  if((lp == NULL) || !lp->bb_trace)
    return;

  va_start(ap, format);
  if(lp->debuginfo != NULL) {
    char buff[DEF_STRBUFSIZE + 1];
    vsnprintf(buff, DEF_STRBUFSIZE, format, ap);
    lp->debuginfo(lp, lp->loghandle, buff);
  }
  va_end(ap);
}

void debug_print_solution(lprec *lp)
{
  (void) lp;
}

void debug_print_bounds(lprec *lp, REAL *upbo, REAL *lowbo)
{
  (void) lp;
  (void) upbo;
  (void) lowbo;
}

void blockWriteLREAL(FILE *output, char *label, LREAL *vector, int first, int last)
{
  (void) output;
  (void) label;
  (void) vector;
  (void) first;
  (void) last;
}

void blockWriteAMAT(FILE *output, const char *label, lprec *lp, int first, int last)
{
  (void) output;
  (void) label;
  (void) lp;
  (void) first;
  (void) last;
}

void blockWriteBMAT(FILE *output, const char *label, lprec *lp, int first, int last)
{
  (void) output;
  (void) label;
  (void) lp;
  (void) first;
  (void) last;
}

MYBOOL REPORT_debugdump(lprec *lp, char *filename, MYBOOL livedata)
{
  (void) lp;
  (void) filename;
  (void) livedata;
  return(FALSE);
}

void REPORT_objective(lprec *lp)
{
  (void) lp;
}

void REPORT_solution(lprec *lp, int columns)
{
  (void) lp;
  (void) columns;
}

void REPORT_constraints(lprec *lp, int columns)
{
  (void) lp;
  (void) columns;
}

void REPORT_duals(lprec *lp)
{
  (void) lp;
}

void REPORT_extended(lprec *lp)
{
  (void) lp;
}

void REPORT_constraintinfo(lprec *lp, char *datainfo)
{
  (void) lp;
  (void) datainfo;
}

void REPORT_modelinfo(lprec *lp, MYBOOL doName, char *datainfo)
{
  (void) lp;
  (void) doName;
  (void) datainfo;
}

void REPORT_lp(lprec *lp)
{
  (void) lp;
}

MYBOOL REPORT_tableau(lprec *lp)
{
  (void) lp;
  return(FALSE);
}

void REPORT_scales(lprec *lp)
{
  (void) lp;
}

MYBOOL REPORT_mat_mmsave(lprec *lp, char *filename, int *colndx, MYBOOL includeOF, char *infotext)
{
  (void) lp;
  (void) filename;
  (void) colndx;
  (void) includeOF;
  (void) infotext;
  return(FALSE);
}
