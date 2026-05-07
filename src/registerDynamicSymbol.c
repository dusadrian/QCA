/*
Copyright (c) 2016 - 2026, Adrian Dusa
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, in whole or in part, are permitted provided that the
following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The names of its contributors may NOT be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ADRIAN DUSA BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "qca_rinternals.h"
#include "qca_rinternals.h"
#include <R_ext/Rdynload.h>
#include "solvePIchart_gurobi.h"
extern SEXP C_findminLpSolveInternal(SEXP chart);
extern SEXP C_findminScpInternal(SEXP chart);
extern SEXP C_findminLagrangianInfo(SEXP chart);
extern SEXP C_getScpProfile(void);
extern SEXP C_resetScpProfile(void);
static const R_CallMethodDef CallEntries[] = {
  {"C_findminLpSolveInternal", (DL_FUNC) &C_findminLpSolveInternal, 1},
  {"C_findminScpInternal", (DL_FUNC) &C_findminScpInternal, 1},
  {"C_findminLagrangianInfo", (DL_FUNC) &C_findminLagrangianInfo, 1},
  {"C_getScpProfile", (DL_FUNC) &C_getScpProfile, 0},
  {"C_resetScpProfile", (DL_FUNC) &C_resetScpProfile, 0},
  {NULL, NULL, 0}
};
void R_init_QCA(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
void R_unload_QCA(DllInfo* info) {
  (void) info;
  gurobi_release_env();
}
