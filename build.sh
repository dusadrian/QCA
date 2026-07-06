echo ${BASH_SOURCE}
BASEDIR=$(dirname "$0")
#echo "$BASEDIR"
cd "$BASEDIR/src" || exit
#cd "${BASH_SOURCE%/src/*}" || exit
R CMD SHLIB QCA.c \
truthTable.c \
CCubes/CCubes.c \
CCubes/consistency.c \
CCubes/consistent_solution.c \
CCubes/consistent_models.c \
CCubes/find_min.c \
CCubes/find_models.c \
CCubes/generate_matrix.c \
CCubes/row_dominance.c \
CCubes/sort_cols.c \
CCubes/sort_matrix.c \
CCubes/super_rows.c \
CCubes/utils.c \
lpSolve/lp_min.c \
lpSolve/colamd.c \
lpSolve/commonlib.c \
lpSolve/ini.c \
lpSolve/isfixedvar.c \
lpSolve/lp_crash.c \
lpSolve/lp_Hash.c \
lpSolve/lp_lib.c \
lpSolve/lp_LUSOL.c \
lpSolve/lp_matrix.c \
lpSolve/lp_MDO.c \
lpSolve/lp_mipbb.c \
lpSolve/lp_MPS.c \
lpSolve/lp_params.c \
lpSolve/lp_presolve.c \
lpSolve/lp_price.c \
lpSolve/lp_pricePSE.c \
lpSolve/lp_report.c \
lpSolve/lp_rlp.c \
lpSolve/lp_scale.c \
lpSolve/lp_simplex.c \
lpSolve/lp_SOS.c \
lpSolve/lp_utils.c \
lpSolve/lp_wlp.c \
lpSolve/lusol.c \
lpSolve/myblas.c \
lpSolve/sparselib.c \
lpSolve/yacc_read.c

find . -type f -name '*.o' -delete
