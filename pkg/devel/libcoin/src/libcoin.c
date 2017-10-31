
/* C Header */

/*
    TO NOT EDIT THIS FILE
    
    Edit `libcoin.w' and run `nuweb -r libcoin.w'
*/

#include "libcoin_internal.h"
#include <R_ext/stats_stubs.h>  /* for S_rcont2 */
#include <mvtnormAPI.h>         /* for calling mvtnorm */
/* Function Definitions */

/* MoreUtils */

/* NROW */

int NROW
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return(INTEGER(a)[0]);
}

/* NCOL */

int NCOL
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return(INTEGER(a)[1]);
}

/* NLEVELS */

int NLEVELS     
(
    SEXP x
) {

    SEXP a;
    int maxlev = 0;

    a = getAttrib(x, R_LevelsSymbol);
    if (a == R_NilValue) {
        Rprintf("no levels!");
        if (TYPEOF(x) != INTSXP)
            error("cannot determine number of levels");
        for (R_xlen_t i = 0; i < XLENGTH(x); i++) {
Rprintf("i %d x %d ", i, INTEGER(x)[i]);
            if (INTEGER(x)[i] > maxlev)
                maxlev = INTEGER(x)[i];
        }
        Rprintf("maxlev %d \n", maxlev);
        return(maxlev);
    }
    return(NROW(a));
}

/* C\_kronecker */

void C_kronecker
(
    const double *A,
    const int m,
    const int n,
    const double *B,
    const int r,
    const int s,
    const int overwrite,
    double *ans
) {

    int i, j, k, l, mr, js, ir;
    double y;

    if (overwrite) {
        for (i = 0; i < m * r * n * s; i++) ans[i] = 0.0;
    }

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++)
                    ans[(js + l) * mr + ir + k] += y * B[l * r + k];
            }
        }
    }
}

/* C\_kronecker\_sym */

void C_kronecker_sym
(
    const double *A,
    const int m,
    const double *B,
    const int r,
    const int overwrite,
    double *ans
) {

    int i, j, k, l, mr, js, ir, s;
    double y;

    mr = m * r;
    s = r;

    if (overwrite) {
        for (i = 0; i < mr * (mr + 1) / 2; i++) ans[i] = 0.0;
    }

    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j <= i; j++) {
            js = j * s;
            y = A[S(i, j, m)];
            for (k = 0; k < r; k++) {
                for (l = 0; l < (j < i ? s : k + 1); l++) {
                    ans[S(ir + k, js + l, mr)] += y * B[S(k, l, r)];
                }
            }
        }
    }
}

/* C\_KronSums\_sym */

/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_
(
    /* C real x Input */
    
        double *x,
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
    
    double *PP_sym_ans
) {

    int pN, qN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] +=  x[qN + i] * x[pN + i];
        }
    }
}

/* C\_MPinv\_sym */

void C_MPinv_sym
(
    const double *x,
    const int n,
    const double tol,
    double *dMP,
    int *rank
) {

    double *val, *vec, dtol, *rx, *work, valinv;
    int valzero = 0, info = 0, kn;

    if (n == 1) {
        if (x[0] > tol) {
            dMP[0] = 1 / x[0];
            rank[0] = 1;
        } else {
            dMP[0] = 0;
            rank[0] = 0;
        }
    } else {
        rx = Calloc(n * (n + 1) / 2, double);
        Memcpy(rx, x, n * (n + 1) / 2);
        work = Calloc(3 * n, double);
        val = Calloc(n, double);
        vec = Calloc(n * n, double);

        F77_CALL(dspev)("V", "L", &n, rx, val, vec, &n, work,
                        &info);

        dtol = val[n - 1] * tol;

        for (int k = 0; k < n; k++)
            valzero += (val[k] < dtol);
        rank[0] = n - valzero;

        for (int k = 0; k < n * (n + 1) / 2; k++) dMP[k] = 0.0;

        for (int k = valzero; k < n; k++) {
            valinv = 1 / val[k];
            kn = k * n;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= i; j++) {
                    /* MP is symmetric */
                    dMP[S(i, j, n)] += valinv * vec[kn + i] * vec[kn + j];
                }
            }
        }
        Free(rx); Free(work); Free(val); Free(vec);
    }
}


/* Memory */

/* C\_get\_P */

int C_get_P
(
/* R LECV Input */

SEXP LECV

) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[0]);
}

/* C\_get\_Q */

int C_get_Q
(
/* R LECV Input */

SEXP LECV

) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[1]);
}

/* C\_get\_varonly */

int C_get_varonly
(
/* R LECV Input */

SEXP LECV

) {

    return(INTEGER(VECTOR_ELT(LECV, varonly_SLOT))[0]);
}

/* C\_get\_Xfactor */

int C_get_Xfactor
(
/* R LECV Input */

SEXP LECV

) {

    return(INTEGER(VECTOR_ELT(LECV, Xfactor_SLOT))[0]);
}

/* C\_get\_LinearStatistic */

double* C_get_LinearStatistic
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, LinearStatistic_SLOT)));
}

/* C\_get\_Expectation */

double* C_get_Expectation
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, Expectation_SLOT)));
}

/* C\_get\_Variance */

double* C_get_Variance
(
/* R LECV Input */

SEXP LECV

) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    double *var, *covar;

    if (isNull(VECTOR_ELT(LECV, Variance_SLOT))) {
        SET_VECTOR_ELT(LECV, Variance_SLOT,
                       allocVector(REALSXP, PQ));
        if (!isNull(VECTOR_ELT(LECV, Covariance_SLOT))) {
            covar = REAL(VECTOR_ELT(LECV, Covariance_SLOT));
            var = REAL(VECTOR_ELT(LECV, Variance_SLOT));
            for (int p = 0; p < PQ; p++)
                var[p] = covar[S(p, p, PQ)];
        }
    }
    return(REAL(VECTOR_ELT(LECV, Variance_SLOT)));
}

/* C\_get\_Covariance */

double* C_get_Covariance
(
/* R LECV Input */

SEXP LECV

) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract covariance from variance only object");
    if (C_get_varonly(LECV) && PQ == 1)
        return(C_get_Variance(LECV));
    return(REAL(VECTOR_ELT(LECV, Covariance_SLOT)));
}

/* C\_get\_MPinv */

double* C_get_MPinv
(
/* R LECV Input */

SEXP LECV

) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract MPinv from variance only object");
    /* allocate memory on as needed basis */
    if (isNull(VECTOR_ELT(LECV, MPinv_SLOT))) {
        SET_VECTOR_ELT(LECV, MPinv_SLOT,
                       allocVector(REALSXP,
                                   PQ * (PQ + 1) / 2));
    }
    return(REAL(VECTOR_ELT(LECV, MPinv_SLOT)));
}

/* C\_get\_ExpectationX */

double* C_get_ExpectationX
(
/* R LECV Input */

SEXP LECV

) {
    return(REAL(VECTOR_ELT(LECV, ExpectationX_SLOT)));
}

/* C\_get\_ExpectationInfluence */

double* C_get_ExpectationInfluence
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, ExpectationInfluence_SLOT)));
}

/* C\_get\_CovarianceInfluence */

double* C_get_CovarianceInfluence
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, CovarianceInfluence_SLOT)));
}

/* C\_get\_VarianceInfluence */

double* C_get_VarianceInfluence
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, VarianceInfluence_SLOT)));
}

/* C\_get\_TableBlock */

double* C_get_TableBlock
(
/* R LECV Input */

SEXP LECV

) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) == R_NilValue)
        error("object does not contain table block slot");
    return(REAL(VECTOR_ELT(LECV, TableBlock_SLOT)));
}

/* C\_get\_Sumweights */

double* C_get_Sumweights
(
/* R LECV Input */

SEXP LECV

) {
    if (VECTOR_ELT(LECV, Sumweights_SLOT) == R_NilValue)
        error("object does not contain sumweights slot");
    return(REAL(VECTOR_ELT(LECV, Sumweights_SLOT)));
}

/* C\_get\_Table */

double* C_get_Table
(
/* R LECV Input */

SEXP LECV

) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(REAL(VECTOR_ELT(LECV, Table_SLOT)));
}

/* C\_get\_dimTable */

int* C_get_dimTable
(
/* R LECV Input */

SEXP LECV

) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(INTEGER(getAttrib(VECTOR_ELT(LECV, Table_SLOT),
                             R_DimSymbol)));
}

/* C\_get\_B */

int C_get_B
(
/* R LECV Input */

SEXP LECV

) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) != R_NilValue)
        return(LENGTH(VECTOR_ELT(LECV, Sumweights_SLOT)));
    return(C_get_dimTable(LECV)[2]);
}

/* C\_get\_nperm */

R_xlen_t C_get_nperm
(
/* R LECV Input */

SEXP LECV

) {
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    return(XLENGTH(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)) / PQ);
}

/* C\_get\_PermutedLinearStatistic */

double* C_get_PermutedLinearStatistic
(
/* R LECV Input */

SEXP LECV

) {

    return(REAL(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)));
}

/* C\_get\_tol */

double C_get_tol
(
/* R LECV Input */

SEXP LECV

) {
    return(REAL(VECTOR_ELT(LECV, tol_SLOT))[0]);
}

/* RC\_init\_LECV\_1d */

SEXP RC_init_LECV_1d
(
    /* C integer P Input */
    
        int P
    ,
    /* C integer Q Input */
    
        int Q
    ,
    int varonly,
    /* C integer B Input */
    
        int B
    ,
    int Xfactor,
    double tol
) {

    SEXP ans;

    /* R\_init\_LECV */
    
        SEXP vo, d, names, tolerance;
        int PQ; 

        /* Memory Input Checks */
        
        if (P <= 0)
            error("P is not positive");

        if (Q <= 0)
            error("Q is not positive");

        if (B <= 0)
            error("B is not positive");

        if (varonly < 0 || varonly > 1)
            error("varonly is not 0 or 1");

        if (Xfactor < 0 || Xfactor > 1)
            error("Xfactor is not 0 or 1");

        if (tol <= DBL_MIN)
            error("tol is not positive");
        

        PQ = P * Q;

        /* Memory Names */
        
        PROTECT(names = allocVector(STRSXP, Table_SLOT + 1));
        SET_STRING_ELT(names, LinearStatistic_SLOT, mkChar("LinearStatistic"));
        SET_STRING_ELT(names, Expectation_SLOT, mkChar("Expectation"));
        SET_STRING_ELT(names, varonly_SLOT, mkChar("varonly"));
        SET_STRING_ELT(names, Variance_SLOT, mkChar("Variance"));
        SET_STRING_ELT(names, Covariance_SLOT, mkChar("Covariance"));
        SET_STRING_ELT(names, MPinv_SLOT, mkChar("MPinv"));
        SET_STRING_ELT(names, ExpectationX_SLOT, mkChar("ExpectationX"));
        SET_STRING_ELT(names, dim_SLOT, mkChar("dimension"));
        SET_STRING_ELT(names, ExpectationInfluence_SLOT,
                       mkChar("ExpectationInfluence"));
        SET_STRING_ELT(names, Xfactor_SLOT, mkChar("Xfactor"));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT,
                       mkChar("CovarianceInfluence"));
        SET_STRING_ELT(names, VarianceInfluence_SLOT,
                       mkChar("VarianceInfluence"));
        SET_STRING_ELT(names, TableBlock_SLOT, mkChar("TableBlock"));
        SET_STRING_ELT(names, Sumweights_SLOT, mkChar("Sumweights"));
        SET_STRING_ELT(names, PermutedLinearStatistic_SLOT,
                       mkChar("PermutedLinearStatistic"));
        SET_STRING_ELT(names, StandardisedPermutedLinearStatistic_SLOT,
                       mkChar("StandardisedPermutedLinearStatistic"));
        SET_STRING_ELT(names, tol_SLOT, mkChar("tol"));
        SET_STRING_ELT(names, Table_SLOT, mkChar("Table"));
        

        /* Table_SLOT is always last and only used in 2d case
           ie omitted here */
        PROTECT(ans = allocVector(VECSXP, Table_SLOT + 1));
        SET_VECTOR_ELT(ans, LinearStatistic_SLOT, allocVector(REALSXP, PQ));
        SET_VECTOR_ELT(ans, Expectation_SLOT, allocVector(REALSXP, PQ));
        SET_VECTOR_ELT(ans, varonly_SLOT, vo = allocVector(INTSXP, 1));
        INTEGER(vo)[0] = varonly;
        if (varonly) {
            SET_VECTOR_ELT(ans, Variance_SLOT, allocVector(REALSXP, PQ));
        } else  {
            SET_VECTOR_ELT(ans, Covariance_SLOT,
                           allocVector(REALSXP, PQ * (PQ + 1) / 2));
        }
        SET_VECTOR_ELT(ans, ExpectationX_SLOT, allocVector(REALSXP, P));
        SET_VECTOR_ELT(ans, dim_SLOT, d = allocVector(INTSXP, 2));
        INTEGER(d)[0] = P;
        INTEGER(d)[1] = Q;
        SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                       allocVector(REALSXP, B * Q));

        /* should always _both_ be there */
        SET_VECTOR_ELT(ans, VarianceInfluence_SLOT,
                       allocVector(REALSXP, B * Q));
        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, B * Q * (Q + 1) / 2));

        SET_VECTOR_ELT(ans, Xfactor_SLOT, allocVector(INTSXP, 1));
        INTEGER(VECTOR_ELT(ans, Xfactor_SLOT))[0] = Xfactor;
        SET_VECTOR_ELT(ans, TableBlock_SLOT, allocVector(REALSXP, B + 1));
        SET_VECTOR_ELT(ans, Sumweights_SLOT, allocVector(REALSXP, B));
        SET_VECTOR_ELT(ans, PermutedLinearStatistic_SLOT,
                       allocMatrix(REALSXP, 0, 0));
        SET_VECTOR_ELT(ans, StandardisedPermutedLinearStatistic_SLOT,
                       allocMatrix(REALSXP, 0, 0));
        SET_VECTOR_ELT(ans, tol_SLOT, tolerance = allocVector(REALSXP, 1));
        REAL(tolerance)[0] = tol;
        namesgets(ans, names);

        /* set inital zeros */
        for (int p = 0; p < PQ; p++) {
            C_get_LinearStatistic(ans)[p] = 0.0;
            C_get_Expectation(ans)[p] = 0.0;
            if (varonly)
                C_get_Variance(ans)[p] = 0.0;
        }
        if (!varonly) {
            for (int p = 0; p < PQ * (PQ + 1) / 2; p++)
                C_get_Covariance(ans)[p] = 0.0;
        }
        for (int q = 0; q < Q; q++) {
            C_get_ExpectationInfluence(ans)[q] = 0.0;
            C_get_VarianceInfluence(ans)[q] = 0.0;
        }
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            C_get_CovarianceInfluence(ans)[q] = 0.0;
    

    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(REALSXP, B + 1));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(REALSXP, B));

    UNPROTECT(2);
    return(ans);
}

/* RC\_init\_LECV\_2d */

SEXP RC_init_LECV_2d
(
    /* C integer P Input */
    
        int P
    ,
    /* C integer Q Input */
    
        int Q
    ,
    int varonly,
    int Lx,
    int Ly,
    /* C integer B Input */
    
        int B
    ,
    int Xfactor,
    double tol
) {
    SEXP ans, tabdim, tab;

    if (Lx <= 0)
        error("Lx is not positive");

    if (Ly <= 0)
        error("Ly is not positive");

    /* R\_init\_LECV */
    
        SEXP vo, d, names, tolerance;
        int PQ; 

        /* Memory Input Checks */
        
        if (P <= 0)
            error("P is not positive");

        if (Q <= 0)
            error("Q is not positive");

        if (B <= 0)
            error("B is not positive");

        if (varonly < 0 || varonly > 1)
            error("varonly is not 0 or 1");

        if (Xfactor < 0 || Xfactor > 1)
            error("Xfactor is not 0 or 1");

        if (tol <= DBL_MIN)
            error("tol is not positive");
        

        PQ = P * Q;

        /* Memory Names */
        
        PROTECT(names = allocVector(STRSXP, Table_SLOT + 1));
        SET_STRING_ELT(names, LinearStatistic_SLOT, mkChar("LinearStatistic"));
        SET_STRING_ELT(names, Expectation_SLOT, mkChar("Expectation"));
        SET_STRING_ELT(names, varonly_SLOT, mkChar("varonly"));
        SET_STRING_ELT(names, Variance_SLOT, mkChar("Variance"));
        SET_STRING_ELT(names, Covariance_SLOT, mkChar("Covariance"));
        SET_STRING_ELT(names, MPinv_SLOT, mkChar("MPinv"));
        SET_STRING_ELT(names, ExpectationX_SLOT, mkChar("ExpectationX"));
        SET_STRING_ELT(names, dim_SLOT, mkChar("dimension"));
        SET_STRING_ELT(names, ExpectationInfluence_SLOT,
                       mkChar("ExpectationInfluence"));
        SET_STRING_ELT(names, Xfactor_SLOT, mkChar("Xfactor"));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT,
                       mkChar("CovarianceInfluence"));
        SET_STRING_ELT(names, VarianceInfluence_SLOT,
                       mkChar("VarianceInfluence"));
        SET_STRING_ELT(names, TableBlock_SLOT, mkChar("TableBlock"));
        SET_STRING_ELT(names, Sumweights_SLOT, mkChar("Sumweights"));
        SET_STRING_ELT(names, PermutedLinearStatistic_SLOT,
                       mkChar("PermutedLinearStatistic"));
        SET_STRING_ELT(names, StandardisedPermutedLinearStatistic_SLOT,
                       mkChar("StandardisedPermutedLinearStatistic"));
        SET_STRING_ELT(names, tol_SLOT, mkChar("tol"));
        SET_STRING_ELT(names, Table_SLOT, mkChar("Table"));
        

        /* Table_SLOT is always last and only used in 2d case
           ie omitted here */
        PROTECT(ans = allocVector(VECSXP, Table_SLOT + 1));
        SET_VECTOR_ELT(ans, LinearStatistic_SLOT, allocVector(REALSXP, PQ));
        SET_VECTOR_ELT(ans, Expectation_SLOT, allocVector(REALSXP, PQ));
        SET_VECTOR_ELT(ans, varonly_SLOT, vo = allocVector(INTSXP, 1));
        INTEGER(vo)[0] = varonly;
        if (varonly) {
            SET_VECTOR_ELT(ans, Variance_SLOT, allocVector(REALSXP, PQ));
        } else  {
            SET_VECTOR_ELT(ans, Covariance_SLOT,
                           allocVector(REALSXP, PQ * (PQ + 1) / 2));
        }
        SET_VECTOR_ELT(ans, ExpectationX_SLOT, allocVector(REALSXP, P));
        SET_VECTOR_ELT(ans, dim_SLOT, d = allocVector(INTSXP, 2));
        INTEGER(d)[0] = P;
        INTEGER(d)[1] = Q;
        SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                       allocVector(REALSXP, B * Q));

        /* should always _both_ be there */
        SET_VECTOR_ELT(ans, VarianceInfluence_SLOT,
                       allocVector(REALSXP, B * Q));
        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, B * Q * (Q + 1) / 2));

        SET_VECTOR_ELT(ans, Xfactor_SLOT, allocVector(INTSXP, 1));
        INTEGER(VECTOR_ELT(ans, Xfactor_SLOT))[0] = Xfactor;
        SET_VECTOR_ELT(ans, TableBlock_SLOT, allocVector(REALSXP, B + 1));
        SET_VECTOR_ELT(ans, Sumweights_SLOT, allocVector(REALSXP, B));
        SET_VECTOR_ELT(ans, PermutedLinearStatistic_SLOT,
                       allocMatrix(REALSXP, 0, 0));
        SET_VECTOR_ELT(ans, StandardisedPermutedLinearStatistic_SLOT,
                       allocMatrix(REALSXP, 0, 0));
        SET_VECTOR_ELT(ans, tol_SLOT, tolerance = allocVector(REALSXP, 1));
        REAL(tolerance)[0] = tol;
        namesgets(ans, names);

        /* set inital zeros */
        for (int p = 0; p < PQ; p++) {
            C_get_LinearStatistic(ans)[p] = 0.0;
            C_get_Expectation(ans)[p] = 0.0;
            if (varonly)
                C_get_Variance(ans)[p] = 0.0;
        }
        if (!varonly) {
            for (int p = 0; p < PQ * (PQ + 1) / 2; p++)
                C_get_Covariance(ans)[p] = 0.0;
        }
        for (int q = 0; q < Q; q++) {
            C_get_ExpectationInfluence(ans)[q] = 0.0;
            C_get_VarianceInfluence(ans)[q] = 0.0;
        }
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            C_get_CovarianceInfluence(ans)[q] = 0.0;
    

    PROTECT(tabdim = allocVector(INTSXP, 3));
    INTEGER(tabdim)[0] = Lx + 1;
    INTEGER(tabdim)[1] = Ly + 1;
    INTEGER(tabdim)[2] = B;
    SET_VECTOR_ELT(ans, Table_SLOT,
                   tab = allocVector(REALSXP,
                       INTEGER(tabdim)[0] *
                       INTEGER(tabdim)[1] *
                       INTEGER(tabdim)[2]));
    dimgets(tab, tabdim);

    UNPROTECT(3);
    return(ans);
}


/* P-Values */

/* C\_chisq\_pvalue */

/* lower = 1 means p-value, lower = 0 means 1 - p-value */
double C_chisq_pvalue
(
    const double stat,
    const int df,
    const int lower,
    const int give_log
) {
    return(pchisq(stat, (double) df, lower, give_log));
}

/* C\_perm\_pvalue */

double C_perm_pvalue
(
    const int greater,
    const double nperm,
    const int lower,
    const int give_log
) {

    double ret;

    if (give_log) {
         if (lower) {
             ret = log1p(- (double) greater / nperm);
         } else {
             ret = log(greater) - log(nperm);
         }
    } else {
        if (lower) {
            ret = 1.0 - (double) greater / nperm;
        } else {
            ret = (double) greater / nperm;
        }
    }
    return(ret);
}

/* C\_norm\_pvalue */

double C_norm_pvalue
(
    const double stat,
    const int alternative,
    const int lower,
    const int give_log
) {

    double ret;

    if (alternative == ALTERNATIVE_less) {
        return(pnorm(stat, 0.0, 1.0, 1 - lower, give_log));
    } else if (alternative == ALTERNATIVE_greater) {
        return(pnorm(stat, 0.0, 1.0, lower, give_log));
    } else if (alternative == ALTERNATIVE_twosided) {
        if (lower) {
            ret = pnorm(fabs(stat)*-1.0, 0.0, 1.0, 1, 0);
            if (give_log) {
                return(log1p(- 2 * ret));
            } else {
                return(1 - 2 * ret);
            }
        } else {
            ret = pnorm(fabs(stat)*-1.0, 0.0, 1.0, 1, give_log);
            if (give_log) {
                return(ret + log(2));
            } else {
                return(2 * ret);
            }
        }
    }
    return(NA_REAL);
}

/* C\_maxtype\_pvalue */

double C_maxtype_pvalue
(
    const double stat,
    const double *Covariance,
    const int n,
    const int alternative,
    const int lower,
    const int give_log,
    int maxpts, /* const? */
    double releps,
    double abseps,
    double tol
) {

    int nu = 0, inform, i, j, sub, nonzero, *infin, *index, rnd = 0;
    double ans, myerror, *lowerbnd, *upperbnd, *delta, *corr, *sd;

    /* univariate problem */
    if (n == 1)
        return(C_norm_pvalue(stat, alternative, lower, give_log));

    if (n == 2)
         corr = Calloc(1, double);
    else
         corr = Calloc(n + ((n - 2) * (n - 1))/2, double);

    sd = Calloc(n, double);
    lowerbnd = Calloc(n, double);
    upperbnd = Calloc(n, double);
    infin = Calloc(n, int);
    delta = Calloc(n, double);
    index = Calloc(n, int);

    /* determine elements with non-zero variance */

    nonzero = 0;
    for (i = 0; i < n; i++) {
        if (Covariance[S(i, i, n)] > tol) {
            index[nonzero] = i;
            nonzero++;
        }
    }

    /* mvtdst assumes the unique elements of the triangular
       covariance matrix to be passes as argument CORREL
    */

    for (int nz = 0; nz < nonzero; nz++) {

        /* handle elements with non-zero variance only */
        i = index[nz];

        /* standard deviations */
        sd[i] = sqrt(Covariance[S(i, i, n)]);

        if (alternative == ALTERNATIVE_less) {
            lowerbnd[nz] = stat;
            upperbnd[nz] = R_PosInf;
            infin[nz] = 1;
        } else if (alternative == ALTERNATIVE_greater) {
            lowerbnd[nz] = R_NegInf;
            upperbnd[nz] = stat;
            infin[nz] = 0;
        } else if (alternative == ALTERNATIVE_twosided) {
            lowerbnd[nz] = fabs(stat) * -1.0;
            upperbnd[nz] = fabs(stat);
            infin[nz] = 2;
        }

        delta[nz] = 0.0;

        /* set up vector of correlations, i.e., the upper
           triangular part of the covariance matrix) */
        for (int jz = 0; jz < nz; jz++) {
            j = index[jz];
            sub = (int) (jz + 1) + (double) ((nz - 1) * nz) / 2 - 1;
            if (sd[i] == 0.0 || sd[j] == 0.0)
                corr[sub] = 0.0;
            else
                corr[sub] = Covariance[S(i, j, n)] / (sd[i] * sd[j]);
        }
    }

    /* call mvtnorm's mvtdst C function defined in mvtnorm/include/mvtnormAPI.h */
    mvtnorm_C_mvtdst(&nonzero, &nu, lowerbnd, upperbnd, infin, corr, delta,
                     &maxpts, &abseps, &releps, &myerror, &ans, &inform, &rnd);

    /* inform == 0 means: everything is OK */
    switch (inform) {
        case 0: break;
        case 1: warning("cmvnorm: completion with ERROR > EPS"); break;
        case 2: warning("cmvnorm: N > 1000 or N < 1");
                ans = 0.0;
                break;
        case 3: warning("cmvnorm: correlation matrix not positive semi-definite");
                ans = 0.0;
                break;
        default: warning("cmvnorm: unknown problem in MVTDST");
                ans = 0.0;
    }
    Free(corr); Free(sd); Free(lowerbnd); Free(upperbnd);
    Free(infin); Free(delta);

    /* ans = 1 - p-value */
    if (lower) {
        if (give_log)
            return(log(ans)); /* log(1 - p-value) */
        return(ans); /* 1 - p-value */
    } else {
        if (give_log)
            return(log1p(ans)); /* log(p-value) */
        return(1 - ans); /* p-value */
    }
}


/* KronSums */

/* C\_KronSums\_dweights\_dsubset */

void C_KronSums_dweights_dsubset
(
    /* C KronSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        const int SYMMETRIC,
        double *centerx,
        double *centery,
        const int CENTER,
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
)
{
    double *s, *w; 
    /* KronSums Body */
    
        double *xx, *yy, cx = 0.0, cy = 0.0, *thisPQ_ans;
        int idx;

        for (int p = 0; p < P; p++) {
            for (int q = (SYMMETRIC ? p : 0); q < Q; q++) {
                /* SYMMETRIC is column-wise, default
                   is row-wise (maybe need to change this) */
                if (SYMMETRIC) {
                    idx = S(p, q, P); 
                } else {
                    idx = q * P + p;  
                }
                PQ_ans[idx] = 0.0;
                thisPQ_ans = PQ_ans + idx;
                yy = y + N * q;
                xx = x + N * p;

                if (CENTER) {
                    cx = centerx[p];
                    cy = centery[q];
                }
                /* init subset loop */
                
                    R_xlen_t diff = 0;
                    s = subset + offset;
                    w = weights;
                    /* subset is R-style index in 1:N */
                    if (Nsubset > 0)
                        diff = (R_xlen_t) s[0] - 1;
                
                /* start subset loop */
                
                    for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
                
                {
                    xx = xx + diff;
                    yy = yy + diff;
                    if (HAS_WEIGHTS) {
                        w = w + diff;
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0] * w[0];
                        }
                    } else {
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0];
                        }
                    }
                    /* continue subset loop */
                    
                        if (Nsubset > 0) {
                            /* NB: diff also works with R style index */
                            diff = (R_xlen_t) s[1] - s[0];
                            if (diff < 0)
                                error("subset not sorted");
                            s++;
                        } else {
                            diff = 1;
                        }
                    
                }
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
            }
        }
    
}

/* C\_KronSums\_iweights\_dsubset */

void C_KronSums_iweights_dsubset
(
    /* C KronSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        const int SYMMETRIC,
        double *centerx,
        double *centery,
        const int CENTER,
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    double *s;
    int *w; 
    /* KronSums Body */
    
        double *xx, *yy, cx = 0.0, cy = 0.0, *thisPQ_ans;
        int idx;

        for (int p = 0; p < P; p++) {
            for (int q = (SYMMETRIC ? p : 0); q < Q; q++) {
                /* SYMMETRIC is column-wise, default
                   is row-wise (maybe need to change this) */
                if (SYMMETRIC) {
                    idx = S(p, q, P); 
                } else {
                    idx = q * P + p;  
                }
                PQ_ans[idx] = 0.0;
                thisPQ_ans = PQ_ans + idx;
                yy = y + N * q;
                xx = x + N * p;

                if (CENTER) {
                    cx = centerx[p];
                    cy = centery[q];
                }
                /* init subset loop */
                
                    R_xlen_t diff = 0;
                    s = subset + offset;
                    w = weights;
                    /* subset is R-style index in 1:N */
                    if (Nsubset > 0)
                        diff = (R_xlen_t) s[0] - 1;
                
                /* start subset loop */
                
                    for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
                
                {
                    xx = xx + diff;
                    yy = yy + diff;
                    if (HAS_WEIGHTS) {
                        w = w + diff;
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0] * w[0];
                        }
                    } else {
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0];
                        }
                    }
                    /* continue subset loop */
                    
                        if (Nsubset > 0) {
                            /* NB: diff also works with R style index */
                            diff = (R_xlen_t) s[1] - s[0];
                            if (diff < 0)
                                error("subset not sorted");
                            s++;
                        } else {
                            diff = 1;
                        }
                    
                }
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
            }
        }
    
}

/* C\_KronSums\_iweights\_isubset */

void C_KronSums_iweights_isubset
(
    /* C KronSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        const int SYMMETRIC,
        double *centerx,
        double *centery,
        const int CENTER,
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    int *s, *w;
    /* KronSums Body */
    
        double *xx, *yy, cx = 0.0, cy = 0.0, *thisPQ_ans;
        int idx;

        for (int p = 0; p < P; p++) {
            for (int q = (SYMMETRIC ? p : 0); q < Q; q++) {
                /* SYMMETRIC is column-wise, default
                   is row-wise (maybe need to change this) */
                if (SYMMETRIC) {
                    idx = S(p, q, P); 
                } else {
                    idx = q * P + p;  
                }
                PQ_ans[idx] = 0.0;
                thisPQ_ans = PQ_ans + idx;
                yy = y + N * q;
                xx = x + N * p;

                if (CENTER) {
                    cx = centerx[p];
                    cy = centery[q];
                }
                /* init subset loop */
                
                    R_xlen_t diff = 0;
                    s = subset + offset;
                    w = weights;
                    /* subset is R-style index in 1:N */
                    if (Nsubset > 0)
                        diff = (R_xlen_t) s[0] - 1;
                
                /* start subset loop */
                
                    for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
                
                {
                    xx = xx + diff;
                    yy = yy + diff;
                    if (HAS_WEIGHTS) {
                        w = w + diff;
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0] * w[0];
                        }
                    } else {
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0];
                        }
                    }
                    /* continue subset loop */
                    
                        if (Nsubset > 0) {
                            /* NB: diff also works with R style index */
                            diff = (R_xlen_t) s[1] - s[0];
                            if (diff < 0)
                                error("subset not sorted");
                            s++;
                        } else {
                            diff = 1;
                        }
                    
                }
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
            }
        }
    
}

/* C\_KronSums\_dweights\_isubset */

void C_KronSums_dweights_isubset
(
    /* C KronSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        const int SYMMETRIC,
        double *centerx,
        double *centery,
        const int CENTER,
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) {
    int *s; 
    double *w;
    /* KronSums Body */
    
        double *xx, *yy, cx = 0.0, cy = 0.0, *thisPQ_ans;
        int idx;

        for (int p = 0; p < P; p++) {
            for (int q = (SYMMETRIC ? p : 0); q < Q; q++) {
                /* SYMMETRIC is column-wise, default
                   is row-wise (maybe need to change this) */
                if (SYMMETRIC) {
                    idx = S(p, q, P); 
                } else {
                    idx = q * P + p;  
                }
                PQ_ans[idx] = 0.0;
                thisPQ_ans = PQ_ans + idx;
                yy = y + N * q;
                xx = x + N * p;

                if (CENTER) {
                    cx = centerx[p];
                    cy = centery[q];
                }
                /* init subset loop */
                
                    R_xlen_t diff = 0;
                    s = subset + offset;
                    w = weights;
                    /* subset is R-style index in 1:N */
                    if (Nsubset > 0)
                        diff = (R_xlen_t) s[0] - 1;
                
                /* start subset loop */
                
                    for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
                
                {
                    xx = xx + diff;
                    yy = yy + diff;
                    if (HAS_WEIGHTS) {
                        w = w + diff;
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0] * w[0];
                        }
                    } else {
                        if (CENTER) {
                            thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                        } else {
                            thisPQ_ans[0] += xx[0] * yy[0];
                        }
                    }
                    /* continue subset loop */
                    
                        if (Nsubset > 0) {
                            /* NB: diff also works with R style index */
                            diff = (R_xlen_t) s[1] - s[0];
                            if (diff < 0)
                                error("subset not sorted");
                            s++;
                        } else {
                            diff = 1;
                        }
                    
                }
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    thisPQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
            }
        }
    
}

/* C\_XfactorKronSums\_dweights\_dsubset */

void C_XfactorKronSums_dweights_dsubset
(
    /* C XfactorKronSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
)
{
    double *s, *w; 
    /* XfactorKronSums Body */
    
        int *xx, ixi;
        double *yy;
     
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            yy = y + N * q;
            xx = x;
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                yy = yy + diff;
                ixi = xx[0] - 1;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0] * w[0];
                } else {
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0];
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            yy = yy + diff;
            ixi = xx[0] - 1;
            if (HAS_WEIGHTS) {
                w = w + diff;
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0] * w[0];
            } else {
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0];
            }
        }
    
}

/* C\_XfactorKronSums\_iweights\_dsubset */

void C_XfactorKronSums_iweights_dsubset
(
    /* C XfactorKronSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    double *s;
    int *w; 
    /* XfactorKronSums Body */
    
        int *xx, ixi;
        double *yy;
     
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            yy = y + N * q;
            xx = x;
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                yy = yy + diff;
                ixi = xx[0] - 1;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0] * w[0];
                } else {
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0];
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            yy = yy + diff;
            ixi = xx[0] - 1;
            if (HAS_WEIGHTS) {
                w = w + diff;
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0] * w[0];
            } else {
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0];
            }
        }
    
}

/* C\_XfactorKronSums\_iweights\_isubset */

void C_XfactorKronSums_iweights_isubset
(
    /* C XfactorKronSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    int *s, *w;
    /* XfactorKronSums Body */
    
        int *xx, ixi;
        double *yy;
     
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            yy = y + N * q;
            xx = x;
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                yy = yy + diff;
                ixi = xx[0] - 1;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0] * w[0];
                } else {
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0];
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            yy = yy + diff;
            ixi = xx[0] - 1;
            if (HAS_WEIGHTS) {
                w = w + diff;
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0] * w[0];
            } else {
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0];
            }
        }
    
}

/* C\_XfactorKronSums\_dweights\_isubset */

void C_XfactorKronSums_dweights_isubset
(
    /* C XfactorKronSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) {
    int *s; 
    double *w;
    /* XfactorKronSums Body */
    
        int *xx, ixi;
        double *yy;
     
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            yy = y + N * q;
            xx = x;
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                yy = yy + diff;
                ixi = xx[0] - 1;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0] * w[0];
                } else {
                    if (ixi >= 0)
                        PQ_ans[ixi + q * P] += yy[0];
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            yy = yy + diff;
            ixi = xx[0] - 1;
            if (HAS_WEIGHTS) {
                w = w + diff;
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0] * w[0];
            } else {
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0];
            }
        }
    
}

/* RC\_KronSums */

/* RC\_KronSums Prototype */

void RC_KronSums
(
    /* RC KronSums Input */
    
        /* R x Input */
        
            SEXP x,
        
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
        /* C real y Input */
        
            double *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        const int SYMMETRIC,
        double *centerx,
        double *centery,
        const int CENTER,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 

{
    if (TYPEOF(x) == INTSXP) {
        if (SYMMETRIC) error("not implemented");
        if (CENTER) error("not implemented");
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_XfactorKronSums_iweights_isubset(INTEGER(x), N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_XfactorKronSums_iweights_dsubset(INTEGER(x), N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_XfactorKronSums_dweights_isubset(INTEGER(x), N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_XfactorKronSums_dweights_dsubset(INTEGER(x), N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
    } else {
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_iweights_isubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_iweights_dsubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_dweights_isubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_dweights_dsubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
    }
}

/* R\_KronSums */

/* R\_KronSums Prototype */

SEXP R_KronSums
(
    /* R x Input */
    
        SEXP x,
    
    SEXP P,
    /* R y Input */
    
        SEXP y,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    SEXP symmetric
) 

{
    SEXP ans;
    /* C integer Q Input */
    
        int Q
    ;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;

    double center;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    if (INTEGER(symmetric)[0]) {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * (INTEGER(P)[0] + 1) / 2));
    } else {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    }
    RC_KronSums(x, N, INTEGER(P)[0], REAL(y), Q, INTEGER(symmetric)[0], &center, &center, 
                !DoCenter, weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

/* C\_KronSums\_Permutation\_isubset */

void C_KronSums_Permutation_isubset
(
    /* C real x Input */
    
        double *x,
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
    
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    int *subsety,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    /* KronSums Permutation Body */
    
        R_xlen_t qP, qN, pN, qPp;

        for (int q = 0; q < Q; q++) {
            qN = q * N;
            qP = q * P;
            for (int p = 0; p < P; p++) {
                qPp = qP + p;
                PQ_ans[qPp] = 0.0;
                pN = p * N;
                for (R_xlen_t i = offset; i < Nsubset; i++)
                    PQ_ans[qPp] += y[qN + (R_xlen_t) subsety[i] - 1] * x[pN + (R_xlen_t) subset[i] - 1];
            }
        }
    
}

/* C\_KronSums\_Permutation\_dsubset */

void C_KronSums_Permutation_dsubset
(
    /* C real x Input */
    
        double *x,
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
    
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    double *subsety,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    /* KronSums Permutation Body */
    
        R_xlen_t qP, qN, pN, qPp;

        for (int q = 0; q < Q; q++) {
            qN = q * N;
            qP = q * P;
            for (int p = 0; p < P; p++) {
                qPp = qP + p;
                PQ_ans[qPp] = 0.0;
                pN = p * N;
                for (R_xlen_t i = offset; i < Nsubset; i++)
                    PQ_ans[qPp] += y[qN + (R_xlen_t) subsety[i] - 1] * x[pN + (R_xlen_t) subset[i] - 1];
            }
        }
    
}

/* C\_XfactorKronSums\_Permutation\_isubset */

void C_XfactorKronSums_Permutation_isubset
(
    /* C integer x Input */
    
        int *x,
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
    
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    int *subsety,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    /* XfactorKronSums Permutation Body */
    
        R_xlen_t qP, qN;

        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            qP = q * P;
            qN = q * N;
            for (R_xlen_t i = offset; i < Nsubset; i++)
                PQ_ans[x[(R_xlen_t) subset[i] - 1] - 1 + qP] += y[qN + (R_xlen_t) subsety[i] - 1];
        }
    
}

/* C\_XfactorKronSums\_Permutation\_dsubset */

void C_XfactorKronSums_Permutation_dsubset
(
    /* C integer x Input */
    
        int *x,
        /* C integer N Input */
        
            R_xlen_t N
        ,
        /* C integer P Input */
        
            int P
        ,
    
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    double *subsety,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 
{
    /* XfactorKronSums Permutation Body */
    
        R_xlen_t qP, qN;

        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

        for (int q = 0; q < Q; q++) {
            qP = q * P;
            qN = q * N;
            for (R_xlen_t i = offset; i < Nsubset; i++)
                PQ_ans[x[(R_xlen_t) subset[i] - 1] - 1 + qP] += y[qN + (R_xlen_t) subsety[i] - 1];
        }
    
}

/* RC\_KronSums\_Permutation */

/* RC\_KronSums\_Permutation Prototype */

void RC_KronSums_Permutation
(
    /* R x Input */
    
        SEXP x,
    
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer P Input */
    
        int P
    ,
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    SEXP subsety,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 

{
    if (TYPEOF(x) == INTSXP) {
    if (TYPEOF(subset) == INTSXP) {
        C_XfactorKronSums_Permutation_isubset(INTEGER(x), N, P, y, Q, 
                                       INTEGER(subset), offset, Nsubset, 
                                       INTEGER(subsety), PQ_ans);
        } else {
        C_XfactorKronSums_Permutation_dsubset(INTEGER(x), N, P, y, Q, 
                                       REAL(subset), offset, Nsubset, 
                                       REAL(subsety), PQ_ans);
    }
    } else {
    if (TYPEOF(subset) == INTSXP) {
        C_KronSums_Permutation_isubset(REAL(x), N, P, y, Q, 
                                       INTEGER(subset), offset, Nsubset, 
                                       INTEGER(subsety), PQ_ans);
        } else {
        C_KronSums_Permutation_dsubset(REAL(x), N, P, y, Q, 
                                       REAL(subset), offset, Nsubset, 
                                       REAL(subsety), PQ_ans);
    }
    }
}

/* R\_KronSums\_Permutation */

/* R\_KronSums\_Permutation Prototype */

SEXP R_KronSums_Permutation
(
    /* R x Input */
    
        SEXP x,
    
    SEXP P,
    /* R y Input */
    
        SEXP y,
    
    /* R subset Input */
    
        SEXP subset
    ,
    SEXP subsety
)

{
    SEXP ans;
    /* C integer Q Input */
    
        int Q
    ;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_KronSums_Permutation(x, N, INTEGER(P)[0], REAL(y), Q, subset, Offset0, Nsubset, 
                            subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/* colSums */

/* C\_colSums\_dweights\_dsubset */

void C_colSums_dweights_dsubset
(
    /* C colSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        const int power,
        double *centerx,
        const int CENTER,
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
        
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C colSums Answer */
    
        double *P_ans
    
) 
{
    double *s, *w; 
    /* colSums Body */
    
        double *xx, cx = 0.0;

        for (int p = 0; p < P; p++) {
            P_ans[0] = 0.0;
            xx = x + N * p;
            if (CENTER) {
                cx = centerx[p];
            }
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    P_ans[0] += pow(xx[0] - cx, power) * w[0];
                } else {
                    P_ans[0] += pow(xx[0] - cx, power);
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            P_ans++;
        }
    
}

/* C\_colSums\_iweights\_dsubset */

void C_colSums_iweights_dsubset
(
    /* C colSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        const int power,
        double *centerx,
        const int CENTER,
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C colSums Answer */
    
        double *P_ans
    
) 
{
    double *s;
    int *w; 
    /* colSums Body */
    
        double *xx, cx = 0.0;

        for (int p = 0; p < P; p++) {
            P_ans[0] = 0.0;
            xx = x + N * p;
            if (CENTER) {
                cx = centerx[p];
            }
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    P_ans[0] += pow(xx[0] - cx, power) * w[0];
                } else {
                    P_ans[0] += pow(xx[0] - cx, power);
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            P_ans++;
        }
    
}

/* C\_colSums\_iweights\_isubset */

void C_colSums_iweights_isubset
(
    /* C colSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        const int power,
        double *centerx,
        const int CENTER,
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C colSums Answer */
    
        double *P_ans
    
) 
{
    int *s, *w;
    /* colSums Body */
    
        double *xx, cx = 0.0;

        for (int p = 0; p < P; p++) {
            P_ans[0] = 0.0;
            xx = x + N * p;
            if (CENTER) {
                cx = centerx[p];
            }
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    P_ans[0] += pow(xx[0] - cx, power) * w[0];
                } else {
                    P_ans[0] += pow(xx[0] - cx, power);
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            P_ans++;
        }
    
}

/* C\_colSums\_dweights\_isubset */

void C_colSums_dweights_isubset
(
    /* C colSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        const int power,
        double *centerx,
        const int CENTER,
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C colSums Answer */
    
        double *P_ans
    
) 
{
    int *s; 
    double *w;
    /* colSums Body */
    
        double *xx, cx = 0.0;

        for (int p = 0; p < P; p++) {
            P_ans[0] = 0.0;
            xx = x + N * p;
            if (CENTER) {
                cx = centerx[p];
            }
            /* init subset loop */
            
                R_xlen_t diff = 0;
                s = subset + offset;
                w = weights;
                /* subset is R-style index in 1:N */
                if (Nsubset > 0)
                    diff = (R_xlen_t) s[0] - 1;
            
            /* start subset loop */
            
                for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
            {
                xx = xx + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    P_ans[0] += pow(xx[0] - cx, power) * w[0];
                } else {
                    P_ans[0] += pow(xx[0] - cx, power);
                }
                /* continue subset loop */
                
                    if (Nsubset > 0) {
                        /* NB: diff also works with R style index */
                        diff = (R_xlen_t) s[1] - s[0];
                        if (diff < 0)
                            error("subset not sorted");
                        s++;
                    } else {
                        diff = 1;
                    }
                
            }
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            P_ans++;
        }
    
}

/* RC\_colSums */

/* RC\_colSums Prototype */

void RC_colSums
(
    /* C colSums Input */
    
        /* C real x Input */
        
            double *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        const int power,
        double *centerx,
        const int CENTER,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C colSums Answer */
    
        double *P_ans
    
) 

{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_colSums_iweights_isubset(x, N, P, power, centerx, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_colSums_iweights_dsubset(x, N, P, power, centerx, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_colSums_dweights_isubset(x, N, P, power, centerx, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_colSums_dweights_dsubset(x, N, P, power, centerx, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    }
}

/* R\_colSums */

/* R\_colSums Prototype */

SEXP R_colSums
(
    /* R x Input */
    
        SEXP x,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{
    SEXP ans;
    int P;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    double center;

    P = NCOL(x);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P));
    RC_colSums(REAL(x), N, P, Power1, &center, !DoCenter, weights, subset, Offset0, 
               Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/* SimpleSums */

/* C\_Sums\_dweights\_dsubset */

double C_Sums_dweights_dsubset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    
) 
{
    double *s, *w; 
    /* Sums Body */
    

        double ans = 0.0;

        if (Nsubset > 0) {
            if (!HAS_WEIGHTS) return((double) Nsubset); 
        } else {
            if (!HAS_WEIGHTS) return((double) N);
        }

        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
            
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
        {
            w = w + diff;
            ans += w[0];
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
                
        }
        w = w + diff;
        ans += w[0];

        return(ans);
    
}

/* C\_Sums\_iweights\_dsubset */

double C_Sums_iweights_dsubset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    
) 
{
    double *s;
    int *w; 
    /* Sums Body */
    

        double ans = 0.0;

        if (Nsubset > 0) {
            if (!HAS_WEIGHTS) return((double) Nsubset); 
        } else {
            if (!HAS_WEIGHTS) return((double) N);
        }

        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
            
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
        {
            w = w + diff;
            ans += w[0];
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
                
        }
        w = w + diff;
        ans += w[0];

        return(ans);
    
}

/* C\_Sums\_iweights\_isubset */

double C_Sums_iweights_isubset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    
) 
{
    int *s, *w;
    /* Sums Body */
    

        double ans = 0.0;

        if (Nsubset > 0) {
            if (!HAS_WEIGHTS) return((double) Nsubset); 
        } else {
            if (!HAS_WEIGHTS) return((double) N);
        }

        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
            
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
        {
            w = w + diff;
            ans += w[0];
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
                
        }
        w = w + diff;
        ans += w[0];

        return(ans);
    
}

/* C\_Sums\_dweights\_isubset */

double C_Sums_dweights_isubset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    
) 
{
    int *s; 
    double *w;
    /* Sums Body */
    

        double ans = 0.0;

        if (Nsubset > 0) {
            if (!HAS_WEIGHTS) return((double) Nsubset); 
        } else {
            if (!HAS_WEIGHTS) return((double) N);
        }

        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
            
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
            
        {
            w = w + diff;
            ans += w[0];
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
                
        }
        w = w + diff;
        ans += w[0];

        return(ans);
    
}

/* RC\_Sums */

/* RC\_Sums Prototype */

double RC_Sums
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    
) 

{
    if (XLENGTH(weights) == 0) {
        if (XLENGTH(subset) == 0) {
            return((double) N);
        } else {
            return((double) Nsubset);
        }
    } 
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_iweights_isubset(N, INTEGER(weights), XLENGTH(weights),
                                           INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_iweights_dsubset(N, INTEGER(weights), XLENGTH(weights), 
                                           REAL(subset), offset, Nsubset));
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_dweights_isubset(N, REAL(weights), XLENGTH(weights),
                                           INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_dweights_dsubset(N, REAL(weights), XLENGTH(weights),
                                           REAL(subset), offset, Nsubset));
        }
    }
}

/* R\_Sums */

/* R\_Sums Prototype */

SEXP R_Sums
(
    /* R N Input */
    
        SEXP N,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{
    SEXP ans;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;

    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = RC_Sums(INTEGER(N)[0], weights, subset, Offset0, Nsubset);
    UNPROTECT(1);

    return(ans);
}


/* Tables */

/* C\_OneTableSums\_dweights\_dsubset */

void C_OneTableSums_dweights_dsubset
(
    /* C OneTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
)
{
    double *s, *w; 
    /* OneTableSums Body */
    
        int *xx;

        for (int p = 0; p < P; p++) P_ans[p] = 0.0;

        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[xx[0]] += (double) w[0];
            } else {
                P_ans[xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[xx[0]] += w[0];
        } else {
            P_ans[xx[0]]++;
        }
    
}

/* C\_OneTableSums\_iweights\_dsubset */

void C_OneTableSums_iweights_dsubset
(
    /* C OneTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
)
{
    double *s;
    int *w; 
    /* OneTableSums Body */
    
        int *xx;

        for (int p = 0; p < P; p++) P_ans[p] = 0.0;

        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[xx[0]] += (double) w[0];
            } else {
                P_ans[xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[xx[0]] += w[0];
        } else {
            P_ans[xx[0]]++;
        }
    
}

/* C\_OneTableSums\_iweights\_isubset */

void C_OneTableSums_iweights_isubset
(
    /* C OneTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
)
{
    int *s, *w;
    /* OneTableSums Body */
    
        int *xx;

        for (int p = 0; p < P; p++) P_ans[p] = 0.0;

        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[xx[0]] += (double) w[0];
            } else {
                P_ans[xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[xx[0]] += w[0];
        } else {
            P_ans[xx[0]]++;
        }
    
}

/* C\_OneTableSums\_dweights\_isubset */

void C_OneTableSums_dweights_isubset
(
    /* C OneTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
)
{
    int *s; 
    double *w;
    /* OneTableSums Body */
    
        int *xx;

        for (int p = 0; p < P; p++) P_ans[p] = 0.0;

        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[xx[0]] += (double) w[0];
            } else {
                P_ans[xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[xx[0]] += w[0];
        } else {
            P_ans[xx[0]]++;
        }
    
}

/* RC\_OneTableSums */

/* RC\_OneTableSums Prototype */

void RC_OneTableSums
(
    /* C OneTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
) 

{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_OneTableSums_iweights_isubset(x, N, P, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_OneTableSums_iweights_dsubset(x, N, P, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_OneTableSums_dweights_isubset(x, N, P, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_OneTableSums_dweights_dsubset(x, N, P, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    }
}

/* R\_OneTableSums */

/* R\_OneTableSums Prototype */

SEXP R_OneTableSums
(
    /* R x Input */
    
        SEXP x,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{

    SEXP ans;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    int P;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    P = NLEVELS(x) + 1;
    
    PROTECT(ans = allocVector(REALSXP, P));
    RC_OneTableSums(INTEGER(x), N, P, weights, subset, 
                    Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

/* C\_TwoTableSums\_dweights\_dsubset */

void C_TwoTableSums_dweights_dsubset
(
    /* C TwoTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C TwoTableSums Answer */
    
        double *PQ_ans
    
)
{
    double *s, *w; 
    /* TwoTableSums Body */
    
        int *xx, *yy;

        for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

        yy = y;
        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQ_ans[yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
    
}

/* C\_TwoTableSums\_iweights\_dsubset */

void C_TwoTableSums_iweights_dsubset
(
    /* C TwoTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C TwoTableSums Answer */
    
        double *PQ_ans
    
)
{
    double *s;
    int *w; 
    /* TwoTableSums Body */
    
        int *xx, *yy;

        for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

        yy = y;
        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQ_ans[yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
    
}

/* C\_TwoTableSums\_iweights\_isubset */

void C_TwoTableSums_iweights_isubset
(
    /* C TwoTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C TwoTableSums Answer */
    
        double *PQ_ans
    
)
{
    int *s, *w;
    /* TwoTableSums Body */
    
        int *xx, *yy;

        for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

        yy = y;
        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQ_ans[yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
    
}

/* C\_TwoTableSums\_dweights\_isubset */

void C_TwoTableSums_dweights_isubset
(
    /* C TwoTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C TwoTableSums Answer */
    
        double *PQ_ans
    
)
{
    int *s; 
    double *w;
    /* TwoTableSums Body */
    
        int *xx, *yy;

        for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

        yy = y;
        xx = x;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQ_ans[yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
    
}

/* RC\_TwoTableSums */

/* RC\_TwoTableSums Prototype */

void RC_TwoTableSums
(
    /* C TwoTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C TwoTableSums Answer */
    
        double *PQ_ans
    
) 

{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_TwoTableSums_iweights_isubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TwoTableSums_iweights_dsubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_TwoTableSums_dweights_isubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TwoTableSums_dweights_dsubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
}

/* R\_TwoTableSums */

/* R\_TwoTableSums Prototype */

SEXP R_TwoTableSums
(
    /* R x Input */
    
        SEXP x,
    
    /* R y Input */
    
        SEXP y,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{

    SEXP ans, dim;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    int P, Q;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    P = NLEVELS(x) + 1;
    Q = NLEVELS(y) + 1;
    
    PROTECT(ans = allocVector(REALSXP, P * Q));
    PROTECT(dim = allocVector(INTSXP, 2));
    INTEGER(dim)[0] = P;
    INTEGER(dim)[1] = Q;
    dimgets(ans, dim);
    RC_TwoTableSums(INTEGER(x), N, P, INTEGER(y), Q, 
                    weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(2);
    return(ans);
}

/* C\_ThreeTableSums\_dweights\_dsubset */

void C_ThreeTableSums_dweights_dsubset
(
    /* C ThreeTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        /* C integer block Input */
        
            int *block,
            /* C integer B Input */
            
                int B
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C ThreeTableSums Answer */
    
        double *PQL_ans
    
)
{
    double *s, *w; 
    /* ThreeTableSums Body */
    
        int *xx, *yy, *bb, PQ = P * Q;

        for (int p = 0; p < PQ * B; p++) PQL_ans[p] = 0.0;

        yy = y;
        xx = x;
        bb = block;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            bb = bb + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        bb = bb + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += w[0];
        } else {
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
        }
    
}

/* C\_ThreeTableSums\_iweights\_dsubset */

void C_ThreeTableSums_iweights_dsubset
(
    /* C ThreeTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        /* C integer block Input */
        
            int *block,
            /* C integer B Input */
            
                int B
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
    
    /* C real subset Input */
    
        double *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C ThreeTableSums Answer */
    
        double *PQL_ans
    
)
{
    double *s;
    int *w; 
    /* ThreeTableSums Body */
    
        int *xx, *yy, *bb, PQ = P * Q;

        for (int p = 0; p < PQ * B; p++) PQL_ans[p] = 0.0;

        yy = y;
        xx = x;
        bb = block;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            bb = bb + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        bb = bb + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += w[0];
        } else {
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
        }
    
}

/* C\_ThreeTableSums\_iweights\_isubset */

void C_ThreeTableSums_iweights_isubset
(
    /* C ThreeTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        /* C integer block Input */
        
            int *block,
            /* C integer B Input */
            
                int B
            ,
        
    
    /* C integer weights Input */
    
        int *weights,
        int HAS_WEIGHTS,
        
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C ThreeTableSums Answer */
    
        double *PQL_ans
    
)
{
    int *s, *w;
    /* ThreeTableSums Body */
    
        int *xx, *yy, *bb, PQ = P * Q;

        for (int p = 0; p < PQ * B; p++) PQL_ans[p] = 0.0;

        yy = y;
        xx = x;
        bb = block;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            bb = bb + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        bb = bb + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += w[0];
        } else {
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
        }
    
}

/* C\_ThreeTableSums\_dweights\_isubset */

void C_ThreeTableSums_dweights_isubset
(
    /* C ThreeTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        /* C integer block Input */
        
            int *block,
            /* C integer B Input */
            
                int B
            ,
        
    
    /* C real weights Input */
    
        double *weights,
        int HAS_WEIGHTS,
    
    /* C integer subset Input */
    
        int *subset,
        /* C subset range Input */
        
            R_xlen_t offset,
            /* C integer Nsubset Input */
            
                R_xlen_t Nsubset
            
        
    ,
    /* C ThreeTableSums Answer */
    
        double *PQL_ans
    
)
{
    int *s; 
    double *w;
    /* ThreeTableSums Body */
    
        int *xx, *yy, *bb, PQ = P * Q;

        for (int p = 0; p < PQ * B; p++) PQL_ans[p] = 0.0;

        yy = y;
        xx = x;
        bb = block;
        /* init subset loop */
        
            R_xlen_t diff = 0;
            s = subset + offset;
            w = weights;
            /* subset is R-style index in 1:N */
            if (Nsubset > 0)
                diff = (R_xlen_t) s[0] - 1;
        
        /* start subset loop */
        
            for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
        
        {
            xx = xx + diff;
            yy = yy + diff;
            bb = bb + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += (double) w[0];
            } else {
                PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
            }
            /* continue subset loop */
            
                if (Nsubset > 0) {
                    /* NB: diff also works with R style index */
                    diff = (R_xlen_t) s[1] - s[0];
                    if (diff < 0)
                        error("subset not sorted");
                    s++;
                } else {
                    diff = 1;
                }
            
        }
        xx = xx + diff;
        yy = yy + diff;
        bb = bb + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += w[0];
        } else {
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
        }
    
}

/* RC\_ThreeTableSums */

/* RC\_ThreeTableSums Prototype */

void RC_ThreeTableSums
(
    /* C ThreeTableSums Input */
    
        /* C integer x Input */
        
            int *x,
            /* C integer N Input */
            
                R_xlen_t N
            ,
            /* C integer P Input */
            
                int P
            ,
        
        /* C integer y Input */
        
            int *y,
            /* C integer Q Input */
            
                int Q
            ,
        
        /* C integer block Input */
        
            int *block,
            /* C integer B Input */
            
                int B
            ,
        
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C ThreeTableSums Answer */
    
        double *PQL_ans
    
) 

{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_iweights_isubset(x, N, P, y, Q, block, B, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_iweights_dsubset(x, N, P, y, Q, block, B,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQL_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_dweights_isubset(x, N, P, y, Q, block, B,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_dweights_dsubset(x, N, P, y, Q, block, B,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQL_ans);
        }
    }
}

/* R\_ThreeTableSums */

/* R\_ThreeTableSums Prototype */

SEXP R_ThreeTableSums
(
    /* R x Input */
    
        SEXP x,
    
    /* R y Input */
    
        SEXP y,
    
    /* R block Input */
    
        SEXP block
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{
    SEXP ans, dim;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    int P, Q, B;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    P = NLEVELS(x) + 1;
    Q = NLEVELS(y) + 1;
    B = NLEVELS(block);
    
    PROTECT(ans = allocVector(REALSXP, P * Q * B));
    PROTECT(dim = allocVector(INTSXP, 3));
    INTEGER(dim)[0] = P;
    INTEGER(dim)[1] = Q;
    INTEGER(dim)[2] = B;
    dimgets(ans, dim);
    RC_ThreeTableSums(INTEGER(x), N, P, INTEGER(y), Q, 
                      INTEGER(block), B,
                      weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(2);
    return(ans);
}


/* Utils */

/* C\_setup\_subset */

void C_setup_subset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    SEXP ans
)
{
    for (R_xlen_t i = 0; i < N; i++) {
        /* ans is R style index in 1:N */
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[i] = i + 1;
        } else {
            REAL(ans)[i] = (double) i + 1;
        }
    }
}

/* C\_setup\_subset\_block */

void C_setup_subset_block
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R block Input */
    
        SEXP block
    ,
    /* R blockTable Input */
    
        SEXP blockTable
    ,
    SEXP ans
)
{
    double *cumtable;
    int Nlevels = LENGTH(blockTable);

    cumtable = Calloc(Nlevels, double);
    for (int k = 0; k < Nlevels; k++) cumtable[k] = 0.0;

    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++)
        cumtable[k] = cumtable[k - 1] + REAL(blockTable)[k - 1];

    for (R_xlen_t i = 0; i < N; i++) {
        /* ans is R style index in 1:N */
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[(int) cumtable[INTEGER(block)[i]]++] = i + 1;
        } else {
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[i]]++] = (double) i + 1;
        }
    }

    Free(cumtable);
}

/* C\_order\_subset\_wrt\_block */

void C_order_subset_wrt_block
(
    /* R subset Input */
    
        SEXP subset
    ,
    /* R block Input */
    
        SEXP block
    ,
    /* R blockTable Input */
    
        SEXP blockTable
    ,
    SEXP ans
)
{
    double *cumtable;    
    int Nlevels = LENGTH(blockTable);

    cumtable = Calloc(Nlevels, double);
    for (int k = 0; k < Nlevels; k++) cumtable[k] = 0.0;

    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++)
        cumtable[k] = cumtable[k - 1] + REAL(blockTable)[k - 1];

    /* subset is R style index in 1:N */
    if (TYPEOF(subset) == INTSXP) {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            INTEGER(ans)[(int) cumtable[INTEGER(block)[INTEGER(subset)[i] - 1]]++] = INTEGER(subset)[i];
    } else {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[(R_xlen_t) REAL(subset)[i] - 1]]++] = REAL(subset)[i];
    }

    Free(cumtable); 
}

/* RC\_order\_subset\_wrt\_block */

/* RC\_order\_subset\_wrt\_block Prototype */

SEXP RC_order_subset_wrt_block
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* R block Input */
    
        SEXP block
    ,
    /* R blockTable Input */
    
        SEXP blockTable
    
)

{
    SEXP ans;
    int NOBLOCK = (XLENGTH(block) == 0 || XLENGTH(blockTable) == 2);

    if (XLENGTH(subset) > 0) {
        if (NOBLOCK) {
            return(subset);
        } else {
            PROTECT(ans = allocVector(TYPEOF(subset), XLENGTH(subset)));
            C_order_subset_wrt_block(subset, block, blockTable, ans);
            UNPROTECT(1);
            return(ans);
        }
    } else {
        PROTECT(ans = allocVector(TYPEOF(subset), N));
        if (NOBLOCK) {
            C_setup_subset(N, ans);
        } else {
            C_setup_subset_block(N, block, blockTable, ans);
        }
        UNPROTECT(1);
        return(ans);
    }
}

/* R\_order\_subset\_wrt\_block */

/* R\_order\_subset\_wrt\_block Prototype */

SEXP R_order_subset_wrt_block
(
    /* R y Input */
    
        SEXP y,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* R block Input */
    
        SEXP block
    
)

{

    /* C integer N Input */
    
        R_xlen_t N
    ;
    SEXP blockTable, ans;

    N = XLENGTH(y) / NCOL(y);

    if (XLENGTH(weights) > 0)
        error("cannot deal with weights here");

    if (NLEVELS(block) > 1) {
        PROTECT(blockTable = R_OneTableSums(block, weights, subset));
    } else {
        PROTECT(blockTable = allocVector(REALSXP, 2));
        REAL(blockTable)[0] = 0.0;
        REAL(blockTable)[1] = RC_Sums(N, weights, subset, Offset0, XLENGTH(subset));
    }
    
    PROTECT(ans = RC_order_subset_wrt_block(N, subset, block, blockTable));

    UNPROTECT(2);
    return(ans);
}


/* LinearStatistics */

/* RC\_LinearStatistic */

/* RC\_LinearStatistic Prototype */

void RC_LinearStatistic
(
    /* R x Input */
    
        SEXP x,
    
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer P Input */
    
        int P
    ,
    /* C real y Input */
    
        double *y,
        /* C integer Q Input */
        
            int Q
        ,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 

{
    double center;

    RC_KronSums(x, N, P, y, Q, !DoSymmetric, &center, &center, !DoCenter, weights, 
                subset, offset, Nsubset, PQ_ans);
}


/* Permutations */

/* RC\_setup\_subset */

/* RC\_setup\_subset Prototype */

SEXP RC_setup_subset
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{
    SEXP ans, mysubset;
    R_xlen_t sumweights;

    if (XLENGTH(weights) == 0 && XLENGTH(subset) > 0)
        return(subset);

    if (XLENGTH(subset) == 0) {
        PROTECT(mysubset = allocVector(REALSXP, N));
        C_setup_subset(N, mysubset);
    } else {
        PROTECT(mysubset = coerceVector(subset, REALSXP));
    }

    if (XLENGTH(weights) == 0) {
        UNPROTECT(1);
        return(mysubset);
    }
        
    sumweights = (R_xlen_t) RC_Sums(N, weights, mysubset, Offset0, XLENGTH(subset));
    PROTECT(ans = allocVector(REALSXP, sumweights));

    R_xlen_t itmp = 0;
    for (R_xlen_t i = 0; i < XLENGTH(mysubset); i++) {
        if (TYPEOF(weights) == REALSXP) {
            for (R_xlen_t j = 0; j < REAL(weights)[(R_xlen_t) REAL(mysubset)[i] - 1]; j++)
                REAL(ans)[itmp++] = REAL(mysubset)[i];
        } else {
            for (R_xlen_t j = 0; j < INTEGER(weights)[(R_xlen_t) REAL(mysubset)[i] - 1]; j++)
                REAL(ans)[itmp++] = REAL(mysubset)[i];
        }
    }
    UNPROTECT(2);
    return(ans);
}

/* C\_Permute */

void C_Permute
(
    double *subset,
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ,
    double *ans
) {

    R_xlen_t n = Nsubset, j;

    for (R_xlen_t i = 0; i < Nsubset; i++) {
        j = n * unif_rand();
        ans[i] = subset[j];
        subset[j] = subset[--n];
    }
}

/* C\_doPermute */

void C_doPermute
(
    double *subset,
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ,
    double *Nsubset_tmp,
    double *perm
) {
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_Permute(Nsubset_tmp, Nsubset, perm);
}

/* C\_PermuteBlock */

void C_PermuteBlock
(
    double *subset,
    double *table,
    int Nlevels,
    double *ans
) {

    double *px, *pans;

    px = subset;
    pans = ans;

    for (R_xlen_t j = 0; j < Nlevels; j++) {
        if (table[j] > 0) {
            C_Permute(px, (R_xlen_t) table[j], pans);
            px += (R_xlen_t) table[j];
            pans += (R_xlen_t) table[j];
        }
    }
}

/* C\_doPermuteBlock */

void C_doPermuteBlock
(
    double *subset,
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ,
    double *table,
    int Nlevels,
    double *Nsubset_tmp,
    double *perm
) {
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_PermuteBlock(Nsubset_tmp, table, Nlevels, perm);
}


/* ExpectationCovariances */

/* RC\_ExpectationInfluence */

/* RC\_ExpectationInfluence Prototype */

void RC_ExpectationInfluence
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R y Input */
    
        SEXP y,
    
    /* C integer Q Input */
    
        int Q
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C sumweights Input */
    
        double sumweights
    ,    
    /* C colSums Answer */
    
        double *P_ans
    
) 

{
    double center;

    RC_colSums(REAL(y), N, Q, Power1, &center, !DoCenter, weights, 
               subset, offset, Nsubset, P_ans);
    for (int q = 0; q < Q; q++) 
        P_ans[q] = P_ans[q] / sumweights;
}

/* R\_ExpectationInfluence */

/* R\_ExpectationInfluence Prototype */

SEXP R_ExpectationInfluence
(
    /* R y Input */
    
        SEXP y,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
) 

{
    SEXP ans;
    /* C integer Q Input */
    
        int Q
    ;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    double sumweights;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    sumweights = RC_Sums(N, weights, subset, Offset0, Nsubset);

    PROTECT(ans = allocVector(REALSXP, Q));
    RC_ExpectationInfluence(N, y, Q, weights, subset, Offset0, Nsubset, sumweights, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

/* RC\_CovarianceInfluence */

/* RC\_CovarianceInfluence Prototype */

void RC_CovarianceInfluence
(
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* R y Input */
    
        SEXP y,
    
    /* C integer Q Input */
    
        int Q
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    double *ExpInf,
    /* C sumweights Input */
    
        double sumweights
    ,
    int VARONLY,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 

{
    if (VARONLY) {
        RC_colSums(REAL(y), N, Q, Power2, ExpInf, DoCenter, weights, 
                   subset, offset, Nsubset, PQ_ans);
        for (int q = 0; q < Q; q++) 
            PQ_ans[q] = PQ_ans[q] / sumweights;
    } else {
        RC_KronSums(y, N, Q, REAL(y), Q, DoSymmetric, ExpInf, ExpInf, DoCenter, weights, 
                    subset, offset, Nsubset, PQ_ans);
        for (int q = 0; q < Q * (Q + 1) / 2; q++) 
            PQ_ans[q] = PQ_ans[q] / sumweights;
    }
}

/* R\_CovarianceInfluence */

/* R\_CovarianceInfluence Prototype */

SEXP R_CovarianceInfluence
(
    /* R y Input */
    
        SEXP y,
    
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    SEXP varonly
)

{
    SEXP ans;
    SEXP ExpInf;
    /* C integer Q Input */
    
        int Q
    ;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    double sumweights;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ExpInf = R_ExpectationInfluence(y, weights, subset));

    sumweights = RC_Sums(N, weights, subset, Offset0, Nsubset);

    if (INTEGER(varonly)[0]) {
        PROTECT(ans = allocVector(REALSXP, Q));
    } else {
        PROTECT(ans = allocVector(REALSXP, Q * (Q + 1) / 2));
    }
    RC_CovarianceInfluence(N, y, Q, weights, subset, Offset0, Nsubset, REAL(ExpInf), sumweights, 
                           INTEGER(varonly)[0], REAL(ans));
    UNPROTECT(2);
    return(ans);
}

/* RC\_ExpectationX */

/* RC\_ExpectationX Prototype */

void RC_ExpectationX
(
    /* R x Input */
    
        SEXP x,
    
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer P Input */
    
        int P
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    /* C OneTableSums Answer */
    
        double *P_ans
    
) 

{
    double center;

    if (TYPEOF(x) == INTSXP) {
        double* Pp1tmp = Calloc(P + 1, double);
        RC_OneTableSums(INTEGER(x), N, P + 1, weights, subset, offset, Nsubset, Pp1tmp);
        for (int p = 0; p < P; p++) P_ans[p] = Pp1tmp[p + 1];
        Free(Pp1tmp);
    } else {
        RC_colSums(REAL(x), N, P, Power1, &center, !DoCenter, weights, subset, offset, Nsubset, P_ans);
    }
}

/* R\_ExpectationX */

/* R\_ExpectationX Prototype */

SEXP R_ExpectationX
(
    /* R x Input */
    
        SEXP x,
    
    SEXP P,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    
)

{
    SEXP ans;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;

    N = XLENGTH(x) / INTEGER(P)[0];
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    RC_ExpectationX(x, N, INTEGER(P)[0], weights, subset, 
                    Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

/* RC\_CovarianceX */

/* RC\_CovarianceX Prototype */

void RC_CovarianceX
(
    /* R x Input */
    
        SEXP x,
    
    /* C integer N Input */
    
        R_xlen_t N
    ,
    /* C integer P Input */
    
        int P
    ,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    /* C subset range Input */
    
        R_xlen_t offset,
        /* C integer Nsubset Input */
        
            R_xlen_t Nsubset
        
    ,
    double *ExpX,
    int VARONLY,
    /* C KronSums Answer */
    
        double *PQ_ans
    
) 

{
    double center;

    if (TYPEOF(x) == INTSXP) {
        if (VARONLY) {
            for (int p = 0; p < P; p++) PQ_ans[p] = ExpX[p];
        } else {
            for (int p = 0; p < P * (P + 1) / 2; p++) 
                PQ_ans[p] = 0.0;
            for (int p = 0; p < P; p++)
                PQ_ans[S(p, p, P)] = ExpX[p];
        }
    } else {
        if (VARONLY) {
            RC_colSums(REAL(x), N, P, Power2, &center, !DoCenter, weights, 
                       subset, offset, Nsubset, PQ_ans);
        } else {
            RC_KronSums(x, N, P, REAL(x), P, DoSymmetric, &center, &center, !DoCenter, weights, 
                        subset, offset, Nsubset, PQ_ans);
        }
    }
}

/* R\_CovarianceX */

/* R\_CovarianceX Prototype */

SEXP R_CovarianceX
(
    /* R x Input */
    
        SEXP x,
    
    SEXP P,
    /* R weights Input */
    
        SEXP weights
    ,
    /* R subset Input */
    
        SEXP subset
    ,
    SEXP varonly
)

{
    SEXP ans;
    SEXP ExpX;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;

    N = XLENGTH(x) / INTEGER(P)[0];
    Nsubset = XLENGTH(subset);

    PROTECT(ExpX = R_ExpectationX(x, P, weights, subset));

    if (INTEGER(varonly)[0]) {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    } else {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * (INTEGER(P)[0] + 1) / 2));
    }
    RC_CovarianceX(x, N, INTEGER(P)[0], weights, subset, Offset0, Nsubset, REAL(ExpX), 
                   INTEGER(varonly)[0], REAL(ans));
    UNPROTECT(2);
    return(ans);
}

/* C\_ExpectationLinearStatistic */

void C_ExpectationLinearStatistic
(
    /* C integer P Input */
    
        int P
    ,
    /* C integer Q Input */
    
        int Q
    ,
    double *ExpInf,
    double *ExpX,
    const int add,
    double *PQ_ans
) {

    if (!add)
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] += ExpX[p] * ExpInf[q];
    }
}

/* C\_CovarianceLinearStatistic */

void C_CovarianceLinearStatistic
(
    /* C integer P Input */
    
        int P
    ,
    /* C integer Q Input */
    
        int Q
    ,
    double *CovInf,
    double *ExpX,
    double *CovX,
    /* C sumweights Input */
    
        double sumweights
    ,
    const int add,
    double *PQPQ_sym_ans
) {

    double f1 = sumweights / (sumweights - 1);
    double f2 = 1.0 / (sumweights - 1);
    double tmp, *PP_sym_tmp;


    if (P * Q == 1) {
        tmp = f1 * CovInf[0] * CovX[0];
        tmp -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
        if (add) {
            PQPQ_sym_ans[0] += tmp;
        } else {
            PQPQ_sym_ans[0] = tmp;
        }
    } else {
        PP_sym_tmp = Calloc(P * (P + 1) / 2, double);
        C_KronSums_sym_(ExpX, 1, P,
                        PP_sym_tmp);
        for (int p = 0; p < P * (P + 1) / 2; p++)
            PP_sym_tmp[p] = f1 * CovX[p] - f2 * PP_sym_tmp[p];
        C_kronecker_sym(CovInf, Q, PP_sym_tmp, P, 1 - (add >= 1),
                        PQPQ_sym_ans);
        Free(PP_sym_tmp);
    }
}

/* C\_VarianceLinearStatistic */

void C_VarianceLinearStatistic
(
    /* C integer P Input */
    
        int P
    ,
    /* C integer Q Input */
    
        int Q
    ,
    double *VarInf,
    double *ExpX,
    double *VarX,
    /* C sumweights Input */
    
        double sumweights
    ,    
    const int add,
    double *PQ_ans
) {


    if (P * Q == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX,
                                    sumweights, (add >= 1),
                                    PQ_ans);
    } else {
        double *P_tmp;
        P_tmp = Calloc(P, double);
        double f1 = sumweights / (sumweights - 1);
        double f2 = 1.0 / (sumweights - 1);
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, 1 - (add >= 1),
                    PQ_ans);
        Free(P_tmp);
    }
}


/* Test Statistics */

/* C\_maxstand\_Covariance */

double C_maxstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

/* C\_maxstand\_Variance */

double C_maxstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

/* C\_minstand\_Covariance */

double C_minstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_PosInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

/* C\_minstand\_Variance */

double C_minstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_PosInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

/* C\_maxabsstand\_Covariance */

double C_maxabsstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = fabs((linstat[p] - expect[p]) /
                  sqrt(covar_sym[S(p, p, PQ)]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

/* C\_maxabsstand\_Variance */

double C_maxabsstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(var[p]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

/* C\_quadform */

double C_quadform
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *MPinv_sym
) {

    double ans = 0.0, tmp = 0.0;

    for (int q = 0; q < PQ; q++) {
        tmp = 0.0;
        for (int p = 0; p < PQ; p++)
            tmp += (linstat[p] - expect[p]) * MPinv_sym[S(p, q, PQ)];
        ans += tmp * (linstat[q] - expect[q]);
    }

    return(ans);
}

/* C\_maxtype */

double C_maxtype
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol,
    const int alternative
) {

    double ret = 0.0;

    if (varonly) {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Variance(PQ, linstat, expect, covar, tol);
        }
    } else {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Covariance(PQ, linstat, expect, covar, tol);
        }
    }
    return(ret);
}

/* C\_standardise */

void C_standardise
(
    const int PQ,
    double *linstat,            /* in place standardisation */
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol
) {

    double var;

    for (int p = 0; p < PQ; p++) {
        if (varonly) {
            var = covar[p];
        } else {
            var = covar[S(p, p, PQ)];
        }
        if (var > tol) {
            linstat[p] = (linstat[p] - expect[p]) / sqrt(var);
        } else {
            linstat[p] = NAN;
        }
    }
}

/* C\_ordered\_Xfactor */

void C_ordered_Xfactor
(
/* maxstat Xfactor Variables */

SEXP LECV,
const int minbucket,
const int teststat,
int *wmax,
double *maxstat,
double *bmaxstat,
double *pval,
const int lower,
const int give_log

) {

    /* Setup maxstat Variables */
    
    double *linstat, *expect, *covar, *varinf, *covinf, *ExpX, *blinstat, tol, *ls;
    int P, Q, B;
    R_xlen_t nperm;

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, 
           tmp, sumleft, sumright, sumweights;
    int rank, PQ, greater;

    Q = C_get_Q(LECV);
    P = C_get_P(LECV);
    PQ = P * Q;
    B = C_get_B(LECV);
    if (B > 1) {
        if (C_get_varonly(LECV))
            error("need covarinance for maximally statistics with blocks");
        covar = C_get_Covariance(LECV);
    }
    linstat = C_get_LinearStatistic(LECV);
    expect = C_get_Expectation(LECV);
    ExpX = C_get_ExpectationX(LECV);
    /* both need to be there */
    varinf = C_get_VarianceInfluence(LECV);
    covinf = C_get_CovarianceInfluence(LECV);
    nperm = C_get_nperm(LECV);
    if (nperm > 0)
        blinstat = C_get_PermutedLinearStatistic(LECV);
    tol = C_get_tol(LECV);
     

    /* Setup maxstat Memory */
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    if (teststat == TESTSTAT_maximum) {
       mvar = Calloc(Q, double);
       /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
       mcovar = Calloc(1, double);
       mMPinv = Calloc(1, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
       /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
       mvar = Calloc(1, double);
    }
    if (nperm > 0) {
        mblinstat = Calloc(Q * nperm, double);
    } else { /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
        mblinstat = Calloc(1, double);
    }

    maxstat[0] = 0.0;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        if (teststat == TESTSTAT_maximum)
            mvar[q] = 0.0;
        for (R_xlen_t np = 0; np < nperm; np++) mblinstat[q + np * Q] = 0.0;
    }
    if (teststat == TESTSTAT_quadratic) {
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            mcovar[q] = 0.0;
    }

    sumleft = 0.0;
    sumright = 0.0;
    for (int p = 0; p < P; p++)
        sumright += ExpX[p];
    sumweights = sumright;
     

    wmax[0] = NA_INTEGER;

    for (int p = 0; p < P; p++) {
        sumleft += ExpX[p];
        sumright -= ExpX[p];

        for (int q = 0; q < Q; q++) {
            mlinstat[q] += linstat[q * P + p];
            for (R_xlen_t np = 0; np < nperm; np++)
                mblinstat[q + np * Q] += blinstat[q * P + p + np * PQ];
            mexpect[q] += expect[q * P + p];
            if (B == 1) {
                /* Compute maxstat Variance / Covariance Directly */
                
                /* does not work with blocks! */
                if (teststat == TESTSTAT_maximum) {
                    C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                                              sumweights, 0, mvar);
                } else {
                    C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                                sumweights, 0, mcovar);
                }
                
            } else {
                /* Compute maxstat Variance / Covariance from Total Covariance */
                
                if (teststat == TESTSTAT_maximum) {
                    for (int pp = 0; pp < p; pp++)
                        mvar[q] += 2 * covar[S(pp + q * P, p + P * q, P * Q)];
                     mvar[q] += covar[S(p + q * P, p + P * q, P * Q)];
                } else {
                     for (int qq = 0; qq <= q; qq++) {
                         for (int pp = 0; pp < p; pp++)
                             mcovar[S(q, qq, Q)] += 2 * covar[S(pp + q * P, p + P * qq, P * Q)];
                         mcovar[S(q, qq, Q)] += covar[S(p + q * P, p + P * qq, P * Q)];
                     }
                }
                
            }
        }

        if ((sumleft >= minbucket) && (sumright >= minbucket) && (ExpX[p] > 0)) {

            ls = mlinstat; 
            /* Compute maxstat Test Statistic */
            
            if (teststat == TESTSTAT_maximum) {
                tmp = C_maxtype(Q, ls, mexpect, mvar, 1, tol,
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, ls, mexpect, mMPinv);
            }
            
            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }

            for (R_xlen_t np = 0; np < nperm; np++) {
                ls = mblinstat + np * Q;
                /* Compute maxstat Test Statistic */
                
                if (teststat == TESTSTAT_maximum) {
                    tmp = C_maxtype(Q, ls, mexpect, mvar, 1, tol,
                                    ALTERNATIVE_twosided);
                } else {
                    C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                    tmp = C_quadform(Q, ls, mexpect, mMPinv);
                }
                
                if (tmp > bmaxstat[np])
                    bmaxstat[np] = tmp;
            }
        }
    }
    /* Compute maxstat Permutation P-Value */
    
    if (nperm > 0) {
        greater = 0;
        for (R_xlen_t np = 0; np < nperm; np++) {
            if (bmaxstat[np] > maxstat[0]) greater++;
        }
        pval[0] = C_perm_pvalue(greater, nperm, lower, give_log);
    }
    
    Free(mlinstat); Free(mexpect); Free(mblinstat); 
    Free(mvar); Free(mcovar); Free(mMPinv);
}

/* C\_unordered\_Xfactor */

void C_unordered_Xfactor
(
/* maxstat Xfactor Variables */

SEXP LECV,
const int minbucket,
const int teststat,
int *wmax,
double *maxstat,
double *bmaxstat,
double *pval,
const int lower,
const int give_log

) {

    double *mtmp;
    int qPp, nc, *levels, Pnonzero, *indl, *contrast;

    /* Setup maxstat Variables */
    
    double *linstat, *expect, *covar, *varinf, *covinf, *ExpX, *blinstat, tol, *ls;
    int P, Q, B;
    R_xlen_t nperm;

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, 
           tmp, sumleft, sumright, sumweights;
    int rank, PQ, greater;

    Q = C_get_Q(LECV);
    P = C_get_P(LECV);
    PQ = P * Q;
    B = C_get_B(LECV);
    if (B > 1) {
        if (C_get_varonly(LECV))
            error("need covarinance for maximally statistics with blocks");
        covar = C_get_Covariance(LECV);
    }
    linstat = C_get_LinearStatistic(LECV);
    expect = C_get_Expectation(LECV);
    ExpX = C_get_ExpectationX(LECV);
    /* both need to be there */
    varinf = C_get_VarianceInfluence(LECV);
    covinf = C_get_CovarianceInfluence(LECV);
    nperm = C_get_nperm(LECV);
    if (nperm > 0)
        blinstat = C_get_PermutedLinearStatistic(LECV);
    tol = C_get_tol(LECV);
    

    /* Setup maxstat Memory */
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    if (teststat == TESTSTAT_maximum) {
       mvar = Calloc(Q, double);
       /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
       mcovar = Calloc(1, double);
       mMPinv = Calloc(1, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
       /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
       mvar = Calloc(1, double);
    }
    if (nperm > 0) {
        mblinstat = Calloc(Q * nperm, double);
    } else { /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
        mblinstat = Calloc(1, double);
    }

    maxstat[0] = 0.0;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        if (teststat == TESTSTAT_maximum)
            mvar[q] = 0.0;
        for (R_xlen_t np = 0; np < nperm; np++) mblinstat[q + np * Q] = 0.0;
    }
    if (teststat == TESTSTAT_quadratic) {
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            mcovar[q] = 0.0;
    }

    sumleft = 0.0;
    sumright = 0.0;
    for (int p = 0; p < P; p++)
        sumright += ExpX[p];
    sumweights = sumright;
    
    mtmp = Calloc(P, double);

    for (int p = 0; p < P; p++) wmax[p] = NA_INTEGER;

    /* Count Levels */
    
    contrast = Calloc(P, int);
    Pnonzero = 0;
    for (int p = 0; p < P; p++) {
        if (ExpX[p] > 0) Pnonzero++;
    }
    levels = Calloc(Pnonzero, int);
    nc = 0;
    for (int p = 0; p < P; p++) {
        if (ExpX[p] > 0) {
            levels[nc] = p;
            nc++;
        }
    }

    if (Pnonzero >= 31)
        error("cannot search for unordered splits in >= 31 levels");

    int mi = 1;
    for (int l = 1; l < Pnonzero; l++) mi *= 2;
    indl = Calloc(Pnonzero, int);
    for (int p = 0; p < Pnonzero; p++) indl[p] = 0;
    

    for (int j = 1; j < mi; j++) { /* go though all splits */
 
        /* Setup unordered maxstat Contrasts */
        
        /* indl determines if level p is left or right */
        int jj = j;
        for (int l = 1; l < Pnonzero; l++) {
            indl[l] = (jj%2);
            jj /= 2;
        }

        sumleft = 0.0;
        sumright = 0.0;
        for (int p = 0; p < P; p++) contrast[p] = 0;
        for (int p = 0; p < Pnonzero; p++) {
            sumleft += indl[p] * ExpX[levels[p]];
            sumright += (1 - indl[p]) * ExpX[levels[p]];
            contrast[levels[p]] = indl[p];
        }
        

        /* Compute unordered maxstat Linear Statistic and Expectation */
        
        for (int q = 0; q < Q; q++) {
            mlinstat[q] = 0.0;
            mexpect[q] = 0.0;
            for (R_xlen_t np = 0; np < nperm; np++)
                mblinstat[q + np * Q] = 0.0;
            for (int p = 0; p < P; p++) {
                qPp = q * P + p;
                mlinstat[q] += contrast[p] * linstat[qPp];
                mexpect[q] += contrast[p] * expect[qPp];
                for (R_xlen_t np = 0; np < nperm; np++)
                    mblinstat[q + np * Q] += contrast[p] * blinstat[q * P + p + np * PQ];
            }
        }
        

        if (B == 1) {
            /* Compute unordered maxstat Variance / Covariance Directly */
            
            if (teststat == TESTSTAT_maximum) {
                C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                                          sumweights, 0, mvar);
            } else {
                C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                            sumweights, 0, mcovar);
            }
            
        } else {
            /* Compute unordered maxstat Variance / Covariance from Total Covariance */
            
            if (teststat == TESTSTAT_maximum) {
                for (int q = 0; q < Q; q++) {
                    mvar[q] = 0.0;
                    for (int p = 0; p < P; p++) {
                        qPp = q * P + p;
                        mtmp[p] = 0.0;
                        for (int pp = 0; pp < P; pp++)
                            mtmp[p] += contrast[pp] * covar[S(pp + q * P, qPp, PQ)];
                    }
                    for (int p = 0; p < P; p++)
                        mvar[q] += contrast[p] * mtmp[p];
                }
            } else {
                for (int q = 0; q < Q; q++) {
                    for (int qq = 0; qq <= q; qq++)
                        mcovar[S(q, qq, Q)] = 0.0;
                    for (int qq = 0; qq <= q; qq++) {
                        for (int p = 0; p < P; p++) {
                            mtmp[p] = 0.0;
                            for (int pp = 0; pp < P; pp++)
                                mtmp[p] += contrast[pp] * covar[S(pp + q * P, p + P * qq, P * Q)];
                        }
                        for (int p = 0; p < P; p++)
                            mcovar[S(q, qq, Q)] += contrast[p] * mtmp[p];
                    }
                }
            }
            
        }

        if ((sumleft >= minbucket) && (sumright >= minbucket)) {

            ls = mlinstat;
            /* Compute maxstat Test Statistic */
            
            if (teststat == TESTSTAT_maximum) {
                tmp = C_maxtype(Q, ls, mexpect, mvar, 1, tol,
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, ls, mexpect, mMPinv);
            }
            
            if (tmp > maxstat[0]) {
                for (int p = 0; p < Pnonzero; p++)
                    wmax[levels[p]] = contrast[levels[p]];
                maxstat[0] = tmp;
            }

            for (R_xlen_t np = 0; np < nperm; np++) {
                ls = mblinstat + np * Q;
                /* Compute maxstat Test Statistic */
                
                if (teststat == TESTSTAT_maximum) {
                    tmp = C_maxtype(Q, ls, mexpect, mvar, 1, tol,
                                    ALTERNATIVE_twosided);
                } else {
                    C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                    tmp = C_quadform(Q, ls, mexpect, mMPinv);
                }
                
                if (tmp > bmaxstat[np])
                    bmaxstat[np] = tmp;
            }
        }
    }

    /* Compute maxstat Permutation P-Value */
    
    if (nperm > 0) {
        greater = 0;
        for (R_xlen_t np = 0; np < nperm; np++) {
            if (bmaxstat[np] > maxstat[0]) greater++;
        }
        pval[0] = C_perm_pvalue(greater, nperm, lower, give_log);
    }
    

    Free(mlinstat); Free(mexpect); Free(levels); Free(contrast); Free(indl); Free(mtmp);
    Free(mblinstat); Free(mvar); Free(mcovar); Free(mMPinv);
}




/* User Interface */

/* RC\_ExpectationCovarianceStatistic */

void RC_ExpectationCovarianceStatistic
(
/* User Interface Inputs */

/* R x Input */

    SEXP x,

/* R y Input */

    SEXP y,

/* R weights Input */

    SEXP weights
,
/* R subset Input */

    SEXP subset
,
/* R block Input */

    SEXP block
,

SEXP ans
) {

    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer P Input */
    
        int P
    ;
    /* C integer Q Input */
    
        int Q
    ;
    /* C integer B Input */
    
        int B
    ;
    double *sumweights, *table;
    double *ExpInf, *VarInf, *CovInf, *ExpX, *ExpXtotal, *VarX, *CovX;
    double *tmpV, *tmpCV;
    SEXP nullvec, subset_block;

    /* Extract Dimensions */
    
    P = C_get_P(ans);
    Q = C_get_Q(ans);
    N = NROW(x);
    B = C_get_B(ans);
    

    /* Compute Linear Statistic */
    
    RC_LinearStatistic(x, N, P, REAL(y), Q, weights, subset, 
                       Offset0, XLENGTH(subset),
                       C_get_LinearStatistic(ans));
    

    /* Setup Memory and Subsets in Blocks */
    
    ExpInf = C_get_ExpectationInfluence(ans);
    VarInf = C_get_VarianceInfluence(ans);
    CovInf = C_get_CovarianceInfluence(ans);
    ExpXtotal = C_get_ExpectationX(ans);
    for (int p = 0; p < P; p++) ExpXtotal[p] = 0.0;
    ExpX = Calloc(P, double);
    VarX = Calloc(P, double);
    CovX = Calloc(P * (P + 1) / 2, double);
    table = C_get_TableBlock(ans);
    sumweights = C_get_Sumweights(ans);
    PROTECT(nullvec = allocVector(INTSXP, 0));

    if (B == 1) {
        table[0] = 0.0;
        table[1] = RC_Sums(N, nullvec, subset, Offset0, XLENGTH(subset));
    } else {
        RC_OneTableSums(INTEGER(block), N, B + 1, nullvec, subset, Offset0, 
                        XLENGTH(subset), table);
    }
    if (table[0] > 0)
        error("No missing values allowed in block");
    PROTECT(subset_block = RC_order_subset_wrt_block(N, subset, block, 
                                                     VECTOR_ELT(ans, TableBlock_SLOT)));
    

    /* start with subset[0] */
    R_xlen_t offset = (R_xlen_t) table[0];

    for (int b = 0; b < B; b++) {

        /* compute sum of weights in block b of subset */
        sumweights[b] = RC_Sums(N, weights, subset_block, 
                                offset, (R_xlen_t) table[b + 1]);

        /* don't do anything for empty blocks */
        if (sumweights[b] > 0) {

            /* Compute Expectation Linear Statistic */
            
            RC_ExpectationInfluence(N, y, Q, weights, subset_block, offset, 
                                    (R_xlen_t) table[b + 1], sumweights[b], ExpInf + b * Q);
            RC_ExpectationX(x, N, P, weights, subset_block, offset, 
                            (R_xlen_t) table[b + 1], ExpX);
            for (int p = 0; p < P; p++) ExpXtotal[p] += ExpX[p];
            C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b, 
                                         C_get_Expectation(ans));
            

            /* Compute Covariance Influence */
            
            /* C_ordered_Xfactor and C_unordered_Xfactor need both VarInf and CovInf */
            RC_CovarianceInfluence(N, y, Q, weights, subset_block, offset, 
                                  (R_xlen_t) table[b + 1], ExpInf + b * Q, sumweights[b], 
                                  !DoVarOnly, CovInf + b * Q * (Q + 1) / 2);
            /* extract variance from covariance */
            tmpCV = CovInf + b * Q * (Q + 1) / 2;
            tmpV = VarInf + b * Q;
            for (int q = 0; q < Q; q++) tmpV[q] = tmpCV[S(q, q, Q)];
            

            if (C_get_varonly(ans)) {
                /* Compute Variance Linear Statistic */
                
                RC_CovarianceX(x, N, P, weights, subset_block, offset, 
                               (R_xlen_t) table[b + 1], ExpX, DoVarOnly, VarX);
                C_VarianceLinearStatistic(P, Q, VarInf + b * Q, ExpX, VarX, sumweights[b], 
                                          b, C_get_Variance(ans));
                
            } else {
                /* Compute Covariance Linear Statistic */
                
                RC_CovarianceX(x, N, P, weights, subset_block, offset, 
                               (R_xlen_t) table[b + 1], ExpX, !DoVarOnly, CovX);
                C_CovarianceLinearStatistic(P, Q, CovInf + b * Q * (Q + 1) / 2,
                                            ExpX, CovX, sumweights[b], b,
                                            C_get_Covariance(ans));
                
            }
        }

        /* next iteration starts with subset[table[b + 1]] */
        offset += (R_xlen_t) table[b + 1];
    }

    Free(ExpX); Free(VarX); Free(CovX);
    UNPROTECT(2);
}

/* R\_ExpectationCovarianceStatistic */

/* R\_ExpectationCovarianceStatistic Prototype */

SEXP R_ExpectationCovarianceStatistic
(
/* User Interface Inputs */

/* R x Input */

    SEXP x,

/* R y Input */

    SEXP y,

/* R weights Input */

    SEXP weights
,
/* R subset Input */

    SEXP subset
,
/* R block Input */

    SEXP block
,

SEXP varonly,
SEXP tol
)

{
    SEXP ans;

    /* Setup Dimensions */
    
        int P, Q, B;

        if (isInteger(x)) {
            P = NLEVELS(x);
        } else {
            P = NCOL(x);
        }
        Q = NCOL(y);

        B = 1;
        if (LENGTH(block) > 0)
            B = NLEVELS(block);
    

    PROTECT(ans = RC_init_LECV_1d(P, Q, INTEGER(varonly)[0], B, isInteger(x), REAL(tol)[0]));

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, ans);

    UNPROTECT(1);
    return(ans);
}

/* R\_PermutedLinearStatistic */

/* R\_PermutedLinearStatistic Prototype */

SEXP R_PermutedLinearStatistic
(
    /* User Interface Inputs */
    
    /* R x Input */

        SEXP x,
    
    /* R y Input */

        SEXP y,
    
    /* R weights Input */

        SEXP weights
    ,
    /* R subset Input */

        SEXP subset
    ,
    /* R block Input */

        SEXP block
    ,
    
    SEXP nperm
)

{
    SEXP ans, expand_subset, block_subset, perm, tmp, blockTable;
    double *linstat;
    int PQ;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    R_xlen_t inperm;

    /* Setup Dimensions */
    
        int P, Q, B;

        if (isInteger(x)) {
            P = NLEVELS(x);
        } else {
            P = NCOL(x);
        }
        Q = NCOL(y);

        B = 1;
        if (LENGTH(block) > 0)
            B = NLEVELS(block);
    
    PQ = P * Q;
    N = NROW(y);
    inperm = (R_xlen_t) REAL(nperm)[0];

    PROTECT(ans = allocMatrix(REALSXP, PQ, inperm));
    PROTECT(expand_subset = RC_setup_subset(N, weights, subset));
    Nsubset = XLENGTH(expand_subset);
    PROTECT(tmp = allocVector(REALSXP, Nsubset));
    PROTECT(perm = allocVector(REALSXP, Nsubset));

    GetRNGstate();
    if (B == 1) {
        for (R_xlen_t np = 0; np < inperm; np++) {
            /* Setup Linear Statistic */
            
            if (np % 256 == 0) R_CheckUserInterrupt();
            linstat = REAL(ans) + PQ * np;
            for (int p = 0; p < PQ; p++)
                linstat[p] = 0.0;
            
            C_doPermute(REAL(expand_subset), Nsubset, REAL(tmp), REAL(perm));
            RC_KronSums_Permutation(x, NROW(x), P, REAL(y), Q, expand_subset, 
                                    Offset0, Nsubset, perm, linstat);
        }
    } else {
        PROTECT(blockTable = allocVector(REALSXP, B + 1));
        /* same as RC_OneTableSums(block, noweights, expand_subset) */
        RC_OneTableSums(INTEGER(block), XLENGTH(block), B + 1, weights, subset, Offset0,
                        XLENGTH(subset), REAL(blockTable));
        PROTECT(block_subset = RC_order_subset_wrt_block(XLENGTH(block), expand_subset, 
                                                         block, blockTable));

        for (R_xlen_t np = 0; np < inperm; np++) {
            /* Setup Linear Statistic */
            
            if (np % 256 == 0) R_CheckUserInterrupt();
            linstat = REAL(ans) + PQ * np;
            for (int p = 0; p < PQ; p++)
                linstat[p] = 0.0;
            
            C_doPermuteBlock(REAL(block_subset), Nsubset, REAL(blockTable), 
                             B + 1, REAL(tmp), REAL(perm));
            RC_KronSums_Permutation(x, NROW(x), P, REAL(y), Q, block_subset, 
                                    Offset0, Nsubset, perm, linstat);
        }
        UNPROTECT(2);
    }
    PutRNGstate();

    UNPROTECT(4);
    return(ans);
}

/* R\_StandardisePermutedLinearStatistic */

/* R\_StandardisePermutedLinearStatistic Prototype */

SEXP R_StandardisePermutedLinearStatistic
(
    SEXP LECV
)

{
    SEXP ans;
    R_xlen_t nperm = C_get_nperm(LECV);
    double *ls;
    if (!nperm) return(R_NilValue);
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    
    PROTECT(ans = allocMatrix(REALSXP, PQ, nperm));

    for (R_xlen_t np = 0; np < nperm; np++) {
        ls = REAL(ans) + PQ * np;
        /* copy first; standarisation is in place */
        for (int p = 0; p < PQ; p++) 
            ls[p] = C_get_PermutedLinearStatistic(LECV)[p + PQ * np];
        if (C_get_varonly(LECV)) {
            C_standardise(PQ, ls, C_get_Expectation(LECV),
                          C_get_Variance(LECV), 1, C_get_tol(LECV));
        } else {
            C_standardise(PQ, ls, C_get_Expectation(LECV),
                          C_get_Covariance(LECV), 0, C_get_tol(LECV));
        }
    }
    UNPROTECT(1);
    return(ans);
}


/* 2d User Interface */

/* RC\_ExpectationCovarianceStatistic\_2d */

void RC_ExpectationCovarianceStatistic_2d
(
/* 2d User Interface Inputs */

/* R x Input */

    SEXP x,

SEXP ix,
/* R y Input */

    SEXP y,

SEXP iy,
/* R weights Input */

    SEXP weights
,
/* R subset Input */

    SEXP subset
,
/* R block Input */

    SEXP block
,

SEXP ans
) {

    SEXP Rcsum, Rrsum;
    int P, Q, Lxp1, Lyp1, B, Xfactor;
    double *ExpInf, *ExpX, *CovX;
    double *table, *table2d, *csum, *rsum, *sumweights, *btab, *linstat;

    P = C_get_P(ans);
    Q = C_get_Q(ans);

    ExpInf = C_get_ExpectationInfluence(ans);
    ExpX = C_get_ExpectationX(ans);
    table = C_get_Table(ans);
    sumweights = C_get_Sumweights(ans);

    Lxp1 = C_get_dimTable(ans)[0];
    Lyp1 = C_get_dimTable(ans)[1];
    B = C_get_B(ans);
    Xfactor = C_get_Xfactor(ans);

    if (C_get_varonly(ans)) {
        CovX = Calloc(P, double);
    } else {
        CovX = Calloc(P * (P + 1) / 2, double);
    }

    table2d = Calloc(Lxp1 * Lyp1, double);
    PROTECT(Rcsum = allocVector(REALSXP, Lyp1));
    csum = REAL(Rcsum);
    PROTECT(Rrsum = allocVector(REALSXP, Lxp1));
    rsum = REAL(Rrsum);

    /* 2d Total Table */
    
    for (int i = 0; i < Lxp1 * Lyp1; i++)
        table2d[i] = 0.0;
    for (int b = 0; b < B; b++) {
        for (int i = 0; i < Lxp1; i++) {
            for (int j = 0; j < Lyp1; j++)
                table2d[j * Lxp1 + i] += table[b * Lxp1 * Lyp1 + j * Lxp1 + i];
        }
    }
    

    linstat = C_get_LinearStatistic(ans);
    for (int p = 0; p < P * Q; p++)
        linstat[p] = 0.0;

    for (int b = 0; b < B; b++) {
        btab = table + Lxp1 * Lyp1 * b;

        /* Linear Statistic 2d */
        
        if (Xfactor) {
            for (int j = 1; j < Lyp1; j++) { /* j = 0 means NA */
                for (int i = 1; i < Lxp1; i++) { /* i = 0 means NA */
                    for (int q = 0; q < Q; q++)
                        linstat[q * (Lxp1 - 1) + (i - 1)] +=
                            btab[j * Lxp1 + i] * REAL(y)[q * Lyp1 + j];
                }
            }
        } else {
            for (int p = 0; p < P; p++) {
                for (int q = 0; q < Q; q++) {
                    int qPp = q * P + p;
                    int qLy = q * Lyp1;
                    for (int i = 0; i < Lxp1; i++) {
                        int pLxi = p * Lxp1 + i;
                        for (int j = 0; j < Lyp1; j++)
                            linstat[qPp] += REAL(y)[qLy + j] * REAL(x)[pLxi] * btab[j * Lxp1 + i];
                    }
                }
            }
        }
        

        /* Col Row Total Sums */
        
        /* Remember: first row / column count NAs */
        /* column sums */
        for (int q = 1; q < Lyp1; q++) {
            csum[q] = 0;
            for (int p = 1; p < Lxp1; p++)
                csum[q] += btab[q * Lxp1 + p];
        }
        csum[0] = 0; /* NA */
        /* row sums */
        for (int p = 1; p < Lxp1; p++)  {
            rsum[p] = 0;
            for (int q = 1; q < Lyp1; q++)
                rsum[p] += btab[q * Lxp1 + p];
        }
        rsum[0] = 0; /* NA */
        /* total sum */
        sumweights[b] = 0;
        for (int i = 1; i < Lxp1; i++) sumweights[b] += rsum[i];
        

        /* 2d Expectation */
        
        RC_ExpectationInfluence(NROW(y), y, Q, Rcsum, subset, Offset0, 0, sumweights[b], ExpInf);

        if (LENGTH(x) == 0) {
            for (int p = 0; p < P; p++)
                ExpX[p] = rsum[p + 1];
            } else {
                RC_ExpectationX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX);
        }

        C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, C_get_Expectation(ans));
        

        /* 2d Covariance */
        
        /* C_ordered_Xfactor needs both VarInf and CovInf */
        RC_CovarianceInfluence(NROW(y), y, Q, Rcsum, subset, Offset0, 0, ExpInf, sumweights[b],
                               !DoVarOnly, C_get_CovarianceInfluence(ans));
        for (int q = 0; q < Q; q++) 
            C_get_VarianceInfluence(ans)[q] = C_get_CovarianceInfluence(ans)[S(q, q, Q)];

        if (C_get_varonly(ans)) {
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P; p++) CovX[p] = ExpX[p];
            } else {
                RC_CovarianceX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX, DoVarOnly, CovX);
            }
            C_VarianceLinearStatistic(P, Q, C_get_VarianceInfluence(ans),
                                      ExpX, CovX, sumweights[b], b,
                                      C_get_Variance(ans));
        } else {
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P * (P + 1) / 2; p++) CovX[p] = 0.0;
                for (int p = 0; p < P; p++) CovX[S(p, p, P)] = ExpX[p];
            } else {
                RC_CovarianceX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX, !DoVarOnly, CovX);
            }
            C_CovarianceLinearStatistic(P, Q, C_get_CovarianceInfluence(ans),
                                        ExpX, CovX, sumweights[b], b,
                                        C_get_Covariance(ans));
        }
        

    }
    Free(table2d); 
    UNPROTECT(2);
}

/* R\_ExpectationCovarianceStatistic\_2d */

/* R\_ExpectationCovarianceStatistic\_2d Prototype */

SEXP R_ExpectationCovarianceStatistic_2d
(
/* 2d User Interface Inputs */

/* R x Input */

    SEXP x,

SEXP ix,
/* R y Input */

    SEXP y,

SEXP iy,
/* R weights Input */

    SEXP weights
,
/* R subset Input */

    SEXP subset
,
/* R block Input */

    SEXP block
,

SEXP varonly,
SEXP tol
)

{
    SEXP ans;
    /* C integer N Input */
    
        R_xlen_t N
    ;
    /* C integer Nsubset Input */
    
        R_xlen_t Nsubset
    ;
    int Xfactor;

    N = XLENGTH(ix);
    Nsubset = XLENGTH(subset);
    Xfactor = XLENGTH(x) == 0;

    /* Setup Dimensions 2d */
    
    int P, Q, B, Lx, Ly;

    if (XLENGTH(x) == 0) {
        P = NLEVELS(ix);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);

    B = 1;
    if (XLENGTH(block) > 0)
        B = NLEVELS(block);

    Lx = NLEVELS(ix);
    Ly = NLEVELS(iy);
    

    PROTECT(ans = RC_init_LECV_2d(P, Q, INTEGER(varonly)[0], 
                                  Lx, Ly, B, Xfactor, REAL(tol)[0]));

    if (B == 1) {
        RC_TwoTableSums(INTEGER(ix), N, Lx + 1, INTEGER(iy), Ly + 1, 
                        weights, subset, Offset0, Nsubset, 
                        C_get_Table(ans));
    } else {
        RC_ThreeTableSums(INTEGER(ix), N, Lx + 1, INTEGER(iy), Ly + 1, 
                          INTEGER(block), B, weights, subset, Offset0, Nsubset, 
                          C_get_Table(ans));
    }
    RC_ExpectationCovarianceStatistic_2d(x, ix, y, iy, weights,
                                         subset, block, ans);

    UNPROTECT(1);
    return(ans);
}

/* R\_PermutedLinearStatistic\_2d */

/* R\_PermutedLinearStatistic\_2d Prototype */

SEXP R_PermutedLinearStatistic_2d
(
    /* 2d User Interface Inputs */
    
    /* R x Input */

        SEXP x,
    
    SEXP ix,
    /* R y Input */

        SEXP y,
    
    SEXP iy,
    /* R weights Input */

        SEXP weights
    ,
    /* R subset Input */

        SEXP subset
    ,
    /* R block Input */

        SEXP block
    ,
    
    SEXP nperm,
    SEXP itable
)

{
    SEXP ans, Ritable;
    int *csum, *rsum, *sumweights, *jwork, *table, *rtable2, maxn = 0, Lxp1, Lyp1, *btab, PQ, Xfactor;
    R_xlen_t inperm;
    double *fact, *linstat;

    /* Setup Dimensions 2d */
    
    int P, Q, B, Lx, Ly;

    if (XLENGTH(x) == 0) {
        P = NLEVELS(ix);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);

    B = 1;
    if (XLENGTH(block) > 0)
        B = NLEVELS(block);

    Lx = NLEVELS(ix);
    Ly = NLEVELS(iy);
    

    PQ = P * Q;
    Xfactor = XLENGTH(x) == 0;
    Lxp1 = Lx + 1;
    Lyp1 = Ly + 1;
    inperm = (R_xlen_t) REAL(nperm)[0];

    PROTECT(ans = allocMatrix(REALSXP, PQ, inperm));

    /* Setup Working Memory */
    
    csum = Calloc(Lyp1 * B, int);
    rsum = Calloc(Lxp1 * B, int);
    sumweights = Calloc(B, int);
    table = Calloc(Lxp1 * Lyp1, int);
    rtable2 = Calloc(Lx * Ly , int);
    jwork = Calloc(Lyp1, int);
    

    /* Convert Table to Integer */
    
    PROTECT(Ritable = allocVector(INTSXP, LENGTH(itable)));
    for (int i = 0; i < LENGTH(itable); i++) {
        if (REAL(itable)[i] > INT_MAX)
            error("cannot deal with weights larger INT_MAX in R_PermutedLinearStatistic_2d");
        INTEGER(Ritable)[i] = (int) REAL(itable)[i];
    }
    

    for (int b = 0; b < B; b++) {
        btab = INTEGER(Ritable) + Lxp1 * Lyp1 * b;
        /* Col Row Total Sums */
        
        /* Remember: first row / column count NAs */
        /* column sums */
        for (int q = 1; q < Lyp1; q++) {
            csum[q] = 0;
            for (int p = 1; p < Lxp1; p++)
                csum[q] += btab[q * Lxp1 + p];
        }
        csum[0] = 0; /* NA */
        /* row sums */
        for (int p = 1; p < Lxp1; p++)  {
            rsum[p] = 0;
            for (int q = 1; q < Lyp1; q++)
                rsum[p] += btab[q * Lxp1 + p];
        }
        rsum[0] = 0; /* NA */
        /* total sum */
        sumweights[b] = 0;
        for (int i = 1; i < Lxp1; i++) sumweights[b] += rsum[i];
        
        if (sumweights[b] > maxn) maxn = sumweights[b];
    }

    /* Setup Log-Factorials */
    
    fact = Calloc(maxn + 1, double);
    /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
    fact[0] = fact[1] = 0.;
    for(int j = 2; j <= maxn; j++)
        fact[j] = fact[j - 1] + log(j);
    

    GetRNGstate();

    for (R_xlen_t np = 0; np < inperm; np++) {

        /* Setup Linear Statistic */
        
        if (np % 256 == 0) R_CheckUserInterrupt();
        linstat = REAL(ans) + PQ * np;
        for (int p = 0; p < PQ; p++)
            linstat[p] = 0.0;
        

        for (int p = 0; p < Lxp1 * Lyp1; p++)
            table[p] = 0;

        for (int b = 0; b < B; b++) {
            /* Compute Permuted Linear Statistic 2d */
            
            S_rcont2(&Lx, &Ly, rsum + Lxp1 * b + 1,
                     csum + Lyp1 *b + 1, sumweights + b, fact, jwork, rtable2);

            for (int j1 = 1; j1 <= Lx; j1++) {
                for (int j2 = 1; j2 <= Ly; j2++)
                    table[j2 * Lxp1 + j1] = rtable2[(j2 - 1) * Lx + (j1 - 1)];
            }
            btab = table;
            /* Linear Statistic 2d */

            if (Xfactor) {
                for (int j = 1; j < Lyp1; j++) { /* j = 0 means NA */
                    for (int i = 1; i < Lxp1; i++) { /* i = 0 means NA */
                        for (int q = 0; q < Q; q++)
                            linstat[q * (Lxp1 - 1) + (i - 1)] +=
                                btab[j * Lxp1 + i] * REAL(y)[q * Lyp1 + j];
                    }
                }
            } else {
                for (int p = 0; p < P; p++) {
                    for (int q = 0; q < Q; q++) {
                        int qPp = q * P + p;
                        int qLy = q * Lyp1;
                        for (int i = 0; i < Lxp1; i++) {
                            int pLxi = p * Lxp1 + i;
                            for (int j = 0; j < Lyp1; j++)
                                linstat[qPp] += REAL(y)[qLy + j] * REAL(x)[pLxi] * btab[j * Lxp1 + i];
                        }
                    }
                }
            }
            
            
        }
    }

    PutRNGstate();

    Free(csum); Free(rsum); Free(sumweights); Free(rtable2);
    Free(jwork); Free(fact);
    UNPROTECT(2);
    return(ans);
}


/* Tests */

/* R\_QuadraticTest */

/* R\_QuadraticTest Prototype */

SEXP R_QuadraticTest
(
    /* R LECV Input */
    
    SEXP LECV
    ,
    SEXP pvalue,
    SEXP lower,
    SEXP give_log,
    SEXP PermutedStatistics
)

{

    SEXP ans, stat, pval, names, permstat;
    double *MPinv, *ls, st, pst, *ex;
    int rank, P, Q, PQ, greater = 0;
    R_xlen_t nperm;

    /* Setup Test Memory */
    
    P = C_get_P(LECV);
    Q = C_get_Q(LECV);
    PQ = P * Q;

    if (C_get_varonly(LECV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    if (C_get_nperm(LECV) > 0 && INTEGER(PermutedStatistics)[0]) {
        PROTECT(ans = allocVector(VECSXP, 3));
        PROTECT(names = allocVector(STRSXP, 3));
        SET_VECTOR_ELT(ans, 2, permstat = allocVector(REALSXP, C_get_nperm(LECV)));
        SET_STRING_ELT(names, 2, mkChar("PermutedStatistics"));
    } else {
        PROTECT(ans = allocVector(VECSXP, 2));
        PROTECT(names = allocVector(STRSXP, 2));
    }
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("TestStatistic"));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 1, mkChar("p.value"));
    namesgets(ans, names);
    REAL(pval)[0] = NA_REAL;
    int LOWER = INTEGER(lower)[0];
    int GIVELOG = INTEGER(give_log)[0];
    int PVALUE = INTEGER(pvalue)[0];
    int PSTAT = INTEGER(PermutedStatistics)[0];
    

    MPinv = C_get_MPinv(LECV);
    C_MPinv_sym(C_get_Covariance(LECV), PQ, C_get_tol(LECV), MPinv, &rank);

    REAL(stat)[0] = C_quadform(PQ, C_get_LinearStatistic(LECV),
                               C_get_Expectation(LECV), MPinv);

    if (!PVALUE) {
        UNPROTECT(2);
        return(ans);
    }

    if (C_get_nperm(LECV) == 0) {
        REAL(pval)[0] = C_chisq_pvalue(REAL(stat)[0], rank, LOWER, GIVELOG);
    } else {
        nperm = C_get_nperm(LECV);
        ls = C_get_PermutedLinearStatistic(LECV);
        st = REAL(stat)[0];
        ex = C_get_Expectation(LECV);
        greater = 0;
        for (R_xlen_t np = 0; np < nperm; np++) {
            pst = C_quadform(PQ, ls + PQ * np, ex, MPinv);
            if (GE(pst, st, C_get_tol(LECV)))
                greater++;
            if (PSTAT) REAL(permstat)[np] = pst;
        }
        REAL(pval)[0] = C_perm_pvalue(greater, nperm, LOWER, GIVELOG);
    }

    UNPROTECT(2);
    return(ans);
}

/* R\_MaximumTest */

/* R\_MaximumTest Prototype */

SEXP R_MaximumTest
(
    /* R LECV Input */
    
    SEXP LECV
    ,
    SEXP alternative,
    SEXP pvalue,
    SEXP lower,
    SEXP give_log,
    SEXP PermutedStatistics,
    SEXP maxpts,
    SEXP releps,
    SEXP abseps
)

{
    SEXP ans, stat, pval, names, permstat;
    double st, pst, *ex, *cv, *ls, tl;
    int P, Q, PQ, vo, alt, greater;
    R_xlen_t nperm;

    /* Setup Test Memory */
    
    P = C_get_P(LECV);
    Q = C_get_Q(LECV);
    PQ = P * Q;

    if (C_get_varonly(LECV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    if (C_get_nperm(LECV) > 0 && INTEGER(PermutedStatistics)[0]) {
        PROTECT(ans = allocVector(VECSXP, 3));
        PROTECT(names = allocVector(STRSXP, 3));
        SET_VECTOR_ELT(ans, 2, permstat = allocVector(REALSXP, C_get_nperm(LECV)));
        SET_STRING_ELT(names, 2, mkChar("PermutedStatistics"));
    } else {
        PROTECT(ans = allocVector(VECSXP, 2));
        PROTECT(names = allocVector(STRSXP, 2));
    }
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("TestStatistic"));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 1, mkChar("p.value"));
    namesgets(ans, names);
    REAL(pval)[0] = NA_REAL;
    int LOWER = INTEGER(lower)[0];
    int GIVELOG = INTEGER(give_log)[0];
    int PVALUE = INTEGER(pvalue)[0];
    int PSTAT = INTEGER(PermutedStatistics)[0];
    

    if (C_get_varonly(LECV)) {
        cv = C_get_Variance(LECV);
    } else {
        cv = C_get_Covariance(LECV);
    }

    REAL(stat)[0] =  C_maxtype(PQ, C_get_LinearStatistic(LECV),
                               C_get_Expectation(LECV),
                               cv, /* C_get_Covariance(LECV), */
                               C_get_varonly(LECV),
                               C_get_tol(LECV),
                               INTEGER(alternative)[0]);

    if (!PVALUE) {
        UNPROTECT(2);
        return(ans);
    }

    if (C_get_nperm(LECV) == 0) {
        if (C_get_varonly(LECV) && PQ > 1) {
            REAL(pval)[0] = NA_REAL;
            UNPROTECT(2);
            return(ans);
        }

        REAL(pval)[0] = C_maxtype_pvalue(REAL(stat)[0], cv, 
                                         PQ, INTEGER(alternative)[0], LOWER, GIVELOG, 
                                         INTEGER(maxpts)[0], REAL(releps)[0],
                                         REAL(abseps)[0], C_get_tol(LECV));
    } else {
        nperm = C_get_nperm(LECV);
        ls = C_get_PermutedLinearStatistic(LECV);
        ex = C_get_Expectation(LECV);
/*        cv = C_get_Covariance(LECV); */
        vo = C_get_varonly(LECV);
        alt = INTEGER(alternative)[0];
        st = REAL(stat)[0];
        tl = C_get_tol(LECV);
        greater = 0;
        for (R_xlen_t np = 0; np < nperm; np++) {
            pst = C_maxtype(PQ, ls + PQ * np, ex, cv, vo, tl, alt);
            if (alt == ALTERNATIVE_less) {
                if (LE(pst, st, tl))
                    greater++;
            } else {
                if (GE(pst, st, tl))
                    greater++;
            }
            if (PSTAT) REAL(permstat)[np] = pst;
        }
        REAL(pval)[0] = C_perm_pvalue(greater, nperm, LOWER, GIVELOG);
    }

    UNPROTECT(2);
    return(ans);
}

/* R\_MaximallySelectedTest */

/* R\_MaximallySelectedTest Prototype */

SEXP R_MaximallySelectedTest
(
    SEXP LECV,
    SEXP ordered,
    SEXP teststat,
    SEXP minbucket,
    SEXP lower,
    SEXP give_log
)

{

    SEXP ans, index, stat, pval, names, permstat;
    int P, mb;

    P = C_get_P(LECV);
    mb = INTEGER(minbucket)[0];

    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(names = allocVector(STRSXP, 4));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("TestStatistic"));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, 1, mkChar("p.value"));
    SET_VECTOR_ELT(ans, 3, permstat = allocVector(REALSXP, C_get_nperm(LECV)));
    SET_STRING_ELT(names, 3, mkChar("PermutedStatistics"));
    REAL(pval)[0] = NA_REAL;

    if (INTEGER(ordered)[0]) {
        SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, 1));
        C_ordered_Xfactor(LECV, mb, INTEGER(teststat)[0],
                          INTEGER(index), REAL(stat), REAL(permstat),
                          REAL(pval), INTEGER(lower)[0],
                          INTEGER(give_log)[0]);
        if (REAL(stat)[0] > 0)
            INTEGER(index)[0]++; /* R style indexing */
    } else {
        SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, P));
        C_unordered_Xfactor(LECV, mb, INTEGER(teststat)[0],
                            INTEGER(index), REAL(stat), REAL(permstat),
                            REAL(pval), INTEGER(lower)[0],
                            INTEGER(give_log)[0]);
    }

    SET_STRING_ELT(names, 2, mkChar("index"));
    namesgets(ans, names);

    UNPROTECT(2);
    return(ans);
}



