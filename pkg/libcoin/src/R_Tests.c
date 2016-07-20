
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "Distributions.h"
#include "MemoryAccess.h"
#include "MaxSelect.h"

SEXP R_ChisqTest
(
    SEXP LEV,
    SEXP linstat, 
    SEXP tol, 
    SEXP lower, 
    SEXP give_log
) {

    SEXP ans, stat, pval;
    double *MPinv, *ls, st, *ex;
    int rank, P, Q, PQ, B, greater = 0;
    
    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;

    if (C_get_varonly(LEV))
        error("cannot compute quadratic form based on variances only");

    MPinv = C_get_MPinv(LEV);
    C_MPinv_sym(C_get_Covariance(LEV), PQ, REAL(tol)[0], MPinv, &rank);
        
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    
    REAL(stat)[0] = C_quadform(PQ, C_get_LinearStatistic(LEV),
                               C_get_Expectation(LEV), MPinv);
    if (LENGTH(linstat) == 0) {
        REAL(pval)[0] = C_chisq_pvalue(REAL(stat)[0], rank, INTEGER(lower)[0],
                                       INTEGER(give_log)[0]);
    } else {
        B = NCOL(linstat);
        ls = REAL(linstat);
        st = REAL(stat)[0];
        ex = C_get_Expectation(LEV);
        greater = 0;
        for (int i = 0; i < B; i++) {
            if (GE(C_quadform(PQ, ls + PQ * i, ex, MPinv), st, REAL(tol)[0]))
                greater++;
        }
        REAL(pval)[0] = C_perm_pvalue(greater, B, INTEGER(lower)[0], INTEGER(give_log)[0]);
    }

    UNPROTECT(1);
    return(ans);
}

SEXP R_MaxtypeTest
(
    SEXP LEV, 
    SEXP linstat, 
    SEXP tol, 
    SEXP alternative, 
    SEXP lower, 
    SEXP give_log, 
    SEXP maxpts, 
    SEXP releps, 
    SEXP abseps
) {

    SEXP ans, stat, pval;
    double st, *ex, *cv, *ls, tl;
    int P, Q, PQ, B, vo, alt, greater;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
            
    if (C_get_varonly(LEV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));

    REAL(stat)[0] =  C_maxtype(PQ, C_get_LinearStatistic(LEV), 
                                   C_get_Expectation(LEV), 
                                   C_get_Covariance(LEV), 
                                   C_get_varonly(LEV),
                                   REAL(tol)[0],
                                   INTEGER(alternative)[0]);

    if (LENGTH(linstat) == 0) {
        if (C_get_varonly(LEV) && PQ > 1) {
            REAL(pval)[0] = NA_REAL;
            UNPROTECT(1);
            return(ans);
        }
        REAL(pval)[0] = C_maxtype_pvalue(REAL(stat)[0], C_get_Covariance(LEV),
                                         PQ, INTEGER(alternative)[0], INTEGER(lower)[0],
                                         INTEGER(give_log)[0],
                                         INTEGER(maxpts)[0], REAL(releps)[0], 
                                         REAL(abseps)[0], REAL(tol)[0]);
    } else {
        B = NCOL(linstat);
        st = REAL(stat)[0];
        ls = REAL(linstat);
        ex = C_get_Expectation(LEV);
        cv = C_get_Covariance(LEV);
        vo = C_get_varonly(LEV);
        alt = INTEGER(alternative)[0];
        tl = REAL(tol)[0];
        greater = 0;
        for (int i = 0; i < B; i++) {
            if (alt == ALTERNATIVE_less) {
                if (LE(C_maxtype(PQ, ls + PQ * i, ex, cv, vo, tl, alt), st, tl))
                    greater++;
            } else {
                if (GE(C_maxtype(PQ, ls + PQ * i, ex, cv, vo, tl, alt), st, tl))
                    greater++;
            }
        }
        REAL(pval)[0] = C_perm_pvalue(greater, B, INTEGER(lower)[0], INTEGER(give_log)[0]);
    }
    UNPROTECT(1);
    return(ans);
}

SEXP R_MaxSelectTest
(
    SEXP LEV, 
    SEXP ordered,
    SEXP linstat, 
    SEXP teststat, 
    SEXP tol, 
    SEXP minbucket, 
    SEXP lower, 
    SEXP give_log
) {

    SEXP ans, index, stat, pval;
    int P, Q, B, mb;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    B = NCOL(linstat);
    mb = INTEGER(minbucket)[0];

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));

    if (INTEGER(ordered)[0]) {
        SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, 1));
        if (C_get_Lb(LEV) == 1) {
            C_ordered_Xfactor(C_get_LinearStatistic(LEV),
                              C_get_Expectation(LEV),
                              C_get_VarianceInfluence(LEV),
                              C_get_CovarianceInfluence(LEV),
                              P, Q, 
                              C_get_ExpectationX(LEV),
                              B,
                              REAL(linstat),
                              mb,
                              REAL(tol)[0], 
                              INTEGER(teststat)[0],
                              INTEGER(index), REAL(stat),
                              REAL(pval), INTEGER(lower)[0], 
                              INTEGER(give_log)[0]);
        } else {
            C_ordered_Xfactor_block(C_get_LinearStatistic(LEV),
                              C_get_Expectation(LEV),
                              C_get_Covariance(LEV),
                              P, Q, 
                              C_get_ExpectationX(LEV),
                              B,
                              REAL(linstat),
                              mb,
                              REAL(tol)[0], 
                              INTEGER(teststat)[0],
                              INTEGER(index), REAL(stat),
                              REAL(pval), INTEGER(lower)[0], 
                              INTEGER(give_log)[0]);
        }
        if (REAL(stat)[0] > 0)
            INTEGER(index)[0]++; /* R style indexing */
    } else {
        SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, P));
        if (C_get_Lb(LEV) == 1) {
            C_unordered_Xfactor(C_get_LinearStatistic(LEV),
                              C_get_Expectation(LEV),
                              C_get_VarianceInfluence(LEV),
                              C_get_CovarianceInfluence(LEV),
                              P, Q, 
                              C_get_ExpectationX(LEV),
                              B,
                              REAL(linstat),
                              mb,
                              REAL(tol)[0], 
                              INTEGER(teststat)[0],
                              INTEGER(index), REAL(stat),
                              REAL(pval), INTEGER(lower)[0], 
                              INTEGER(give_log)[0]);
        } else {
            C_unordered_Xfactor_block(C_get_LinearStatistic(LEV),
                              C_get_Expectation(LEV),
                              C_get_Covariance(LEV),
                              P, Q, 
                              C_get_ExpectationX(LEV),
                              B,
                              REAL(linstat),
                              mb,
                              REAL(tol)[0], 
                              INTEGER(teststat)[0],
                              INTEGER(index), REAL(stat),
                              REAL(pval), INTEGER(lower)[0], 
                              INTEGER(give_log)[0]);
        }
    }

    UNPROTECT(1);
    return(ans);
}                                      
