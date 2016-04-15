
#include "libcoin.h"

double C_quadform(int PQ, double *linstat, double *expect, double *MPinv)
{
    int qPQ;
    double ans = 0.0, tmp = 0.0;
    
    for (int q = 0; q < PQ; q++) {
        qPQ = q * PQ;
        tmp = 0.0;
        for (int p = 0; p < PQ; p++)
            tmp += (linstat[p] - expect[p]) * MPinv[qPQ + p];
        ans += tmp * (linstat[q] - expect[q]);
    }
    return(ans);
}

double CR_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv)
{
    return(C_quadform(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                      REAL(Expectation), REAL(MPinv)));
}

double C_maxstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (covar[p*PQ + p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar[p*PQ + p]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_maxstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_minstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (covar[p*PQ + p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar[p*PQ + p]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

double C_minstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

double C_maxabsstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (covar[p*PQ + p] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(covar[p*PQ + p]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_maxabsstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(var[p]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double CR_maxstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol)
{
    if (isMatrix(CoVariance))
        return(C_maxstat_Covariance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                                    REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
    return(C_maxstat_Variance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                              REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
}

double CR_minstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol)
{
    if (isMatrix(CoVariance))
        return(C_minstat_Covariance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                                    REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
    return(C_minstat_Variance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                              REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
}

double CR_maxabsstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol)
{
    if (isMatrix(CoVariance))
        return(C_maxabsstat_Covariance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                                       REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
    return(C_maxabsstat_Variance(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                                 REAL(Expectation), REAL(CoVariance), REAL(tol)[0]));
}
