

### Class for raw data: a set of `x' variables and a set of `y' variables,
### possibly blocked and with weights
setClass(Class = "IndependenceProblem",
    representation = representation(
        x       = "data.frame",
        y       = "data.frame",
        block   = "factor",
	weights = "numeric"
    ),
    validity = function(object) {
        dims <- ((nrow(object@x) == nrow(object@y)) && 
                 (nrow(object@x) == length(object@block))) 
        dims <- dims && (length(object@block) == length(object@weights))
        Wint <- max(abs(object@weights - floor(object@weights))) < eps()
        block <- all(table(object@block) > 1)
        NAs <- all(complete.cases(object@x) & complete.cases(object@y))
        return((dims && block) && (NAs && Wint))
    }
)


### Class for transformed data, the `x' variables are transformed 
### to a (n x p) matrix `xtrans' and the `y' variables to `ytrans' (n x q).
### `scores' is a matrix of scores
setClass(Class = "IndependenceTestProblem",
    representation = representation(
        xtrans     = "matrix",
        ytrans     = "matrix",
        xtrafo     = "function",
        ytrafo     = "function"
    ),
    contains = "IndependenceProblem",
    validity = function(object) 
        (storage.mode(object@xtrans) == "double" && 
         storage.mode(object@ytrans) == "double")
)

### Covariance matrix
setClass(Class = "CovarianceMatrix",
    representation = representation(
        covariance = "matrix"
    )
)

### Variances only
setClass(Class = "Variance",
    representation = representation(
        variance = "numeric"
    )
)

setClassUnion("VarCovar", c("CovarianceMatrix", "Variance"))

### Linear statistic, expectation and covariance according to 
### Strasser & Weber (1999)
setClass(Class = "IndependenceLinearStatistic",
    representation = representation(
        linearstatistic             = "numeric",
        expectation                 = "numeric",
        covariance                  = "VarCovar"
    ),
    contains = "IndependenceTestProblem",
)

### Tests based on linear statistics
setClass(Class = "IndependenceTestStatistic",
    representation = representation(
        teststatistic               = "numeric",
        standardizedlinearstatistic = "numeric"

    ),
    contains = "IndependenceLinearStatistic",
)

### teststatistic = standardizedlinearstatistic
setClass(Class = "ScalarIndependenceTestStatistic",
    representation = representation(
        alternative   = "character"
    ),
    contains = "IndependenceTestStatistic",
    validity = function(object) 
        object@alternative %in% c("two.sided", "less", "greater")
)

### teststatistic = max(abs(standardizedlinearstatistic))
setClass(Class = "MaxTypeIndependenceTestStatistic",
    representation = representation(
        alternative                 = "character"
    ),
    contains = "IndependenceTestStatistic",
    validity = function(object)
        object@alternative %in% c("two.sided", "less", "greater")
)

### teststatistic = quadform(linearstatistic)
setClass(Class = "QuadTypeIndependenceTestStatistic",
    representation = representation(
        covarianceplus              = "matrix",
        df                          = "numeric"
    ),
    contains = "IndependenceTestStatistic"
)

### p values
setClass(Class = "PValue",
    representation = representation(
        pvalue = "function",
        p      = "function",
        name   = "character"
    )
)

### Null distribution
setClass(Class = "NullDistribution",
    representation = representation(
        q          = "function",
        d          = "function",
        support    = "function",
        parameters = "list"
    ),
    contains = "PValue"
)

### There are essentially three types of computing null distributions:
setClass(Class = "AsymptNullDistribution", contains = "NullDistribution")
setClass(Class = "ApproxNullDistribution", contains = "NullDistribution")
setClass(Class = "ExactNullDistribution", contains = "NullDistribution")

### the "fitted" test including data and everything
setClass(Class = "IndependenceTest",
    representation = representation(
        distribution = "PValue", ### was: "NullDistribution",
        statistic    = "IndependenceTestStatistic",
        estimates    = "list",
        method       = "character"
    ),
    prototype = list(method = "General Independence Test")
)

### the "fitted" test for scalar linear statistics
setClass(Class = "ScalarIndependenceTest",
    representation = representation(
        nullvalue    = "numeric"
    ),
    contains = "IndependenceTest",
    validity = function(object)
        extends(class(object@statistic), 
                "ScalarIndependenceTestStatistic")
)

### possibly with confidence intervals
setClass(Class = "ScalarIndependenceTestConfint",
    representation = representation(
        confint    = "function",
        conf.level = "numeric"
    ),
    contains = "ScalarIndependenceTest"
)

### max type test statistics
setClass(Class = "MaxTypeIndependenceTest",
    contains = "IndependenceTest",
    validity = function(object)
        extends(class(object@statistic), 
                "MaxTypeIndependenceTestStatistic")
)

### quad form test statistics
setClass(Class = "QuadTypeIndependenceTest",
    contains = "IndependenceTest",
    validity = function(object)
        extends(class(object@statistic), 
                "QuadTypeIndependenceTestStatistic")
)

### SymmetryProblems
setClass(Class = "SymmetryProblem",
    contains = "IndependenceProblem",
    validity = function(object) {
        if (ncol(object@x) != 1 || !is.factor(object@x[[1]]))
            stop(sQuote("x"), " slot does not contain a single factor")
        if (!is_completeblock(object))
            stop(sQuote("object"), 
                 " is not a an unreplicated complete block design")
        TRUE
    }
)
