cap program drop reg_copy 
program reg_copy, eclass
version 13.0

syntax varlist(min = 2) [if] [in] [, vce(string)]
token `varlist'

tempvar select
mark `select' `if' `in'

if length("`if'") > 0 | length("`in'") > 0 {
    preserve
    qui keep if `select'
}
local N = _N

local regvar = 0
foreach v in `varlist' {

    qui sum `v'
    if r(N) == 0 {
        dis as error "nonnumeric found where numeric required"
        exit 3251
    }

    local regvar = `regvar' + 1
}

mata: mata clear
mata: st_view(W = ., ., "`varlist'")
mata: y =  W[.,1]
mata: X = W[., 2..cols(W)], J(rows(W), 1, 1)

mata: beta = invsym(cross(X,X)) * cross(X, y)
mata: sqerr = (y - X*beta):^2
mata: ssr = sum((X*beta - J(rows(W), 1, mean(y))):^2)
mata: sse = sum((y - X*beta):^2)
mata: sst = ssr + sse
mata: msr = ssr/(`regvar' - 1)
mata: mse = sse/(`N' - `regvar')
mata: mst = sst/(`N' - 1)
foreach v in sse ssr sst mse msr mst {
    mata: st_numscalar("`v'", `v')
}


mata: std_err = sqrt(diagonal((sum(sqerr)/(rows(W) - cols(W))) * invsym(cross(X,X))))
mata: t_stat = J(rows(std_err), 1, 1)
mata: p = J(rows(std_err), 1, 1)
mata: ci_lb = J(rows(std_err), 1, 1)
mata: ci_ub = J(rows(std_err), 1, 1)
forv i = 1/`regvar' {
    mata: t_stat[`i', 1] = beta[`i', 1]/std_err[`i', 1]
}

mata: res = beta, std_err, t_stat, p, ci_lb, ci_ub
mata: st_matrix("res", res)
forv j = 1/`regvar' {
    matrix res[`j', 4] = tprob(`N' - `regvar', res[`j', 3])
    matrix res[`j', 5] = res[`j', 1] + invt(`N' - `regvar', 0.025) * res[`j', 2]
    matrix res[`j', 6] = res[`j', 1] + invt(`N' - `regvar', 0.975) * res[`j', 2]
}

local df1 = `regvar' - 1
local df2 = `N' - `regvar'
local df = `N' - 1
local rsq = ssr/sst
local adj_rsq = 1 - ((1 -`rsq') * (`N' - 1) / ((`N' - 1) - (`regvar' - 1)))
local rmse = sqrt(sse/((`N' - 1) - (`regvar' - 1)))
local j = 1
local fstat = F(`df1', `df2', ssr/sse)
di as text ""
di as text %12s "  Source" " {c |}       SS           df       MS  " "     Number of obs   ="  as result %9.4g `N'
di as text "{hline 13}{c +}{hline 34}""    F(`df1',`df2')         =" as result %9.2f msr/mse
di as text %12s "  Model" " {c |}"   as result %12.1g ssr "         " as result `df1'  as result %12.2g msr  "    Prob > F        =" as result %9.4f 1 - F(`df1', `df2', msr/mse)
di as text %12s "  Residual" " {c |}" as result %12.1g sse "        " as result `df2' as result %12.2g mse  "    R-squared       ="  as result %9.4f `rsq'    
di as text "{hline 13}{c +}{hline 34}"  "    Adj R-squared   =" as result %9.4g `adj_rsq'
di as text %12s "  Total" " {c |}" as result %12.1g sst "        " as result `df' as result %12.2g mst "    Root MSE        ="    as result %9.1g `rmse'
di as text ""
di as text "{hline 13}{c TT}{hline 65}"
di as text %12s abbrev("`1'",12) " {c |} Coefficient   Std. Err.    t     P>|t|      [95% Conf. Interval]" 
di as text "{hline 13}{c +}{hline 65}"
foreach v in `varlist' {
    if "`v'" == "`1'" {
        continue
    }
    di as text %12s abbrev("`v'",12) " {c |} "  as result %9.5g res[`j', 1] "  " %9.5g res[`j', 2] %9.2f res[`j', 3] %9.3f res[`j', 4] "     " %9.3f res[`j', 5] "   " %9.3f res[`j', 6]
    local j = `j' + 1
}
di as text %12s abbrev("_cons",12) " {c |} "  as result %9.5g res[`regvar', 1] "  " %9.5g res[`regvar', 2] %9.2f res[`regvar', 3] %9.3f res[`regvar', 4] "     " %9.3f res[`regvar', 5] "   " %9.3f res[`regvar', 6]
di as text "{hline 13}{c BT}{hline 65}"
end