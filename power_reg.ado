*! version 0.2 3Jan2017 Mauricio Caceres, caceres@nber.org
*! Power based on regression specification

capture program drop power_reg
program power_reg, rclass
    version 13.1
	syntax varlist(numeric ts fv) /// dependent_var covariates
           [if] [in] ,            /// subset
	[                             ///
		cluster(varname)          /// Grouping variable
        rho(real 0)               /// ICC [0 = compute it from data]
        strata(varlist)           /// Stratify by varlist
                                  /// - Continuous: specify # of quantiles
                                  /// - Categorical/dummy: specify '0' to use
                                  ///   the variable's categories.
        nstrata(numlist)          /// Number of strata for each stratum
                                  ///
        Ptreat(real 0.5)          /// Proportion treated
        alpha(real 0.05)          /// Confidence level
        kappa(real 0.8)           /// Power level
                                  ///
        binary                    /// Binary variable adjustment
        direction(str)            /// Direction of MDE for binary variable
        icccontrol(str)           /// Variable control for icc
        pcteffect(numlist)        /// % effect for sample size
        abseffect(numlist)        /// Absolute effect for sample size
                                  ///
        usestata                  /// Use Stata's built-in power command
	]

    * Parse varlist and sample to use
    * -------------------------------

    tempvar notouse
    gettoken depvar controls: varlist
    marksample touse
    markout `touse' `strata', strok
    markout `touse' `cluster', strok
    markout `touse' `varlist'
	_rmcoll `controls' if `touse', expand
    local controls `r(varlist)'
    gen byte `notouse' = !`touse'

    * Check options are sane
    * ----------------------

    if ("`usestata'" != "") {
        if ("`cluster'" != "") {
            di as err "cluster adjustment using Stata not implemented"
            exit
        }
        if ("`controls'" != "") & ("`binary'" != "") {
            di as err "{p}covariate adjustment for a binary variable" ///
                      " using Stata is not implemented{p_end}"
            exit
        }
    }

    * Check strata options are sane
    * -----------------------------

    if ("`strata'" != "") {
        qui ds `strata'
        local strata `r(varlist)'

        * if !`:list strata in controls' {
        *     di as err "{p}stratifying variable(s) '`strata'' " ///
        *               "should be in the controls{p_end}"
        *     exit 198
        * }

        if `:list sizeof strata' != `:list sizeof nstrata' {
            di as err "{p}specify # of strata for each " ///
                      "stratifying variable{p_end}"
            exit
        }

        local nstrata2 ""
        local i = 0
        foreach ns of local nstrata {
            local ++i
            local sv: word `i' of `strata'
            if (`ns' < 2) & (`ns' != 0) {
                di as err "{p}asked for `ns' strata for `sv'" ///
                          "; specify n > 1 or 0{p_end}"
                exit
            }
            if (`ns' == 0) {
                if ("`cluster'" != "") {
                    di "{p}{it:warning:} `sv' should be at the cluster" ///
                       " level or there may be unexpected results. {p_end}"
                }
                qui duplicates report `sv'
                local nstrata2 `nstrata2' `r(unique_value)'
            }
            else local nstrata2 `nstrata2' `ns'
        }
        local st  = `:di subinstr("`nstrata2'", " ", "*", .)'
    }
    else local st = 1

    * Set up strata
    * -------------

    if ("`strata'" != "") {
        if ("`cluster'" != "") {
            tempvar groups
            preserve
                keep `cluster' `touse' `strata'
                qui keep if `touse'
                qui collapse (sum) `touse' (mean) `strata', by(`cluster')
                local n2 = floor(_N / 2)

                if (`st' > `n2') {
                    di as err "{p}asked for `st' strata with `=_N'" ///
                              " clusters; should be <= `n2'{p_end}"
                    restore
                    exit
                }

                tempvar groups group
                qui gen `group'  = 1
                qui gen `groups' = 1
                qui forvalues i = 1 / `:list sizeof strata' {
                    local ns:  word `i' of `nstrata'
                    local ns2: word `i' of `nstrata2'
                    local sv:  word `i' of `strata'
                    sort `groups' `sv'
                    if (`ns' == 0) {
                        by `groups' `sv': replace `group' = (_n == 1)
                        by `groups': replace `group' = sum(`group')
                        local ns = `ns2'
                    }
                    else {
                        by `groups': replace `group' = ceil(`ns' * (_n / _N))
                    }
                    replace `groups' = (`groups' - 1) * `ns' + `group'
                }

                mata: get_info(st_local("groups"), st_local("touse"), 1, `st')
                sort `cluster' `groups'
                keep `cluster' `groups'
                tempfile stratafile
                qui save `stratafile'
            restore
            qui merge m:1 `cluster' using `stratafile'
            qui assert (_merge == 3) | `notouse'
            drop _merge
        }
        else {
            qui count if `touse'
            local N = `r(N)'
            local n2 = floor(`N' / 2)

            if (`st' > `n2') { // Stratified individual-level randomization
                di as err "{p}asked for `st' strata with `N' obs; " ///
                          "should be <= `n2'{p_end}"
                exit
            }

            tempvar groups group
            qui gen `group'  = 1 if `touse'
            qui gen `groups' = 1 if `touse'
            qui forvalues i = 1 / `:list sizeof strata' {
                local ns:  word `i' of `nstrata'
                local ns2: word `i' of `nstrata2'
                local sv:  word `i' of `strata'
                sort `notouse' `groups' `sv'
                if (`ns' == 0) {
                    by `notouse' `groups' `sv': replace `group' = (_n == 1)
                    by `notouse' `groups': replace `group' = sum(`group')
                    local ns = `ns2'
                }
                else {
                    by `notouse' `groups': replace `group' = ceil(`ns' * (_n / _N))
                }
                replace `groups' = (`groups' - 1) * `ns' + `group'
            }
            mata: get_info(st_local("groups"), st_local("touse"), 0, `st')
        }
        local controls `controls' i.`groups'
    }

    * Get stats based on depvar
    * -------------------------

    qui sum `depvar' if `touse'
    local var  = `r(Var)'
    local mean = `r(mean)'
    local N    = `r(N)'

    * Get regression adjustment
    * -------------------------

    qui xi: reg `depvar' `controls' if `touse'
    * local adjust = 1 - e(r2_a)
    local adjust = 1 - e(r2)
    local rmse   = e(rmse)

    * Adjust for clusters. We estimate ICC using ANOVA.
    * -------------------------------------------------

    if ("`cluster'" != "") {
        local vif_opts rho(`rho') cluster(`cluster') control(`icccontrol')
        get_vif `depvar' if `touse', `vif_opts'
        local vif = `r(vif)' * `adjust'
        if (`rho' == 0) local rho = `r(rho)'
        return scalar J   = `r(J)'
        return scalar mn  = `r(mn)'
        return scalar cv  = `r(cv)'
        return scalar sn  = `r(sn)'

        local J   = `r(J)'
        local mn  = `r(mn)'
        local cv  = `r(cv)'
        local sn  = `r(sn)'
        local rowname = "clusters"
    }
    else {
        local vif = 1 * `adjust'
        local rowname = "obs"
    }


    if ("`usestata'" != "") {
        * Compute MDE
             if ("`direction'" == "neg") local direction direction(lower)
        else if ("`direction'" == "pos") local direction direction(upper)
        else if ("`direction'" == "")    local direction ""
        local nratio = `ptreat' / (1 - `ptreat')
        if ("`binary'" != "") {
            local opts alpha(`alpha') power(`kappa') nratio(`nratio')
            local power power twoproportions
        }
        else {
            local opts alpha(`alpha') power(`kappa') nratio(`nratio') sd(`rmse')
            local power power twomeans
        }
        qui `power' `mean', n(`N') `opts' `direction'
        local mde = `r(diff)'

        * If asked, compute sampsi for % deviations of the mean
        if ("`pcteffect'" != "") {
            local pctnames ""
            tempname sampsi_pct
            foreach delta of local pcteffect {
                local dmde = `mean' + `delta' * `mean'
                qui `power' `mean' `dmde', `opts'
                matrix `sampsi_pct' = nullmat(`sampsi_pct'), `:di %15.0f r(N)'
                local pctnames `pctnames' `:di %9.1fc 100 * `delta''%
            }
            matrix colnames `sampsi_pct' = `pctnames'
            matrix rownames `sampsi_pct' = `rowname'
            return matrix sampsi_pct = `sampsi_pct'
        }

        * If asked, compute sampsi for absolute deviations of the mean
        if ("`abseffect'" != "") {
            local absnames ""
            tempname sampsi_abs
            foreach delta of local pcteffect {
                local dmde = `mean' + `delta'
                qui `power' `mean' `dmde', `opts'
                matrix `sampsi_abs' = nullmat(`sampsi_abs'), `:di %15.0f r(N)'
                local absnames `absnames' `:di %9.3fc `delta''%
            }
            matrix colnames `sampsi_abs' = `absnames'
            matrix rownames `sampsi_abs' = `rowname'
            return matrix sampsi_abs = `sampsi_abs'
        }
    }
    else {
        * Compute MDE
        local opts alpha(`alpha') kappa(`kappa') p(`ptreat') `binary'
        get_mde `var' `N' `mean', `opts' vif(`vif') direction(`direction')
        if ("`direction'" != "") {
                 if ("`direction'" == "pos") local mde = abs(`r(mde)')
            else if ("`direction'" == "neg") local mde = -abs(`r(mde)')
        }
        else local mde = `r(mde)'

        * If asked, compute sampsi for % deviations of the mean
        if ("`pcteffect'" != "") {
            local pctnames ""
            tempname sampsi_pct
            foreach delta of local pcteffect {
                local dmde = `mean' + `delta' * `mean'
                get_sampsi `var' `dmde' `mean', `opts' vif(`vif')
                matrix `sampsi_pct' = nullmat(`sampsi_pct'), `:di %15.0f r(sampsi)'
                local pctnames `pctnames' `:di %9.1fc 100 * `delta''%
            }
            matrix colnames `sampsi_pct' = `pctnames'
            matrix rownames `sampsi_pct' = `rowname'
        }

        * If asked, compute sampsi for absolute deviations of the mean
        if ("`abseffect'" != "") {
            local absnames ""
            tempname sampsi_abs
            foreach delta of local pcteffect {
                local dmde = `mean' + `delta'
                get_sampsi `var' `dmde' `mean', `opts' vif(`vif')
                matrix `sampsi_abs' = nullmat(`sampsi_abs'), `:di %15.0f r(sampsi)'
                local absnames `absnames' `:di %9.3fc `delta''%
            }
            matrix colnames `sampsi_abs' = `absnames'
            matrix rownames `sampsi_abs' = `rowname'
        }
    }

    * Pretty printing
    * ---------------

    if ("`binary'" == "") {
        local sd  "sd = sd1 = sd2"
        local bin ""
    }
    else {
        local sd  "var = m1 * (1 - m1) * P + m2 * (1 - m2) * (1 - P)"
        local bin " (binary outcome)"
    }

    if ("`controls'" == "") {
        local power "Power calculations for two-sample means test."
        local null  "Ho: m2 = m1  versus  Ha: m2 != m1"
    }
    else {
        local power "Power calculations for linear regression`bin'."
        local null  "Y = a + b T + X + e; Ho: b = 0 versus Ha: b != 0"
    }

    di ""
    di "{p}`power'{p_end}"
    di "{p}`null'{p_end}"
    di "{p}`sd'{p_end}"
    di ""
    di "Study parameters:"
    di "        P = `:di %9.4f `ptreat''"
    if ("`cluster'" == "") {
        di "       PN = `:di %9.0gc round(`ptreat' * `N')'"
    }
    di "    alpha = `:di %9.4f `alpha''"
    di "    power = `:di %9.4f `kappa''"
    di "       m1 = `:di %9.4f `mean''"
    if ("`cluster'" != "") {
        di ""
        di "Cluster parameters:"
        di "      ICC = `:di %9.4f `rho''"
        di "    nclus = `:di %9.0gc `J''"
        di "  treated = `:di %9.0gc round(`ptreat' * `J')'"
        di "      m_n = `:di %9.4fc `mn''"
        di "       cv = `:di %9.4f `cv''"
    }
    di ""
    di "Results:"
    di "      MDE = `:di %9.4f `mde''"
    di " MDE / m1 = `:di %9.4f `mde' / `mean''"
    if ("`pcteffect'" != "") | ("`abseffect'" != "")  {
        di ""
        di "Sample size:"
        if ("`pcteffect'" != "") {
            matrix list `sampsi_pct', nob noh f(%15.0fc)
        }
        if ("`abseffect'" != "") {
            matrix list `sampsi_abs', nob noh f(%15.0fc)
        }
    }

    * Report
    * ------

    qui count if `touse'
    return scalar N          = `r(N)'
    return scalar mean       = `mean'
    return scalar std        = `:di sqrt(`var')'
    return scalar rmse       = `rmse'
    return scalar rho        = `rho'
    return scalar adjust     = `adjust'
    return scalar vif        = `:di `vif' / `adjust''
    return scalar mde        = `mde'
    return scalar mde_pct    = `:di `mde' / `mean''
    if ("`sampsi_pct'" != "") {
        return matrix sampsi_pct = `sampsi_pct'
    }
    if ("`sampsi_abs'" != "") {
        return matrix sampsi_abs = `sampsi_abs'
    }
end

* Get variance inflation factor (VIF)
* -----------------------------------

capture program drop get_vif
program get_vif, rclass sortpreserve
    syntax varname [if] [in], cluster(str) [rho(real 0) control(str)]

    marksample touse
    markout `touse' `cluster', strok
    tempvar nclus nj vif

    * Compute mean and variance of cluster sizes
    get_rho `varlist' if `touse', cluster(`cluster') control(`control')
    if (`rho' == 0) {
        if (`r(rho)' > 0) local rho = `r(rho)'
        else {
            local rho = 0
            di "{p}{cmd:ICC (rho) truncated at 0.}{p_end}"
        }
    }
    local mn  = `r(mn)'
    local cv  = `r(cv)'
    local sn  = `r(sn)'
    local J   = `r(J)'

    * Compute VIF and adjust N, MDE accordingly
    qui {
        egen `nclus' = tag(`cluster') if `touse'
        bys `cluster' `touse': gen `nj' = _N  if `nclus'
        gen `vif' = `mn' * `J' * `nj' * (1 + (`nj' - 1) * `rho') / `sn'^2
        sum `vif'
        local vif = `r(sum)'
    }

    return scalar vif = `vif'
    return scalar rho = `rho'
    return scalar J   = `J'
    return scalar mn  = `mn'
    return scalar cv  = `cv'
    return scalar sn  = `sn'
end

* Get intra-cluster correlation (ICC, rho)
* ----------------------------------------

capture program drop get_rho
program get_rho, rclass sortpreserve
    syntax varname [if] [in], cluster(str) [control(str)]
    marksample touse
    markout `touse' `cluster', strok

    * Some cluster stats
    tempvar nclus nj nj2
    qui {
        egen `nclus' = tag(`cluster') if `touse'
        bys `cluster' `touse': gen `nj' = _N if `nclus'
        gen `nj2'    = `nj'^2
    }

    qui count if `nclus'
    local J = `r(N)'

    qui sum `nj'
    local mn = `r(mean)'
    local cv = `r(sd)' / `r(mean)'
    local sn = `r(sum)'

    qui sum `nj2'
    local sn2 = `r(sum)'

    * Compute ANOVA, adjusting for a covariate if requested. See:
    * Stanish, W. M. and Taylor, N. (1983). Estimation of the Intraclass
    * Correlation Coefficient for the Analysis of Covariance Model. _The
    * American Statistician_, 37(3):221–224.
    qui if ("`control'" != "") {
        tempvar xmean xerr xerrc
        qui sum `control' if `touse'
        local mean   = `r(mean)'
        egen `xmean' = mean(`varlist') if `touse', by(`cluster')
        gen  `xerr'  = `nj2' * (`mean'  - `xmean')^2 if `touse'
        gen  `xerrc' = (`mean' - `control')^2        if `touse'

        qui sum `xerr'
        local xssb = r(sum)
        qui sum `xerrc'
        local xssw = r(sum)

        local k = 1 + (`xssb' / `xssw') / (`J' - 1)
    }
    else {
        local k = 1
    }

    qui anova `varlist' `cluster' `control' if `touse'
    local msw = e(rss) / e(df_r)
    local msb = e(ss_1) / e(df_1)
    local n0  = (e(N)  - `sn2' / e(N)) / e(df_1)
    local rho = (`msb' - `msw') / (`msb' + (`n0' - `k') * `msw')

    return scalar J   = `J'
    return scalar sn  = `sn'
    return scalar mn  = `mn'
    return scalar cv  = `cv'
    return scalar rho = `rho'
end

* Get minimum detectable effect (MDE)
* -----------------------------------

cap program drop get_mde
program get_mde, rclass
    syntax anything, [   /// var mean_treat mean_control
        Ptreat(real 0.5) /// Proportion treated
        alpha(real 0.05) /// Confidence level
        kappa(real 0.8)  /// Power level
        rho(real 0)      /// ICC
        cv(real 0)       /// Coefficient of variation, sigma_n / mu_n
        mn(real 1)       /// mu_n
        vif(real 1)      /// Variance Inflation Factor (or design effect)
        cluster          /// Cluster adjustment using rho, cv, mn
        binary           /// Whether to treat the outcome as binary
        unequalclus      /// Unequal-cluster adjustment
        UPPERbound       /// Assume variance is 0.25 for binary variable
        direction(str)   /// How to adjust MDE
    ]

    tokenize `anything'
    local var  = `1'
    local N    = `2'
    local mean = `3'

    /*
     * # Formulas
     *
     * For a continuous outcome, power is given by
     *
     *     MDE = (z_α + z_(1 - κ)) * √σy / (N * P * (1 - P))
     *
     * where z are the quantiles of the normal distribution corresponding
     * to α, 1 - κ the probability of Type I and Type II error,
     * respectively. If Y is binary, σy becomes
     *
     *     σy = μy * (1 - μy) * P + μt * (1 - μt) * (1 - P)
     *     μt = μy + MDE
     *
     * MDE then solves the quadratic equation that is defined by replacing
     * σy in the first squation. Note that in this case MDE is asymetric,
     * as it has an uneven impact on the variance determined by whether MDE
     * is positive or negative. `direction` specifies the direction of the
     * effect. Last, with clustering
     *
     *     vif     = 1 + ρ * ((cv^2 + 1) * μn - 1)
     *     MDEclus = MDE * √vif
     *
     * Where cv = 0 if eqsize = true. The user can specify vif and μn with
     * cluster = false and the same adjustment will be applied.
     *
     * # Sources
     * - Kong, S.-H., Ahn, C. W., and Jung, S.-H. (2003). Sample Size
     *   Calculation for Dichotomous Outcomes in Cluster Randomization
     *   Trials with Varying Cluster Size. _Drug Information Journal_,
     *   37(1):109–114.
     * - Manatunga, A. K., Hudgens, M. G., and Chen, S. (2001). Sample Size
     *   Estimation in Cluster Randomized Studies with Varying Cluster
     *   Size. _Biometrical Journal_, 43(1):75–86.
     */

     local t = abs(invnormal(`alpha' / 2) + invnormal(1 - `kappa'))

    * Adjust for clusters
    if ("`cluster'" != "") {
        if ("`unequalclus'" == "") local vif = (1 + `rho' * (`mn' - 1))
        else local vif = (1 + `rho' * ((`cv'^2 + 1) * `mn' - 1))
    }

    if ("`binary'" != "") & ("`upperbound'" != "") {
        local var = 0.25
    }

    * Adjust for binary outcome (NOTE: I'm not sure this is standard; so
    * while the algebra works out fine, it makes me uneasy to use).
    if ("`binary'" != "") & ("`upperbound'" == "") {
        local G   = `vif' * (`t')^2 / (`ptreat' * (1 - `ptreat') * `N')
        local a   = 1 + `G' * (1 - `ptreat')
        local b   = `G' * (2 * `mean' - 1) * (1 - `ptreat')
        local c   = - `G' * `mean' * (1 - `mean')
        local t1  = - `b' / (2 * `a')
        local t2  = sqrt(`b'^2 - 4 * `a' * `c') / (2 * `a')

        if ("`direction'" == "neg")  local MDE = `t1' - `t2'
        else local MDE = `t1' + `t2'
    }
    else {
        local MDE = sqrt(`vif') * `t' * sqrt(`var') / sqrt(`ptreat' * (1 - `ptreat') * `N')
    }

    return scalar mde = `MDE'
end

* Get sample size (N)
* -------------------

cap program drop get_sampsi
program get_sampsi, rclass
    syntax anything, [   /// var mean_treat mean_control
        Ptreat(real 0.5) /// Proportion treated
        alpha(real 0.05) /// Confidence level
        kappa(real 0.8)  /// Power level
        rho(real 0)      /// ICC
        cv(real 0)       /// Proportion treated
        mn(real 1)       /// Proportion treated
        vif(real 1)      /// Variance Inflation Factor
        cluster          /// Cluster adjustment using rho, cv, mn
        binary           /// Whether to treat the outcome as binary
        unequalclus      /// Unequal-cluster adjustment
    ]

    tokenize `anything'
    local var = `1'
    local mt  = `2'
    local mc  = `3'

    /*
     * # Formulas
     *
     * For a continuous outcome, power is given by
     *
     *     N = σy * ((z_α + z_(1 - κ)) / (μt - μc))^2 / (P * (1 - P))
     *
     * where z are the quantiles of the normal distribution corresponding
     * to α, 1 - κ the probability of Type I and Type II error,
     * respectively. If Y is binary, σy becomes
     *
     *     σy = μc * (1 - μc) * P + μt * (1 - μt) * (1 - P)
     *
     * If there is clustering, then
     *
     *     vif   = 1 + ρ * ((cv^2 + 1) * μn - 1)
     *     Nclus = N * vif / μn
     *
     * Where cv = 0 if eqsize = true. The user can specify vif and μn with
     * cluster = false and the same adjustment will be applied.
     *
     * # Sources
     * - Kong, S.-H., Ahn, C. W., and Jung, S.-H. (2003). Sample Size
     *   Calculation for Dichotomous Outcomes in Cluster Randomization
     *   Trials with Varying Cluster Size. _Drug Information Journal_,
     *   37(1):109–114.
     * - Manatunga, A. K., Hudgens, M. G., and Chen, S. (2001). Sample Size
     *   Estimation in Cluster Randomized Studies with Varying Cluster
     *   Size. _Biometrical Journal_, 43(1):75–86.
     */

     local t = invnormal(`alpha' / 2) + invnormal(1 - `kappa')

    * Adjust for proportions model
    if ("`binary'" != "") {
        local var = (`mc' * (1 - `mc') * `ptreat') + (`mt' * (1 - `mt') * (1 - `ptreat'))
    }

    * Adjust for clusters
    if ("`cluster'" != "") {
        if ("`unequalclus'" == "") local vif = (1 + `rho' * (`mn' - 1))
        else local vif = (1 + `rho' * ((`cv'^2 + 1) * `mn' - 1))
    }

    local N = (`vif' / `mn') * `var' * (`t' / (`mt' - `mc'))^2 / (`ptreat' * (1 - `ptreat'))
    return scalar sampsi = `N'
end

***********************************************************************
*                            Mata helpers                             *
***********************************************************************

capture mata: mata drop get_info()
capture mata: mata drop panel_info()
mata:
    void function get_info(string scalar groups,
                           string scalar touse,
                           real scalar cluster,
                           real scalar nstrata)
    {

        if (cluster == 1) {
            nj = st_data(., touse)
            panel_info(rows(nj), min(nj), max(nj))
        }

        info  = panelsetup(st_data(., groups, touse), 1)
        panel = panelstats(info)
        panel_info(panel[1], panel[3], panel[4], "strata")
    }

    // Print number of clusters/strata and number of obs per cluster/strata
    void function panel_info(real scalar N,
                             real scalar pmin,
                             real scalar pmax,
                             | string scalar what)
    {
        if (args() == 3) {
            what = "panel"
        }
        if (pmin == pmax) {
            nstr = strtrim(sprintf("%21.0gc", N))
            jstr = strtrim(sprintf("%21.0gc", pmin))
            printf("Balanced " + what + ". J = " + nstr + ", n_j = " + jstr + "\n")
        }
        else {
            nstr = strtrim(sprintf("%21.0gc", N))
            jstr = strtrim(sprintf("%21.0gc", pmin))
            jstr = jstr + " to " + strtrim(sprintf("%21.0gc", pmax))
            printf("Unbalanced " + what + ". J = " + nstr + ", n_j = " + jstr + "\n")
        }
    }
end
