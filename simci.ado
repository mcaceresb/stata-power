*! version 0.5 19Oct016 Mauricio Caceres, caceres@nber.org
*! Simulated CI, MDE, and power for regression specification (accepts
*! clusters and arbitrary number of stratifying/blocking variables)

* For debugging
* -------------

capture program drop simci

capture mata: mata drop simci()
capture mata: mata drop simci_continuous()
capture mata: mata drop simci_binary()

capture mata: mata drop power_search()
capture mata: mata drop parse_simci()
capture mata: mata drop parse_power()

capture mata: mata drop simci_setup_basic()
capture mata: mata drop simci_setup_strata()
capture mata: mata drop simci_setup_cluster()
capture mata: mata drop simci_panel_info()

capture mata: mata drop shuffle_basic()
capture mata: mata drop shuffle_strata()
capture mata: mata drop shuffle_cluster()
capture mata: mata drop shuffle_strata_cluster()

capture mata: mata drop expand_rows()
capture mata: mata drop pctile()
capture mata: mata drop bc_var()

***********************************************************************
*                            Main Program                             *
***********************************************************************

* TODO: "autostratify"? // 2016-09-23 14:43 EDT
program simci, rclass sortpreserve
	syntax varlist(numeric ts fv) /// dependent_var covariates
           [if] [in] ,            /// subset
	[                             ///
        Ptreat(real 0.5)          /// Proportion treated
        alpha(real 0.05)          /// Confidence level
        reps(int 100)             /// Non-parametric repetitions
                                  ///
		cluster(varname)          /// Grouping variable
        forcestrata               /// Force stratification
        nstrata(numlist)          /// Number of strata for each var in varlist
        strata(varlist)           /// Stratify by varlist
                                  /// - Continuous: specify # of quantiles
                                  /// - Categorical/dummy: specify '0' to use
                                  ///   the variable's categories.
        effect(str)               /// Induce artificial effect
        power(str)                /// Search for power
	]
    local savelist `varlist'

    * Figure out what the function will output
	if ("`effect'" != "") & ("`power'" != "") {
		di as err "{p}options effect and power are mutually exclusive{p_end}"
		exit
	}
    else if ("`effect'" != "") & ("`power'"  == "") local compute effect
    else if ("`effect'" == "") & ("`power'"  != "") local compute power
    else if ("`effect'" == "") & ("`power'"  == "") local compute ci

    * If asked for an effect, parse
    if ("`compute'" == "effect") {
        local 0 , `effect'
        syntax, effect(real) [binary bounds(str)]
        if ("`bounds'" != "auto") {
            local 0 , bounds(`bounds')
            syntax, [bounds(numlist)]
        }
    }
    else local effect 0

    * If asked for power search, parse
    if ("`compute'" == "power") {
        local 0 , `power'
        syntax, DIRection(str) [kappa(real 0.8) binary tol(real 0) startat(real 0)]
    }

    local varlist `savelist'

    /*
     * # Notes
     *
     * The formulas and citations below only justify the logic for when
     * -compute- is set to 'ci'; inducing an effect and searching for
     * power is in beta. Use with caution.
     *
     * # Simulated CI
     *
     * The idea is to simulate a non-parametric CI based on placebo
     * assignments of a treatment variable. The program assigns
     * treatment at random, hence a null effect, to individuals or
     * clusters, optionally stratifying by any number of variables (or
     * the means thereof, in the case of clusters). Consider
     *
     *     Y_ij = a + b T_j + g X_ij + e_ij
     *
     * There are C = J choose PJ ways to treat the clusters (or C =
     * P choose PN in the case of individuals). If we computed b_ols
     * for c = 1, ..., C we would know the exact distribution of our
     * estimator, conditional on the data being representative of the
     * study data. C is typically intractably large, hence we simulate
     * K draws with sum(T_jk = PJ) and run
     *
     *     Y_ij = a + b_k T_jk + g X_ij + e_ij
     *
     * Let Fhat be the empirical cdf of b_k; a valid 1 - alpha CI for
     * b is given by
     *
     *     CI(1 - a) = [Fhat^-1(a / 2), Fhat^-1(1 - a / 2)]
     *
     * Duflo et al. (2007) recommends this for clusters as the
     * regression will naturally take into account the correlation
     * structure of the data, but the logic is easy to extend to
     * individual-level randomization and stratified randomization.
     *
     * # Simulated power
     *
     * The idea here is simple: Instead of inducing placebo assignments,
     * introduce an effect at each run of the randomization. Power is
     * then the average number of times the null is rejected (in this
     * case, we use the simulated CI from the previous section as our
     * rejection criterion). We conduct a search over various MDEs until
     * power is approximately the level we want.
     *
     * # Sources
     *
     * - Duflo, E., Glennerster, R., and Kremer, M. (2007). Chapter
     *   61 Using Randomization in Development Economics Research: A
     *   Toolkit. In Handbook of Development Economics, volume 4, pages
     *   3895–3962.
     */

    * Parse varlist and sample to use
    * -------------------------------

    tempvar notouse
    gettoken depvar controls: varlist
    marksample touse
    markout `touse' `strata', strok
    markout `touse' `cluster', strok
	_rmcoll `controls' if `touse', expand
    local controls `r(varlist)'
    gen byte `notouse' = !`touse'

    * Check options are sane
    * ----------------------

    * Effect opts
    if ("`compute'" == "effect") {
        if (`effect' == 0) & ("`binary'" != "") {
            di "{p}{it:warning:} ignoring option -binary- with a 0 effect{p_end}"
        }

        if !(inlist(`:list sizeof bounds', 0, 1, 2)) {
            di as err "-bounds- must be empty, 'auto', or an upper and lower bound"
            exit
        }
    }

    * Power opts
    if ("`compute'" == "power") {
        if !(inlist("`direction'"), "pos", "neg") {
            di as err "please specify `direction' = 'pos', 'neg'"
            exit
        }

        if ((`kappa' <= 0) | (`kappa' >= 1)) {
            di as err "can't search for power outside (0, 1)"
            exit
        }

        if (`tol' == 0) {
            local tol = min(1e-2, 10 / `reps')
        }

        local mintol = min(`kappa', 1 - `kappa')
        if (`tol' >= `mintol') {
            di as err "with power `kappa', tolerance should smaller than `mintol'"
            exit
        }

        local maxtol = 1 / `reps'
        if (`tol' < `maxtol') {
            di as err "with `reps' repetitions, tolerance should greater than `maxtol'"
            exit
        }
    }

    * Strata opts
    if ("`strata'" != "") {
        qui ds `strata'
        local strata `r(varlist)'

        * if !`:list strata in controls' {
        *     di as err "{p}stratifying variable '`strata''" ///
        *               " must all be in the controls{p_end}"
        *     exit 198
        * }

        if `:list sizeof strata' != `:list sizeof nstrata' {
            di as err "{p}specify # of strata for each stratifying variable{p_end}"
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

    * Set up randomization for Mata
    * -----------------------------

    * Clusters and strata require sorting! They need to appear
    * sequentially in the data and in the same order across simulations.
    * Hence we force clusters and strata to appear in some arbitrary
    * order and we subsequently fix that order for the simulation.
    *   - Clusters: Encode cluster if not numeric so order does not vary.
    *   - Strata: Generate strata (possibly by cluster) as numeric (1 to
    *             prod(`nstrata'))
    * Then simply create a Mata pointer to the correct function to
    * shuffle the treatment indicator.

    tempname shufflefun atreat nt
    if ("`cluster'" != "") {
        tempvar clusorder
        qui if regexm("`:type `cluster''", "str") {
            encode `cluster' if `touse', gen(`clusorder')
        }
        else local clusorder `cluster'

        if ("`strata'" == "") {
            * Cluster-level randomization: We shuffle a (0, 1) indicator with PJ
            * 1s and (1 - P)J 0s and then expand each element by the number of
            * observations in each cluster.
            sort `notouse' `clusorder'
            mata: `atreat' = simci_setup_cluster(st_local("clusorder"), ///
                                                 st_local("touse"), `ptreat')
            mata: st_local("nclus", strofreal(rows(asarray(`atreat', "2"))))
            mata: `nt' = asarray(`atreat', "0")
            mata: `shufflefun' = &shuffle_cluster()
        }
        else {
            * Stratified cluster-level randomization: We shuffle prod(`nstrata')
            * (0, 1) indicators with P(J / prod(`nstrata')) 1s and the rest 0s
            * at each strata. We stack the indicators and then expand each element
            * by the number of observations in each cluster.
            tempvar groups
            preserve
                keep `clusorder' `touse' `strata'
                qui keep if `touse'
                qui collapse (sum) `touse' (mean) `strata', by(`clusorder')
                local n2 = floor(_N / 2)

                if (`st' > `n2') {
                    di as err "{p}asked for `st' strata with `=_N' clusters; " ///
                              "should be <= `n2'{p_end}"
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

                qui duplicates report `groups'
                cap assert `r(unique_value)' == `st'
                if (_rc != 0) & ("`forcestrata'" == "") {
                    di as err "asked for `st' strata but computed"                ///
                              " `r(unique_value)'." _n(1) "check the stratifying" ///
                              " covariates or run again with -forcestrata-"
                    restore
                    exit
                }

                sort `groups' `clusorder'
                mata: `atreat' = simci_setup_strata(st_local("groups"), ///
                                                    st_local("touse"),  ///
                                                    1, `ptreat', `st')

                mata: `nt' = asarray(`atreat', "0")
                keep `clusorder' `groups'
                sort `clusorder' `groups'
                tempfile stratafile
                qui save `stratafile'
                local nclus = _N
            restore
            qui merge m:1 `clusorder' using `stratafile'
            qui assert (_merge == 3) | `notouse'
            drop _merge
            mata: `shufflefun' = &shuffle_strata_cluster()
        }
    }
    else {
        qui count if `touse'
        local N = `r(N)'
        qui if ("`strata'" == "") {
            * Individual-level randomization: We shuffle a (0, 1) indicator with
            * PN 1s and the rest 0s.
            scalar ntreat   = `:di %15.0f ceil(`ptreat' * `r(N)')'
            scalar ncontrol = `N' - scalar(ntreat)

            mata: `atreat' = simci_setup_basic(st_numscalar("ntreat"), ///
                                               st_numscalar("ncontrol"))
            mata: `nt' = st_numscalar("ntreat")
            mata: `shufflefun' = &shuffle_basic()
        }
        else {
            * Stratified individual-level randomization: We shuffle prod(`nstrata')
            * (0, 1) indicators with P(N / prod(`nstrata')) 1s and the rest 0s
            * at each strata; we then stack the indicators.
            local n2 = floor(`N' / 2)
            if (`st' > `n2') { // Stratified individual-level randomization
                di as err "{p}asked for `st' strata with `N' obs; " ///
                          "should be <= `n2'{p_end}"
                exit
            }

            tempvar groups group
            qui gen `group'  = 1
            qui gen `groups' = 1
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

            mata: `atreat' = simci_setup_strata(st_local("groups"), ///
                                                st_local("touse"),  ///
                                                0, `ptreat', `st')
            mata: `nt' = asarray(`atreat', "0")
            mata: `shufflefun' = &shuffle_strata()
        }
    }

    * Run the simulation
    * ------------------

    sort `notouse' `groups' `clusorder'

    if ("`compute'" == "ci")   local binary ""
    else if ("`binary'" == "") local binary _continuous
    else if ("`binary'" != "") local binary _binary

    tempname results output
    if ("`compute'" == "ci") | ("`compute'" == "effect") {
        mata: `results' = simci`binary'(st_local("depvar"),   /// depvar
                                        st_local("controls"), /// controls
                                        st_local("touse"),    /// touse
                                        `atreat',             /// treat
                                        `reps',               /// reps
                                        `effect',             /// effect
                                        `shufflefun')         //  shuffle
        if ("`bounds'" == "") {
            mata: `output' = parse_simci(`results', `alpha', `nt')
        }
        else if ("`bounds'" == "auto") {
            tempname resb cib
            mata: `resb' = simci(st_local("depvar"),   /// depvar
                                 st_local("controls"), /// controls
                                 st_local("touse"),    /// touse
                                 `atreat',             /// treat
                                 `reps',               /// reps
                                 0,                    /// effectret
                                 `shufflefun')         //  shuffle
            mata: `cib' = pctile(`resb'[, 1], (`alpha' / 2, 1 - `alpha' / 2))
            mata: `output' = parse_simci(`results', `alpha', `nt', `cib'[1], `cib'[2])
            mata: mata drop `cib' `resb'
        }
        else {
            gettoken lower upper: bounds
            mata: `output' = parse_simci(`results', `alpha', `nt', `lower', `upper')
        }
    }
    else if ("`compute'" == "power") {
        mata: `results' = power_search(st_local("depvar"),    /// depvar
                                       st_local("controls"),  /// controls
                                       st_local("touse"),     /// touse
                                       `atreat',              /// treat
                                       `reps',                /// reps
                                       `nt',                  /// nt
                                       0,                     /// effect
                                       `shufflefun',          /// shuffle
                                       `alpha',               /// alpha
                                       `kappa',               /// kappa
                                       "`direction'",         /// direction
                                       `tol',                 /// tol
                                       `startat',             /// startat
                                       &simci`binary'())      //  effectfun
            mata: `output' = parse_power(`results', `nt', `tol')
    }
    mata: mata drop `atreat' `nt' `results' `shufflefun'

    * Pretty printing
    * ---------------

    if ("`compute'" == "ci") | ("`compute'" == "effect") {
        local pow  "CI simulation for linear regression"
        local eq   "Y = a + b T + X + e"
        local null "Ho: b = `effect' versus Ha: b != `effect'""

        local trimci     = "(`:di trim("`:di %9.4f `r(lower)''")'"
        local trimci     = "`trimci', `:di trim("`:di %9.4f `r(upper)''")')"
        local trimci_pct = "(`:di trim("`:di %9.4f `r(lower_pct)''")'"
        local trimci_pct = "`trimci_pct', `:di trim("`:di %9.4f `r(upper_pct)''")')"

        di ""
        di "`pow': `eq'"
        di "`null'"
        di ""
        di "Study parameters:"
        di "    alpha = `:di %9.4f `alpha''"
        di "        P = `:di %9.4f `ptreat''"
        di "       PN = `:di %9.0gc `r(nt)''"
        di ""
        di "Simulated results:"
        di "       m1 = `:di %9.4f `r(mu)''"
        di "        b = `:di %9.4f `r(b)''"
        di "     sd_b = `:di %9.4f `r(sd)''"
        di "       ci = `trimci'"
        di "  ci / m1 = `trimci_pct'"
        if ("`r(power)'" != "") {
            di ""
            di "Simulated power for b = `effect':"
            di "    power = `:di %9.4f `r(power)''"
            di "    lower = `:di %9.4f `r(lower)''"
            di "    upper = `:di %9.4f `r(upper)''"
        }
    }
    else {
        local pow  "Power search for linear regression"
        local eq   "Y = a + b T + X + e"
        di ""
        di "Study parameters:"
        di "    alpha = `:di %9.4f `alpha''"
        di "        P = `:di %9.4f `ptreat''"
        if ("`cluster'" != "") di "     clus = `:di %9.0gc `nclus''"
        else di "      obs = `:di %9.0gc `N''"
        di "  treated = `:di %9.0gc `r(nt)''"
        di "    kappa = `:di %9.4f `kappa''"
        di "      tol = `:di %9.4f `tol''"
        di ""
        di "Simulated results:"
        di "       m1 = `:di %9.4f `r(mu)''"
        di "      MDE = `:di %9.4f `r(mde)''"
        di "    power = `:di %9.4f `r(power)''"
        di " MDE / m1 = `:di %9.4f `r(mde)' / `r(mu)''"
        di ""
        di "`stopped'"
        if ("`increase'" != "") di "`increase'"
    }
    di ""
end

***********************************************************************
*                           Mata functions                            *
***********************************************************************

mata:

// Simulate CI
// -----------

// Speed considerations:
//   - Randomization, regression, _and_ loop in pure Stata: Very slow.
//     Unless you want to give up precision you'll need to sort by a
//     random variable at least once per loop (or maybe crazy I/O if you
//     pre-generate all the treatment indicators and merge at each run).
//     Then AFAIK there's no way to compute just the OLS parameters in
//     Stata (it always computes the standard errors, etc.)
//   - Randomization in Mata, regression and loop in Stata: Slow. Mata
//     can shuffle the indicator and you can use `getmata` to put it in
//     memory. This saves you the sorting and the merging, but `getmata`
//     is not fast and `reg` is still a problem.
//   - Randomization and regression in Mata, loop in Stata: Fast. Mata
//     can shuffle the indicator, which is faster than sort and merge,
//     and then run just the OLS matrix algebra.
//   - Randomization, regression, and loop in Mata: Fastest. While the
//     speed improvements from the prior step are marginal, the loop in
//     pure Mata does run faster than the Stata loop doing Mata matrix
//     algebra at each turn.

real matrix function simci(string scalar depvar,
                           string scalar controls,
                           string scalar touse,
                           transmorphic treat,
                           real scalar reps,
                           real scalar effect,
                           pointer(real colvector function) shuffle)
{
    y = X = .
    st_view(y, ., depvar,   touse)
    st_view(X, ., controls, touse)
    results = J(reps, 3, 0)
    for (r = 1; r <= reps; r++) {
        st = (*shuffle)(treat)
        mu = mean(y[selectindex(!st)])
        XX = (st, X)
        b  = invsym(cross(XX, 1, XX, 1)) * cross(XX, 1, y, 0)
        results[r, 1::2] = b[1], mu
    }
    return(results)
}

real matrix function simci_continuous(string scalar depvar,
                                      string scalar controls,
                                      string scalar touse,
                                      transmorphic treat,
                                      real scalar reps,
                                      real scalar effect,
                                      pointer(real colvector function) shuffle)
{
    y = X = .
    st_view(y, ., depvar,   touse)
    st_view(X, ., controls, touse)
    results = J(reps, 3, 0)
    for (r = 1; r <= reps; r++) {
        st = (*shuffle)(treat)
        mu = mean(y[selectindex(!st)])
        XX = (st, X)
        yy = y :+ st * effect
        b  = invsym(cross(XX, 1, XX, 1)) * cross(XX, 1, yy, 0)
        results[r, 1::2] = b[1], mu
    }
    return(results)
}

real matrix function simci_binary(string scalar depvar,
                                  string scalar controls,
                                  string scalar touse,
                                  transmorphic treat,
                                  real scalar reps,
                                  real scalar effect,
                                  pointer(real colvector function) shuffle)
{
    y = X = .
    st_view(y, ., depvar,   touse)
    st_view(X, ., controls, touse)

    status  = 0
    sel     = (effect < 0)
    sign    = (effect > 0) - (effect < 0)
    results = J(reps, 3, 0)
    realeff = effect
    for (r = 1; r <= reps; r++) {
        st = (*shuffle)(treat)
        mu = mean(y[selectindex(!st)])
        tk = sum(st)
        sk = y' * st
        t0 = tk - sk

        // Check effect is OK
        if (effect > t0 / tk) {
            // printf("warning: effect larger than the portion of 0s in treatment\n")
            // printf("    effect = %9.4f\n", effect)
            // printf("      mu_t = %9.4f\n", mean(y[selectindex(st)]))
            // printf("      mu_c = %9.4f\n", mu)
            // printf("      sy_t = %9.0fc\n", sk)
            // printf("        nt = %9.0fc\n", tk)
            // printf("    bounds = (%9.4f, %9.4f)\n", - sk / tk, t0 / tk)
            results[r, 3] = 1
            realeff = t0 / tk
        }
        else if (effect < - sk / tk)  {
            // printf("warning: effect smaller than the portion of 1s in treatment\n")
            // printf("    effect = %9.4f\n", effect)
            // printf("      mu_t = %9.4f\n", mean(y[selectindex(st)]))
            // printf("      mu_c = %9.4f\n", mu)
            // printf("      sy_t = %9.0fc\n", sk)
            // printf("        nt = %9.0fc\n", tk)
            // printf("    bounds = (%9.4f, %9.4f)\n", - sk / tk, t0 / tk)
            results[r, 3] = 1
            realeff = - sk / tk
        }

        // Swap 0s or 1s, as requested
        yy  = y
        sk1 = abs(round(realeff * tk))
        if (realeff > 0) {
            sk2 = t0 - sk1
        }
        else {
            sk2 = sk - sk1
        }
        subs     = selectindex((y :== sel) :* st)
        yy[subs] = yy[subs] :+ jumble(J(sk1, 1, sign) \ J(sk2, 1, 0))

        // Run the regression
        XX = (st, X)
        b  = invsym(cross(XX, 1, XX, 1)) * cross(XX, 1, yy, 0)
        results[r, 1::2] = b[1], mu
        realeff = effect
    }
    return(results)
}

// Set up the randomization
// ------------------------

// Set up a vector with PN 1s and (1 - P)N 0s
transmorphic function simci_setup_basic(real scalar nt,
                                        real scalar nc)
{
    return(J(nt, 1, 1) \ J(nc, 1, 0))
}

// Set up a vector with PJ 1s and (1 - P)J 0s, recording nj the
// number of obs per cluster to correctly expand the indicator.
transmorphic function simci_setup_cluster(string scalar cluster,
                                          string scalar touse,
                                          real scalar ptreat)
{
    info  = panelsetup(st_data(., cluster, touse), 1)
    panel = panelstats(info)
    simci_panel_info(panel[1], panel[3], panel[4])

    nj     = info[, 2] - info[, 1] :+ 1
    nt     = round(rows(info) * ptreat)
    mtreat = J(nt, 1, 1) \ J(rows(info) - nt, 1, 0)
    atreat = asarray_create()
    asarray(atreat, "0", nt)
    asarray(atreat, "1", mtreat)
    asarray(atreat, "2", nj)
    return(atreat)
}

// Set ns = prod(`nstrata') vectors with PJ / ns 1s and (1 - P)J / ns 0s;
// if applicable, record nj the number of obs per cluster to correctly
// expand the indicator.
transmorphic function simci_setup_strata(string scalar groups,
                                         string scalar touse,
                                         real scalar cluster,
                                         real scalar ptreat,
                                         real scalar nstrata)
{
    if (cluster == 1) {
        nj = st_data(., touse)
        simci_panel_info(rows(nj), min(nj), max(nj))
    }

    info  = panelsetup(st_data(., groups, touse), 1)
    panel = panelstats(info)
    simci_panel_info(panel[1], panel[3], panel[4], "strata")

    ng     = info[, 2] - info[, 1] :+ 1
    tmat   = round(ng * ptreat)
    nt     = sum(tmat)
    tmat   = tmat, ng - tmat
    mtreat = J(sum(ng), 1, missingof(tmat))
    streat = asarray_create()
    for (i = 1; i <= nstrata; i++) {
        asarray(streat, strofreal(i), J(tmat[i, 1], 1, 1) \ J(tmat[i, 2], 1, 0))
    }
    atreat = asarray_create()
    asarray(atreat, "0", nt)
    asarray(atreat, "1", mtreat)
    asarray(atreat, "2", streat)
    asarray(atreat, "3", info)
    if (cluster == 1) {
        asarray(atreat, "4", nj)
    }
    return(atreat)
}

// Print number of clusters/strata and number of obs per cluster/strata
void function simci_panel_info(real scalar N,
                               real scalar pmin,
                               real scalar pmax,
                               | string scalar what)
{
    if (args() == 3) {
        what = "panel"
    }
    if (pmin == pmax) {
        nstr = strtrim(sprintf("%21.0fc", N))
        jstr = strtrim(sprintf("%21.0fc", pmin))
        printf("Balanced " + what + ". J = " + nstr + ", n_j = " + jstr + "\n")
    }
    else {
        nstr = strtrim(sprintf("%21.0fc", N))
        jstr = strtrim(sprintf("%21.0fc", pmin))
        jstr = jstr + " to " + strtrim(sprintf("%21.0fc", pmax))
        printf("Unbalanced " + what + ". J = " + nstr + ", n_j = " + jstr + "\n")
    }
}

// Shuffle the treatment, optionally by cluster, strata, cluster-strata
// --------------------------------------------------------------------

// Just shuffle the treatment vector
real colvector function shuffle_basic(real colvector treat)
{
    return(jumble(treat))
}

// Shuffle clusters, then expand using the number of individuals per cluster
real colvector function shuffle_cluster(transmorphic T)
{
    return(expand_rows(jumble(asarray(T, "1")), asarray(T, "2")))
}

// Shuffle each strata and stack
real colvector function shuffle_strata(transmorphic T)
{
    treat = asarray(T, "1")
    A     = asarray(T, "2")
    info  = asarray(T, "3")
    for (i = 1; i <= rows(info); i++) {
        treat[info[i, 1]::info[i, 2]] = jumble(asarray(A, strofreal(i)))
    }
    return(treat)
}

// Shuffle each strata, stack, and expand by cluster
real colvector function shuffle_strata_cluster(transmorphic T)
{
    treat = asarray(T, "1")
    A     = asarray(T, "2")
    info  = asarray(T, "3")
    nj    = asarray(T, "4")
    for (i = 1; i <= rows(info); i++) {
        treat[info[i, 1]::info[i, 2]] = jumble(asarray(A, strofreal(i)))
    }
    return(expand_rows(treat, nj))
}

// Expand rows (for cluster randomization)
// ---------------------------------------

function expand_rows(matrix X, real vector nj)
{
    real scalar i, b, e, n, add

    f = trunc(nj) :* (nj :> 0)
    _editmissing(f, 0)
    n = rows(X)

    if (n != length(f)) _error(3200)
    add = sum(f)
    Y   = J(add, cols(X), missingof(X))

    if (add > 0) {
        b = 1
        for (i = 1; i <= n; i++) {
            if (f[i] < 1) continue
            e = b + f[i] - 1
            Y[|b, 1 \ e, .|] = X[J(f[i], 1, i), .]
            b = b + f[i]
        }
    }
    return(Y)
}

// Search for power
// ----------------

real rowvector function power_search(string scalar depvar,
                                     string scalar controls,
                                     string scalar touse,
                                     transmorphic treat,
                                     real scalar reps,
                                     real scalar nt,
                                     real scalar effect,
                                     pointer(real colvector function) shuffle,
                                     real scalar alpha,
                                     real scalar kappa,
                                     string scalar direction,
                                     real scalar tol,
                                     real scalar startat,
                                     pointer(real matrix function) effectfun)
{
    // Baseline CI
    // -----------

    // We'll search for the power of this test (i.e. the rejection
    // criterion is the CI from the simulation)
    res = simci(depvar, controls, touse, treat, reps, effect, shuffle)
    ci  = pctile(res[, 1], (alpha / 2, 1 - alpha / 2))
    if (direction == "neg") {
        sign  = -1
        bound = ci[1]
    }
    else {
        sign  = 1
        bound = ci[2]
    }
    parse_simci(res, alpha, nt)

    // Search for power
    // ----------------

    // The strategy is to find a lower and upper bound for power
    // and hone in by taking the mid-point between the two. The
    // simulation stops when the estimated power is within tolerance
    // of the power level requested or when, because of randomness,
    // we can no longer improve upon the closest approximation.
    lower_val  = upper_val = upper_pow = lower_pow = .
    upper_diff = 1 - kappa
    lower_diff = kappa

    // Roughly, power should be ~ 0.5 at the CI bound. So we start
    // with 2 * bound and go from there (or at the start requested).
    factor   = (startat == 0)? 2: 1.25
    effect   = (kappa > 0.5)? factor * bound: bound / factor
    effect   = (startat == 0)? effect: startat
    stop     = 0
    iter     = 0
    printh   = 0
    printed  = 0

    // Pretty printing of CI info and the power estimate definition
    print_bound = strofreal(bound, "%9.4fc")
    print_alpha = strofreal(1 - alpha, "%9.3f")
    print_cil   = strofreal(ci[1], "%9.4f")
    print_ciu   = strofreal(ci[2], "%9.4f")
    print_ci    = "\nCI(" + print_alpha + ") = "
    print_ci    = print_ci + "(" + print_cil + ", " + print_ciu + ")"
    if (sign == -1) {
        printf(print_ci + "; power = mean(betas < " + print_bound + ")\n")
    }
    else {
        printf(print_ci + "; power = mean(betas > " + print_bound + ")\n")
    }

    // Will report on each iteration
    while (!stop) {
        ++iter
        if (printh & !printed) {
            printf("Iteration     Lower (power)       Upper (power)   Candidate (power)\n")
            printed = 1
        }

        // Get the distribution of betas after inducing -effect-
        eres   = (*effectfun)(depvar, controls, touse, treat,
                              reps, effect, shuffle)
        coefs  = eres[, 1]
        status = eres[, 3]
        ntrunc = sum(status)

        // Get the correct power estimate, based on the direction of the effect
        if (sign == -1) {
            power = mean(coefs :< bound)
        }
        else {
            power = mean(coefs :> bound)
        }

        // Update estimated power
        if (abs(power - kappa) < tol) {
            // Report on each iteration, once we have our bounds
            if (!missing(upper_val) & !missing(lower_val)) {
                printf("%9.0f %9.4f (%6.4f)  %9.4f (%6.4f)  %9.4f (%6.4f)\n",
                       iter,
                       lower_val,
                       lower_pow,
                       upper_val,
                       upper_pow,
                       effect,
                       power)
            }

            // If within tolerance, stop
            stop     = 1
            addprint = "\nFound power within tolerance"
            sim_mde  = effect
            sim_pow  = power
        }
        else if (power > kappa) {
            // If larger than power, check if it's an improvement
            power_diff = power - kappa
            // printf("%9.4f < %9.4f\n", power_diff, upper_diff)
            if (missing(upper_val) | (power_diff < upper_diff)) {
                // Report on each iteration, once we have our bounds
                if (!missing(upper_val) & !missing(lower_val)) {
                    printf("%9.0f %9.4f (%6.4f)  %9.4f (%6.4f)  %9.4f (%6.4f)\n",
                           iter,
                           lower_val,
                           lower_pow,
                           upper_val,
                           upper_pow,
                           effect,
                           power)
                }

                // If power improved, this is the new upper bound
                upper_val  = effect
                upper_diff = power_diff
                upper_pow  = power
                if (missing(lower_val)) {
                    effect = (effect + bound) / 2
                }
                else {
                    effect = (upper_val + lower_val) / 2
                }
            }
            else {
                // If not an improvement, stop
                stop     = 1
                addprint = "\nThe simulation stopped improving"
                sim_mde  = (upper_diff <  lower_diff) * upper_val +
                           (upper_diff >= lower_diff) * lower_val
                sim_pow  = (upper_diff <  lower_diff) * upper_pow +
                           (upper_diff >= lower_diff) * lower_pow
            }
        }
        else if (power < kappa) {
            // If smaller than power, check if it's an improvement
            power_diff = kappa - power
            // printf("%9.4f < %9.4f\n", power_diff, lower_diff)
            if (missing(lower_val) | (power_diff < lower_diff)) {
                // Report on each iteration, once we have our bounds
                if (!missing(upper_val) & !missing(lower_val)) {
                    printf("%9.0f %9.4f (%6.4f)  %9.4f (%6.4f)  %9.4f (%6.4f)\n",
                           iter,
                           lower_val,
                           lower_pow,
                           upper_val,
                           upper_pow,
                           effect,
                           power)
                }

                // If power improved, this is the new lower bound
                lower_val  = effect
                lower_diff = power_diff
                lower_pow  = power

                if (missing(upper_val) & (ntrunc == 0)) {
                    // If no upper bound, search for one unless the current
                    // lower bound was an upper bound in disguise.
                    effect = effect * factor
                }
                else if (!missing(upper_val)) {
                    // If there is an upper bound, continue
                    effect = (upper_val + lower_val) / 2
                }
                else if (ntrunc > 0) {
                    // If no upper bound and this was an upper bound in
                    // disguise, then you can't reach the requisite power
                    // level. End program with a warning.
                    stop       = 1
                    sim_mde    = effect
                    sim_pow    = power
                    ntrunc_str = strtrim(sprintf("%9.0fc times.", ntrunc))
                    addprint   = "\nWARNING: MDE was truncated " + ntrunc_str +
                                 "\nWARNING: Simulated MDE and power are" +
                                 " upper bounds; stopped"
                    addprint   = sprintf(addprint, ntrunc)
                }
            }
            else {
                // If not an improvement, stop
                stop     = 1
                addprint = "\nThe simulation stopped improving"
                sim_mde  = (upper_diff <= lower_diff) * upper_val +
                           (upper_diff >  lower_diff) * lower_val
                sim_pow  = (upper_diff <= lower_diff) * upper_pow +
                           (upper_diff >  lower_diff) * lower_pow
            }
        }
        printh = (!missing(upper_val) & !missing(lower_val))
    }
    st_numscalar("r(iter)", iter)
    stata(sprintf("return scalar iter = %15.0f", iter))

    // Pretty printing of results
    iter = strofreal(iter, "%9.0fc")
    sim_diff   = abs(sim_pow - kappa)
    print_pow  = strofreal(sim_pow,  "%9.4f")
    print_mde  = strofreal(sim_mde,  "%9.4fc")
    print_diff = strofreal(sim_diff, "%9.4f")
    print_abs  = "|" + strofreal(kappa, "%9.4f") + " - " + print_pow + "|"
    print_abs  = print_abs + " = " + print_diff
    st_local("power_diffstr", print_abs)
    printf(addprint + " after " + iter + " iterations.\n")
    printf(invtokens(("MDE = ", print_mde, ", ",
                      "power = ", print_pow, ", ",
                      print_abs, "\n")))
    return((sim_mde, sim_pow, sim_diff))
}

// Misc aux functions
// ------------------

real colvector function pctile(real vector x, real vector pctiles)
{
    _sort(x, 1)
    quantiles = J(length(pctiles), 1, missingof(pctiles))
    len = length(x)
    qq  = (1::len) / len
    for (j = 1; j <= length(pctiles); j++) {
        i = sum(qq :< pctiles[j]) + 1
        if (qq[i] == pctiles[j]) quantiles[j] = mean(x[i::(i + 1)])
        else  quantiles[j] = x[i]
    }
    return(quantiles)
}

real scalar bc_var(real vector x)
{
    n   = length(x)
    mux = mean(x)
    sse = sum((x :- mux):^2)
    return(sse / (n - 1))
}

// Parse output from simulations
// -----------------------------

void function parse_simci(real matrix results,
                          real scalar alpha,
                          real scalar nt,
                          | real scalar lower,
                            real scalar upper)
{
    // Check there were no problems (should only matter for binary outcomes)
    status = results[, 3]
    ntrunc = sum(status)
    if (ntrunc > 0) {
        addprint = "\nWARNING: MDE was truncated %9.0fc times." +
                   "\nWARNING: Simulated effect and power are upper bounds."
        printf(addprint, ntrunc)
    }

    // Parse results
    coefs = results[, 1]
    ci = pctile(coefs, (alpha / 2, 1 - alpha / 2))
    l  = ci[1]
    u  = ci[2]

    mres = mean(results)
    sdmu = sqrt(variance(colshape(results[, 2], 1)))
    b    = mres[1]
    mu   = mres[2]
    sd   = sqrt(bc_var(coefs))
    lpct = l / mu
    upct = u / mu

    // Set Stata's return values
    st_numscalar("r(nt)",    nt)
    st_numscalar("r(mu)",    mu)
    st_numscalar("r(mu_sd)", sdmu)
    st_numscalar("r(b)",     b)
    st_numscalar("r(sd)",    sd)

    st_numscalar("r(lower)",     l)
    st_numscalar("r(lower_pct)", lpct)
    st_numscalar("r(upper)",     u)
    st_numscalar("r(upper_pct)", upct)

    stata("return scalar nt    = r(nt)")
    stata("return scalar mu    = r(mu)")
    stata("return scalar mu_sd = r(mu_sd)")
    stata("return scalar b     = r(b)")
    stata("return scalar sd    = r(sd)")

    stata("return scalar lower     = r(lower)")
    stata("return scalar lower_pct = r(lower_pct)")
    stata("return scalar upper     = r(upper)")
    stata("return scalar upper_pct = r(upper_pct)")

    // If asked for bounds check
    simci = (nt, mu, b, sd, l, lpct, u, upct)
    if (args() > 3) {
        power = mean(!((coefs :< upper) :* (coefs :> lower)))
        st_numscalar("r(power)", power)
        stata(sprintf("return scalar power = %15.9f", power))
        simci = simci, power
    }
    st_matrix("simci", simci)
    stata("return matrix simci = simci")
}

void function parse_power(real rowvector results,
                          real scalar nt,
                          real scalar tol)
{
    st_numscalar("r(nt)", nt)
    st_numscalar("r(mde)", results[1])
    st_numscalar("r(power)", results[2])
    st_numscalar("r(power_diff)", results[3])

    stata("return scalar nt    = r(nt)")
    stata("return scalar mde   = r(mde)")
    stata("return scalar power = r(power)")
    stata("return scalar power_diff = r(power_diff)")

    stopped = st_local("power_diffstr")
    if (abs(results[3]) < tol) {
        stopped = stopped + " < " + strofreal(tol, "%9.4f")
    }
    else {
        stopped = stopped + "; no further improvements after "
        stopped = stopped + strofreal(st_numscalar("r(iter)"), "%9.0f") + " iterations."
        st_local("increase", "Consider increasing the # of repetitions!")
    }
    st_local("stopped", stopped)
}
end
