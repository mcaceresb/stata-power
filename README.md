Description
-----------

Various notes and Stata programs to aid in power analysis.

I write the notes mainly to explain to myself how to write `power_reg` and `simci`. The former can compute parametric power with clustering and stratification. The latter can simulate it. I wrote them mainly because I could not figure out how to do either for the simple (OLS) case with built-in Stata commands.

Requirements
------------

I only have access to Stata 13.1, so I impose that to be the minimum. The command is really simple, however, so I would not be surprised if it worked with earlier versions.

Installation
------------

```stata
net install power_tools, from(https://raw.githubusercontent.com/mcaceresb/stata-power/master/)
```

Examples
---------

```stata
sysuse auto, clear
local depvar      price
local controls    mpg rep78
local cluster     rep78
local stratum     gear_ratio
local categorical foreign

* Parametric power
power_reg `depvar' `controls'
power_reg `depvar' `controls', cluster(`cluster')
power_reg `depvar' `controls', cluster(`cluster') strata(`stratum') nstrata(2)

* Simulate a CI
simci `depvar' `controls', reps(1000)
simci `depvar' `controls', reps(1000) cluster(`cluster')
simci `depvar' `controls' `stratum', reps(1000) cluster(`cluster') ///
    strata(`stratum') nstrata(2)
simci `depvar' `controls' i.`categorical', reps(1000) ///
    strata(`categorical') nstrata(0)

* Simulate MDE given power
simci `depvar' `controls', reps(1000) power(kappa(0.8) direction(neg))
simci `depvar' `controls', reps(1000) cluster(`cluster') ///
    power(kappa(0.8) direction(neg))
simci foreign  `controls' `stratum', reps(1000) strata(`stratum') nstrata(2) ///
    power(binary dir(pos))
simci `depvar' `controls', reps(1000) cluster(`cluster') strata(`stratum') ///
    nstrata(2) power(dir(pos))

* Simulate power given MDE
simci foreign  `controls', reps(1000) effect(effect(-0.5) ///
    bounds(-0.2 0.2) binary)
simci `depvar' `controls', reps(1000) cluster(`cluster') ///
    effect(effect(-0.5) bounds(-0.2 0.2))
simci `depvar' `controls', reps(1000) strata(`stratum') nstrata(2) ///
    effect(effect(-0.5) bounds(-0.2 0.2))
```

TODO
----

- [ ] Add documentation for `simci`
- [ ] Add documentation for `power_reg`
- [ ] Add examples for `power_reg`
