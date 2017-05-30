Description
-----------

Various notes and Stata programs to aid in power analysis.

I write the notes mainly to explain to myself how to write `power_reg` and `simci`. The former can compute [parametric power with clustering and stratification](https://github.com/mcaceresb/stata-power/blob/master/notes/power-clustered-notes.pdf). The latter can [simulate it](https://github.com/mcaceresb/stata-power/blob/master/notes/power-simulation-notes.pdf). I wrote them mainly because I could not figure out how to do either for the simple (OLS) case with built-in Stata commands.

Requirements
------------

I only have access to Stata 13.1, so I impose that to be the minimum.
The command is really simple, however, so I would not be surprised if it
worked with earlier versions. The exception would be the `fast` option,
which was compiled using C and v2.0 of the Stata Plugin Interface (SPI).
This might be tied to Stata 13.1. See how to recompile below.

Installation
------------

```stata
net install power_tools, from(https://raw.githubusercontent.com/mcaceresb/stata-power/master/)
```

To update, run
```stata
adoupdate, update
```

To uninstall, run
```stata
ado uninstall power_tools
```

Examples
---------

```stata
sysuse auto, clear
tempfile auto
save `auto'
qui forvalues i = 1 / 20 {
    append using `auto'
}
replace price = price + runiform()

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

* Simulate a CI using C plugin (only tested under Linux)
simci `depvar' `controls', reps(1000) fast

* You can check it's actually faster
net install benchmark, from(https://raw.githubusercontent.com/mcaceresb/stata-benchmark/master/)
benchmark, disp reps(10): qui simci `depvar' `controls', reps(1000)
benchmark, disp reps(10): qui simci `depvar' `controls', reps(1000) fast
```

Note that the `fast` option depends on the [GNU Scientific Library
(GSL)](https://www.gnu.org/software/gsl). If your system's `libgsl*.so`
and `libgslcblas*.so` files are not in `/usr/lib`, you should point
to them _**before**_ starting Stata by setting `LD_LIBRARY_PATH`.
I regularly `ssh` into a RedHat server, and the files were in
`/usr/local/lib`, so I ran
```bash
LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH
```

before starting Stata. You can add those lines to `~/.bashrc` to avoid
having to do that every time you log into a session.

Compiling
---------

The `fast` option uses a Stata plugin (compiled in C). To compile in Linux/Unix:
```bash
git clone https://github.com/mcaceresb/stata-power
cd stata-power
make SPI=3.0 # SPI v3.0, Stata 14 and later
make SPI=2.0 # SPI v2.0, Stata 13 and earlier
```

The advantage is twofold

1. First, C runs much faster than mata, which is how the function is implemented.
2. Second, C allows parallel loop execution. Since the simulation
   computes regression coefficients `reps` times, using N threads
   should result in an approximately Nx speed improvement. This works
   even with Stata/IC.

Note Mata runs faster than Stata's reg largely because this simulation
uses just the regression coefficients; reg computes a lot of additional
elements that the program does not use.

Dependencies
------------

- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl)
- [OpenMP](http://www.openmp.org)
- [Stata Plugin Interface](http://www.stata.com/plugins) (SPI version 2.0 for Stata < 14; version 3.0 for Stata >= 14)

I have only tested this in Linux so far. See [Stata's plugin documentation](http://www.stata.com/plugins/) for instructions on how to build the plugin in other platforms.

TODO
----

- [ ] Add documentation for `simci`
- [ ] Add documentation for `power_reg`
- [ ] Add examples for `power_reg`
- [ ] Compile `fast` option for Windows and OSX.
- [ ] Finish writing `fast` plugin so it works for clustering and stratification.
- [ ] Finish writing `fast` plugin so it works with `effect()` and `power()`.

License
-------

[MIT](https://github.com/mcaceresb/stata-power/blob/master/LICENSE)
