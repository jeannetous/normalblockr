
This release hopefully fixes the ERROR currently shown on
`r-devel-linux-x86_64-fedora-gcc`, where re-building
`breast-cancer-proteomics.Rmd` aborts with

    malloc(): unsorted double linked list corrupted

The abort signature has varied from run to run -- an earlier one on the same
check was

    *** caught segfault ***
    address 0x580, cause 'memory not mapped'

which is itself consistent with the diagnosis below: heap corruption, whose
symptom depends on what the allocator happens to touch next.

It is submitted sooner than the usual interval for that reason; please let
us know if you would rather we wait.

## The fix

The crash is in the vignette's graphical-lasso section. Until now that step called back into R -- `glassoFast::glassoFast()` -- from inside the C++
(V)EM loop, once per M-step, i.e. thousands of R re-entries from within a
single `.Call()`. One memory-safety bug on that boundary had been found, and the crash signature matched what our own CI showed intermittently on Linux
across several R versions.

This new release removes that boundary: the graphical lasso is now implemented in the package's own C++ (`src/graphical_lasso.h`). There is no longer any call from C++ back into R inside the recursion.

We were never able to reproduce the crash locally (it appeared on roughly one CI run in ten, across R versions and platforms, which is the expected profile of heap corruption), so we cannot demonstrate the fix directly on a failing case. The argument is that the suspected mechanism is gone. As supporting evidence, the vignette's full  pipeline has been run 40 times under `MALLOC_CHECK_=3` on the new code, without incident.

## Other changes

The release also adds a second model family (variables clustered by their
regression response to the covariates rather than by their covariance). See
NEWS.md.

## Tested environments

* locally on Ubuntu Linux 24.04, R-release, GCC (`R CMD check --as-cran`)

* remotely with github-actions:
  - Linux ubuntu 24.04, R-release
  - Linux ubuntu 24.04, R-devel
  - Linux ubuntu 24.04, R-oldrel-1
  - macOS, R-release
  - Windows Server, R-release

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

The local `--as-cran` run reports notes that are artifacts of this machine
rather than of the package: the non-portable `-mno-omit-leaf-frame-pointer`
that Debian/Ubuntu's own `r-base` puts in `Makeconf` (the package's
`src/Makevars` sets only `CXX_STD`, `-DARMA_WARN_LEVEL=1` and the BLAS/LAPACK
libs), HTML Tidy not being installed here, and -- intermittently, when GitHub
rate-limits us -- a 429 on the README's last-commit badge URL, which resolves
normally otherwise.

Should the incoming checks flag "Days since last update" (0.2.1 was published
on 2026-09-03): this submission is early on purpose, to answer the check
ERROR on r-devel-linux-x86_64-fedora-gcc described above. We are happy to
wait if you would rather we did.
