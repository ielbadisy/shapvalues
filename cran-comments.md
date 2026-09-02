## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.
* The second note ("Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'") is produced by flags in the local R
  installation's `CFLAGS`, not by anything the package sets; it does not
  appear on CRAN's build machines.

## Test environments

* local: R 4.5.1 on Ubuntu 24.04
