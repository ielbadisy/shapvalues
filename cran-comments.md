## Submission

This is a new release.

The source package is built with

    R CMD build --compact-vignettes=both

so the PDF vignette under `inst/doc` is already GS+qpdf compacted.

## R CMD check results

0 errors | 0 warnings | 2 notes

* One note is the standard "New submission" note.
* The other note ("Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'") comes from flags in the local R
  installation's `CFLAGS`, not from anything the package sets; it does not
  appear on CRAN's build machines.

## Test environments

* local: R 4.5.1 on Ubuntu 24.04
