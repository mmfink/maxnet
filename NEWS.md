---
editor_options: 
  markdown: 
    wrap: 72
---

# maxnet 0.01-6

Fork of <https://github.com/mrmaxent/maxnet> starting from commit
[271b85e](https://github.com/mrmaxent/maxnet/commit/271b85ec1c1476ac1af54da891ad8f885abb79b5)

**New functions:**

-   *maxnet.varImp* Produces a vector and plot of variable importance
-   *maxnet.coefficients* A table of 'features' and lambdas similar to
    what was produced in the java Maxent

**Improvements:**

-   Improved and expanded documentation
-   *plot.maxnet* Tweaks to the appearance of response plots

**Fixes:**

-   Coerce p to vector if possible (upstream issue #16)
-   Remove ambiguous call to "..." (upstream issue #7)

# maxnet 0.1.4

**Note that upstream has made changes since without bumping version**

Two small changes in version 0.1.4:

addsamplestobackground is now true by default, which fixes convergence
problems with some datasets character levels for factors are now
supported
