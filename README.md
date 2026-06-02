# PSinference

[![R-CMD-check](https://github.com/ricardomourarpm/PSinference/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ricardomourarpm/PSinference/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/PSinference)](https://cran.r-project.org/package=PSinference)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

**PSinference** provides exact finite-sample inferential procedures
for singly and **multiply** released plug-in sampling (PS) synthetic
datasets under a multivariate normal model.

The key insight is simple: an analyst who receives $M$ independent
synthetic datasets $V_1, \ldots, V_M$ of size $n$ naturally treats
all released data as a whole by stacking them into a single dataset
of size $Mn$. This stacking is statistically justified — the $Mn$
rows are conditionally i.i.d. given the original data — and
immediately extends the exact procedures of Klein et al. (2021) to
arbitrary $M \geq 1$ via the substitution $n \to Mn$.

## Features

- **Four exact inferential procedures**: generalised variance,
  sphericity, independence, and regression
- **Multiple releases**: supports $M \geq 1$ synthetic datasets via
  the stacked-data representation
- **S3 class system**: `ps_test` objects with `print`, `summary`,
  and `plot` methods
- **Comparison with Reiter's combining rules**: via `compare_reiter()`
- **Full vignette** with real-data illustrations
- **Comprehensive input validation** and informative error messages

## Installation

```r
# Stable version from CRAN
install.packages("PSinference")

# Development version
remotes::install_github("mnrzrad/PSinference2")
```

## Quick Start

```r
library(PSinference)
data(nhanes_ps_clean)

# Generate 5 synthetic releases (stacked)
set.seed(42)
V <- simSynthDataMultiple(nhanes_ps_clean, M = 5)

# Sphericity test
res <- sphericity_test(V, M = 5)
print(res)
plot(res)

# Or use the unified wrapper
ps_test(V, M = 5, test = "independence", part = 4L)

# M = 1 recovers Klein et al. (2021)
V1 <- simSynthData(nhanes_ps_clean)
ps_test(V1, M = 1, test = "sphericity")
```

## References

Klein, M., Moura, R., and Sinha, B. (2021). Multivariate normal
inference based on singly imputed synthetic data under plug-in
sampling. *Sankhya B*, 83, 273--287.

Raghunathan, T. E., Reiter, J. P., and Rubin, D. B. (2003).
Multiple imputation for statistical disclosure limitation.
*Journal of Official Statistics*, 19, 1--16.

## Funding

This work was supported by the Fundação para a Ciência e a Tecnologia
(FCT, Portugal) under projects UIDB/00297/2020 and UIDP/00297/2020
(NOVAMath).

## License

GPL-3.
