#' Brittany Loess-Soil Physicochemical Properties
#'
#' @description
#' Standardised physicochemical measurements from 37 cultivated fields on
#' loess-silt parent material in Brittany, France (Morvan et al., 2023).
#' The dataset provides a realistic statistical disclosure control (SDC)
#' scenario with \eqn{p = 6} variables in two natural blocks of
#' \eqn{p_1 = p_2 = 3}.
#'
#' \strong{Block 1 -- public indicators} (\eqn{p_1 = 3}): standard
#' agronomic measurements already published in public soil maps, released
#' as observed data: \code{pH_water}, \code{pH_KCl},
#' \code{log_CEC_Metson}.
#'
#' \strong{Block 2 -- sensitive indicators} (\eqn{p_2 = 3}): commercially
#' sensitive farm-management measurements, released as PS synthetic draws:
#' \code{log_Organic_C}, \code{log_Total_N}, \code{log_P_Olsen}.
#'
#' @details
#' \subsection{SDC scenario}{
#' A regional agency releases standardised soil data from loess farms.
#' The public block (pH, log CEC) is released as observed. The sensitive
#' block (log OC, log N, log P) is released as plug-in sampling (PS)
#' synthetic draws to protect individual farm values.
#' }
#'
#' \subsection{Preprocessing}{
#' The full dataset (Morvan et al., 2023) spans 137 fields across three
#' parent-material types. This subset retains only the most-frequent
#' parent-material group (loess-silt). Variables with log-normal
#' distributions (CEC, OC, N, P) were log-transformed. Mahalanobis
#' outliers at the \eqn{\chi^2_{6, 0.975}} quantile were removed.
#' All variables were standardised to mean 0 and standard deviation 1.
#' All five tests in \code{\link{mvn_test}} fail to reject multivariate
#' normality at the 5\% level.
#' }
#'
#' @format A numeric matrix with 37 rows and 6 columns (all standardised,
#'   mean 0, sd 1):
#' \describe{
#'   \item{pH_water}{Soil pH measured in water suspension.}
#'   \item{pH_KCl}{Soil pH measured in 1 M KCl suspension.}
#'   \item{log_CEC_Metson}{Log cation exchange capacity, Metson method
#'     (log meq per 100 g soil).}
#'   \item{log_Organic_C}{Log soil organic carbon (log g/kg).}
#'   \item{log_Total_N}{Log total soil nitrogen (log g/kg).}
#'   \item{log_P_Olsen}{Log Olsen-P available phosphorus
#'     (log g P2O5/kg).}
#' }
#'
#' @section Multivariate Normality:
#' \tabular{lrr}{
#'   \strong{Test} \tab \strong{Statistic} \tab \strong{p-value} \cr
#'   SW: pH_water       \tab W = 0.968 \tab 0.355 \cr
#'   SW: pH_KCl         \tab W = 0.977 \tab 0.614 \cr
#'   SW: log_CEC_Metson \tab W = 0.955 \tab 0.139 \cr
#'   SW: log_Organic_C  \tab W = 0.957 \tab 0.165 \cr
#'   SW: log_Total_N    \tab W = 0.977 \tab 0.638 \cr
#'   SW: log_P_Olsen    \tab W = 0.975 \tab 0.568 \cr
#'   Mardia skewness    \tab chi-sq = 46.95 \tab 0.800 \cr
#'   Mardia kurtosis    \tab z = -1.39      \tab 0.166 \cr
#'   Henze-Zirkler      \tab HZ = 0.026     \tab 0.336 \cr
#'   Royston H          \tab H = 1.07       \tab 0.710 \cr
#' }
#'
#' @source
#' Morvan, T., Lambert, Y., Germain, P., Lemercier, B., Moreira, M.
#' and Beff, L. (2023). A dataset of physico-chemical properties,
#' extractable organic N, N mineralisation and physical organic matter
#' fractionation of soils. \emph{Data in Brief}, \strong{51}, 109776.
#' \doi{10.1016/j.dib.2023.109776}
#'
#' Data repository (CC BY 4.0): \doi{10.57745/DGIPGR}
#'
#' @seealso \code{\link{mvn_test}}, \code{\link{simSynthData}},
#'   \code{\link{ps_test}}
#'
#' @examples
#' data(brittany_soil_ps)
#' dim(brittany_soil_ps)
#' colnames(brittany_soil_ps)
#'
#' # Verify multivariate normality
#' mvn_test(brittany_soil_ps, hz_nsim = 500, plot = FALSE)
#'
#' # Generate 3 PS synthetic releases
#' set.seed(1)
#' V3 <- simSynthData(brittany_soil_ps, M = 3)
#'
#' # Test independence: public vs sensitive block
#' independence_test(V3, M = 3,
#'   group_a = c("pH_water", "pH_KCl", "log_CEC_Metson"),
#'   group_b = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
#'   iterations = 2000L)
"brittany_soil_ps"


#' @title Employee Attitude Survey Data (Multivariate Normal)
#'
#' @description
#' A \eqn{30 \times 4} numeric matrix derived from the \code{attitude}
#' dataset in base R, retaining the four variables that most closely
#' satisfy multivariate normality. The data represent employee attitude
#' survey scores from clerical employees in a large financial
#' organisation. Each row corresponds to a work group; the four columns
#' are percentage scores for different aspects of attitude.
#'
#' This dataset is included in \pkg{PSinference} as a real-data example
#' for the multiple-release plug-in sampling procedures. It has been
#' verified to satisfy multivariate normality using both univariate
#' Shapiro-Wilk tests (all p-values > 0.25) and Mardia's multivariate
#' skewness and kurtosis tests (skewness p = 0.11, kurtosis p = 0.99).
#'
#' @format A numeric matrix with 30 rows and 4 columns:
#' \describe{
#'   \item{rating}{Overall rating of the work group (percentage).}
#'   \item{complaints}{Score for handling of employee complaints
#'     (percentage).}
#'   \item{privileges}{Score for employee privileges (percentage).}
#'   \item{learning}{Score for learning opportunities (percentage).}
#' }
#'
#' @details
#' The full \code{attitude} dataset has 7 variables and 30 observations.
#' The variables \code{raises}, \code{critical}, and \code{advance}
#' were excluded because \code{critical} showed marginal normality
#' (Shapiro-Wilk p = 0.035) and the full 7-variable dataset showed
#' borderline Mardia skewness (p = 0.025). The retained four variables
#' give a clean multivariate normal dataset suitable for illustrating
#' the exact PS inference procedures with \code{p_1 = p_2 = 2}.
#'
#' @section Multivariate Normality Verification:
#' \tabular{lrr}{
#'   \strong{Test} \tab \strong{Statistic} \tab \strong{p-value} \cr
#'   SW: rating    \tab W = 0.955 \tab 0.254 \cr
#'   SW: complaints\tab W = 0.966 \tab 0.552 \cr
#'   SW: privileges\tab W = 0.971 \tab 0.643 \cr
#'   SW: learning  \tab W = 0.958 \tab 0.304 \cr
#'   Mardia skewness\tab chi-sq = 19.56 \tab 0.111 \cr
#'   Mardia kurtosis\tab z = 0.013 \tab 0.989 \cr
#' }
#'
#' @source
#' Derived from the \code{attitude} dataset in base R.
#' Chatterjee, S. and Price, B. (1977). \emph{Regression Analysis by
#' Example}. New York: Wiley.
#'
#' @seealso \code{\link{ps_mtcars}}, \code{\link{mvn_test}},
#'   \code{\link{ps_test}}
#'
#' @examples
#' data(ps_attitude)
#'
#' # Check multivariate normality
#' mvn_test(ps_attitude)
#'
#' # Generate 5 synthetic releases and run sphericity test
#' set.seed(1)
#' V <- simSynthDataMultiple(ps_attitude, M = 5)
#' sphericity_test(V, M = 5)
"ps_attitude"


#' @title Motor Trend Car Road Tests Data (Log-Transformed,
#'   Multivariate Normal)
#'
#' @description
#' A \eqn{32 \times 4} numeric matrix derived from the \code{mtcars}
#' dataset in base R, after log-transforming the right-skewed variables
#' to achieve approximate multivariate normality. The four columns are
#' \code{log(mpg)}, \code{log(disp)}, \code{drat}, and \code{log(wt)},
#' representing fuel efficiency, engine displacement, rear axle ratio,
#' and vehicle weight for 32 automobile models from the 1974 Motor
#' Trend US magazine.
#'
#' This dataset illustrates that continuous survey or administrative
#' data often require a variance-stabilising transformation before
#' applying multivariate normal inference. After log-transformation,
#' the dataset passes all standard multivariate normality tests
#' (Mardia skewness p = 0.44, kurtosis p = 0.61), making it suitable
#' for illustrating the plug-in sampling procedures.
#'
#' @format A numeric matrix with 32 rows and 4 columns:
#' \describe{
#'   \item{log_mpg}{Natural logarithm of miles per US gallon.}
#'   \item{log_disp}{Natural logarithm of engine displacement
#'     (cubic inches).}
#'   \item{drat}{Rear axle ratio (not transformed; already
#'     approximately normal).}
#'   \item{log_wt}{Natural logarithm of weight (1000 lbs).}
#' }
#'
#' @details
#' The raw \code{mtcars} variables \code{mpg}, \code{disp}, and
#' \code{wt} are right-skewed; log-transformation gives approximately
#' symmetric, bell-shaped marginal distributions. The variable
#' \code{drat} is already approximately normal and is included without
#' transformation. This dataset provides a natural two-block structure:
#' Block 1 (\code{log_mpg}, \code{log_disp}) represents engine
#' performance, and Block 2 (\code{drat}, \code{log_wt}) represents
#' vehicle dynamics, enabling meaningful independence and regression
#' tests between the two blocks.
#'
#' @section Multivariate Normality Verification:
#' \tabular{lrr}{
#'   \strong{Test} \tab \strong{Statistic} \tab \strong{p-value} \cr
#'   SW: log_mpg   \tab W = 0.982 \tab 0.699 \cr
#'   SW: log_disp  \tab W = 0.954 \tab 0.053 \cr
#'   SW: drat      \tab W = 0.956 \tab 0.110 \cr
#'   SW: log_wt    \tab W = 0.959 \tab 0.156 \cr
#'   Mardia skewness\tab chi-sq = 13.95 \tab 0.443 \cr
#'   Mardia kurtosis\tab z = 0.508 \tab 0.611 \cr
#' }
#'
#' @source
#' Derived from the \code{mtcars} dataset in base R.
#' Henderson, H. V. and Velleman, P. F. (1981). Building multiple
#' regression models interactively.
#' \emph{Biometrics}, 37, 391--411.
#'
#' @seealso \code{\link{ps_attitude}}, \code{\link{mvn_test}},
#'   \code{\link{ps_test}}
#'
#' @examples
#' data(ps_mtcars)
#'
#' # Check multivariate normality
#' mvn_test(ps_mtcars)
#'
#' # Generate 3 synthetic releases
#' set.seed(1)
#' V <- simSynthDataMultiple(ps_mtcars, M = 3)
#'
#' # Test independence: engine performance vs vehicle dynamics
#' independence_test(V, M = 3, part = 2L)
#'
#' # Test regression: does engine performance predict vehicle dynamics?
#' regression_test(V, M = 3, part = 2L)
"ps_mtcars"
