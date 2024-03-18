\name{calibrate}

\alias{calibrate}

\title{Calibrate raw data to crisp or fuzzy sets}

\description{
This function transforms (calibrates) the raw data to either crisp or fuzzy sets
values, using both the \dfn{direct} and the \dfn{indirect} methods of calibration.
}

\usage{
calibrate(x, type = "fuzzy", method = "direct", thresholds = NA,
          logistic = TRUE, idm = 0.95, ecdf = FALSE, below = 1, above = 1, ...)
}

\arguments{
  \item{x}{A numerical causal condition.}
  \item{type}{Calibration type, either \bold{\code{"crisp"}} or \bold{\code{"fuzzy"}}.}
  \item{method}{Calibration method, either \bold{\code{"direct"}}, \bold{\code{"indirect"}}
        or \bold{\code{"TFR"}}.}
  \item{thresholds}{A vector of (named) thresholds.}
  \item{logistic}{Calibrate to fuzzy sets using the logistic function.}
  \item{idm}{The set inclusion degree of membership for the logistic function.}
  \item{ecdf}{Calibrate to fuzzy sets using the empirical cumulative distribution 
        function of the raw data.}
  \item{below}{Numeric (non-negative), determines the shape below crossover.}
  \item{above}{Numeric (non-negative), determines the shape above crossover.}
  \item{...}{Additional parameters, mainly for backwards compatibility.}
}

\details{
Calibration is a transformational process from raw numerical data (interval or
ratio level of measurement) to set membership scores, based on a certain number of
qualitative anchors.

When \bold{\code{type = "crisp"}}, the process is similar to recoding the original
values to a number of categories defined by the number of thresholds. For one
threshold, the calibration produces two categories (intervals): 0 if below, 1 if above.
For two thresholds, the calibration produces three categories: 0 if below the first threshold,
1 if in the interval between the thresholds and 2 if above the second threshold etc.

When \bold{\code{type = "fuzzy"}}, calibration produces fuzzy set membership scores, using
three anchors for the increasing or decreasing \dfn{s-shaped} distributions (including
the logistic function), and six anchors for the increasing or decreasing \dfn{bell-shaped}
distributions.

The argument \bold{\code{thresholds}} can be specified either as a simple numeric vector, or as a
named numeric vector. If used as a named vector, for the first category of \dfn{s-shaped}
distributions, the names of the thresholds should be:\cr
\tabular{rl}{
    \bold{\code{"e"}} \tab for the full set exclusion\cr
    \bold{\code{"c"}} \tab for the set crossover\cr
    \bold{\code{"i"}} \tab for the full set inclusion\cr
}
For the second category of \dfn{bell-shaped} distributions, the names of the thresholds
should be:\cr
\tabular{rl}{
    \bold{\code{"e1"}} \tab for the first (left) threshold for full set exclusion\cr
    \bold{\code{"c1"}} \tab for the first (left) threshold for set crossover\cr
    \bold{\code{"i1"}} \tab for the first (left) threshold for full set inclusion\cr
    \bold{\code{"i2"}} \tab for the second (right) threshold for full set inclusion\cr
    \bold{\code{"c2"}} \tab for the second (right) threshold for set crossover\cr
    \bold{\code{"e2"}} \tab for the second (right) threshold for full set exclusion\cr
}
If used as a simple numerical vector, the order of the values matter.

If \bold{\code{e}} \eqn{<} \bold{\code{c}} \eqn{<} \bold{\code{i}}, then the membership 
function is increasing from \bold{\code{e}} to \bold{\code{i}}. If \bold{\code{i}} \eqn{<} 
\bold{\code{c}} \eqn{<} \bold{\code{e}}, then the membership function is decreasing from 
\bold{\code{i}} to \bold{\code{e}}.

Same for the \dfn{bell-shaped} distribution, if \bold{\code{e1}} \eqn{<} \bold{\code{c1}}
\eqn{<} \bold{\code{i1}} \eqn{\le} \bold{\code{i2}} \eqn{<} \bold{\code{c2}} \eqn{<}
\bold{\code{e2}}, then the membership function is first increasing from \bold{\code{e1}}
to \bold{\code{i1}}, then flat between \bold{\code{i1}} and  \bold{\code{i2}}, and then
decreasing from \bold{\code{i2}} to \bold{\code{e2}}. In contrast, if \bold{\code{i1}}
\eqn{<} \bold{\code{c1}} \eqn{<} \bold{\code{e1}} \eqn{\le} \bold{\code{e2}} \eqn{<}
\bold{\code{c2}} \eqn{<} \bold{\code{i1}}, then the membership function is first decreasing
from \bold{\code{i1}} to \bold{\code{e1}}, then flat between \bold{\code{e1}} and 
\bold{\code{e2}}, and finally increasing from \bold{\code{e2}} to \bold{\code{i2}}.

When \bold{\code{logistic = TRUE}} (the default), the argument \bold{\code{idm}} specifies the
inclusion degree of membership for the logistic function. If \bold{\code{logistic = FALSE}}, the
function returns linear \dfn{s-shaped} or \dfn{bell-shaped} distributions (curved using the
arguments \bold{\code{below}} and \bold{\code{above}}), unless activating the argument
\bold{\code{ecdf}}.

If there is no prior knowledge on the shape of the distribution, the argument \bold{\code{ecdf}}
asks the computer to determine the underlying distribution of the empirical, observed points,
and the calibrated measures are found along that distribution.

Both \bold{\code{logistic}} and \bold{\code{ecdf}} arguments can be used only for \dfn{s-shaped}
distributions (using 3 thresholds), and they are mutually exclusive.

The parameters \bold{\code{below}} and \bold{\code{above}} (active only when both
\bold{\code{logistic}} and \bold{\code{ecdf}} are deactivated, establish the degree of
concentration and  dilation (convex or concave shape) between the threshold and crossover:\cr
\tabular{ll}{
    \bold{\code{0 < below < 1}} \tab dilates in a concave shape below the crossover\cr
    \bold{\code{below = 1}}     \tab produces a linear shape (neither convex, nor concave)\cr
    \bold{\code{below > 1}}     \tab concentrates in a convex shape below the crossover\cr
    \bold{\code{0 < above < 1}} \tab dilates in a concave shape above the crossover\cr
    \bold{\code{above = 1}}     \tab produces a linear shape (neither convex, nor concave)\cr
    \bold{\code{above > 1}}     \tab concentrates in a convex shape above the crossover\cr
}
Usually, \bold{\code{below}} and \bold{\code{above}} have equal values, unless specific reasons
exist to make them different.

For the \bold{\code{type = "fuzzy"}} it is also possible to use the \bold{\code{"indirect"}}
method to calibrate the data, using a procedure first introduced by Ragin (2008). The indirect method
assumes a vector of thresholds to cut the original data into equal intervals, then it applies
a (quasi)binomial logistic regression with a fractional polynomial equation.

The results are also fuzzy between 0 and 1, but the method is entirely different: it has no anchors
(specific to the direct method), and it doesn't need to specify a calibration function to calculate
the scores with.

The third method applied to fuzzy calibrations is called \bold{\code{"TFR"}} and calibrates categorical
data (such as Likert type response scales) to fuzzy values using the Totally Fuzzy and Relative
method (Chelli and Lemmi, 1995).
}

\value{
A numeric vector of set membership scores, either crisp (starting from 0 with increments of 1),
or fuzzy numeric values between 0 and 1.
}

\author{
Adrian Dusa
}

\references{
Cheli, B.; Lemmi, A. (1995) \dQuote{A 'Totally' Fuzzy and Relative
Approach to the Multidimensional Analysis of Poverty}. In \emph{Economic Notes},
vol.1, pp.115-134.

Dusa, A. (2019) \emph{QCA with R. A Comprehensive Resource}. 
Springer International Publishing, \doi{10.1007/978-3-319-75668-4}.

Ragin, C. (2008) \dQuote{Fuzzy Sets: Calibration Versus Measurement.}
In \emph{The Oxford Handbook of Political Methodology}, edited by Janet
Box-Steffensmeier, Henry E. Brady, and David Collier, pp.87-121.
Oxford: Oxford University Press.
}

\examples{

# generate heights for 100 people
# with an average of 175cm and a standard deviation of 10cm
set.seed(12345)
x <- rnorm(n = 100, mean = 175, sd = 10)


cx <- calibrate(x, type = "crisp", thresholds = 175)
plot(x, cx, main="Binary crisp set using 1 threshold",
     xlab = "Raw data", ylab = "Calibrated data", yaxt="n")
axis(2, at = 0:1)


cx <- calibrate(x, type = "crisp", thresholds = c(170, 180))
plot(x, cx, main="3 value crisp set using 2 thresholds",
     xlab = "Raw data", ylab = "Calibrated data", yaxt="n")
axis(2, at = 0:2)


# calibrate to a increasing, s-shaped fuzzy-set
cx <- calibrate(x, thresholds = "e=165, c=175, i=185")
plot(x, cx, main = "Membership scores in the set of tall people", 
     xlab = "Raw data", ylab = "Calibrated data")

     
# calibrate to an decreasing, s-shaped fuzzy-set
cx <- calibrate(x, thresholds = "i=165, c=175, e=185")
plot(x, cx, main = "Membership scores in the set of short people", 
     xlab = "Raw data", ylab = "Calibrated data")


# when not using the logistic function, linear increase
cx <- calibrate(x, thresholds = "e=165, c=175, i=185", logistic = FALSE)
plot(x, cx, main = "Membership scores in the set of tall people", 
     xlab = "Raw data", ylab = "Calibrated data")


# tweaking the parameters "below" and "above" the crossover,
# at value 3.5 approximates a logistic distribution, when e=155 and i=195
cx <- calibrate(x, thresholds = "e=155, c=175, i=195", logistic = FALSE,
      below = 3.5, above = 3.5)
plot(x, cx, main = "Membership scores in the set of tall people", 
     xlab = "Raw data", ylab = "Calibrated data")


# calibrate to a bell-shaped fuzzy set
cx <- calibrate(x, thresholds = "e1=155, c1=165, i1=175, i2=175, c2=185, e2=195",
      below = 3, above = 3)
plot(x, cx, main = "Membership scores in the set of average height",
     xlab = "Raw data", ylab = "Calibrated data")


# calibrate to an inverse bell-shaped fuzzy set
cx <- calibrate(x, thresholds = "i1=155, c1=165, e1=175, e2=175, c2=185, i2=195",
      below = 3, above = 3)
plot(x, cx, main = "Membership scores in the set of non-average height",
     xlab = "Raw data", ylab = "Calibrated data")


# the default values of "below" and "above" will produce a triangular shape
cx <- calibrate(x, thresholds = "e1=155, c1=165, i1=175, i2=175, c2=185, e2=195")
plot(x, cx, main = "Membership scores in the set of average height",
     xlab = "Raw data", ylab = "Calibrated data")


# different thresholds to produce a linear trapezoidal shape
cx <- calibrate(x, thresholds = "e1=155, c1=165, i1=172, i2=179, c2=187, e2=195")
plot(x, cx, main = "Membership scores in the set of average height",
     xlab = "Raw data", ylab = "Calibrated data")


# larger values of above and below will increase membership in or out of the set
cx <- calibrate(x, thresholds = "e1=155, c1=165, i1=175, i2=175, c2=185, e2=195",
      below = 10, above = 10)
plot(x, cx, main = "Membership scores in the set of average height",
     xlab = "Raw data", ylab = "Calibrated data")


# while extremely large values will produce virtually crisp results
cx <- calibrate(x, thresholds = "e1=155, c1=165, i1=175, i2=175, c2=185, e2=195",
      below = 10000, above = 10000)
plot(x, cx, main = "Binary crisp scores in the set of average height",
     xlab = "Raw data", ylab = "Calibrated data", yaxt="n")
axis(2, at = 0:1)
abline(v = c(165, 185), col = "red", lty = 2)

# check if crisp
round(cx, 0)


# using the empirical cumulative distribution function
# require manually setting logistic to FALSE
cx <- calibrate(x, thresholds = "e=155, c=175, i=195", logistic = FALSE,
      ecdf = TRUE)
plot(x, cx, main = "Membership scores in the set of tall people", 
     xlab = "Raw data", ylab = "Calibrated data")


## the indirect method, per capita income data from Ragin (2008)
inc <- c(40110, 34400, 25200, 24920, 20060, 17090, 15320, 13680, 11720,
         11290, 10940, 9800, 7470, 4670, 4100, 4070, 3740, 3690, 3590,
         2980, 1000, 650, 450, 110)

cinc <- calibrate(inc, method = "indirect",
        thresholds = "1000, 4000, 5000, 10000, 20000")

plot(inc, cinc, main = "Membership scores in the set of high income", 
     xlab = "Raw data", ylab = "Calibrated data")
     

# calibrating categorical data
set.seed(12345)
values <- sample(1:7, 100, replace = TRUE)

TFR <- calibrate(values, method = "TFR")

table(round(TFR, 3))
}
\keyword{functions}
