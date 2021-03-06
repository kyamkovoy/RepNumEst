% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_rts.R
\name{find_rts2}
\alias{find_rts2}
\title{Finds Reproductive Numbers.}
\usage{
find_rts2(cases, si_shape = 1.47, si_rate = 0.04529, k = 100)
}
\arguments{
\item{cases}{List of cases to be used for estimating reproductive numbers.}

\item{si_shape}{Shape parameter for serial interval.}

\item{si_rate}{Rate parameter for serial interval.}

\item{k}{Maximum length of serial interval.}
}
\value{
List of effective reproductive numbers and corresponding week
  numbers.
}
\description{
Using data and specified serial interval, finds the effective reproductive
  number. The default serial interval uses an estimate for the United
  States on the weekly time scale. Only requires a list of cases and not
  the entire dataframe for the data.
}
