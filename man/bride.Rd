\name{bride}
\alias{bride}
\alias{bride.default}
\alias{print.bride}
\alias{summary.bride}
\title{Brier Score Decomposition}
\description{Decomposes the empirical Brier score into reliability, resolution and uncertainty. Two different estimators for the components are provided: The original estimators proposed by Murphy (1974), and the bias-corrected estimators proposed by Ferro and Fricker (2012). Variances of all the components are estimated. Applies only to probabilistic predictions of binary events.}
\usage{
  bride(p, y, n.bins=10, \dots)
}
\arguments{
  \item{p}{ a vector of forecast probabilities. No default. }
  \item{y}{ a vector of binary event indicators. No default. }
  \item{n.bins}{ number of bins used for calculating the calibration function. Default \code{n.bins = 10}.  }
  \item{...}{ possible additional arguments. Not used at the moment.  }
}
\details{
  The values of the forecast probabilities in \code{p} are binned into \code{n.bins} bins of equal length on the unit interval. All probabilities are replaced with their in-bin average. Based on this binning, the calibration function \eqn{P(y=1|p)} is estimated, which is required to estimate the components of the Brier Score decomposition of \code{p}.  }
\value{
An object of class \code{bride}, essentially a list containing:
  \item{p,y}{ the objects of the original request.}
  \item{n.bins}{ number of equidistant, exhaustive bins used to calculate the cross table.}
  \item{rel,res,unc}{ reliability, resolution, uncertainty estimates derived by Murphy (1970).}
  \item{rel2,res2,unc2}{ bias-corrected reliability, resolution, uncertainty estimates derived by Ferro and Fricker (2012).}
  \item{rel.var,res.var,unc.var,rel2.var,res2.var,unc2.var}{ variance estimators derived by Siegert (2013).}
}
\examples{
## number of forecasts and verifications
N <- 100
## produce some forecasts
p <- runif(n=N)
## produce events y[i] that occur with probability p[i]
y <- rbinom(n=N, size=1, prob=p)
## perform brier score decomposition
b <- bride(p=p,y=y)
## print
print(b)
}

\author{Stefan Siegert}
\keyword{forecast verification}
\references{Murphy, AH (1974) \emph{A new vector partition of the probability score}, Journal of Applied Meteorology, 12:595-600\cr\cr
Ferro CAT, Fricker TE (2012) \emph{A bias-corrected decomposition of the Brier Score}, Quarterly Journal of the Royal Meteorological Society, 138(668): 1954-1960\cr\cr
Siegert, S. (2013) \emph{Variance estimation for Brier Score decomposition}, http://arxiv.org/abs/1303.6182 \cr
}


