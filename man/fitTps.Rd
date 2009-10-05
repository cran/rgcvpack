\name{fitTps}
\alias{fitTps}
\title{ Fitting Thin Plate Smoothing Spline }
\description{
  Fit thin plate splines of any order with user specified knots
}
\usage{
fitTps(x, y, m = 2, knots = NULL, scale.type = "range", method = "v",
    lambda = NULL, cost = 1, nstep.cv = 80, verbose = FALSE, tau = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the design data points }
  \item{y}{ the observation vector }
  \item{m}{ the order of the spline }
  \item{knots}{ the placement the thin plate spline basis }
  \item{scale.type}{ "range" (default), the x and knots
    will be rescaled with respect to x; "none", nothing is done on
    x and knots }
  \item{method}{ "v", GCV is used for choosing lambda; "d", user
    specified lambda }
  \item{lambda}{ only used when method="d" }
  \item{cost}{ the fudge factor for inflating the model degrees of
    freedom, default to be 1 }
  \item{nstep.cv}{ the number of initial steps for GCV grid search }
  \item{verbose}{ whether some computational details should be outputed }
  \item{tau}{ the truncation ratio used in SVD when knots is specified by
    the user, some possible values are 1, 10, 100, ... }
}
\details{
  The minimization problem for this function is
  \deqn{ \sum_{i=1}^{n}{(y_i - f(x_i))^2} + \lambda*J_m(f), }{%
         sum((y_i - f(x_i))^2) + lambda*J_m(f), }
  where \eqn{J_m(.)} is the m-the order thin plate spline penalty functional.
  
  If scale.type="range", each column of x is rescaled to [0 1] in the
  following way x' = (x - min(x))/range(x), and the knots is rescaled
  w.r.t. min(x) and range(x) in the same way.
  
  When the cost argument is used, the GCV score is computed as
  \deqn{ \mbox{GCV}(\lambda) =
    \frac{n*\mbox{RSS}(\lambda)}{(n - cost*\mbox{tr}(A))^2}. }{%
      GCV(lambda) = n*RSS(\lambda)/(n - cost*tr(A))^2. }
   
}
\value{
  A Tps object of the following components
  \item{x }{same as input}
  \item{y }{same as input}
  \item{m }{same as input}
  \item{knots}{same as input}
  \item{scale.type}{same as input}
  \item{method}{same as input}
  \item{lambda}{same as input}
  \item{cost}{same as input}
  \item{nstep.cv}{same as input}
  \item{tau}{same as input}
  \item{df}{model degrees of freedom}
  \item{gcv}{gcv score of the model adjusted for the fudge factor}
  \item{xs}{scaled design points}
  \item{ks}{scaled knots design}
  \item{c}{coefficient c}
  \item{d}{coefficient d}
  \item{yhat}{predicted values at the data points}
  \item{svals}{singular values of the matrix decomposition}
  \item{gcv.grid}{gcv grid table, number of rows=nstep.cv}
  \item{call}{the call to this function}
}
\references{
  D. Bates, M. Lindstrom, G. Wahba, B. Yandell (1987), GCVPACK -- routines
  for generalized cross-validation. Commun. Statist.-Simula., 16(1), 263-297. 
}
\author{ Xianhong Xie }
\note{
  This function uses GCVPACK fortran code with some addition and
  modification by the author.  
}
\seealso{ \code{\link{predict.Tps}} }
\examples{

#define the test function
f <- function(x, y) { .75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
                      .75*exp(-((9*x+1)^2/49 + (9*y+1)^2/10)) +
                      .50*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
                      .20*exp(-((9*x-4)^2 + (9*y-7)^2)) }

#generate a data set with the test function
set.seed(200)
N <- 13; xr <- (2*(1:N) - 1)/(2*N); yr <- xr
zr <- outer(xr, yr, f); zrmax <- max(abs(zr))

noise <- rnorm(N^2, 0, 0.07*zrmax)
zr <- zr + noise #this is the noisy data we will use

#convert the data into column form
xc <- rep(xr, N)
yc <- rep(yr, rep(N,N))
zc <- as.vector(zr)

#fit the thin plate spline with all the data points as knots
tpsfit1 <- fitTps(cbind(xc,yc), zc, m=2, scale.type="none")
persp(xr, yr, matrix(predict(tpsfit1),N,N), theta=130, phi=20,
      expand=0.45, xlab="x1", ylab="x2", zlab="y", xlim=c(0,1),
      ylim=c(0,1),zlim=range(zc), ticktype="detailed", scale=FALSE,
      main="GCV Smooth I")

#fit the thin plate spline with subset of data points as knots
grid.list  <- list(xc=seq(2/13,11/13,len=10),
                   yc=seq(2/13,11/13,len=10))
knots.grid <- expand.grid(grid.list)

tpsfit2 <- fitTps(cbind(xc,yc), zc, m=2, knots=knots.grid)
persp(xr, yr, matrix(predict(tpsfit2),N,N), theta=130, phi=20,
      expand=0.45, xlab="x1", ylab="x2", zlab="y", xlim=c(0,1),
      ylim=c(0,1),zlim=range(zc), ticktype="detailed", scale=FALSE,
      main="GCV Smooth II")
}
\keyword{smooth}

