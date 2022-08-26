\name{AutoWeatherIndices}
\alias{AutoWeatherIndices}
\title{Calculating Weather Indices Automatically Using Weather Variables}
\usage{
AutoWeatherIndices(x,nw)
}
\arguments{
  \item{x}{
Matrix of weather variables recorded weekly, along with the response variable (eg. Yiled) recorded over years in the first column.
}
  \item{nw}{
Number of weeks for which data has been recorded.
}
}
\description{
Weather indices are computed from weather variables. The users can input any number of weather varibles recorded over any number of weeks. This package has no restriction on the number of weeks and wetaher variables to be taken as inout.The details of the method can be seen (i)'Joint effects of weather variables on rice yields' by  R. Agrawal, R. C. Jain  and M. P. Jha in Mausam, vol. 34, pp. 189-194, 1983,<doi:10.54302/mausam.v34i2.239>, (ii)'Improved weather indices based Bayesian regression model for forecasting crop yield' by  M. Yeasin, K. N. Singh, A. Lama and B. Gurung in Mausam, vol. 72, pp.879-886, 2021,<doi:10.54302/mausam.v72i4.670>.
}
\value{
  \item{return_list}{The matrix of weather indices along with two correlation matrices between the weather variables and the response variable.}
}
\author{
Achal Lama, Kamlesh N Singh and Bishal Gurung
}
\references{
Agrawal, R. Jain, R.C. and Jha, M.P.(1983). Joint effects of weather variables on rice yields.Mausam,34(2):189-194. (doi:10.54302/mausam.v34i2.239).

Yeasin, Md. Singh, K.N. Lama, A. and Gurung, B. (2021). Improved weather indices based Bayesian regression model for forecasting crop yield. Mausam, 72(4):879-886. (doi:10.54302/mausam.v72i4.670).

}

\examples{
x=matrix(c(2011,31.66,33.23,32.94,32.91,33.37,30.59,25.51,26.30,26.44,
26.19,28.10,26.34,21.64,22.14,23.23,1978,32.83,32.69,30.10,30.71,29.79,28.87,
28.81,27.69,25.01,26.34,26.31,25.31,26.53,27.59,22.31)
,nrow=2,ncol=16,byrow=TRUE)
AutoWeatherIndices(x,5)
}
