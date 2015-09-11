\name{qpcrImpute}
\alias{qpcrImpute}
\title{Impute Non-detects in qPCR data}
\description{
  This function models the missing data mechanism and uses an EM
  algorithm to impute the non-detect values in qPCR data. 
}
\usage{
qpcrImpute(object, dj=NULL, pyfit=NULL, groupVars=NULL, tol=1, iterMax=100, 
          outform=NULL, formula=NULL, numsam=5, 
          vary_fit=TRUE, vary_model=TRUE, add_noise=TRUE)
}
\arguments{
  \item{object}{a qPCRset}
  \item{dj}{normalization values. If NULL, features with "control"
    in featureType(object) are used to normalize the data. If no control
    features are found, the data are not normalized.}
  \item{pyfit}{initial estimate of the relationship between the
    probability of a non-detect and average expression. If NULL, this
    relationship is estimated from the data.} 
  \item{groupVars}{which columns in pData(object) should be used to
    determine replicate samples. If NULL, all columns are used.} 
  \item{tol}{likelihood convergence criterion of the EM algorithm.}
  \item{iterMax}{maximimum number of iterations of the EM algorithm.}
  
  \item{outform}{the form of the output requested.If "Single" performes a single imputation of missing values. If "Param" returnes estimated model parameters: mean and variance. If "Multy" performes a multiple imputation of missing values, and creats multiple data sets with imputed values.}
  \item{formula}{specifies the model.}
  \item{numsam}{number of the datasets to be created if outform="Multy". The default value is 5.}
  \item{vary_fit}{if outform="Multy", includes the model uncertainty due to the logit of the probability of being missing. The default value is "TRUE".}
  \item{vary_model}{if outform="Multy", includes the model uncertainty due to the estimating mean of the data. The default value is "TRUE".}
  \item{vary_noise}{if outform="Multy", introduses the variance component due to the random noise. The default value is "TRUE".}
}
\value{
  The function returns a qPCRset object with non-detects replaced by
  their imputed values. 
}
\author{Matthew N. McCall}
\examples{
  data(sagmb2011)
  tst <- qpcrImpute(sagmb2011,groupVars="sampleType")
}
\keyword{manip}

