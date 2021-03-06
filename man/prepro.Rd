% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{prepro}
\alias{prepro}
\title{Preprocess CGM data.}
\usage{
prepro(
  inputdir = "",
  outputdir = "",
  outlierdet = TRUE,
  interval = 15,
  imputation = FALSE,
  immethod = "linear",
  maxgap = 60,
  compeleteday = TRUE,
  removeday = FALSE,
  device = 0,
  transunits = FALSE,
  removeflday = TRUE
)
}
\arguments{
\item{inputdir}{Path of input directory containing CGM files.}

\item{outputdir}{Path of output directory where preprocessed CGM data will be stored.}

\item{outlierdet}{Logical. If TRUE the outliers will be detected.}

\item{interval}{The interval of CGM data.}

\item{imputation}{Logical. If TRUE the missing data will be imputed.}

\item{maxgap}{If the missing gap is greater than max gap, the missing gap will be ignore when imputed.}

\item{compeleteday}{Logical. If TRUE the day with missing data will be filtered.}

\item{removeday}{Logical. If TRUE the day with missing gap greater than maxgap will be filtered.}

\item{device}{Device type: 0 (manual format); 1 (Abbott libre freestyle); 2 (Medtronic ipro2); 3 (Dexcom G6), default 0.}

\item{transunits}{Logical. If TURE the glucose values will be divided by 18.}

\item{removeflday}{Logical. If TRUE the data of first and last day will be filter.}

\item{imethod}{The method that will be used to impute missing data.}
}
\description{
Preprocess CGM data.
}
