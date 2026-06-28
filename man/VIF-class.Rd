\name{VIF-class}

\docType{class}

\alias{VIF-class}

\alias{show,VIF-method}


\title{ VIF class}

\description{

An object of the \code{VIF} class contains information about collinearity in relavant variables. The object can be created with the following functions:  \code{\link[usdm]{vifcor}} and \code{\link[usdm]{vifstep}}.

}

\section{Slots}{
Slots for VIF object
	\describe{
    \item{\code{variables}:}{Character vector of all input variable names}
    \item{\code{excluded}:}{character vector of excluded variable names}
    \item{\code{corMatrix}:}{a correlation \code{matrix} among the retained variables}
    \item{\code{results}:}{\code{data.frame} including VIF values for the remained (not excluded) variables}
    \item{\code{exclusionLog}:}{\code{data.frame} recording every exclusion step with columns: \code{step} (integer, order of exclusion), \code{excluded} (variable dropped), \code{correlated_with} (the variable it was paired/compared with), \code{correlation} (absolute correlation between the pair), \code{vif_excluded} (VIF of the dropped variable at that step), \code{vif_retained} (VIF of the retained variable; \code{NA} for \code{vifstep})}
    \item{\code{chains}:}{named \code{list} giving the transitive exclusion chain for each final retained variable that "won" at least one pairing. Each element is a character vector ordered from the most recently excluded variable to the most distally excluded one. E.g., \code{list(Bio7 = c("Bio10", "Bio5"))} means the chain reads \code{Bio7 <- Bio10 <- Bio5}.}
    }
}


\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}

\examples{
showClass("VIF")
}

\keyword{classes}
\keyword{spatial}
