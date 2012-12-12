\name{KappaV}
\alias{KappaV}
\title{
Calculates KappaV
}
\description{
This function calculates the vectorial Kappa, the associated standard deviation and also returns the confusion matrix calculated between two vectorial landscapes.
}

\usage{
KappaV(shp1.path, shp2.path,
       shp1.fieldID = "ID", shp2.fieldID = shp1.fieldID,
       shp1.fieldOS = "OS", shp2.fieldOS = shp1.fieldOS,
       plot = FALSE)
}

\arguments{
  \item{shp1.path}{\code{character}. The path to the first shape file.}
  \item{shp2.path}{\code{character}. The path to the second shape file.}
  \item{shp1.fieldID}{\code{character}. The column name in the .dbf file to indicate the polygons IDs.}
  \item{shp2.fieldID}{\code{character}. The column name in the .dbf file to indicate the polygons IDs.}
  \item{shp1.fieldOS}{\code{character}. The column name in the .dbf file to indicate the nominal variable of interest.}
  \item{shp2.fieldOS}{\code{character}. The column name in the .dbf file to indicate the nominal variable of interest.}
  \item{plot}{\code{logical}. Whether to plot the two landscapes.}
}

\details{
If nto specified blabla
}

\value{
A list with two components:
\item{confusion.matrix }{Like it say, the confusion matrix with the corresponding areas.}
\item{kappa.v }{a numeric with two values: Kappa (the value of vectorial Kappa) and Kappa.sd (the associated standard deviation)}
}

\references{
Soon.
}


\seealso{
\link{Kappa} in the \code{PresenceAbsence} package that handles the Kappa and its SD calculation from the confusion matrix.
}

\examples{
# have a look to the package's vignette: http://www.vincentbonhomme.fr/KappaV
}

