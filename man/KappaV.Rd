\docType{package}
\name{KappaV}
\alias{KappaV}
\alias{KappaV-package}
\title{Calculates vectorial kappa between two shapes}
\usage{
KappaV(shp1.path, shp2.path, shp1.fieldID = "ID",
  shp2.fieldID = shp1.fieldID, shp1.fieldOS = "OS",
  shp2.fieldOS = shp1.fieldOS, plot = FALSE)
}
\arguments{
  \item{shp1.path}{\code{character}. The path to the first
  shape file.}

  \item{shp2.path}{\code{character}. The path to the second
  shape file.}

  \item{shp1.fieldID}{\code{character}. The column name in
  the .dbf file to indicate the polygons IDs.}

  \item{shp2.fieldID}{\code{character}. The column name in
  the .dbf file to indicate the polygons IDs.}

  \item{shp1.fieldOS}{\code{character}. The column name in
  the .dbf file to indicate the nominal variable of
  interest.}

  \item{shp2.fieldOS}{\code{character}. The column name in
  the .dbf file to indicate the nominal variable of
  interest.}

  \item{plot}{\code{logical}. Whether to plot the two
  landscapes.}
}
\value{
A list with two components: \code{confusion.matrix}: the
confusion matrix with the corresponding areas. and
\code{$kappa.v}: a numeric with two values: \code{$Kappa}
(the value of vectorial Kappa) and \code{$Kappa.sd }(the
associated standard deviation)
}
\description{
This function calculates the vectorial Kappa, the
associated standard deviation and also returns the
confusion matrix calculated between two vectorial
landscapes.

Calculates vectorial Kappa, i.e. congruence between two
vectorial landscapes.
}
\details{
If not specified the default parameters are
\code{shp1.fieldID = "ID"}, \code{shp2.fieldID =
shp1.fieldID}, \code{shp1.fieldOS =
"OS"},\code{shp2.fieldOS = shp1.fieldOS} and \code{plot =
FALSE}
}
\examples{
# Have a look to package's vignette above.
}
\references{
Paper submitted.

Paper submitted.
}
\seealso{
\link{Kappa} in the \code{PresenceAbsence} package that
handles the Kappa and its SD calculation from the confusion
matrix.

Have a look to the package's vignette:
\url{http://www.vincentbonhomme.fr/KappaV}

KappaV' homepage :
\url{http://www.vincentbonhomme.fr/KappaV} with tutorials
and hotline.

KappaV' GitHub repo :
\url{https://github.com/vbonhomme/KappaV} to contribute,
among other things.
}
\keyword{Abtract}

