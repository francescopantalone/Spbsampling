#' The income of municipalities of "Emilia Romagna".
#'
#' The dataset contains the total income of the municipalities in the region "Emilia Romagna", in Italy. Each municipality is defined by their own ISTAT (Istituto nazionale di statistica, Italy) code, a name and its geographical positions (through a pair of coordinates).
#' For each of them we have the following auxiliary variables: province, number of taxpayers and total income of the municipality.
#'
#' @format A data frame with 334 rows and 7 variables:
#' \describe{
#'   \item{municipality_code}{identification municipality code}
#'   \item{municipality}{name of the municipality}
#'   \item{province}{province of the municipality}
#'   \item{numtaxpay}{number of taxpayers in the municipality}
#'   \item{tot_inc}{average income of the municipality}
#'   \item{x_coord}{coordinate x of the municipality}
#'   \item{y_coord}{coordinate y of the municipality}
#'          }
#'
#' @source
#' The dataset is a rearrangement from the data released by MEF - Dipartimento delle Finanze (Italy).
"income_emilia"

#' Simulated Population 1.
#'
#' The dataset contains a simulated georeferenced population of dimension \eqn{N=1000}.
#' The coordinates are generated in the range \eqn{[0,1]} as a simulated realization of a particular random point pattern: the Neyman-Scott process with Cauchy cluster kernel.
#' The nine values for each unit are generated according to the outcome of a Gaussian stochastic process, with an intensity dependence parameter \eqn{\rho=0.001} (that means low dependence) and with no spatial trend.
#'
#' @format A data frame with  1000 rows and 11 variables:
#' \describe{
#'   \item{x}{coordinate x}
#'   \item{y}{coordinate y}
#'   \item{z11}{first value of the unit}
#'   \item{z12}{second value of the unit}
#'   \item{z13}{third value of the unit}
#'   \item{z14}{fourth value of the unit}
#'   \item{z15}{fifth value of the unit}
#'   \item{z16}{sixth value of the unit}
#'   \item{z17}{seventh value of the unit}
#'   \item{z18}{eighth value of the unit}
#'   \item{z19}{ninth value of the unit}
#'          }
#'
#' @source
#' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
"simul1"

#' Simulated Population 2.
#'
#' The dataset contains a simulated georeferenced population of dimension \eqn{N=1000}.
#' The coordinates are generated in the range \eqn{[0,1]} as a simulated realization of a particular random point pattern: the Neyman-Scott process with Cauchy cluster kernel.
#' The nine values for each unit are generated according to the outcome of a Gaussian stochastic process, with an intensity dependence parameter \eqn{\rho=0.01} (that means medium dependence) and with a spatial trend \eqn{x_{1}+x_{2}+\epsilon}.
#'
#' @format A data frame with  1000 rows and 11 variables:
#' \describe{
#'   \item{x}{coordinate x}
#'   \item{y}{coordinate y}
#'   \item{z21}{first value of the unit}
#'   \item{z22}{second value of the unit}
#'   \item{z23}{third value of the unit}
#'   \item{z24}{fourth value of the unit}
#'   \item{z25}{fifth value of the unit}
#'   \item{z26}{sixth value of the unit}
#'   \item{z27}{seventh value of the unit}
#'   \item{z28}{eighth value of the unit}
#'   \item{z29}{ninth value of the unit}
#'          }
#'
#' @source
#' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
"simul2"

#' Simulated Population 3.
#'
#' The dataset contains a simulated georeferenced population of dimension \eqn{N=1000}.
#' The coordinates are generated in the range \eqn{[0,1]} as a simulated realization of a particular random point pattern: the Neyman-Scott process with Cauchy cluster kernel.
#' The nine values for each unit are generated according to the outcome of a Gaussian stochastic process, with an intensity dependence parameter \eqn{\rho=0.1} (that means high dependence) and with a spatial trend \eqn{x_{1}+x_{2}+\epsilon}.
#'
#' @format A data frame with  1000 rows and 11 variables:
#' \describe{
#'   \item{x}{coordinate x}
#'   \item{y}{coordinate y}
#'   \item{z31}{first value of the unit}
#'   \item{z32}{second value of the unit}
#'   \item{z33}{third value of the unit}
#'   \item{z34}{fourth value of the unit}
#'   \item{z35}{fifth value of the unit}
#'   \item{z36}{sixth value of the unit}
#'   \item{z37}{seventh value of the unit}
#'   \item{z38}{eighth value of the unit}
#'   \item{z39}{ninth value of the unit}
#'         }
#'
#' @source
#' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
"simul3"

#' LUCAS data for the region "Abruzzo", Italy.
#'
#' The land use/cover area frame statistical survey, abbreviated as LUCAS, is a European field survey program funded and executed by Eurostat.
#' Its objective is to set up area frame surveys for the provision of coherent and harmonised statistics on land use and land cover in the European Union (EU).
#' Note that in LUCAS survey the concept of land is extended to inland water areas (lakes, river, coastal areas, etc.) and does not embrace uses below the earth's surface (mine deposits, subways, etc.).
#' The LUCAS survey is a point survey, in particular the basic unit of observation is a circle with a radius of 1.5m (corresponding to an identifiable point on an orthophoto).
#' In the classification there is a clear distinction between land cover and land use: land cover means physical cover ("material") observed at the earth's surface; land use means socio-economic function of the observed earth's surface.
#' For each of both we assign a code to identified which type the point is.
#' Land cover has 8 main categories, which are indicated by letter:
#' \describe{
#'     \item{A}{artificial land}
#'     \item{B}{cropland}
#'     \item{C}{woodland}
#'     \item{D}{shrubland}
#'     \item{E}{grassland}
#'     \item{F}{bareland}
#'     \item{G}{water}
#'     \item{H}{wetlands}
#' }
#' Every main category has subclasses, which are indicated by the combination of the letter of the category and digits. Altogether there are 84 classes.
#' Land use has 14 main categories. It has altogether 33 classes, which are indicated by the combination of the letter “U” and three digits.
#'
#' @format A data frame with 2699 rows and 7 variables:
#' \describe{
#'   \item{id}{identified code for the unit spatial point}
#'   \item{prov}{province}
#'   \item{elev}{elevation of the unit spatial point, meant as the height above or below sea level}
#'   \item{x}{coordinate x}
#'   \item{y}{coordinate y}
#'   \item{lc}{land cover code}
#'   \item{lu}{land use code}
#'         }
#'
#' @source
#' The dataset contains the data from LUCAS 2012 for the region "Abruzzo", Italy.
"lucas_abruzzo"


