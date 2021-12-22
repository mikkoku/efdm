#' Example dataset.
#'
#' A list containing the necessary datasets for EFDM forest scenario example.
#'
#' @format A list of data frames:
#' \describe{
#'   \item{actprob}{Activity probabilities}
#'   \item{noman_pairs}{Pair data for growth without management activities}
#'   \item{thin_pairs}{Pair data for thinning}
#'   \item{initial_state}{Initial forest state}
#'   \item{drain_coef}{Coefficient to transform harvested areas into harvest accumulation by timber assortments}
#'   \item{vol_coef}{Coefficients to transform volume classes into volumes m3/ha}
#'   \item{income_coef}{Coefficients to transform harvest accumulation into income}
#' }
"example"

#' Finnish bio-geographical regions
#'
#' A low resolution copy of the finnish bio-geographical regions.
#' The original shapefile was provided by Finnish Environment Institute.
#' See \url{https://ckan.ymparisto.fi/fi/dataset/metsakasvillisuusvyohykkeet}.
#'
#' @format An sf-object:
#' \describe{
#'   \item{region}{Factor with three levels: North, Middle, South}
#'   \item{geometry}{Polygon, each region is composed of multiple polygons}
#' }
"MetsaKasvVyoh"
