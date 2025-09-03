#' Random seed generator for analysis reproducibility
#'
#' This function generates 8-digit random seeds
#' @param n.seeds Number of random seeds
#' @return Random seeds
#' @export
#' @examples
#' scprep_seeds(n.seeds=1000)
#

#
scprep_seeds <- function(
	n.seeds=n.seeds) {
	#
	seeds <- round(sample(1:1e8, n.seeds, replace=FALSE))
	#
	return(seeds)
	#
}