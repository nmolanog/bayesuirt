#' nlf_2pl Function
#'
#' this functions is the non-linear function associated with the 2pl uirt model as estated in molano thesis.
#' @param X  design matrix indexing each observation to corresponding item
#' @param psi vector of item parameters
#' @param z design matrix indexing each observation to corresponding individual
#' @param theta vector of latent traits (i.e. random effects associated to individual level)
#' @param nitems integer especifying the number of items. It is used to separate parameters as psi[1:nitems] and psi[(nitems+1):(2*nitems)]
#' @return add style to sheet in wb
#' @keywords a vector with reals
#' @export
#' @examples
#' ##none for now

nlf_2pl<-function(X,psi,z,theta,nitems){
  return(drop(X%*%psi[1:nitems]+exp(X%*%psi[(nitems+1):(2*nitems)])*(z%*%theta)))
}
