#' uirt_DM Function
#'
#' this functions converts a test dataframe to a vectoriced form following molano's thesis.
#' @param dat a data frame containing test info. rows for individuals and columns for items. aditionally, a variable which identify individuals must be supied.
#' @param idname atomic character giving the name of the id variable
#' @return a list with the following objects:
#' \itemize{
#'   \item long.test: a data frame with id, item and y (0 or 1). Long format for the test.
#'   \item Y: Vector containing individual answers (0 or 1).
#'   \item X: desing matrix indexing each value of Y to the corresponding item.
#'   \item Z: desing matrix indexing each value of Y to the corresponding individual.
#' }
#' @keywords a vector with reals
#' @export
#' @examples
#' data("test_data")
#' temp<-uirt_DM(test_data,"id")

uirt_DM<-function(dat,idname){
  nitems<-ncol(test_data)-1
  nind<-nrow(test_data)
  long.test<-tidyr::gather(test_data,"item","y",-idname)
  Y<-long.test$y
  X<-model.matrix(y ~-1+item,data=long.test)
  Z<-model.matrix(formula(paste0("y ~-1+as.factor(",idname,")")),data=long.test)
  res<-list(long.test=long.test,Y=Y,X=X,Z=Z)
  return(res)
}



#' nlf_2pl Function
#'
#' this functions is the non-linear function associated with the 2pl uirt model as estated in molano thesis.
#' @param X  design matrix indexing each observation to corresponding item
#' @param psi vector of item parameters
#' @param z design matrix indexing each observation to corresponding individual
#' @param theta vector of latent traits (i.e. random effects associated to individual level)
#' @param nitems integer especifying the number of items. It is used to separate parameters as psi[1:nitems] and psi[(nitems+1):(2*nitems)]
#' @return a double vector
#' @keywords a vector with reals
#' @export
#' @examples
#' data("test_data")
#' temp<-uirt_DM(test_data,"id")
#' nitems<-ncol(temp$X)
#' nind<-ncol(temp$Z)
#' alpha<-runif(nitems,0.6,1.5) ##simulation of discrimination parameters
#' beta<-runif(nitems,-2.4,2.4) ##simulation of dificulty parameters
#' theta<-rnorm(nind) ##simulation of random effects
#' psi<-c(-alpha*beta,log(alpha)) ###vector of parameters based on alpha and beta
#' res<-nlf_2pl(temp$X,psi,temp$Z,theta,nitems) ###nonlinear function evaluated at given values

nlf_2pl<-function(X,psi,z,theta,nitems){
  return(drop(X%*%psi[1:nitems]+exp(X%*%psi[(nitems+1):(2*nitems)])*(z%*%theta)))
}
