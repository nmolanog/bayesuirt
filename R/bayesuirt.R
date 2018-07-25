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
#' @keywords a double vector
#' @export
#' @examples
#' data("test_data")
#' temp<-uirt_DM(test_data,"id")

uirt_DM<-function(dat,idname){
  nitems<-ncol(dat)-1
  nind<-nrow(dat)
  long.test<-tidyr::gather(dat,"item","y",-idname)
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

#' link_2pl Function
#'
#' Inverse logit function, used as link in the 2pl model
#' @param x a double vector
#' @return a double vector
#' @keywords logit
#' @export
#' @examples
#' link_2pl(rnorm(10))
link_2pl<-function(x){1/(1+exp(-x))}


#' grad.fix_2pl Function
#'
#' derivative of nlf_2pl with respect to psi
#' @param X  design matrix indexing each observation to corresponding item
#' @param psi vector of item parameters
#' @param z design matrix indexing each observation to corresponding individual
#' @param theta vector of latent traits (i.e. random effects associated to individual level)
#' @param nitems integer especifying the number of items. It is used to separate parameters as psi[1:nitems] and psi[(nitems+1):(2*nitems)]
#' @return a double vector
#' @keywords derivative
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
#' res<-grad.fix_2pl(temp$X,psi,temp$Z,theta,nitems) ###derivative of nlf_2pl with respect to psi evaluated at given values

grad.fix_2pl<-function(X,psi,z,theta,nitems){
  return(cbind(X,drop(exp(X%*%psi[(nitems+1):(2*nitems)])*(z%*%theta))*X))
}

#' grad.mix_2pl Function
#'
#' derivative of nlf_2pl with respect to theta
#' @param X  design matrix indexing each observation to corresponding item
#' @param psi vector of item parameters
#' @param z design matrix indexing each observation to corresponding individual
#' @param nitems integer especifying the number of items. It is used to separate parameters as psi[1:nitems] and psi[(nitems+1):(2*nitems)]
#' @return a double vector
#' @keywords derivative
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
#' res<-grad.mix_2pl(temp$X,psi,temp$Z,nitems) ###derivative of nlf_2pl with respect to theta evaluated at given values

grad.mix_2pl<-function(X,psi,z,nitems){
  return(drop(exp(X%*%psi[(nitems+1):(2*nitems)]))*z)
}


#' log_lik_2pl Function
#'
#' log likelihood function for 2pl
#' @param y  vectorised test
#' @param X  design matrix indexing each observation to corresponding item
#' @param psi vector of item parameters
#' @param z design matrix indexing each observation to corresponding individual
#' @param theta vector of latent traits (i.e. random effects associated to individual level)
#' @param nitems integer especifying the number of items. It is used to separate parameters as psi[1:nitems] and psi[(nitems+1):(2*nitems)]
#' @return a double atomic vector
#' @keywords log likelihood
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
#' res<-log_lik_2pl(temp$Y,temp$X,psi,temp$Z,theta,nitems) ###2pl log likelihood evaluated at given values

log_lik_2pl<-function(y,X,psi,z,theta,nitems){
  eta_y<-nlf_2pl(X,psi,z,theta,nitems)
  mu<-link_2pl(eta_y)
  return(sum(dbinom(y,1,mu,log=T)))
}


#' mh_gibbs_2pl Function
#'
#' this functions implements the metropolis-hastings within gibbs algorithm proposed in molano thesis.
#' @param dat a data frame containing test info. rows for individuals and columns for items. aditionally, a variable which identify individuals must be supied.
#' @param idname atomic character giving the name of the id variable in dat.
#' @param psi0 a vector of initial values for parameter psi0. See nlf_2pl documentation.
#' @param theta0  vector of initial values for latent traits. See nlf_2pl documentation.
#' @param B Prior covariance matrix for psi.
#' @param b Prior expected values for psi.
#' @param Nc number of mcmc iterations.
#' @return a list with the following objects:
#' \itemize{
#'   \item fixp: a matrix where each row corresponds to an mcmc simulation of psi parameters.
#'   \item theta: a matrix where each row corresponds to an mcmc simulation of latent traits.
#'   \item Dev: a vector where deviance is calculated for each mcmc iteration.
#'   \item B: Prior covariance matrix used for psi.
#'   \item b: Prior Expected values used for psi.
#' }
#' @keywords a double vector
#' @export
#' @examples
#' nitems<-7
#'nind<-100
#'sort(alpha_o<-runif(nitems,0.6,1.5))
#'sort(beta_o<-runif(nitems,-2.4,2.4))
#'d_o<- -alpha_o*beta_o
#'
#'summary(theta_o<-rnorm(nind))
#'temp.b<-matrix(beta_o,ncol=nitems,nrow=nind, byrow =T)
#'temp.a<-matrix(alpha_o,ncol=nitems,nrow=nind, byrow =T)
#'temp.t<-matrix(theta_o,ncol=nitems,nrow=nind)
#'eta<-temp.a*(theta_o-temp.b)
#'pmod<-exp(eta)/(1+exp(eta))
#'test<-data.frame(ifelse(runif(nitems*nind)<pmod,1,0))
#'colnames(test)<-paste("t",1:nitems,sep="")
#'sort(apply(test,2,function(i)sum(i)*100/length(i)))
#'test<-data.frame(test,id=1:nind)

#'bres<-mh_gibbs_2pl(test,"id",psi0=c(d_o,log(alpha_o)),theta_o,Nc=1000)

mh_gibbs_2pl<-function(dat,idname,psi0,theta0,B=diag(c(rep(9,dim(X)[2]),rep(1,dim(X)[2]))),b=c(rep(0,length(psi0))),Nc=10000){
  vectorized<-uirt_DM(dat,idname)
  Y<-vectorized$Y
  X<-vectorized$X
  Z<-vectorized$Z
  id<-vectorized$long.test[,idname]
  itm<-vectorized$long.test[,"item"]
  rm(list=c("vectorized"))
  nind<-dim(Z)[2]
  nitems<-dim(X)[2]
  B.inv<-solve(B)
  B.inv_b<-B.inv%*%b
  G<-diag(1,nind)
  #####################################
  ####### objects to store chains
  #####################################
  psi.c<-matrix(NA,ncol=length(psi0),nrow=Nc)
  theta.c<-matrix(NA,ncol=nind,nrow=Nc)
  Dev<-rep(NA,Nc)
  #########################################
  ######### storing inits
  #########################################
  psi.c[1,]<-psi0
  theta.c[1,]<-theta0
  Dev[1]<- -2*log_lik_2pl(Y,X,psi0,Z,theta0,nitems)
  ########################
  #####gibbs
  ########################
  for(i in 2:Nc){
    #########
    ##mh for bta
    ########
    X_hat<-grad.fix_2pl(X,psi.c[i-1,],Z,theta.c[i-1,],nitems)
    eta<-nlf_2pl(X,psi.c[i-1,],Z,theta.c[i-1,],nitems)
    mu<-link_2pl(eta)
    Vy<-((1+exp(eta))^2)*exp(-eta)
    y1_hat<-drop(X_hat%*%c(psi.c[i-1,]))+Vy*(Y-mu)

    B.<-solve(B.inv+t(X_hat)%*%(X_hat/Vy))
    b.<-drop(B.%*%(B.inv_b+t(X_hat)%*%(y1_hat/Vy)))
    psi.p<-drop(mvtnorm::rmvnorm(1,b.,B.))#####proposal

    X_hat.p<-grad.fix_2pl(X,psi.p,Z,theta.c[i-1,],nitems)
    eta.p<-nlf_2pl(X,psi.p,Z,theta.c[i-1,],nitems)
    mu.p<-link_2pl(eta.p)
    Vy.p<-((1+exp(eta.p))^2)*exp(-eta.p)
    y1_hat.p<-drop(X_hat.p%*%psi.p)+Vy.p*(Y-mu.p)

    B.p<-solve(B.inv+t(X_hat.p)%*%(X_hat.p/Vy.p))
    b.p<-drop(B.p%*%(B.inv_b+t(X_hat.p)%*%(y1_hat.p/Vy.p)))

    r1<-exp(sapply(1:nitems,FUN=function(j){mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)
      -mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)}) +
        tapply(dbinom(Y,1,mu.p,log =T)-dbinom(Y,1,mu,log =T),itm,sum)+
        sapply(1:nitems,FUN=function(j){mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b.p[c(j,nitems+j)],B.p[c(j,nitems+j),c(j,nitems+j)],log=T)
          -mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b.[c(j,nitems+j)],B.[c(j,nitems+j),c(j,nitems+j)],log=T)}))

    acpt.psi <- ifelse(runif(nitems)<r1,TRUE,FALSE)
    psi.c[i,1:nitems]<-sapply(1:nitems,function(x){if(acpt.psi[x]){psi.p[x]}else{psi.c[i-1,x]}})
    psi.c[i,(nitems+1):(2*nitems)]<-sapply(1:nitems,function(x){if(acpt.psi[x]){psi.p[nitems+x]}else{psi.c[i-1,nitems+x]}})

    #########
    ##mh for theta
    ########
    Z_hat<-grad.mix_2pl(X,psi.c[i,],Z,nitems)
    eta_t<-nlf_2pl(X,psi.c[i,],Z,theta.c[i-1,],nitems)
    mu_t<-link_2pl(eta_t)
    Vy_t<-((1+exp(eta_t))^2)*exp(-eta_t)
    yt_hat<-drop(Z_hat%*%theta.c[i-1,]+Vy_t*(Y-mu_t))

    G.<-1/diag((G+t(Z_hat)%*%(Z_hat/Vy_t)))
    g.<-drop(G.*(t(Z_hat)%*%(yt_hat/Vy_t)))
    theta.p<-drop(rnorm(nind,g.,sqrt(G.)))#####proposal

    eta_t.p<-nlf_2pl(X,psi.c[i,],Z,theta.p,nitems)
    mu_t.p<-link_2pl(eta_t.p)
    Vy_t.p<-((1+exp(eta_t.p))^2)*exp(-eta_t.p)
    yt_hat.p<-drop(Z_hat%*%theta.p+Vy_t.p*(Y-mu_t.p))

    G.p<-1/(diag(G+t(Z_hat)%*%(Z_hat/Vy_t.p)))
    g.p<-drop(G.p*(t(Z_hat)%*%(yt_hat.p/Vy_t.p)))

    r2<-exp(dnorm(theta.p,0,sqrt(diag(G)),log=T)-dnorm(theta.c[i-1,],0,sqrt(diag(G)),log=T)+
              tapply(dbinom(Y,1,mu_t.p,log =T)-dbinom(Y,1,mu_t,log =T),id,sum)+
              dnorm(theta.c[i-1,],g.p,sqrt(G.p),log=T)-dnorm(theta.p,g.,sqrt(G.),log=T))
    acpt <- ifelse(runif(nind)<r2,TRUE,FALSE)
    theta.c[i,]<-sapply(1:nind,function(x){if(acpt[x]){theta.p[x]}else{theta.c[i-1,x]}})
    Dev[i]<- -2*log_lik_2pl(Y,X,psi.c[i,],Z,theta.c[i,],nitems)
    cat(i," 2pl: accpt.psi:",sum(acpt.psi)/nitems,"--"," accpt.theta:",sum(acpt)/nind,"\n")
  }
  return(list(fixp=psi.c,theta=theta.c,Dev=Dev,B=B,b=b))
}

