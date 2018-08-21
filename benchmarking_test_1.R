rm(list=ls())
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

nitems<-7
nind<-100
sort(alpha_o<-runif(nitems,0.6,1.5))
sort(beta_o<-runif(nitems,-2.4,2.4))
d_o<- -alpha_o*beta_o

summary(theta_o<-rnorm(nind))
temp.b<-matrix(beta_o,ncol=nitems,nrow=nind, byrow =T)
temp.a<-matrix(alpha_o,ncol=nitems,nrow=nind, byrow =T)
temp.t<-matrix(theta_o,ncol=nitems,nrow=nind)
eta<-temp.a*(theta_o-temp.b)
pmod<-exp(eta)/(1+exp(eta))
test<-data.frame(ifelse(runif(nitems*nind)<pmod,1,0))
colnames(test)<-paste("t",1:nitems,sep="")
sort(apply(test,2,function(i)sum(i)*100/length(i)))
test<-data.frame(test,id=1:nind)

nlf_2pl<-function(X,psi,z,theta,nitems){
  return(drop(X%*%psi[1:nitems]+exp(X%*%psi[(nitems+1):(2*nitems)])*(z%*%theta)))
}
inv_link_2pl<-function(x){1/(1+exp(-x))}

grad.fix_2pl<-function(X,psi,z,theta,nitems){
  return(cbind(X,drop(exp(X%*%psi[(nitems+1):(2*nitems)])*(z%*%theta))*X))
}

grad.mix_2pl<-function(X,psi,z,theta,nitems){
  return(drop(exp(X%*%psi[(nitems+1):(2*nitems)]))*z)
}
hprim_2pl<-function(x){1/(x*(1-x))}



log_lik<-function(y,X,psi,z,theta,nitems){
  eta_y<-nlf(X,psi,z,theta,nitems)
  mu<-inv_link(eta_y)
  return(sum(dbinom(y,1,mu,log=T)))
}


dat<-test;idname<-"id";psi0<-c(d_o,log(alpha_o));theta0<-theta_o;B<-NULL;b<-NULL;
grad.fix<-grad.fix_2pl;grad.mix<-grad.mix_2pl;nlf<-nlf_2pl;inv_link<-inv_link_2pl;hprim<-hprim_2pl;Nc<-1000

vectorized<-uirt_DM(dat,idname)
Y<-vectorized$Y
X<-vectorized$X
Z<-vectorized$Z
if(is.null(B)){
  B<-diag(c(rep(9,dim(X)[2]),rep(1,dim(X)[2])))
}
if(is.null(b)){
  b<-c(rep(0,length(psi0)))
}
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
Dev[1]<- -2*log_lik(Y,X,psi0,Z,theta0,nitems)
########################
#####gibbs
########################
i<-2
#########
##mh for bta
########
X_hat<-grad.fix(X,psi.c[i-1,],Z,theta.c[i-1,],nitems)
eta<-nlf(X,psi.c[i-1,],Z,theta.c[i-1,],nitems)
mu<-inv_link(eta)
hprim_hinv<-hprim(mu)
y1_hat<-drop(X_hat%*%c(psi.c[i-1,]))+hprim_hinv*(Y-mu)
Vy1_hat<-mu*(1-mu)*hprim_hinv^2

B.<-solve(B.inv+t(X_hat)%*%(X_hat/Vy1_hat))
b.<-drop(B.%*%(B.inv_b+t(X_hat)%*%(y1_hat/Vy1_hat)))
psi.p<-drop(mvtnorm::rmvnorm(1,b.,B.))#####proposal

X_hat.p<-grad.fix(X,psi.p,Z,theta.c[i-1,],nitems)
eta.p<-nlf(X,psi.p,Z,theta.c[i-1,],nitems)
mu.p<-inv_link(eta.p)
hprim_hinv.p<-hprim(mu.p)
y1_hat.p<-drop(X_hat.p%*%psi.p)+hprim_hinv.p*(Y-mu.p)
Vy1_hat.p<-mu.p*(1-mu.p)*hprim_hinv.p^2

B.p<-solve(B.inv+t(X_hat.p)%*%(X_hat.p/Vy1_hat.p))
b.p<-drop(B.p%*%(B.inv_b+t(X_hat.p)%*%(y1_hat.p/Vy1_hat.p)))

library(parallel)
num_cores<-detectCores()-1
cl<-makeCluster(num_cores)
sol_1<-function(y){sapply(y,FUN=function(j){mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)
  -mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)})}
sol_2<-function(y){purrr::map_dbl(y,function(j){mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)
    -mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)})}
sol_3<-function(y){vapply(y,FUN=function(j){mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)
  -mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)},FUN.VALUE=numeric(1))}
sol_4<-function(y,psi.p,nitems,psi.c,b,B,i){
  num_cores<-detectCores()-1
  cl<-makeCluster(num_cores)
  temp<-parSapply(cl,y,FUN=function(j,psi.p,nitems,psi.c,b,B,i){mvtnorm::dmvnorm(psi.p[c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)
    -mvtnorm::dmvnorm(psi.c[i-1,c(j,nitems+j)],b[c(j,nitems+j)],B[c(j,nitems+j),c(j,nitems+j)],log=T)},psi.p=psi.p,nitems=nitems,psi.c=psi.c,b=b,B=B,i=i)
  stopCluster(cl)
  return(temp)
  }

library(microbenchmark)
sol_1(1:nitems)
sol_2(1:nitems)
sol_3(1:nitems)
sol_4(1:nitems,psi.p,nitems,psi.c,b,B,i)
microbenchmark(sol_1(1:nitems),sol_2(1:nitems),sol_3(1:nitems),sol_4(1:nitems,psi.p,nitems,psi.c,b,B,i),times =100)
