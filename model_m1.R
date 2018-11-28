rm(list=ls())
library(tidyverse)
library(bayesuirt)
library(pacman)
p_load(car)
p_load(coda)
###Simulate a data set with model m1 of Yuli
nitems<-7
nind<-100
sort(alpha_o<-rnorm(nitems))
sort(beta_o<-rnorm(nitems))
summary(theta_o<-rnorm(nind))

temp.b<-matrix(beta_o,ncol=nitems,nrow=nind, byrow =T)
temp.a<-matrix(alpha_o,ncol=nitems,nrow=nind, byrow =T)
temp.t<-matrix(theta_o,ncol=nitems,nrow=nind)
logit_expt_theta_M<-exp(temp.t)/(1+exp(temp.t))
exp_alpha_M<-exp(temp.a)
exp_beta_M<-exp(temp.b)
### m1 model of yuli
pmod<-1-(1-(logit_expt_theta_M)^(exp_alpha_M))^(exp_beta_M)

summary(pmod)
test<-data.frame(ifelse(runif(nitems*nind)<pmod,1,0))
colnames(test)<-paste("t",1:nitems,sep="")
test%>%{apply(., 1,sum)/nitems} %>%summary
test%>%apply( 2,sum)/nind
test<-data.frame(test,id=1:nind)

###non linear function
nlf_m1<-function(X,psi,z,theta,nitems){
  logit_zt<-exp(z%*%theta)/(1+exp(z%*%theta))
  exp_xa<-exp(X%*%psi[(nitems+1):(2*nitems)])
  exp_xb<-exp(X%*%psi[1:nitems])
  return(drop(1-(1-(logit_zt)^(exp_xa))^(exp_xb)))
}
###inverse link
inv_link_m1<-function(x){x}
###fixed effects gradient
grad.fix_m1<-function(X,psi,z,theta,nitems){
  logit_zt<-exp(z%*%theta)/(1+exp(z%*%theta))
  exp_xa<-exp(X%*%psi[(nitems+1):(2*nitems)])
  exp_xb<-exp(X%*%psi[1:nitems])
  logitzt_high_expa<-logit_zt^(exp_xa)
  a_grad<-drop(exp_xa*exp_xb*(logitzt_high_expa)*log(logit_zt)*((1-logitzt_high_expa)^(exp_xb-1)))*X
  b_grad<- drop(-exp_xb*log(1-(logitzt_high_expa))*((1-(logit_zt^exp_xa))^(exp_xb)))*X
  return(cbind(b_grad,a_grad))
}
###random effects gradient
grad.mix_m1<-function(X,psi,z,theta,nitems){
  logit_zt<-exp(z%*%theta)/(1+exp(z%*%theta))
  exp_xa<-exp(X%*%psi[(nitems+1):(2*nitems)])
  exp_xb<-exp(X%*%psi[1:nitems])
  logitzt_high_expa<-logit_zt^(exp_xa)
  return(drop(exp_xa*exp_xb*logitzt_high_expa*((1-logitzt_high_expa)^(exp_xb-1))*(1-logit_zt))*z)
}

### derivative of link
hprim_m1<-function(x){rep(1,length(x))}

bres<-mh_gibbs_modular(test,"id",psi0=c(beta_o,alpha_o),theta0=theta_o,B=NULL,b=NULL,
                       grad.fix=grad.fix_m1,grad.mix=grad.mix_m1,nlf=nlf_m1,inv_link=inv_link_m1,
                       hprim=hprim_m1,Nc=50000)

brn<-1
chains<-bres$fixp
chns.sum<-sapply(1:dim(chains)[2],function(x){c(mean=mean(chains[-brn,x],na.rm = T),
                                                quantile(chains[-brn,x], probs = c(0.025,.5,.975), na.rm = T),
                                                sd=sd(chains[-brn,x],na.rm = T))})
oldpar2<-par(no.readonly =TRUE)
pdf(file="chains_m1.pdf",width=10,height=12)
par(mfrow=c(4,3))
for(i in 1:ncol(chains)){
  plot(as.ts(chains[,i]),cex.main=0.9,ylab="",xlab="iterations",main=paste("mcmc ",i))
  autocorr.plot(mcmc(chains[,i]),cex.main=0.9, col="red", lwd=2, auto.layout=FALSE ,lag.max=100,
                main=paste("acf. ",i))
  plot(density(chains[,i]),cex.main=0.9, col="red", lwd=2,main=paste("kds. ",i))

}
dev.off()
par(oldpar2)


