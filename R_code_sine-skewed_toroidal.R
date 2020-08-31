############################################################################################################
######       Sine-skewed toroidal distributions and their application in protein bioinformatics       ######
######                         Jose Ameijeiras-Alonso   and   Christophe Ley                          ######
######                               KU Leuven               Ghent University                         ######
######												      ######
######      Biostatistics, to appear (available at https://doi.org/10.1093/biostatistics/kxaa039)     ######
######           arXiv preprint: 1910.13293 (available at https://arxiv.org/pdf/1910.13293)           ######
############################################################################################################

## Supplementary R code for the paper "Sine-skewed toroidal distributions and their application in protein bioinformatics" 
## To cite this material in publications use
## Ameijeiras-Alonso, J., Ley, C. (2020). Sine-skewed toroidal distributions and their application in protein bioinformatics.  To appear in Biostatistics. DOI: 10.1093/biostatistics/kxaa039.

# install.packages("Rsolnp")



## 2.1 The general expression of sine-skewed densities. 
## Equation (2.1) for the bivariate case
## This function can be employed over any predefined toroidal symmetric density: d"model"=function(x,y,mu1,mu2,kappa1,kappa2,rho)

dsinesk=function(x,y,mu1,mu2,kappa1,kappa2,rho,lambda1,lambda2,model=NULL){ 
  
  if (!is.numeric(lambda1)|!is.numeric(lambda2)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}  
  if ((length(lambda1)!=1)|(length(lambda2)!=1)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}  
  if ((abs(lambda1)+abs(lambda2))>1) {stop("The sum of the absolute values of 'lambda1' and 'lambda2' must be lower or equal than one")}
  if (is.null(model)){
	warning("Argument 'model' is missing. By default, the sine-skewed Sine density is employed, i.e., 'model='sine''.")
	model="sine"
  }

  
  lamdamod=function(x,y,mu1,mu2,lambda1,lambda2) (1+lambda1*sin(x-mu1)+lambda2*sin(y-mu2)) 
  devalmodel=eval(parse(text=paste("d",model,"(x,y,mu1,mu2,kappa1,kappa2,rho)",sep="")))
  
  
  return(devalmodel*lamdamod(x,y,mu1,mu2,lambda1,lambda2))
}


## Sine distribution (Section 3.2)

dsine=function(x,y,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){ 
  
  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}  
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  
  
  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}

  indexes=0:99
  sCvec=0
  
  Cvec=choose(2*indexes,indexes) *(rho^2/(4*kappa1*kappa2))^indexes*besselI(kappa1,indexes)*besselI(kappa2,indexes)
  sCvec=sCvec+sum( Cvec )
  while(Cvec[length(Cvec)]>tol){
    indexes=indexes+100
    Cvec=choose(2*indexes,indexes) *(rho^2/(4*kappa1*kappa2))^indexes*besselI(kappa1,indexes)*besselI(kappa2,indexes)
    sCvec=sCvec+sum( Cvec )
  }
  
  C=4*pi^2*sCvec
  val=exp(kappa1*cos(x-mu1)+kappa2*cos(y-mu2)+rho*sin(x-mu1)*sin(y-mu2))/C
  return(val)
}


## Cosine distribution (Section 3.3)

dcosine=function(x,y,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){ 

  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}  
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  
  
  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}
  
  indexes=1:100
  sCvec=0
  if(rho<0){
    Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(abs(rho),indexes)*(-1)^(indexes)
  }else{
    Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(rho,indexes)
  }
  sCvec=sCvec+sum( Cvec )
  while(abs(Cvec[length(Cvec)])>tol){
    indexes=indexes+100
    if(rho<0){
       Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(abs(rho),indexes)*(-1)^(indexes)
    }else{
       Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(rho,indexes)
    }    
    sCvec=sCvec+sum( Cvec )
  }
  
  
  
  C=4*pi^2*(besselI(kappa1,0)*besselI(kappa2,0)*besselI(abs(rho),0)+2*sCvec)
  val=exp(kappa1*cos(x-mu1)+kappa2*cos(y-mu2)+rho*cos(x-mu1-y+mu2))/C
  return(val)
}




## Bivariate Wrapped Cauchy (Section 3.4)

dbwc=function(x,y,mu1,mu2,kappa1,kappa2,rho){
  
  
  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}  
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  
  
  if (kappa1>=1|kappa2>=1|kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be in the interval [0,1)")}
  if ((rho>=1)|(rho<= (-1))) {stop("Argument 'rho' must be in the interval (-1,1)")}
  
  
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  val= (1-rho^2)*(1-kappa1^2)*(1-kappa2^2)/(4*pi^2*(c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)))
  return(val)
}



## 2.3 Generating mechanism
## Random generation (Sine-skewed bivariate distribution)
## This function can be employed on any predefined random generation mechanism associated to a toroidal symmetric density: r"model"=function(x,y,mu1,mu2,kappa1,kappa2,rho)

rsinesk=function(n,mu1,mu2,kappa1,kappa2,rho,lambda1,lambda2,model=NULL){ 
  
  if (!is.numeric(lambda1)|!is.numeric(lambda2)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}  
  if ((length(lambda1)!=1)|(length(lambda2)!=1)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}  
  if ((abs(lambda1)+abs(lambda2))>1) {stop("The sum of the absolute values of 'lambda1' and 'lambda2' must be lower or equal than one")}
  if (is.null(model)){
	warning("Argument 'model' is missing. By default, the sine-skewed Sine distribution is employed, i.e., 'model='sine''.")
	model="sine"
  }
  
  newsample=eval(parse(text=paste("r",model,"(n,mu1,mu2,kappa1,kappa2,rho)",sep="")))
  
  ssnewsample=newsample
  un=runif(n)
  cons=un>(1+lambda1*sin(newsample[,1]-mu1)+lambda2*sin(newsample[,2]-mu2))/2
  ssnewsample[cons,1]=-ssnewsample[cons,1]+2*mu1
  ssnewsample[cons,2]=-ssnewsample[cons,2]+2*mu2

  return(ssnewsample)

}


## Sine distribution (acceptance-rejection method with enveloping constant function)

rsine=function(n,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){ 
  
  if (!is.numeric(n)) {stop("Argument 'n' must be a positive integer number")}
  if ((length(n) != 1) | (n%%1 != 0) | (n <= 0)){stop("Argument 'n' must be a positive integer number")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  

  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}
  
  newsample=matrix(0,ncol=2,nrow=n)
  if(kappa1*kappa2>rho^2) {
    Mf=dsine(mu1,mu2,mu1,mu2,kappa1,kappa2,rho,tol)
  }else{
    seqx=seq(-pi,pi,len=1000)
    Mf=max(outer(seqx,seqx,function(x,y) dsine(x,y,mu1,mu2,kappa1,kappa2,rho,tol)))  
    Mf=1.1*Mf
  }
  
  i=1
  while(i<=n){
    xu=runif(1,-pi,pi)
    yu=runif(1,-pi,pi)
    u=runif(1)
    if(u<dsine(xu,yu,mu1,mu2,kappa1,kappa2,rho,tol)/Mf){
      newsample[i,]=c(xu,yu)
      i=i+1
    }
  }
  return(newsample)
}


## Cosine distribution (acceptance-rejection method)

rcosine=function(n,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){ 
  
  if (!is.numeric(n)) {stop("Argument 'n' must be a positive integer number")}
  if ((length(n) != 1) | (n%%1 != 0) | (n <= 0)){stop("Argument 'n' must be a positive integer number")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  

  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}
  
  newsample=matrix(0,ncol=2,nrow=n)
  if(-rho<(kappa1*kappa2/(kappa1+kappa2))) {
    Mf=dcosine(mu1,mu2,mu1,mu2,kappa1,kappa2,rho,tol)
  }else{
    seqx=seq(-pi,pi,len=1000)
    Mf=max(outer(seqx,seqx,function(x,y) dcosine(x,y,mu1,mu2,kappa1,kappa2,rho,tol)))  
    Mf=1.1*Mf
  }
  
  i=1
  while(i<=n){
    xu=runif(1,-pi,pi)
    yu=runif(1,-pi,pi)
    u=runif(1)
    if(u<dcosine(xu,yu,mu1,mu2,kappa1,kappa2,rho,tol)/Mf){
      newsample[i,]=c(xu,yu)
      i=i+1
    }
  }
  return(newsample)
}


## Bivariate wrapped Cauchy distribution (acceptance-rejection method)

rbwc=function(n,mu1,mu2,kappa1,kappa2,rho){ 
  
  if (!is.numeric(n)) {stop("Argument 'n' must be a positive integer number")}
  if ((length(n) != 1) | (n%%1 != 0) | (n <= 0)){stop("Argument 'n' must be a positive integer number")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}  
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}  
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}  
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}  

  if (kappa1>=1|kappa2>=1|kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be in the interval [0,1)")}
  if ((rho>=1)|(rho<= (-1))) {stop("Argument 'rho' must be in the interval (-1,1)")}
  
  newsample=matrix(0,ncol=2,nrow=n)
  if((kappa1>0)&(kappa2>0)) {
    Mf=dbwc(mu1,mu2,mu1,mu2,kappa1,kappa2,rho)
  }else{
    seqx=seq(-pi,pi,len=1000)
    Mf=max(outer(seqx,seqx,function(x,y) dbwc(x,y,mu1,mu2,kappa1,kappa2,rho)))  
    Mf=1.1*Mf
  }
  
  i=1
  while(i<=n){
    xu=runif(1,-pi,pi)
    yu=runif(1,-pi,pi)
    u=runif(1)
    if(u<dbwc(xu,yu,mu1,mu2,kappa1,kappa2,rho)/Mf){
      newsample[i,]=c(xu,yu)
      i=i+1
    }
  }
  return(newsample)
}


## 4.1 Maximum likelihood estimation
## Obtained with the Rsolnp library 
## For the symmetric submodels set symmetric=TRUE 
## ninipar is the number of initial starting points employed in the optimization algorithm
## vecmax is a vector of 3 elements containing the largest values employed for the parameters (kappa1, kappa2, rho) in the random initial points 
## vecmin is a vector of 3 elements containing the smallest values employed for the parameters (kappa1, kappa2, rho)  in the random initial points 
## lowercons and uppercons are the lower and upper paremeter constrains

mle.sinesk=function(x,symmetric=FALSE,model="sine",ninipar=NULL,vecmax=NULL,vecmin=NULL,lowercons=NULL,uppercons=NULL){ 

  if(!is.matrix(x)) { stop("Argument 'x' must be a numeric matrix")}
  if(dim(x)[2]!=2) { stop("Argument 'x' must be a matrix with two columns")}
  if(symmetric != T & symmetric != F) { 
    warning("Argument 'symmetric' must be T or F. Default value of 'symmetric' was used")
    symmetric=FALSE
  }
  
  if(is.null(ninipar)){
    ninipar=100
    if(model=="cosine"){ninipar=1000}
  }
  if (!is.numeric(ninipar)) {
    warning("Argument 'ninipar' must be a positive integer number. Default value of 'ninipar' was used")
    ninipar=100
    if(model=="cosine"){ninipar=1000}
  }
  if ((length(ninipar) != 1) | (ninipar%%1 != 0) | (ninipar <= 0)){
    warning("Argument 'ninipar' must be a positive integer number Default value of 'ninipar' was used")
    ninipar=100
    if(model=="cosine"){ninipar=1000}
  }
  
  if(is.null(vecmax)){
    if(model=="sine"|model=="cosine"){valmax=10}
    if(model=="bwc"){valmax=1}
    vecmax=rep(valmax,3)
  }
  if (!is.numeric(vecmax)) {stop("Argument 'vecmax' must be a vector of three real numbers")}
  if (length(vecmax) != 3) {stop("Argument 'vecmax' must be a vector of three real numbers")} 
  
  if(is.null(vecmin)){
    if(model=="sine"|model=="cosine"){valmax=10}
    if(model=="bwc"){valmax=1}
    vecmin=c(0,0,-valmax)
  }  
  if (!is.numeric(vecmin)) {stop("Argument 'vecmin' must be a vector of three real numbers")}
  if (length(vecmin) != 3) {stop("Argument 'vecmin' must be a vector of three real numbers")} 
  
  if(is.null(lowercons)){
    lowercons=c(-pi,-pi,vecmin)
    if(model=="sine"|model=="cosine"){
      lowercons[5]=-Inf
    }
    if(symmetric==F){
      lowercons=c(lowercons,-1,-1)
    }
  }

  if (!is.numeric(lowercons)) {stop("Argument 'lowercons' must be a vector of three real numbers")}
  if(symmetric==T){
  if (length(lowercons) != 5) {stop("Argument 'lowercons' must be a vector of three real numbers")} 
  }else{
    if (length(lowercons) != 7) {stop("Argument 'lowercons' must be a vector of three real numbers")}
  }
  
  if(is.null(uppercons)){  
    uppercons=c(pi,pi,vecmax)
    if(model=="sine"|model=="cosine"){
      uppercons[3:5]=Inf
    }
    if(symmetric==F){
      uppercons=c(uppercons,1,1)
    }  
  }

  
  if (!is.numeric(uppercons)) {stop("Argument 'uppercons' must be a vector of three real numbers")}
  if(symmetric==T){
    if (length(uppercons) != 5) {stop("Argument 'uppercons' must be a vector of three real numbers")} 
  }else{
    if (length(uppercons) != 7) {stop("Argument 'uppercons' must be a vector of three real numbers")}
  }  
  
  
  
  require(Rsolnp)
  
  if(symmetric==F){
    llpar=function(par){
      -sum(log(dsinesk(x[,1],x[,2],par[1],par[2],par[3],par[4],par[5],par[6],par[7],model)))
    }
  }else{
    llpar=function(par){
      -sum(log(dsinesk(x[,1],x[,2],par[1],par[2],par[3],par[4],par[5],0,0,model)))
    }
  }

  suml1=function(par){abs(par[6])+abs(par[7])}
  valllfin=Inf
  paramfin=numeric()
  
  for(iterip in 1:ninipar){
  inipar=numeric()
  inipar[1]=runif(1,-pi,pi)
  inipar[2]=runif(1,-pi,pi)
  inipar[3]=runif(1,vecmin[1],vecmax[1])
  inipar[4]=runif(1,vecmin[2],vecmax[2])
  inipar[5]=runif(1,vecmin[3],vecmax[3])
  if(symmetric==F){
    lambt=runif(1,-1,1)
    lambt=c(lambt,runif(1,-1+abs(lambt),1-abs(lambt)))
    ilambt=sample(1:2,2)        
    inipar[6]=lambt[ilambt[1]]
    inipar[7]=lambt[ilambt[2]]
  }

  if(symmetric==F){
    paramtot=try(solnp(inipar,llpar,ineqfun=suml1,ineqLB=0,ineqUB=1,LB=lowercons,UB=uppercons,control = list(trace=0)),silent=T)
  }else{
    paramtot=try(solnp(inipar,llpar,ineqLB=0,ineqUB=1,LB=lowercons,UB=uppercons,control = list(trace=0)),silent=T)
  }
  
  if(class(paramtot)=="try-error"){
    valll=Inf
  }else{
    
    if(symmetric==F){
      if(suml1(paramtot$pars)<=1){
        valll=paramtot$values[length(paramtot$values)]
      }else{
        valll=Inf
      }
    }else{
      valll=paramtot$values[length(paramtot$values)]
    }
    
    if(valll<valllfin){
      paramfin=paramtot$pars
      valllfin=valll
    }
  }
  }
  
  resultmle=list()
  resultmle$mu1=paramfin[1]
  resultmle$mu2=paramfin[2]
  resultmle$kappa1=paramfin[3]
  resultmle$kappa2=paramfin[4]
  resultmle$rho=paramfin[5]
  if(symmetric==F){
    resultmle$lambda1=paramfin[6]
    resultmle$lambda2=paramfin[7]    
  }  
  
  #Estimated LL value
  resultmle$LL=-valllfin
  
  return(resultmle)
  
}  

## 4.3 Testing for symmetry

symmetry.test=function(x,model=NULL,mu1=NULL,mu2=NULL,kappa1=NULL,kappa2=NULL,rho=NULL,lambda1=NULL,lambda2=NULL,ninipar=NULL){
  
  if(!is.matrix(x)) { stop("Argument 'x' must be a numeric matrix")}
  if(dim(x)[2]!=2) { stop("Argument 'x' must be a matrix with two columns")}
  
  ndata=length(x[,1])
  if(!is.null(mu1)&!is.null(mu2)&!is.null(kappa1)&!is.null(kappa2)&!is.null(rho)&!is.null(lambda1)&!is.null(lambda2)){
    statval=-2*(sum(log(dsinesk(x[,1],x[,2],mu1,mu2,kappa1,kappa2,rho,0,0,model)))-sum(log(dsinesk(x[,1],x[,2],mu1,mu2,kappa1,kappa2,rho,lambda1,lambda2,model))))
  }else{
    msym=mle.sinesk(x,symmetric=TRUE,model=model,ninipar=ninipar)
    msinesk=mle.sinesk(x,symmetric=FALSE,model=model,ninipar=ninipar)
    LL0=msym$LL
    LL1=msinesk$LL
    statval=-2*(LL0-LL1)
  }
  pval=1-pchisq(statval,2)
  names(statval)="LR"
  lamb0=0
  names(lamb0)="Value of the vector lambda"
  totval=list(statistic = statval, p.value = pval, null.value = lamb0, alternative = "two.sided", method = "Symmetry Likelihood Ratio Test", sample.size = ndata,data.name=deparse(substitute(x)))
  class(totval) <- "htest"
  return(totval)
}  
