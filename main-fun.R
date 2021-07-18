################################################################################
## R CODE FOR CATTANEO-JANSSON (2017, ECMA)
## DATE: 19-DEC-2018
## Main file
################################################################################
## Functions
################################################################################

###############################
## DGPS
###############################
dgps = function(n, m, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    ## out=dgps(n=1000000,m=2); mean(out$y>=dnorm(out$x[1:n])*dnorm(out$x[(n+1):(2*n)]))

    dgps.par = read.csv("dgps.csv",stringsAsFactors=FALSE);
    eval.frac = function(x) eval(parse(text=x))

    muy   = eval.frac(dgps.par[dgps.par[,"M"]==m,"muy"])
    sdy   = eval.frac(dgps.par[dgps.par[,"M"]==m,"sdy"])
    mux   = eval.frac(dgps.par[dgps.par[,"M"]==m,"mux"])
    sdx   = eval.frac(dgps.par[dgps.par[,"M"]==m,"sdx"])
    d     = dgps.par[dgps.par[,"M"]==m,"d"]
    P     = dgps.par[dgps.par[,"M"]==m,"P"]
    theta = dgps.par[dgps.par[,"M"]==m,"theta0"]
    B     = dgps.par[dgps.par[,"M"]==m,"B0"]
    S     = dgps.par[dgps.par[,"M"]==m,"S0"]

    y = rnorm(n,mean=muy,sd=sdy);
    x = matrix(rnorm(d*n,mean=mux,sd=sdx),ncol=d);
    return(list(d=d, P=P, theta=theta, S=S, B=B, n=n, m=m, y=y, x=x))
}

###############################
## ROT Bandwidth Selectors
###############################
scale = function(x) {
    q = quantile(x, c(.25,.75))
    return(min(sd(x), (q[2]-q[1])/1.349))
}

ROT.h.1d = function(S, B, P, d, n) {
    if (sign(B) == sign(S)) {
        h.opt = (d * abs(B) / (P * abs(S) * n))^(1/(P+d))
    } else {
        h.opt = (abs(B) / (abs(S) * n))^(1/(P+d))
    }
    return(h.opt)
}

###############################
## Bootstrap Statistics
###############################
B.stats = function(x) c(quantile(x, probs=c(.025,.975,.5), na.rm=TRUE, type=1), mean(x), var(x))

###############################
## C Code Wrapper
###############################
# source("main-fun.R"); if (is.loaded("phi", PACKAGE="C_Kernel_Fns")) dyn.unload("C_Kernel_Fns.so"); dyn.load("C_Kernel_Fns.so")

C.Kest = function(y,x,d,P,n,B,h,ids) {
    tmp11=rep.int(0,B+1);
    tmp21=rep.int(0,n*(n-1)/2);
    tmp31=rep.int(0,n);

    out.C = .C("Kest",theta.hat=as.double(tmp11),
                      K.ij=as.double(tmp21),
                      f_i=as.double(tmp31),
                      y=as.double(y),x=as.double(x),d=as.integer(d),
                      P=as.integer(P),n=as.integer(n),B=as.integer(B),
                      h=as.double(h),ids=as.integer(ids))

    return(list(theta.hat=out.C$theta.hat[B+1],
                theta.hat.B=out.C$theta.hat[1:B],
                f_i=out.C$f_i))
}

## TEST: C Code
if (FALSE){
    d=1; n=100; B=2; h=.07; ids = c(sample(0:(n-1),B*n,replace=TRUE),0:(n-1))
    z=cbind(rexp(n),rexp(n));#mvrnorm(n,mu=rep(0,2),Sigma=matrix(c(1,1/2,1/2,1),2,2));

    out.C = C.Kest(y=z[,1],x=z[,2],d=d,P=2,n=n,B=B,h=h,ids=ids)

    ## 1) Testing density estimators & point estimators
    Kij = dnorm(outer(z[,2],z[,2],FUN="-"),0,h); Kii = diag(Kij)
    f.hat = colSums(Kij)/n;
    theta.hat = sum(z[,1]>=f.hat)/n

    cbind(out.C$f_i/n,f.hat,round(out.C$f_i/n-f.hat,6))
    cbind(out.C$theta.hat,theta.hat,round(out.C$theta.hat-theta.hat,6))

    ## 2) Testing bootstrap density estimators & point estimators
    k=1; ids.B=ids[(n*k+1):(n*k+n)]+1
    Kij = dnorm(outer(z[ids.B,2],z[ids.B,2],FUN="-"),0,h); Kii = diag(Kij)
    f.hat = colSums(Kij)/n;
    theta.hat.B = sum(z[ids.B,1]>=f.hat)/n

    cbind(out.C$f_i/n,f.hat,round(out.C$f_i/n-f.hat,6))
    cbind(out.C$theta.hat.B[k+1],theta.hat.B,round(out.C$theta.hat.B[k+1]-theta.hat.B,6))
}











