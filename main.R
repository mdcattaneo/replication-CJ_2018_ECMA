################################################################################
## R CODE FOR CATTANEO-JANSSON (2017, ECMA)
## DATE: 19-DEC-2018
## Main file
################################################################################
## Compile C code:
## gcc -shared -fPIC -o C_Kernel_Fns.so C_Kernel_Fns.c -lm
################################################################################
rm(list=ls(all=TRUE))
path=""; source("main-fun.R")
if (is.loaded("phi", PACKAGE="C_Kernel_Fns")) dyn.unload("C_Kernel_Fns.so"); dyn.load("C_Kernel_Fns.so")

################################################################################
## MONTECARLO SETUP
################################################################################
## Parallelization
s.part=0
s.block=10
S = seq.int(1,10000,s.block); s.start = S[s.part]; s.end = s.start + (s.block-1)

## Sample size & Bootstrap replications
n=3000; B=1000
## Bandwidth Grids
c.grd = seq(.5,1.3,.1)
#c.hat = seq(.7,1.3,.1)

## Design
#m=1
dgps.par=read.csv("dgps.csv"); M=nrow(dgps.par)
for (m in 1:M){
## Output Tables
thetas  = c("theta.hat","theta.hat_theta")
stats.B = c("B.q025","B.q975","B.median","B.mean","B.S2")

col.names = c("s","m","n","d","P","c","h","B","theta0",thetas)
for (i in c(thetas)) for (j in stats.B) col.names = c(col.names, paste(i,j,sep="."))

table.grd = matrix(NA, nrow=s.block*length(c.grd), ncol=length(col.names), dimnames=list(NULL,col.names))
#table.hat = matrix(NA, nrow=s.block*length(c.hat), ncol=length(col.names), dimnames=list(NULL,col.names))

## Load Design
out = dgps(n=1,m=m);
d = out$d; P = out$P; theta0 = out$theta;
h.opt = ROT.h.1d(S=out$S, B=out$B, P=P, d=d, n=n)

################################################################################
# Monte Carlo Experiment
################################################################################
# s.start=1; s=s.start; c.i=floor(length(c.grd)/2); n=100; B=10
showwhen = s.start; showevery=1; time.old.all = Sys.time()
cat(paste("\nSimulations began:",Sys.time()),"\n")

for (s in s.start:s.end) {
    ## Generate data
    data = dgps(n=n,m=m,seed=s*666); 
    ids = matrix(c(sample(0:(n-1),B*n,replace=TRUE),0:(n-1)), ncol=B+1)

    ## h.grd
    for (c.i in 1:length(c.grd)) {
        row.table = (s-s.start)*length(c.grd) + c.i
        h = c.grd[c.i]*h.opt
        table.grd[row.table,1:9] = c(s,m,n,d,P,c.grd[c.i],h,B,theta0)

        out.C = C.Kest(y=data$y,x=data$x,d=d,P=P,n=n,B=B,h=h,ids=ids)

        ## theta.hat Estimate
        table.grd[row.table, 10:11] = c(out.C$theta.hat,out.C$theta.hat-theta0)

        ## Bootstrap stats
        table.grd[row.table, 12:16] = B.stats(out.C$theta.hat.B)
        table.grd[row.table, 17:21] = B.stats(out.C$theta.hat.B-out.C$theta.hat)
    }

    ## h.hat
    ## Estimates for ROT bandwidths

    if (s==showwhen) {
        cat(paste("\nSimulations Completed:",s,"- From",s.start,"to",s.end,"- Diff:", round(difftime(Sys.time(),time.old.all,units="mins"),2),"mins"));
        time.old.all = Sys.time(); showwhen=showwhen+showevery
    }

};
## Save table
write.csv(table.grd, file=paste(path,"table_m",m,"_n",n,"_p",s.part,"_grd.csv",sep=""))
cat('\n Design ',m,' completed.\n')
}


