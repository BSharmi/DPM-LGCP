### R code from vignette source 'scale-model.Rnw'

###################################################
### code chunk number 1: scale-model.Rnw:68-82
####################################################

if (Sys.getenv("USER") == "sigrunnsrbye") {
setwd("/Users/sigrunnsrbye/Dropbox/Tutorial_igmrfs")
} 

library(spatstat)
library(mvtnorm)
library(lattice)
library(mgcv)
library(pixmap)
library(numDeriv)
library(fields)
library(INLA)
#inla.upgrade(testing=TRUE)


###################################################
### code chunk number 2: scale-model.Rnw:106-107 (eval = FALSE)
###################################################
## formula = y ~  f(idx, model = "rw2", scale.model = TRUE, hyper=..)


###################################################
### code chunk number 3: scale-model.Rnw:113-114 (eval = FALSE)
###################################################
## f(idx, model="..", group=g, control.group = list(model="rw2", scale.model=TRUE))


###################################################
### code chunk number 4: scale-model.Rnw:186-241
###################################################
ss.rw1.sdref = function(values,tau) { 
    n = length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3) 
    y = NA
    
    formula = y ~ -1 + f(values, 
            model="rw1", 
            constr =TRUE, 
            values = values, 
            diagonal = 1e-10, 
            hyper = list(prec = list(
                                 initial = log(tau), 
                                 fixed=TRUE))) 
    r =inla(formula,data=data.frame(y,values),family = "gaussian",
            control.family = list(fixed =
                    TRUE), 
            control.compute=list(return.marginals=F),verbose=F)
    sd=r$summary.random$values$sd 
    return (sd) 
    }

ss.rw2.sdref= function(values,tau) { 
    n=length(values)
    stopifnot(!missing(n)) 
    stopifnot(n > 3)
    y = NA #idx =1:n
    
    A1 = matrix(1, n, 1) 
    A2 = matrix(0, n, 1) 
    for(i in 1:n) { 
        A2[i, ]= i 
    }
    
    A = rbind(c(A1), c(A2)) 
    e = rep(0, 2) 
    extraconstr = list(A = A, e= e)
    
    formula = y ~ -1 + f( values, model="rw2", constr = FALSE,
            values=values, diagonal = 1e-10, extraconstr = extraconstr, hyper
            = list(prec = list(initial = log(tau), fixed=TRUE)))
    
    r = inla(formula, data = data.frame(y, values), family =
            "gaussian", control.family = list(fixed = TRUE),
            control.compute=list(return.marginals=F),verbose=F)
    
    sd=r$summary.random$values$sd 
    return (sd) 
}

n=100 
sd.rw1=ss.rw1.sdref(1:n,1) 
sd.rw2=ss.rw2.sdref(1:n,1)
##sd.rw1.scaled=ss.rw1.sdref(1:n,1) 
##sd.rw2.scaled=ss.rw2.sdref(1:n,1)


###################################################
### code chunk number 5: scale-model.Rnw:248-249
###################################################
plot(sd.rw1,xlab="Node i", ylab="Marginal standard deviation") 


###################################################
### code chunk number 6: scale-model.Rnw:254-255
###################################################
plot(sd.rw2,xlab="Node i",ylab="Marginal standard deviation") 


###################################################
### code chunk number 7: scale-model.Rnw:310-311 (eval = FALSE)
###################################################
## U = sqrt(b/qgamma(alpha,a,1)) 


###################################################
### code chunk number 8: scale-model.Rnw:313-321
###################################################
func.u = function(a,b,alpha,sigma.ref) { 
    upper.limit = sqrt(b)*sigma.ref/sqrt(qgamma(alpha,a,1)) 
    return(upper.limit) 
}
a=1 
b=5*10^{-5} 
alpha=0.001 
upper.limit=func.u(a,b,alpha,1) 


###################################################
### code chunk number 9: scale-model.Rnw:340-353
###################################################
func<-function(x,t)
{
    f.x=x/t+sin(2*pi* x/t)-0.5
    return(f.x)
}

f.est <- function(y,x,a,b,option)
{
    formula = y~1+f(x,model="rw2",hyper=list(prec=list(prior="loggamma",param=c(a,b))),scale.model=option)
    result = inla(formula,family="gaussian",data=data.frame(y,x))
    f.est=result$summary.random$x$mean
    return(f.est)
}


###################################################
### code chunk number 10: scale-model.Rnw:356-381
###################################################
### Same values of the function, same observations, but  defined on different intervals ###
set.seed(89236)
x2=0:n
x1=x2/100
x3=x2*10

sigma=0.5
f.x=func(x2,length(x2))
y=f.x+rnorm(length(f.x),0,sigma)

a=1
b=0.00005

mat1=cbind(f.x,f.est(y,x1,a,b,T),f.est(y,x1,a,b,F))
mat2=cbind(f.x,f.est(y,x2,a,b,T),f.est(y,x2,a,b,F))
mat3=cbind(f.x,f.est(y,x3,a,b,T),f.est(y,x3,a,b,F))

# Generalized marginal variances  
v1=exp(mean(log(ss.rw2.sdref(x1,1)^2)))
v2=exp(mean(log(ss.rw2.sdref(x2,1)^2)))
v3=exp(mean(log(ss.rw2.sdref(x3,1)^2)))

u1=func.u(a,b,alpha,sqrt(v1))
u2=func.u(a,b,alpha,sqrt(v2))
u3=func.u(a,b,alpha,sqrt(v3))


###################################################
### code chunk number 11: scale-model.Rnw:396-397 (eval = FALSE)
###################################################
## formula = y ~  f(x, model = "rw2", scale.model = TRUE, hyper=...)


###################################################
### code chunk number 12: scale-model.Rnw:408-409 (eval = FALSE)
###################################################
## result = inla(formula, family = "gaussian", data = data.frame(y,x))


###################################################
### code chunk number 13: scale-model.Rnw:412-413 (eval = FALSE)
###################################################
## result$summary.random$x$mean


###################################################
### code chunk number 14: scale-model.Rnw:423-428
###################################################
# Estimate of f using unscaled and scaled model
r=range(c(mat1,mat2,mat3))
matplot(x1,mat1,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3) 
#legend("topright",c("true f","scale.model=F","scale.model=T"),lty=c(1,2,4),col=c(1,2,4))
#points(x,y)


###################################################
### code chunk number 15: scale-model.Rnw:433-434
###################################################
matplot(x2,mat2,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)


###################################################
### code chunk number 16: scale-model.Rnw:439-440
###################################################
matplot(x3,mat3,type="lll",lty=c(1,2,4),col=c(1,4,2),xlab="x",ylab="f(x)",ylim=r,lwd=3)


###################################################
### code chunk number 17: scale-model.Rnw:504-505 (eval = FALSE)
###################################################
## hyper.prec = list(prec = list(param = c(1, 0.01)))


###################################################
### code chunk number 18: scale-model.Rnw:508-523 (eval = FALSE)
###################################################
## data(Munich)
## g = system.file("demodata/munich.graph", package="INLA")
## 
## ## Note that here we what to have an estimator of the effect of year
## ## also the for years where we have no observation, therefore we give a
## ## vector with all possible values assumed by the covariate year, that
## ## is seq(1918,2001)
## 
## formula = rent ~ 
##    f(location, model = "besag", graph = g, initial = 1, hyper = hyper.prec) +
##    f(year, model = "rw2", values = seq(1918,2001), hyper = hyper.prec) +
##    f(floor.size, model = "rw2", hyper = hyper.prec) +
##    Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
##    Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
##    zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1


###################################################
### code chunk number 19: scale-model.Rnw:526-527 (eval = FALSE)
###################################################
## mod  =  inla(formula, data = Munich)


###################################################
### code chunk number 20: scale-model.Rnw:530-545
###################################################
data(Munich)
g = system.file("demodata/munich.graph", package="INLA")

func.munich <- function(a,b,option)
{
    hyper.prec=list(prec=list(param=c(a,b)))
    formula= rent ~ f(location,model="besag",graph=g,initial=1, hyper=hyper.prec, scale.model=option) +
        f(year,model="rw2",values=seq(1918,2001), hyper=hyper.prec, scale.model=option) +
            f(floor.size,model="rw2",param=c(1,0.01), hyper=hyper.prec, scale.model=option) +
                Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh +
                    Kein.Badkach  + Besond.Bad + Gehobene.Kueche +
                        zim1 + zim2 + zim3 + zim4 + zim5 + zim6 -1
    
    result =  inla(formula,data=Munich)
}


###################################################
### code chunk number 21: scale-model.Rnw:548-553
###################################################
a=1
b=0.01
u.munich=func.u(a,b,alpha,1)
mod=func.munich(a,b,FALSE)
mod.scaled=func.munich(a,b,TRUE)


###################################################
### code chunk number 22: scale-model.Rnw:571-573
###################################################
mat=cbind(mod$summary.random$year$mean,mod.scaled$summary.random$year$mean)
matplot(mod$summary.random$year$ID,mat,type="ll",lty=c(2,1),col=c(2,1),xlab="Year of construction",ylab="",lwd=2.5)


###################################################
### code chunk number 23: scale-model.Rnw:578-580
###################################################
mat=cbind(mod$summary.random$floor.size$mean,mod.scaled$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID,mat,type="ll",lty=c(2,1),col=c(2,1),xlab="Floor size",ylab="",lwd=2.5)


###################################################
### code chunk number 24: scale-model.Rnw:595-606
###################################################
alpha=0.001
#u.vec=c(0.01,1,10,100)
u.vec=c(0.001,0.1,10,30)
a.vec=c(1,1,1,1)
b.vec=u.vec^2*qgamma(alpha,a.vec,1)

option=TRUE
c1=func.munich(a.vec[1],b.vec[1],option)
c2=func.munich(a.vec[2],b.vec[2],option)
c3=func.munich(a.vec[3],b.vec[3],option)
c4=func.munich(a.vec[4],b.vec[4],option)


###################################################
### code chunk number 25: scale-model.Rnw:629-632
###################################################
mat=cbind(c1$summary.random$year$mean,c2$summary.random$year$mean,c3$summary.random$year$mean,c4$summary.random$year$mean)
matplot(mod$summary.random$year$ID,mat,type="llll",lty=c(1,2,3,4),col=c(1,2,3,4),xlab="Year of construction",ylab="",lwd=2.5)
legend("topleft",paste("U= ",u.vec),lty=1:4,col=1:4,lwd=2.5)


###################################################
### code chunk number 26: scale-model.Rnw:637-640
###################################################
mat=cbind(c1$summary.random$floor.size$mean,c2$summary.random$floor.size$mean,c3$summary.random$floor.size$mean,c4$summary.random$floor.size$mean)
matplot(mod$summary.random$floor.size$ID,mat,type="llll",lty=c(1,2,3,4),col=c(1,2,3,4),xlab="Floor size",ylab="",lwd=2.5)
legend("topright",paste("U= ",u.vec),lty=1:4,col=1:4,lwd=2.5)


###################################################
### code chunk number 27: scale-model.Rnw:666-667 (eval = FALSE)
###################################################
## inla.setOption(scale.model.default = TRUE)


