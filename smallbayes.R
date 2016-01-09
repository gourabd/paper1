infor<-function(o,f,m)
{
 d1=which(f+m==1);
 d2=which(f+m==3);
 d3=which(f==1 & m==1);
 infor=sort(c(d1,d2,d3))
 info=o[infor]
 n0I=sum(info==0)
 n1I=sum(info==1)
 n2I=sum(info==2)
 c(n0I,n1I,n2I)
}

child=function(f)
{
 c=rep(0,length(f))
 c[f==2]=1
 c[f==1]=sample(c(0,1),length(which(f==1)),replace=TRUE,c(0.5,0.5))
 c
}
llh<-function(logalpha,o,f,m)
{
 alpha1=exp(logalpha[1])
 alpha2=exp(logalpha[2])
 n01P=sum(f+m==1);d1=which(f+m==1);
 n12P=sum(f+m==3);d2=which(f+m==3);
 n11P=sum(f==1 & m==1);d3=which(f==1 & m==1);
 infor=sort(c(d1,d2,d3))
 info=o[infor]
 n0I=sum(info==0)
 n1I=sum(info==1)
 n2I=sum(info==2)
 -((n1I*log(alpha1))+(n2I*log(alpha2))-(n01P*log(1+alpha1))
    -(n12P*log(alpha1+alpha2))-(n11P*log(1+(2*alpha1)+alpha2)))
}  
gradllh<-function(logalpha,o,f,m)
{
 alpha1=exp(logalpha[1])
 alpha2=exp(logalpha[2])
 n01P=sum(f+m==1);d1=which(f+m==1);
 n12P=sum(f+m==3);d2=which(f+m==3);
 n11P=sum(f==1 & m==1);d3=which(f==1 & m==1);
 infor=sort(c(d1,d2,d3))
 info=o[infor]
 n0I=sum(info==0)
 n1I=sum(info==1)
 n2I=sum(info==2)
 r=0;
 r[1]=alpha1*((n1I/alpha1)-(n01P/(1+alpha1))-(n12P/(alpha1+alpha2))
      -(2*n11P/(1+(2*alpha1)+alpha2)));
 r[2]=alpha2*((n2I/alpha2)-(n12P/(alpha1+alpha2))
      -(n11P/(1+(2*alpha1)+alpha2)));
  -r
}

sims <- function(beta,mu,n)
{
 p=0.5
 prob=rbind(p^2,2*p*(1-p),(1-p)^2)
#parent genotype
 f=apply(prob,2,sample,x=c(0,1,2),size=n,replace=TRUE)
 m=apply(prob,2,sample,x=c(0,1,2),size=n,replace=TRUE)
#offspring genotype
 o=child(f)+child(m)
 eta=0
 eta[o==0]=mu
 eta[o==1]=mu+beta[1]
 eta[o==2]=mu+beta[2]
 pro=exp(-eta)
 r=runif(n,0,1)
 y=as.numeric(pro>r)
 sampo=o[y==1]
 sampf=f[y==1]
 sampm=m[y==1]
 init1=sum(sampo==1)/sum(sampo==0)
 init2=sum(sampo==2)/sum(sampo==0)
 init=log(c(init1,init2))
 logalpha=(((optim(init,fn=llh,gr=gradllh,method="CG",o=sampo,f=sampf,m=sampm)))$par)
 nullt=(-llh(c(0,0),sampo,sampf,sampm))
 samplogalpha=matrix(rnorm(2*n,0,5),n,2)
 samps=log(mean(exp(-apply(samplogalpha,1,llh,o=sampo,f=sampf,m=sampm))))
 sampt=(-llh(logalpha,sampo,sampf,sampm))
 t1=2*(sampt-nullt)
 t2=2*(samps-nullt)
 p1=as.numeric(t1>qchisq(0.95,2))
 p2=as.numeric(t2>qchisq(0.95,2))
 c(beta,mu,n,logalpha,length(sampo),nullt,sampt,samps,t1,t2,p1,p2)
}
param=rep(3000,1000)
res=t(sapply(param,sims,beta=c(-1,-2),mu=3))
resmean=apply(res,2,mean)
resvar=apply(res,2,var)
write(beta
write(resmean,"resu")
write(resvar,"resu",append=TRUE)

fulllh<-function(logalpha,o,f,m,me,sd)
{
 part1=exp(-llh(logalpha,o,f,m))
 part2=prod(pnorm(logalpha,0,1000))
 part1*part2
}


