library(AlgDesign)
ptrans<-function(dp,q)
{
 o=dp[1:q]
 f=dp[(q+1):(2*q)]
 m=dp[((2*q)+1):(3*q)]
 # o,f,m each vector corresponds to multiple markers for one subject
 pr=rep(1,q)
 for(i in 1:q)
 {
 if(o[i]==0 & f[i]==1 & m[i]==0) pr[i]=0.5
 if(o[i]==1 & f[i]==1 & m[i]==0) pr[i]=0.5
 if(o[i]==0 & f[i]==0 & m[i]==1) pr[i]=0.5
 if(o[i]==1 & f[i]==0 & m[i]==1) pr[i]=0.5
 if(o[i]==1 & f[i]==1 & m[i]==2) pr[i]=0.5
 if(o[i]==2 & f[i]==1 & m[i]==2) pr[i]=0.5
 if(o[i]==1 & f[i]==2 & m[i]==1) pr[i]=0.5
 if(o[i]==2 & f[i]==2 & m[i]==1) pr[i]=0.5
 if(o[i]==0 & f[i]==1 & m[i]==1) pr[i]=0.25
 if(o[i]==1 & f[i]==1 & m[i]==1) pr[i]=0.5
 if(o[i]==2 & f[i]==1 & m[i]==1) pr[i]=0.25
 }
 prod(pr)
}
 
setoff<-function(par,q)
{
 f=par[1:q];
 m=par[(q+1):(2*q)]
 # f and m are row vectors - multiple markers for one subject
 loff=vector("list",length(f))
 for(i in 1:length(f))
 {
  if (f[i]+m[i]==0) loff[[i]]=0;
  if (f[i]+m[i]==4) loff[[i]]=2;
  if (f[i]==1 & m[i]==1) loff[[i]]=c(0,1,2);
  if (f[i]+m[i]==1) loff[[i]]=c(0,1);
  if (f[i]+m[i]==3) loff[[i]]=c(1,2);
  if (f[i]==0 & m[i]==2) loff[[i]]=1;
  if (f[i]==2 & m[i]==0) loff[[i]]=1;
 }
 loff # list with possible offspring - each entry in the list correspond to a marker
}

setd<-function(par,q)
{
 f=par[1:q];
 m=par[(q+1):(2*q)]
 # f and m are row vectors - levels of multiple markers for any one subject in a strata  
 loff=setoff(par,q);
 len=unlist(lapply(loff,length));
 g=as.matrix(gen.factorial(len,nVars=q,center=FALSE))
 for(i in 1:length(f))
 {
  if (len[i]==1) g[,i]=loff[[i]]
  if (len[i]==2 & f[i]+m[i]==1) g[,i]=g[,i]-1
  if (len[i]==3) g[,i]=g[,i]-1
 }
 g  
#S(f,m) for one specific parental genotype combination - each column for a marker, each row for possible offspring combination
}
child=function(f)
{
 # f is a vector - for a single marker
 c=rep(0,length(f))
 c[f==2]=1
 c[f==1]=sample(c(0,1),length(which(f==1)),replace=TRUE,c(0.5,0.5))
 c
}
stratacount<-function(o,f,m)
{
 #q=length(o[1,])
 #bigmat=data.frame(o);
 #ones=rep(1,length(bigmat[,1]));
 #strata=aggregate(ones,by=as.list(bigmat),sum)
 #strata=data.matrix(strata)
 #strata[,1:(q)]=strata[,1:(q)]
 #strata

 q=length(o[1,])
 bigmat=data.frame(cbind(o,f,m));
 ones=rep(1,length(bigmat[,1]));
 strata=aggregate(ones,by=as.list(bigmat),sum)
 strata=data.matrix(strata)
 strata[,1:(3*q)]=strata[,1:(3*q)]
 strata

}
matingstrata<-function(p)
{
 q=length(p[1,])/2
 bigmat=data.frame(p);
 ones=rep(1,length(bigmat[,1]));
 strata=aggregate(ones,by=as.list(bigmat),sum)
 strata=data.matrix(strata)
 strata[,1:(2*q)]=strata[,1:(2*q)]
 strata
}

init<-function(code,q,gencode)
{
 base=paste(rep(0,q),collapse='')
 beta=1
 if (sum(gencode==base)!= 0) beta=sum(gencode==code)/sum(gencode==base)
 beta
}

nulll<-function(data,no_strat,q)
{
 sum(no_strat*log(data[,(3*q)+1]))
}
 
addllh<-function(logalpha,denomdata,no_mstrat,data,no_strat,q)
{
 workbeta=exp(denomdata[,1:q]%*% logalpha[1:q])*exp(denomdata[,1]*denomdata[,2]*logalpha[q+1])
 denom1=(denomdata[,(3*q)+1]*workbeta)
 pardata=denomdata[,(q+1):(3*q)]
 denom=sum(no_mstrat*log(aggregate(denom1,by=as.list(data.frame(pardata)),sum)$V1))
 offbeta=exp(data[,1:q]%*% logalpha[1:q])*exp(data[,1]*data[,2]*logalpha[q+1])
 numer1=data[,(3*q)+1]*offbeta
 numer=sum(no_strat*log(numer1))
-numer+denom

}




llh<-function(logbeta,n,newn,S,denomdata,no_mstrat,offn,newoffn,offS,data,no_strat,q)
{
 beta=rep(1,length(n))
 beta[match(newn,n)]=exp(logbeta)
 workbeta=rep(1,length(denomdata[,1]))
 for(i in match(newn,n)){ workbeta[S[[i]]]=beta[i]}

 denom1=(denomdata[,(3*q)+1]*workbeta)
 pardata=denomdata[,(q+1):(3*q)]
 denom=sum(no_mstrat*log(aggregate(denom1,by=as.list(data.frame(pardata)),sum)$x))
 
 offbeta=rep(1,length(data[,1]))
 offb=beta[match(offn,n)]
 for(i in match(newoffn,offn)){ offbeta[offS[[i]]]=offb[i]}

 numer1=data[,(3*q)+1]*offbeta
 numer=sum(no_strat*log(numer1))
 

 #numer=sum(no_strat*log(offbeta))
 
-numer+denom
}
 

sim<-function(l,num,mu,betas,intbeta,q,pr)
{

prob=rbind(pr^2,2*pr*(1-pr),(1-pr)^2)

R=matrix(0,l,(3^q)+7+q)
R=data.frame(R)

SS=matrix(0,l,(q)+7+q)
SS=data.frame(SS)
f=0 
mo=0
gencode=0
samppar=0
uniquepar=0
for(t in 1:l)
{
#parent genotype
 f=apply(prob,2,sample,x=c(0,1,2),size=num,replace=TRUE)
 mo=apply(prob,2,sample,x=c(0,1,2),size=num,replace=TRUE)
#offspring genotype
 o=apply(f,2,child)+apply(mo,2,child)
eta=(mu+o%*%betas+(o[,1]*o[,2]*intbeta))
pro=exp(-eta)
r=runif(num,0,1)
y=as.numeric(pro>r)
sampo=as.matrix(o[y==1,])
sampf=f[y==1,]
sampm=mo[y==1,]
samppar=cbind(sampf,sampm)
prop=mean(y==1)


uniquepar=matingstrata(samppar)

possoff=apply(uniquepar[,1:(2*q)],1,setd,q)



possoffdata=do.call(rbind,possoff)
possoffcode=apply(possoffdata,1,paste,collapse='')



gencode=apply(sampo,1,paste,collapse='')

S=(split(seq(nrow(possoffdata)),as.list(data.frame(possoffcode))))
n=names(S)
base=paste(rep(0,q),collapse='')
newn=n[n!=base]
m=mapply(init,n,MoreArgs=list(q,gencode))
initbeta=1
for(i in match(newn,n)){ initbeta[S[[i]]]=m[i]}
initbeta[is.na(initbeta)]=1


length= matrix(unlist(lapply(possoff,dim)),length(possoff),2,byrow=T)[,1]
x.row <- split(seq(nrow(uniquepar[,1:(2*q)])),seq(nrow(uniquepar)))
res=lapply(seq_along(x.row), function(.fact){uniquepar[rep(x.row[[.fact]], length[.fact]),1:(2*q)]})
pardata=do.call(rbind,res)
possdata=cbind(possoffdata,pardata)
obstransdenom=apply(possdata,1,ptrans,q)
denomdata=cbind(possdata,obstransdenom)
no_mstrat=uniquepar[,(2*q)+1]

strata=stratacount(sampo,sampf,sampm)
strcode=apply(as.matrix(strata[,1:q]),1,paste,collapse='')

obstrans=apply(strata[,1:(3*q)],1,ptrans,q)
data=cbind(strata[,1:(3*q)],obstrans)
no_strat=strata[,(3*q)+1]
offS=(split(seq(nrow(strata[,1:(3*q)])),as.list(data.frame(strcode))))


#data=as.matrix(strata[,1:(q)])
#no_strat=strata[,(q)+1]
#offS=(split(seq(nrow(data)),as.list(data.frame(strcode))))
offn=names(offS)
newoffn=offn[offn!=base]

x=which(names(m)==base)
if(length(x)==0) into=log(m); 
if(length(x)==1) into=log(m[-x])
into[into==-Inf]=0


null=rep(0,length(into))

St=optim(into,llh,n=n,newn=newn,S=S,denomdata=denomdata,no_mstrat=no_mstrat,offn=offn,newoffn=newoffn,offS=offS,data=data,no_strat=no_strat,q=q,control=list(maxit=10000))
#ASt=optim(c(0,0,0),addllh,denomdata=denomdata,no_mstrat=no_mstrat,data=data,no_strat=no_strat,q=q,control=list(maxit=10000))
#Anullf=addllh(c(0,0,0),denomdata,no_mstrat,data,no_strat,q)
#SS[t,]=c(prop,ASt$par,ASt$value,Anullf,betas,intbeta,ASt$convergence,ASt$counts[1])


nullf=llh(null,n=n,newn=newn,S=S,denomdata=denomdata,no_mstrat=no_mstrat,offn=offn,newoffn=newoffn,offS=offS,data=data,no_strat=no_strat,q=q)

tst=2*(nullf-St$value)


R[t,]=c(prop,St$par,St$value,nullf,tst,qchisq(0.95,length(into)),as.numeric(tst>qchisq(0.95,length(into))), betas,St$convergence,St$counts[1])
}
R
}




#R=sim(1000,15000,4,c(0),1,c(0.5))
#write(apply(R,2,mean),"resfull2",append=TRUE)
#write(apply(R,2,sd),"resfull2",append=TRUE)
#bigS=sim(100,60000,4,c(-0.5,0),-0.2,2,c(0.7,0.7))

R=sim(100,60000,4,c(0,0),-0.2,2,c(0.7,0.7))
write(apply(SSR,2,mean),"resfull5",append=TRUE)
write(apply(SSR,2,sd),"resfull5",append=TRUE)

R=sim(100,60000,4,c(-0.5,0),0,2,c(0.7,0.7))
newR=as.matrix(R[,2:9])%*% t(rtrans)



write(apply(SSR,2,mean),"resfull6",append=TRUE)
write(apply(SSR,2,sd),"resfull6",append=TRUE)

R=sim(100,60000,4,c(0,0),0,2,c(0.7,0.7))
write(apply(SSR,2,mean),"resfull6",append=TRUE)
write(apply(SSR,2,sd),"resfull6",append=TRUE)


#write.table(R,"resfull4",append=TRUE)
