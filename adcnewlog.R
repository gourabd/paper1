library("Rlab")
options(error=dump.frames)

gencode<-function(compgen,p)
{
 newcode=rep(0,length(compgen[,1]))
 newcode[compgen[,p]+compgen[,p+1] == 0]=NA
 newcode[compgen[,p]+compgen[,p+1] == 1]=NA
 newcode[compgen[,p]+compgen[,p+1] == 3]=1
 newcode[compgen[,p]+compgen[,p+1] == 4]=2
 newcode
}
mating<-function(r)
{
 toss=runif(1,0,1)
 if(r[1]+r[2]==0) o=0;
 if(r[1]+r[2]==1 & toss <0.5) o=0;
 if(r[1]+r[2]==1 & toss >=0.5) o=1;
 if(r[1]+r[2]==2 & r[1]*r[2]==0) o=1;
 if (r[1]==1 & r[2]==1 & toss<0.25) o=0;
 if (r[1]==1 & r[2]==1 & toss>=0.75) o=2;
 if (r[1]==1 & r[2]==1 & toss>=0.25 & toss<0.75) o=1;
 if(r[1]+r[2]==3 & toss <0.5) o=1;
 if(r[1]+r[2]==3 & toss >=0.5) o=2;
 if(r[1]+r[2]==4) o=2;
o
}
counts<-function(compgen)
{
 str11=function(x)(x[2]==0 & x[3]+x[4]==1) 
 str12=function(x)(x[2]==1 & x[3]+x[4]==1) 
 str31=function(x)(x[2]==1 & x[3]+x[4]==3) 
 str32=function(x)(x[2]==2 & x[3]+x[4]==3)
 str21=function(x)(x[2]==0 & x[3]==1 & x[4]==1) 
 str22=function(x)(x[2]==1 & x[3]==1 & x[4]==1) 
 str23=function(x)(x[2]==2 & x[3]==1 & x[4]==1) 
 n11=sum(apply(compgen,1,str11))
 n12=sum(apply(compgen,1,str12))
 n21=sum(apply(compgen,1,str21))
 n22=sum(apply(compgen,1,str22))
 n23=sum(apply(compgen,1,str23))
 n31=sum(apply(compgen,1,str31))
 n32=sum(apply(compgen,1,str32))
 c(n11,n12,n21,n22,n23,n31,n32)
}

noninfo<-function(compgen)
{
 str11=which(compgen[,2]==0 & compgen[,3]+compgen[,4]==0) 
 str12=which(compgen[,2]==1 & compgen[,3]+compgen[,4]==2 & compgen[,3]*compgen[,4]==0) 
 str13=which(compgen[,2]==2 & compgen[,3]+compgen[,4]==4) 
 c(str11,str12,str13)
}

addllh<-function(logbeta,n,m)
{
 p0=exp(logbeta[1])
 p1=exp(logbeta[1]+logbeta[2])
 p2=exp(logbeta[1]+(2*logbeta[2]))
 
 n1= n[1]+n[2]
 n2=n[3]+n[4]+n[5]
 n3=n[6]+n[7]
 
 numer1=((n[1]+n[3])*log(p0))+ ((n[2]+n[4]+n[6])*log(p1))+ ((n[5]+n[7])*log(p2))
 denom1 = (n1*log(p0+p1))+(n2*log(p0+(2*p1)+p2))+(n3*log(p1+p2))
 part1 = -numer1+denom1

 m1= m[1]+m[2]
 m2=m[3]+m[4]+m[5]
 m3=m[6]+m[7]
 numer=((m[1]+m[3])*log(1-p0))+ ((m[2]+m[4]+m[6])*log(1-p1))+ ((m[5]+m[7])*log(1-p2))
 denom = (m1*log(2-p0-p1))+(m2*log(4-p0-(2*p1)-p2))+(m3*log(2-p1-p2))
 res=part1-numer+denom 
 return(res)
 
}

 
grllh<-function(logbeta,n,m)
{
 p0=exp(logbeta[1])
 p1=exp(logbeta[1]+logbeta[2])
 p2=exp(logbeta[1]+(2*logbeta[2]))
  
 delp00=p0; delp01=0;
 delp10=p1;delp11=p1;
 delp20=p2;delp21=2*p2;

 bb=(n[2]+n[4]+(2*n[5])+n[7])
 cc=(n[1]+n[4]+(2*n[3])+n[6])
 
 part1=0
 part1[1]=0
 part1[2]=-bb+((bb+cc)*(exp(logbeta[2])/(1+exp(logbeta[2])))) 
 m1= m[1]+m[2]
 m2=m[3]+m[4]+m[5]
 m3=m[6]+m[7]
 
 numer1=((m[1]+m[3])*delp00/(1-p0))+ ((m[2]+m[4]+m[6])*delp10/(1-p1))+ ((m[5]+m[7])*delp20/(1-p2))
 numer2=((m[1]+m[3])*delp01/(1-p0))+ ((m[2]+m[4]+m[6])*delp11/(1-p1))+ ((m[5]+m[7])*delp21/(1-p2))

 denom1=(m1*(delp00+delp10)/(2- p0-p1))+(m2*(delp00+(2*delp10)+delp20)/(4-p0-(2*p1)-p2))+(m3*(delp10+delp20)/(2-p1-p2))
 denom2=(m1*(delp01+delp11)/(2- p0-p1))+(m2*(delp01+(2*delp11)+delp21)/(4-p0-(2*p1)-p2))+(m3*(delp11+delp21)/(2-p1-p2))
 

 gr=c(part1[1]+numer1-denom1,part1[2]+numer2-denom2) 
 gr
}

part2gr<-function(logbeta1,logbeta0,n,m)
{
 grllh(c(logbeta0,logbeta1),n,m)[2]
} 

part2lh<-function(logbeta1,logbeta0,n,m)
{
 res=addllh(c(logbeta0,logbeta1),n,m)
 attr(res, "gradient") <- part2gr(logbeta1,logbeta0,n,m)
 return(res)
}



subscore<-function(data)
{
 x=data[1]
 p1=data[2]
 p2=data[3]
 y=data[4]
 z=data[5]
 logbeta=c(data[6],data[7])
 P0=exp(logbeta[1])
 P1=exp(logbeta[1]+logbeta[2])
 P2=exp(logbeta[1]+(2*logbeta[2]))

 delp01=0; delp11=P1; delp21=2*P2;
 R=0
 if (x==0 & p1+p2==1 & y==1) R=(delp01/(P0))-((delp01+delp11)/(P0+P1))
 if (x==1 & p1+p2==1 & y==1) R=(delp11/(P1))-((delp01+delp11)/(P0+P1))
 if (x==1 & p1+p2==3 & y==1) R=(delp11/(P1))-((delp11+delp21)/(P1+P2))
 if (x==2 & p1+p2==3 & y==1) R=(delp21/(P2))-((delp11+delp21)/(P1+P2))
 if (x==0 & p1==1 & p2==1 & y==1)R=(delp01/(P0))-((delp01+(2*delp11)+delp21)/(P0+(2*P1)+P2))
 if (x==1 & p1==1 & p2==1 & y==1)R=(delp11/(P1))-((delp01+(2*delp11)+delp21)/(P0+(2*P1)+P2))
 if (x==2 & p1==1 & p2==1 & y==1)R=(delp21/(P2))-((delp01+(2*delp11)+delp21)/(P0+(2*P1)+P2))

 if (x==0 & p1+p2==1 & y==0) R= -(delp01/(1-P0))+((delp01+delp11)/(2- P0-P1))
 if (x==1 & p1+p2==1 & y==0) R=-(delp11/(1-P1))+((delp01+delp11)/(2- P0-P1))
 if (x==1 & p1+p2==3 & y==0) R=-(delp11/(1-P1))+((delp11+delp21)/(2-P1-P2))
 if (x==2 & p1+p2==3 & y==0) R=-(delp21/(1-P2))+((delp11+delp21)/(2-P1-P2))
 if (x==0 & p1==1 & p2==1 & y==0)R=-(delp01/(1-P0))+((delp01+(2*delp11)+delp21)/(4-P0-(2*P1)-P2))
 if (x==1 & p1==1 & p2==1 & y==0)R=-(delp11/(1-P1))+((delp01+(2*delp11)+delp21)/(4-P0-(2*P1)-P2))
 if (x==2 & p1==1 & p2==1 & y==0)R=-(delp21/(1-P2))+((delp01+(2*delp11)+delp21)/(4-P0-(2*P1)-P2))
 
 return(R)
}

q=100
R=matrix(0,q,8)
S=matrix(0,q,4)
grad=matrix(0,q,2)



simalgo<-function(sigma,W1,W2,Exp1,Exp2,base1,base2,basex1,basex2,main,inter)
{
#set.seed(runif(1,50000,99999))
R=0
S=0
grad=0
#W1=0.2
#W2=0.2
w1=c((1-W1)^2,2*W1*(1-W1),W1^2)
w2=c((1-W2)^2,2*W2*(1-W2),W2^2)
#Exp1=Exp2=0.5
expi1=c(1-Exp1,Exp1)
expi2=c(1-Exp2,Exp2)
n=500
expg=as.numeric(runif(n,0,1)<0.5)
n1=length(which(expg==1))
n2=n-n1

simgf1=sample(c(0,1,2),length(which(expg==1)),replace=TRUE,w1)
simgf2=sample(c(0,1,2),n-length(which(expg==1)),replace=TRUE,w2)
simgm1=sample(c(0,1,2),length(which(expg==1)),replace=TRUE,w1)
simgm2=sample(c(0,1,2),n-length(which(expg==1)),replace=TRUE,w2)
simparsibn1=1+rpois(length(which(expg==1)),0.7)
simparsibn2=1+rpois(n-length(which(expg==1)),0.7)

gpar1=cbind((1:n1),simgf1,simgm1)
simgpar1=gpar1[rep(1:n1,simparsibn1),]
gpar2=cbind(((n1+1):n),simgf2,simgm2)
simgpar2=gpar2[rep(1:n2,simparsibn2),]

sibpars1=apply(simgpar1[,2:3],1,mating)

sibpars2=apply(simgpar2[,2:3],1,mating)

simoldpar=cbind(c(simgpar1[,1],simgpar2[,1]),c(sibpars1,sibpars2),rbind(simgpar1[,2:3],simgpar2[,2:3]))



sibparm1=sample(c(0,1,2),length(sibpars1),replace=TRUE,w1)
sibparm2=sample(c(0,1,2),length(sibpars2),replace=TRUE,w2)

simsibn1=1+rpois(length(sibparm1),0.7)
simsibn2=1+rpois(length(sibparm2),0.7)

bigpid1=simgpar1[,1]
bigpid2=simgpar2[,1]

par1=cbind(((n+1):(n+length(sibpars1))),sibpars1,sibparm1,bigpid1)
par2=cbind(((n+1+length(sibpars1)):(n+length(sibpars1)+length(sibpars2))),sibpars2,sibparm2,bigpid2)
simpar1=par1[rep(1:length(sibpars1),simsibn1),]
simpar2=par2[rep(1:length(sibpars2),simsibn2),]

off1=apply(simpar1[,2:3],1,mating)
off2=apply(simpar2[,2:3],1,mating)

id=c(simoldpar[,1],simpar1[,1],simpar2[,1])
bigid=c(simoldpar[,1],simpar1[,4],simpar2[,4])
off=c(simoldpar[,2],off1,off2)
parent=rbind(simoldpar[,3:4],simpar1[,2:3],simpar2[,2:3])
villold=c(rep(1,length(simgpar1[,1])),rep(0,length(simgpar2[,1])))
vill=c(villold,rep(1,length(off1)),rep(0,length(off2)))

simsibn=c(simparsibn1,simparsibn2,simsibn1,simsibn2)

simold=cbind(id,off,parent)

#expo=sample(c(0,1),length(off),replace=TRUE)
expo=rep(0,length(vill))
expo[vill==1]=sample(c(0,1),sum(vill==1),replace=TRUE,prob=expi1)
expo[vill==0]=sample(c(0,1),sum(vill==0),replace=TRUE,prob=expi2)

#noise=rnorm(length(unique(id)),0,sigma)
noise= sample(c(-sigma,sigma),length(unique(id)),prob=c(0.5,0.5),replace=TRUE)
simnoise=exp(rep(noise,simsibn))

prob=rep(0,length(vill))
#base1=base2=-0.555
#main=-0.396
#basex1=basex2=-0.198
#inter=-1.153
prob[vill==1]=pmax(pmin((exp(base1+(main*off[vill==1])+(basex1*expo[vill==1])+(inter*off[vill==1]*expo[vill==1])))*simnoise[vill==1],1),0)
prob[vill==0]=pmax(pmin((exp(base2+(main*off[vill==0])+(basex2*expo[vill==0])+(inter*off[vill==0]*expo[vill==0]))*simnoise[vill==0]),1),0)

uni=runif(length(prob),0,1)
aff=rbern(length(prob),prob)
table(aff,expo)
table(aff,off)

simdata=cbind(simold,aff,expo)

#idold=unique(simoldpar[,1])
#miss=sample(idold,40,replace=FALSE)

#missdata=which(simdata[,1] %in% miss)
#simdata=simdata[-missdata,]

off=simdata[,2]
aff=simdata[,5]
expo=simdata[,6]


v1=(glm(aff~off+expo+off*expo,family=binomial(link="log")))
co=v1$coefficients


off1=off[noninfo(simdata)]
aff1=aff[noninfo(simdata)]
expo1=expo[noninfo(simdata)]
v2=v1
v2=(glm(aff1~off1+expo1+off1*expo1,family=binomial(link="log")))
co1=v2$coefficients

#v1=(glm(aff~off+vill+off*vill,family=binomial(link="log")))
#co=(glm(aff~off+vill+off*vill,family=binomial(link="log")))$coefficients

b00=co1[1]
be=co1[3]
b1=co[2]
bint=co[4]


nsimaff=counts(simdata[simdata[,5]==1,])
nsimunaff=counts(simdata[simdata[,5]==0,])
nsimaff1=counts(simdata[simdata[,5]==1 & simdata[,6]==1,])
nsimunaff1=counts(simdata[simdata[,5]==0 & simdata[,6]==1,])
nsimaff0=counts(simdata[simdata[,5]==1 & simdata[,6]==0,])
nsimunaff0=counts(simdata[simdata[,5]==0 & simdata[,6]==0,])

n=nsimaff1;
bb=(n[2]+n[4]+(2*n[5])+n[7])
cc=(n[1]+n[4]+(2*n[3])+n[6])
r1=log(bb/cc)
var1=(1/bb)+(1/cc)
n=nsimaff0;
bb=(n[2]+n[4]+(2*n[5])+n[7])
cc=(n[1]+n[4]+(2*n[3])+n[6])
r2=log(bb/cc)
var2=(1/bb)+(1/cc)

addx1=nlm(part2lh,r1,b00+be,nsimaff1,nsimunaff1,hessian=TRUE)
beta1e=addx1$estimate


addx0=nlm(part2lh,r2,b00,nsimaff0,nsimunaff0,hessian=TRUE)
beta1=addx0$estimate
betaint=beta1e-beta1

param1=param2=rep(0,length(simdata[,1]))
param1[which(expo==0)]=b00;param1[which(expo==1)]=b00+be;
param2[which(expo==0)]=beta1;param2[which(expo==1)]=beta1e;
out=cbind(simdata,param1,param2,bigid)
score=apply(out[,2:8],1,subscore)
score1=rep(0,length(out[,1]));score1[which(expo==0)]=score[which(expo==0)]
score2=rep(0,length(out[,1]));score2[which(expo==1)]=score[which(expo==1)]
familyscore1=aggregate(data.frame(out[,1],score1),list(pid=out[,1]),sum)$score1
familyscore2=aggregate(data.frame(out[,1],score2),list(pid=out[,1]),sum)$score2

pedscore1=aggregate(data.frame(out[,9],score1),list(pid=out[,9]),sum)$score1
pedscore2=aggregate(data.frame(out[,9],score2),list(pid=out[,9]),sum)$score2

familyscore=cbind(familyscore1,familyscore2)

pedscore=cbind(pedscore1,pedscore2)


G= t(familyscore)%*%familyscore

pedG= t(pedscore)%*%pedscore


D=diag(c(1/addx0$hessian,1/addx1$hessian))
var=D%*%G%*%D
pedvar=D%*%pedG%*%D

affdata=simdata[simdata[,5]==1,]
bigaffid=bigid[simdata[,5]==1]
param1=param2=rep(0,length(affdata[,1]))
param1[which(affdata[,6]==0)]=b00;param1[which(affdata[,6]==1)]=b00+be;
param2[which(affdata[,6]==0)]=r2;param2[which(affdata[,6]==1)]=r1;
out=cbind(affdata,param1,param2,bigaffid)
affscore=apply(out[,2:8],1,subscore)
affscore1=rep(0,length(affdata[,1]));affscore1[which(affdata[,6]==0)]=score[which(affdata[,6]==0)]
affscore2=rep(0,length(affdata[,1]));affscore2[which(affdata[,6]==1)]=score[which(affdata[,6]==1)]
afffamilyscore1=aggregate(data.frame(out[,1],affscore1),list(out[,1]),sum)$affscore1
afffamilyscore2=aggregate(data.frame(out[,1],affscore2),list(out[,1]),sum)$affscore2

affpedscore1=aggregate(data.frame(out[,9],affscore1),list(pid=out[,9]),sum)$affscore1
affpedscore2=aggregate(data.frame(out[,9],affscore2),list(pid=out[,9]),sum)$affscore2


afffamilyscore=cbind(afffamilyscore1,afffamilyscore2)
affpedscore=cbind(affpedscore1,affpedscore2)

affG= t(afffamilyscore)%*%afffamilyscore
affpedG= t(affpedscore)%*%affpedscore


affD=diag(c(var1,var2))
affvar=affD%*%affG%*%affD
affpedvar=affD%*%affpedG%*%affD


R=c(b00,be,b1,bint,beta1,betaint,r2,r1-r2)
S=c((1/addx1$hessian)+(1/addx0$hessian),var1+var2,var[2,2]+var[1,1]-(2*var[1,2]),affvar[2,2]+affvar[1,1]-(2*affvar[1,2]),pedvar[2,2]+pedvar[1,1]-(2*pedvar[1,2]),affpedvar[2,2]+affpedvar[1,1]-(2*affpedvar[1,2]))
grad=c(addx1$gradient,addx0$gradient)
prob=mean(aff)
res=c(R,S,grad,prob)
res
}

rm(res); rm(V); rm(meanV) ;rm(varV);rm(var2V);rm(R1);rm(R2);rm(S1);rm(S2)
#res <- lapply(1:5000, function(i) try(simalgo(0.0,0.1,0.1,0.5,0.5,-0.499,-0.499,-0.1152,-0.1152,-0.202,-0.283), TRUE))

res <- lapply(1:5000, function(i) try(simalgo(0.0,0.1,0.3,0.2,0.5,-0.499,-0.222,-0.1152,-0.333,-0.202,-0.283), TRUE))
#res <- lapply(1:500, function(i) try(simalgo(0.2,0.5,0.5,0.5,0.5,-0.555,-0.555,-0.198,-0.198,-0.396,-1.153), TRUE))

V=(t(data.frame(res[sapply(res, function(x) !inherits(x, "try-error"))])))
meanV=apply(V[,1:8],2,mean)
varV=apply(V[,1:8],2,var)
var2V=apply(V[,9:14],2,mean)
summary(V[,17])



upperV=meanV[c(6,8)]+(1.96*sqrt(varV[c(6,8)]/length(V[,1])))
lowerV=meanV[c(6,8)]-(1.96*sqrt(varV[c(6,8)]/length(V[,1])))
upperV1=meanV[c(6,8)]+(1.96*sqrt(var2V[c(1,2)]/length(V[,1])))
upperV2=meanV[c(6,8)]+(1.96*sqrt(var2V[c(3,4)]/length(V[,1])))
upperV3=meanV[c(6,8)]+(1.96*sqrt(var2V[c(5,6)]/length(V[,1])))

lowerV1=meanV[c(6,8)]-(1.96*sqrt(var2V[c(1,2)]/length(V[,1])))
lowerV2=meanV[c(6,8)]-(1.96*sqrt(var2V[c(3,4)]/length(V[,1])))
lowerV3=meanV[c(6,8)]-(1.96*sqrt(var2V[c(5,6)]/length(V[,1])))


summary(V[,15])
summary(V[,16])

varV
meanV

R1=as.numeric((V[,6]-(1.96*sqrt(V[,9]))) < -0.283&   (V[,6]+(1.96*sqrt(V[,9]))) > -0.283)
R2=as.numeric((V[,6]-(1.96*sqrt(V[,11]))) < -0.283 &   (V[,6]+(1.96*sqrt(V[,11]))) > -0.283)
R3=as.numeric((V[,6]-(1.96*sqrt(V[,13]))) < -0.283 &   (V[,6]+(1.96*sqrt(V[,13]))) > -0.283)

S1=as.numeric((V[,8]-(1.96*sqrt(V[,10]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,10]))) > -0.283)
S2=as.numeric((V[,8]-(1.96*sqrt(V[,12]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,12]))) > -0.283)
S3=as.numeric((V[,8]-(1.96*sqrt(V[,14]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,14]))) > -0.283)


c(abs(meanV[c(8)]+0.416)/0.416,varV[8],var2V[2],var2V[4],(lowerV[2]< -0.416 & upperV[2]>-0.416),(lowerV1[2]< -0.416 & upperV1[2]>-0.416),(lowerV2[2]< -0.416 & upperV2[2]>-0.416),mean(S1),mean(S2))
c(abs(meanV[c(6)]+0.416)/0.416,varV[6],var2V[1],var2V[3],(lowerV[1]< -0.416 & upperV[1]>-0.416),(lowerV1[1]< -0.416 & upperV1[1]>-0.416),(lowerV2[1]< -0.416 & upperV2[1]>-0.416),mean(R1),mean(R2))



write.table(V,"Vmixedlog0af01")

write(c(abs(meanV[c(8)]+0.283)/0.283,varV[8],var2V[2],var2V[4],var2V[6],(lowerV[2]< -0.283 & upperV[2]>-0.283),(lowerV1[2]< -0.283 & upperV1[2]>-0.283),(lowerV2[2]< -0.283 & upperV2[2]>-0.283),(lowerV3[2]< -0.283 & upperV3[2]>-0.283),mean(S1),mean(S2),mean(S3),mean(V[,17])),"adcmixedlog0af01",ncolumns=13)
write(c(abs(meanV[c(6)]+0.283)/0.283,varV[6],var2V[1],var2V[3],var2V[5],(lowerV[1]< -0.283 & upperV[1]>-0.283),(lowerV1[1]< -0.283 & upperV1[1]>-0.283),(lowerV2[1]< -0.283 & upperV2[1]>-0.283),(lowerV3[1]< -0.283 & upperV3[1]>-0.283),mean(R1),mean(R2),mean(R3),mean(V[,17])),"adcmixedlog0af01",ncolumns=13,append=TRUE)


rm(res); rm(V); rm(meanV) ;rm(varV);rm(var2V);rm(R1);rm(R2);rm(S1);rm(S2)

#res <- lapply(1:5000, function(i) try(simalgo(0.2,0.1,0.1,0.5,0.5,-0.499,-0.499,-0.1152,-0.1152,-0.202,-0.283), TRUE))

res <- lapply(1:5000, function(i) try(simalgo(0.2,0.1,0.3,0.2,0.5,-0.499,-0.222,-0.1152,-0.333,-0.202,-0.283), TRUE))


V=(t(data.frame(res[sapply(res, function(x) !inherits(x, "try-error"))])))
meanV=apply(V[,1:8],2,mean)
varV=apply(V[,1:8],2,var)
var2V=apply(V[,9:14],2,mean)
summary(V[,17])



upperV=meanV[c(6,8)]+(1.96*sqrt(varV[c(6,8)]/length(V[,1])))
lowerV=meanV[c(6,8)]-(1.96*sqrt(varV[c(6,8)]/length(V[,1])))
upperV1=meanV[c(6,8)]+(1.96*sqrt(var2V[c(1,2)]/length(V[,1])))
upperV2=meanV[c(6,8)]+(1.96*sqrt(var2V[c(3,4)]/length(V[,1])))
upperV3=meanV[c(6,8)]+(1.96*sqrt(var2V[c(5,6)]/length(V[,1])))

lowerV1=meanV[c(6,8)]-(1.96*sqrt(var2V[c(1,2)]/length(V[,1])))
lowerV2=meanV[c(6,8)]-(1.96*sqrt(var2V[c(3,4)]/length(V[,1])))
lowerV3=meanV[c(6,8)]-(1.96*sqrt(var2V[c(5,6)]/length(V[,1])))


summary(V[,15])
summary(V[,16])

varV
meanV

R1=as.numeric((V[,6]-(1.96*sqrt(V[,9]))) < -0.283&   (V[,6]+(1.96*sqrt(V[,9]))) > -0.283)
R2=as.numeric((V[,6]-(1.96*sqrt(V[,11]))) < -0.283 &   (V[,6]+(1.96*sqrt(V[,11]))) > -0.283)
R3=as.numeric((V[,6]-(1.96*sqrt(V[,13]))) < -0.283 &   (V[,6]+(1.96*sqrt(V[,13]))) > -0.283)

S1=as.numeric((V[,8]-(1.96*sqrt(V[,10]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,10]))) > -0.283)
S2=as.numeric((V[,8]-(1.96*sqrt(V[,12]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,12]))) > -0.283)
S3=as.numeric((V[,8]-(1.96*sqrt(V[,14]))) < -0.283 &   (V[,8]+(1.96*sqrt(V[,14]))) > -0.283)

write.table(V,"Vmixedlog1af01")

write(c(abs(meanV[c(8)]+0.283)/0.283,varV[8],var2V[2],var2V[4],var2V[6],(lowerV[2]< -0.283 & upperV[2]>-0.283),(lowerV1[2]< -0.283 & upperV1[2]>-0.283),(lowerV2[2]< -0.283 & upperV2[2]>-0.283),(lowerV3[2]< -0.283 & upperV3[2]>-0.283),mean(S1),mean(S2),mean(S3),mean(V[,17])),"adcmixedlog1af01",ncolumns=13)
write(c(abs(meanV[c(6)]+0.283)/0.283,varV[6],var2V[1],var2V[3],var2V[5],(lowerV[1]< -0.283 & upperV[1]>-0.283),(lowerV1[1]< -0.283 & upperV1[1]>-0.283),(lowerV2[1]< -0.283 & upperV2[1]>-0.283),(lowerV3[1]< -0.283 & upperV3[1]>-0.283),mean(R1),mean(R2),mean(R3),mean(V[,17])),"adcmixedlog1af01",ncolumns=13,append=TRUE)

