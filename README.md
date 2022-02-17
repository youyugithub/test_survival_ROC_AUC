# test_survival_ROC_AUC

```

nind<-1000
Fi<-rexp(nind,0.5)
Ci<-rexp(nind,0.5)
Ti<-pmin(Fi,Ci)
delta<-Fi<=Ci
prediction<-rnorm(nind)-Fi
time<-1

survtime<-sort(Ti[delta])
hazard<-rep(NA,length(survtime))
survival<-rep(NA,length(survtime))
hS_<-rep(NA,length(survtime))
for(ttidx in 1:length(survtime)){
  tt<-survtime[ttidx]
  hazard[ttidx]<-sum(Ti==tt&delta)/sum(Ti>=tt)
  if(ttidx==1){
    survival[ttidx]<-1-hazard[ttidx]
    hS_[ttidx]<-hazard[ttidx]
  }else{
    survival[ttidx]<-survival[ttidx-1]*(1-hazard[ttidx])
    hS_[ttidx]<-hazard[ttidx]*survival[ttidx-1]
  }
}


K_vec<-c(-Inf,sort(prediction),Inf)

TP<-matrix(NA,length(survtime),length(K_vec))
FN<-matrix(NA,length(survtime),length(K_vec))

for(ttidx in 1:length(survtime)){
  tt<-survtime[ttidx]
  for(Kidx in 1:length(K_vec)){
    K<-K_vec[Kidx]
    TP[ttidx,Kidx]<-sum(Ti==tt&delta&prediction>K)
    FN[ttidx,Kidx]<-sum(Ti==tt&delta&prediction<=K)
  }
}


# check
# sum(hS_[1:20])-(1-survival[20])

ttidx<-sum(survtime<=time)

sen_vec<-rep(NA,length(K_vec))
spe_vec<-rep(NA,length(K_vec))
for(Kidx in 1:length(K_vec)){
  K<-K_vec[Kidx]
  sen_vec[Kidx]<-sum(hS_[1:ttidx]*TP[1:ttidx,Kidx]/(TP[1:ttidx,Kidx]+FN[1:ttidx,Kidx]))/sum(hS_[1:ttidx])
  spe_vec[Kidx]<-(sum(prediction<=K)/nind-sum(hS_[1:ttidx]*FN[1:ttidx,Kidx]/(TP[1:ttidx,Kidx]+FN[1:ttidx,Kidx])))/survival[ttidx]
}

range(sen_vec)
range(spe_vec)
plot(sen_vec,spe_vec)

simple_auc <- function(FPR,TPR){
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

simple_auc(rev(1-spe_vec),rev(sen_vec))
library(timeROC)
library(survival)
temp<-timeROC(Ti,delta,prediction,cause=1,times=time,ROC=TRUE)
plot(temp,time=time)


tempS_<-c(1,survival[1:(ttidx-1)])
hazard[1:ttidx]*(1-hazard[1:ttidx])*tempS_^2
hazard[1:ttidx]*(1-tempS_)*tempS_
(survival[ttidx]*(1-survival[ttidx]))




library(survAUC)

AUC.cd(Surv(Ti,delta),Surv(Ti,delta),lp=prediction,lpnew=prediction,times=time)
AUC.hc(Surv(Ti,delta),Surv(Ti,delta),lpnew=prediction,times=time)
AUC.sh(Surv(Ti,delta),Surv(Ti,delta),lp=prediction,lpnew=prediction,times=time)







(sen_true<-sum(prediction>K&Fi<time)/sum(Fi<time))
(spe_true<-sum(prediction<=K&Fi>=time)/sum(Fi>=time))

print(c(sen,sen_true))
print(c(spe,spe_true))


##########
# Chiang #
##########

survtime<-sort(unique(Ti))
hazard<-rep(NA,length(survtime))
survival<-rep(NA,length(survtime))
hazard_c<-rep(NA,length(survtime))
survival_c<-rep(NA,length(survtime))
for(ttidx in 1:length(survtime)){
  tt<-survtime[ttidx]
  hazard[ttidx]<-sum(Ti==tt&delta)/sum(Ti>=tt)
  hazard_c[ttidx]<-sum(Ti==tt&!delta)/sum(Ti>=tt)
  if(ttidx==1){
    survival[ttidx]<-1-hazard[ttidx]
    survival_c[ttidx]<-1-hazard_c[ttidx]
  }else{
    survival[ttidx]<-survival[ttidx-1]*(1-hazard[ttidx])
    survival_c[ttidx]<-survival_c[ttidx-1]*(1-hazard_c[ttidx])
  }
}

weight<-rep(NA,nind)
for(ii in 1:nind){
  if(Ti[ii]<=time&delta[ii]){
    ttidx<-sum(survtime<=Ti[ii])
    weight[ii]<-1/survival_c[ttidx]
  }else if(Ti[ii]>time){
    ttidx<-sum(survtime<=time)
    weight[ii]<-1/survival_c[ttidx]
  }
}

K_vec<-c(-Inf,sort(prediction),Inf)
sen_vec2<-rep(NA,length(K_vec))
spe_vec2<-rep(NA,length(K_vec))
for(Kidx in 1:length(K_vec)){
  K<-K_vec[Kidx]
  sen_vec2[Kidx]<-sum(weight[prediction>K&Ti<=time&delta])/sum(weight[Ti<=time&delta])
  spe_vec2[Kidx]<-sum(weight[prediction<=K&Ti>time])/sum(weight[Ti>time])
}


simple_auc(rev(FPR_vec),rev(TPR_vec))
temp<-timeROC(Ti,delta,prediction,cause=1,weighting="marginal",times=time,ROC=TRUE)
plot(temp,time=time)
lines(1-spe_vec2,sen_vec2,type="l")

c(TPR,sen,sen_true)
c(1-FPR,spe,spe_true)

```
