library(matrixStats)
sett=out$tnow


ALStochModelV4SetT=function(time){
  tnow=c(0)
  ABival=c(1*10^-6)
  ABoval=c(1*10^-8)
  Tauval=c(1.37*10^(-10))
  Fival=c(3.36*10^(-10))
  Foval=c(3.36*10^(-11))
  Aval=c(.14)
  M1val=c(.02)
  M2val=c(.02)
  Nval=c(.14)
  R0=6
  i=2
  while(tnow[i-1]<time) {
    #if(tnow<=100){
    # R=R0*(tnow/100)
    #}else(R=R0)
    chainlist=c(((3.4*10^-4)/(3.36*10^-10))/Nval[i-1],((((9.51*10^-6)*7)-9.51)/.14)/ABival[i-1],
                (((.05/365)+(8*10^-9)/.14+(8*10^-10)/.14+((2*10^-3*.9))*1/(7*10^-3)))/ABoval[i-1],((((8.1*10^-11)+((1.35*10^-11)*6))+.277)/.14)/Tauval[i-1],
                (((1.662*10^-3)+(2.77*10^-3))/.14)/Fival[i-1],((.05/365)+(2.77*10^-4))/Foval[i-1],(1.793+ .045+((1.2*10^(-3))))/Aval[i-1],
                ((.047*((2*10^(-2))/(2.58*10^(-11))))*.9+(.015))/M1val[i-1],((.047*((2*10^(-2))/((2.58*10^(-11)))))*.09+(.015))/M2val[i-1])
    chainorder=c("Neurons","ABi","ABo","Tau","Fi","Fo","Astro","M1","M2")
    chain=data.frame(((3.4*10^-4)/(3.36*10^-10))/Nval[i-1],((((9.51*10^-6)*7)-9.51)/.14)/ABival[i-1],
                     (((.05/365)+(8*10^-9)/.14+(8*10^-10)/.14+((2*10^-3*.9))*1/(7*10^-3)))/ABoval[i-1],((((8.1*10^-11)+((1.35*10^-11)*6))+.277)/.14)/Tauval[i-1],
                     (((1.662*10^-3)+(2.77*10^-3))/.14)/Fival[i-1],((.05/365)+(2.77*10^-4))/Foval[i-1],(1.793+ .045+((1.2*10^(-3))))/Aval[i-1],
                     ((.047*((2*10^(-2))/(2.58*10^(-11))))*.9+(.015))/M1val[i-1],((.047*((2*10^(-2))/((2.58*10^(-11)))))*.09+(.015))/M2val[i-1])
    colnames(chain)=chainorder
    u=runif(length(chainlist))
    possiblet=sett[i]-tnow[i-1]
    #print(possiblet)
    tchange=max(possiblet)
    chain2=as.data.frame(t(chain))
    chain2=arrange(chain2)
    for (j in 0:length(chain2[,1])) {
      if(isTRUE(chain2[j,]==chain$Neurons)){
        Nval[i]=NeuronChain(Nval[i-1],Fival[i-1],Foval[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$ABi)){
        ABival[i]=ABetaiChain(Nval[i-1],ABival[i-1],tnow[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$ABo)){
        ABoval[i]=ABetaochain(ABoval[i-1],ABival[i-1],Nval[i-1],Aval[i-1],M1val[i-1],M2val[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$Tau)){
        Tauval[i]=TauChainV2(Tauval[i-1],Fival[i-1],Foval[i-1],Nval[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$Fi)){
        Fival[i]=FivalChainV2(Fival[i-1],Tauval[i-1],Nval[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$Fo)){
        Foval[i]=FovalChainV2(Fival[i-1],Foval[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$Astro)){
        Aval[i]=AvalChainV2(ABoval[i-1],M1val[i-1],Aval[i-1],max(possiblet),chainlist)
      } else if(isTRUE(chain2[j,]==chain$M1)){
        M1val[i]=M1valChainV2(M1val[i-1],Foval[i-1],max(possiblet),chainlist)
      } else{
        M2val[i]=M2valChainV2(M2val[i-1],Foval[i-1],max(possiblet),chainlist)
      }
      
    }
    tnow[i]=tnow[i-1]+possiblet
    i=i+1
  }
  return(data.frame(Nval,ABival,ABoval,Tauval,Fival,Foval,M1val,M2val,Aval,tnow))
  
}
outSett=ALStochModelV4SetT(3650)

##Neuron Variance Between Runs
nValVarDF=data.frame(outSett$Nval)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  nValVarDF=cbind(nValVarDF,out$Nval)
}


nValVarDF$rowVar=rowVars(as.matrix(nValVarDF))
nvalvarmean=mean(nValVarDF$rowVar)
plot(nValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Neuron Value Variance",col='brown3')
abline(h=nvalvarmean)

##ABi Variance Between Runs
ABiValVarDF=data.frame(outSett$ABival)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  ABiValVarDF=cbind(ABiValVarDF,out$ABival)
}


ABiValVarDF$rowVar=rowVars(as.matrix(ABiValVarDF))
ABivalvarmean=mean(ABiValVarDF$rowVar)
plot(ABiValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Amyloid Beta Inside Value Variance",col='blue')
##ABo Variance Between Runs
ABoValVarDF=data.frame(outSett$ABoval)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  ABoValVarDF=cbind(ABoValVarDF,out$ABoval)
}


ABoValVarDF$rowVar=rowVars(as.matrix(ABoValVarDF))
ABovalvarmean=mean(ABoValVarDF$rowVar)
plot(ABoValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Amyloid Beta Outside Value Variance",col='chartreuse4')

##Tau Variance Between Runs
TauValVarDF=data.frame(outSett$Tauval)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  TauValVarDF=cbind(TauValVarDF,out$Tauval)
}


TauValVarDF$rowVar=rowVars(as.matrix(TauValVarDF))
Tauvalvarmean=mean(TauValVarDF$rowVar)
plot(TauValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Tau Value Variance",col='deeppink')


##Fi Variance Between Runs
FiValVarDF=data.frame(outSett$Fival)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  FiValVarDF=cbind(FiValVarDF,out$Fival)
}


FiValVarDF$rowVar=rowVars(as.matrix(FiValVarDF))
Fivalvarmean=mean(FiValVarDF$rowVar)
plot(FiValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "NFT Inside Value Variance",col='darkgoldenrod1')


##Fo Variance Between Runs
FoValVarDF=data.frame(outSett$Foval)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  FoValVarDF=cbind(FoValVarDF,out$Foval)
}


FoValVarDF$rowVar=rowVars(as.matrix(FoValVarDF))
Fovalvarmean=mean(FoValVarDF$rowVar)
plot(FoValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "NFT Outside Value Variance",col='darkorange')


##M1 Variance Between Runs
M1ValVarDF=data.frame(outSett$M1val)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  M1ValVarDF=cbind(M1ValVarDF,out$M1val)
}


M1ValVarDF$rowVar=rowVars(as.matrix(M1ValVarDF))
M1valvarmean=mean(M1ValVarDF$rowVar)
plot(M1ValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Microglia Type 1 Value Variance",col='darkmagenta')


##M2 Variance Between Runs
M2ValVarDF=data.frame(outSett$M2val)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  M2ValVarDF=cbind(M2ValVarDF,out$M2val)
}


M2ValVarDF$rowVar=rowVars(as.matrix(M2ValVarDF))
M2valvarmean=mean(M2ValVarDF$rowVar)
plot(M2ValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Microglia Type 2 Value Variance",col='aquamarine')


##Astrocytes Variance Between Runs
AValVarDF=data.frame(outSett$Aval)
#plot(outSett$tnow,outSett$Nval)
for (i in 2:25) {
  out=ALStochModelV4SetT(3650)
  AValVarDF=cbind(AValVarDF,out$Aval)
}


AValVarDF$rowVar=rowVars(as.matrix(AValVarDF))
Avalvarmean=mean(AValVarDF$rowVar)
plot(AValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Astrocyte Value Variance",col='burlywood4')


##All Graphs
par(mfrow=c(3,3))
plot(nValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Neuron Value Variance",col='brown3')
plot(ABiValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Amyloid Beta Inside Value Variance",col='blue')
plot(ABoValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Amyloid Beta Outside Value Variance",col='chartreuse4')
plot(TauValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Tau Value Variance",col='deeppink')
plot(FiValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "NFT Inside Value Variance",col='darkgoldenrod1')
plot(FoValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "NFT Outside Value Variance",col='darkorange')
plot(M1ValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Microglia Type 1 Value Variance",col='darkmagenta')
plot(M2ValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Microglia Type 2 Value Variance",col='aquamarine')
plot(AValVarDF[2:7918,26],xlab = "Time Points",ylab = "Variance",main = "Astrocyte Value Variance",col='burlywood4')

