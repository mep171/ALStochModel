library(dplyr)


NeuronChain=function(Nval,Fival,Foval,tchange,chainl){
  uN=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uN<(rnorm(1,((3.4*10^-4)/(3.36*10^-10))))){#/sum(chainl)))){
    Nval=Nval-((3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*tchange*Nval
  } else{
    Nval=Nval
  } 
  return(Nval)
}

ABetaiChain=function(Nval,ABival,tnow,tchange,chainl){
  uABi1=runif(1)
  if(tnow<=100){
    R=6*tnow/100
  }else{
    R=6
  }
  if((uABi1<rnorm(1,((((9.51)*1+R))-(9.51))/.14))){#/sum(chainl)))){#&(ABichangeboth>=0)){
    ABival=ABival+(((9.51*10^-6*(1+R)-(9.51*ABival))*tchange))*(Nval/.14)
  }
  else{
    ABival=ABival
  } 
  return(ABival)
}

ABetaochain=function(ABoval,ABival,Nval,Aval,M1val,M2val,tchange,chainl){
  uABo=runif(1)
  if(uABo<rnorm(1,((.05/365)+(8*10^-9)/.14+(8*10^-10)/.14-((2*10^-3*.9))*1/(7*10^-3)))){#/sum(chainl))){#((((.05/365)/sum(chainlist)/(8*10^-9)/((8*10^-10)))/sum(chainlist)))/(((2*10^-3)/(7*10^-3))/sum(chainlist))){
    #ABoval=ABoval+(ABival*(.05/365)+((8*10^-9)*(Nval/.14))+((8*10^-10)*(Aval/.14))-(((2*10^-3)*(((M1val+.9*M2val))))*(ABoval/(ABoval+(7*10^-3)))))
    ABoval=ABoval+(ABival*0.0001369863+(8*10^-9)*(Nval/.14)+(8*10^-10)*Aval/.14-((2*10^-3)*(M1val+.9*M2val))*(ABoval/(ABoval+(7*10^-3))))*tchange
  } else{
    ABoval=ABoval
  }
  return(ABoval)
}

TauChainV2=function(Tauval,Fival,Foval,Nval,tchange,chainl){
  uTau=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uTau<rnorm(1,((((8.1*10^-11)+((1.35*10^-11)*6))-.277)/.14))){#/sum(chainl))){
    Tauval=Tauval+(((8.1*10^-11)+(1.35*10^-11)*6)-.277*Tauval)*tchange*(Nval/.14)
  } else{
    Tauval=Tauval
  } 
  return(Tauval)
}

FivalChainV2=function(Fival,Tauval,Nval,tchange,chainl){
  uFi=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uFi<rnorm(1,(((1.662*10^-3)-(2.77*10^-3))/.14))){#/sum(chainl))){
    Fival=Fival+(((1.662*10^-3)*Tauval)*(Nval/.14)-((2.77*10^-3)*Fival)*tchange*(Nval/.14))
  } else{
    Fival=Fival
  } 
  return(Fival)
}

FovalChainV2=function(Fival,Foval,tchange,chainl){
  uFo=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uFo<rnorm(1,((.05/365)-(2.77*10^-4)))){#/sum(chainl))){
    Foval=Foval+(((.05/365)*Fival)-((2.77*10^-4)*Foval))*tchange
  } else{
    Foval=Foval
  } 
  return(Foval)
}


AvalChainV2=function(ABoval,M1val,Aval,tchange,chainl){
  uA=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uA<rnorm(1,(1.793+ .045-((1.2*10^(-3)))))){#/sum(chainl))){
    Aval=Aval+(1.793*ABoval+.045*M1val-((1.2*10^(-3))*Aval))*tchange
  } else{
    Aval=Aval
  } 
  return(Aval)
}

M1valChainV2=function(M1val,Foval,tchange,chainl){
  uM1=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uM1<rnorm(1,(((.047*((2*10^(-2))/(2.58*10^(-11)))))*.9-(.015)))){#/sum(chainl))){
    M1val=M1val+((.047*((2*10^(-2))*(Foval/(Foval+(2.58*10^(-11)))))*.9)-(.015*M1val))*tchange
  } else{
    M1val=M1val
  } 
  return(M1val)
}

M2valChainV2=function(M2val,Foval,tchange,chainl){
  uM2=runif(1)
  #print(3.4*10^-4)*((Fival)/(Fival+3.36*10^-10)))*Nval))
  if(uM2<rnorm(1,((.047*((2*10^(-2))/((2.58*10^(-11)))))*.09-(.015)))){#/sum(chainl))){
    M2val=M2val+((.047*((2*10^(-2))*(Foval/(Foval+(2.58*10^(-11)))))*.09)-(.015*M2val))*tchange
  } else{
    M2val=M2val
  } 
  return(M2val)
}



ALStochModelV5=function(time){
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
  while (tnow[i-1]<time) {
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
    possiblet=c(-log(u)/chainlist)
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
    tnow[i]=tnow[i-1]+tchange
    i=i+1
  }
  return(data.frame(Nval,ABival,ABoval,Tauval,Fival,Foval,M1val,M2val,Aval,tnow))
  
}
out=ALStochModelV5(3650)
out2=ALStochModelV5(3650)
out3=ALStochModelV5(3650)
out4=ALStochModelV5(3650)
out5=ALStochModelV5(3650)




##Stoch Neuron Graph
plot(out$tnow,out$Nval,type = 'l',lty=1,ylab = 'Neuron Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Neuron Concentration Over 10 Years",col="red")
meanN=c(mean(out$Nval))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('red','brown1','brown3','brown4','firebrick','firebrick1','firebrick3','firebrick4','darkred','maroon','orangered3','orangered4','orangered2','red3','red4','indianred4','tomato','tomato2','tomato3','tomato4','brown','brown2','firebrick2','red2','orangered','indianred')
  style=c(1:25)
  lines(out$tnow,out$Nval,lty=style[i],col=cols[i],lwd=1.5)
}
legend('topright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")

##Stoch ABi Graph

plot(out$tnow,out$ABival,type = 'l',lty=1,ylab = 'Amyloid Beta Concentration Inside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Amyloid Beta Concentration Inside",col="blue")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('blue','blue4','darkblue','cyan1','cyan4','cornflowerblue','cadetblue1','dodgerblue','dodgerblue4','deepskyblue4','darkturquoise','lightblue1','lightskyblue','mediumblue','midnightblue','paleturquoise','paleturquoise3','navyblue','skyblue1','skyblue4','royalblue1','royalblue3','royalblue4','powderblue','steelblue2','turquoise')
  style=c(1:25)
  lines(out$tnow,out$ABival,lty=style[i],col=cols[i],lwd=1.5)
}
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")



##Stoch ABo graph

plot(out$tnow,out$ABoval,type = 'l',lty=1,ylab = 'Amyloid Beta Concentration Outside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Amyloid Beta Concentration Outside",col="green")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('green','aquamarine','aquamarine4','darkolivegreen','darkolivegreen3','darkolivegreen4','darkgreen','chartreuse','chartreuse4','forestgreen','darkseagreen3','darkseagreen4','green3','green4','lightgreen','lawngreen','mediumspringgreen','mediumseagreen','limegreen','palegreen','palegreen4','olivedrab','olivedrab3','seagreen','seagreen3','yellowgreen')
  style=c(1:25)
  lines(out$tnow,out$ABoval,lty=style[i],col=cols[i],lwd=1.5)
}
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")


##Stoch Tau Graph

plot(out$tnow,out$Tauval,type = 'l',lty=1,ylab = 'Tau Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Tau Concentration Over 10 Years",col="deeppink")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('deeppink','deeppink3','deeppink4','hotpink','hotpink3','lightcoral','maroon1','maroon2','maroon3','magenta','magenta1','magenta2','palevioletred1','palevioletred2','pink1','orchid1','orchid2','violet','violetred','violetred1','violetred2','violetred3','palevioletred','palevioletred3','hotpink1','hotpink2')
  style=c(1:25)
  lines(out$tnow,out$Tauval,lty=style[i],col=cols[i],lwd=1.5)
}

legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 )


##Stoch Fi Graph

plot(out$tnow,out$Fival,type = 'l',lty=1,ylab = 'NFT Concentration Inside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "NFT Concentration Inside Over 10 Years",col="gold",ylim = c(3.3*10^-10,7.5*10^-10))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('gold','goldenrod','goldenrod1','goldenrod2','goldenrod3','gold1','gold2','khaki','khaki1','lightgoldenrod','lightgoldenrod1','lightgoldenrod2','palegoldenrod','navajowhite2','yellow','yellow1','yellow2','yellow3','darkgoldenrod1','darkgoldenrod2','darkgoldenrod3','khaki2','khaki3','gold3','gold4','darkgoldenrod')
  style=c(1:25)
  lines(out$tnow,out$Fival,lty=style[i],col=cols[i],lwd=1.5)
}
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 5 ,title = "Run #")

##perfect
##Orange

plot(out$tnow,out$Foval,type = 'l',lty=1,ylab = 'NFT Concentration Outside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "NFT Concentration Outside Over 10 Years",col="darkorange")
#meanFo=c(mean(out$Foval))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  #meanFo[i]=c(mean(out$Foval))
  cols=c('darkorange','darkorange1','darkorange2','darkorange3','chocolate','chocolate1','chocolate2','orange','orange1','orange2','orange3','sienna1','sienna2','sienna3','sandybrown','salmon1','salmon2','coral','tan1','tan2','chocolate3','chocolate4','darkorange4','peru','lightsalmon2','lightsalmon3')
  style=c(1:25)
  lines(out$tnow,out$Foval,lty=style[i],col=cols[i],lwd=1.5)
}
legend('topleft',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #",bty='n')

##Stoch M1

plot(out$tnow,out$M1val,type = 'l',lty=1,ylab = 'Type 1 Microglia Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Type 1 Microglia Concentration Over 10 Years",col="darkmagenta")
#meanM1=c(mean(out$M1val))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  #meanM1[i]=c(mean(out$M1val))
  cols=c('darkmagenta','darkviolet','darkslateblue','darkorchid','darkorchid1','darkorchid2','darkorchid3','darkorchid4','mediumslateblue','mediumpurple','mediumpurple1','mediumpurple2','mediumpurple3','mediumpurple4','mediumorchid','mediumorchid1','mediumorchid2','mediumorchid3','mediumorchid4','magenta4','orchid','orchid3','orchid4','plum4','purple','purple3')
  style=c(1:25)
  lines(out$tnow,out$M1val,lty=style[i],col=cols[i],lwd=1.5)
}
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")


##Stoch M2

plot(out$tnow,out$M2val,type = 'l',lty=1,ylab = 'Type 2 Microglia Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Type 2 Microglia Concentration Over 10 Years",col="aquamarine")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('aquamarine','aquamarine1','aquamarine2','aquamarine3','aquamarine4','cyan2','darkturquoise','darkslategray2','darkslategray1','lightblue','lightblue3','lightblue4','lightcyan3','lightcyan4','paleturquoise2','paleturquoise3','paleturquoise4','lightskyblue3','lightskyblue4','lightsteelblue3','lightsteelblue4','lightsteelblue','steelblue','slategray2','slategray3','slategray4')
  style=c(1:25)
  lines(out$tnow,out$M2val,lty=style[i],col=cols[i],lwd=1.5)
}
legend('topright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")


##Perfect
##Stoch A

plot(out$tnow,out$Aval,type = 'l',lty=1,ylab = 'Astrocyte Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Astrocyte Concentration Over 10 Years",col="burlywood4")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('burlywood4','burlywood3','burlywood2','lightgoldenrod4','lightsalmon4','navajowhite2','navajowhite3','navajowhite4','khaki4','orange4','chocolate4','darkorange4','sienna','sienna4','sandybrown','salmon4','saddlebrown','lemonchiffon4','rosybrown4','rosybrown','tan1','tan2','tan3','tan4','peru','goldenrod4')
  style=c(1:25)
  lines(out$tnow,out$Aval,lty=style[i],col=cols[i],lwd=1.5)
  #lines(output$time,output$A)
}
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th'),col =cols,lty = style,ncol = 2 ,title = "Run #")
