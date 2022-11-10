##Stoch & Det Neuron Graph

plot(out$tnow,out$Nval,type = 'l',lty=1,ylab = 'Neuron Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Neuron Concentration Over 10 Years, Stochastic Vs Deterministic",col="red")
#meanN=c(mean(out$Nval))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('red','brown1','brown3','brown4','firebrick','firebrick1','firebrick3','firebrick4','darkred','maroon','orangered3','orangered4','orangered2','red3','red4','indianred4','tomato','tomato2','tomato3','tomato4','brown','brown2','firebrick2','red2','orangered','indianred')
  style=c(1:25)
  lines(out$tnow,out$Nval,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$N,lty=26)
legend('bottomleft',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('red','brown1','brown3','brown4','firebrick','firebrick1','firebrick3','firebrick4','darkred','maroon','orangered3','orangered4','orangered2','red3','red4','indianred4','tomato','tomato2','tomato3','tomato4','brown','brown2','firebrick2','red2','orangered','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #",bty ='n')


##Stoch & Det ABi Graph
plot(out$tnow,out$ABival,type = 'l',lty=1,ylab = 'Amyloid Beta Concentration Inside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Amyloid Beta Concentration Inside Over 10 Years, Stochastic Vs Deterministic",col="blue")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('blue','blue4','darkblue','cyan1','cyan4','cornflowerblue','cadetblue1','dodgerblue','dodgerblue4','deepskyblue4','darkturquoise','lightblue1','lightskyblue','mediumblue','midnightblue','paleturquoise','paleturquoise3','navyblue','skyblue1','skyblue4','royalblue1','royalblue3','royalblue4','powderblue','steelblue2','turquoise','black')
  style=c(1:25)
  lines(out$tnow,out$ABival,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$ABi,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('blue','blue4','darkblue','cyan1','cyan4','cornflowerblue','cadetblue1','dodgerblue','dodgerblue4','deepskyblue4','darkturquoise','lightblue1','lightskyblue','mediumblue','midnightblue','paleturquoise','paleturquoise3','navyblue','skyblue1','skyblue4','royalblue1','royalblue3','royalblue4','powderblue','steelblue2','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #")


##Stoch & Det ABo Graph
plot(out$tnow,out$ABoval,type = 'l',lty=1,ylab = 'Amyloid Beta Concentration Outside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Amyloid Beta Concentration Outside Over 10 Years, Stochastic Vs Deterministic",col="green")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('green','aquamarine','aquamarine4','darkolivegreen','darkolivegreen3','darkolivegreen4','darkgreen','chartreuse','chartreuse4','forestgreen','darkseagreen3','darkseagreen4','green3','green4','lightgreen','lawngreen','mediumspringgreen','mediumseagreen','limegreen','palegreen','palegreen4','olivedrab','olivedrab3','seagreen','seagreen3','yellowgreen')
  style=c(1:25)
  lines(out$tnow,out$ABoval,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$ABo,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('green','aquamarine','aquamarine4','darkolivegreen','darkolivegreen3','darkolivegreen4','darkgreen','chartreuse','chartreuse4','forestgreen','darkseagreen3','darkseagreen4','green3','green4','lightgreen','lawngreen','mediumspringgreen','mediumseagreen','limegreen','palegreen','palegreen4','olivedrab','olivedrab3','seagreen','seagreen3','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #")

##Stoch & Det Tau Graph
plot(out$tnow,out$Tauval,type = 'l',lty=1,ylab = 'Tau Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Tau Concentration Over 10 Years, Stochastic Vs Deterministic",col="deeppink")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('deeppink','deeppink3','deeppink4','hotpink','hotpink3','lightcoral','maroon1','maroon2','maroon3','magenta','magenta1','magenta2','palevioletred1','palevioletred2','pink1','orchid1','orchid2','violet','violetred','violetred1','violetred2','violetred3','palevioletred','palevioletred3','hotpink1','hotpink2')
  style=c(1:25)
  lines(out$tnow,out$Tauval,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$Tau,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('deeppink','deeppink3','deeppink4','hotpink','hotpink3','lightcoral','maroon1','maroon2','maroon3','magenta','magenta1','magenta2','palevioletred1','palevioletred2','pink1','orchid1','orchid2','violet','violetred','violetred1','violetred2','violetred3','palevioletred','palevioletred3','hotpink1','black')
       ,lty = c(1:26),ncol = 2 )


##Stoch & Det Fi Graph
plot(out$tnow,out$Fival,type = 'l',lty=1,ylab = 'NFT Concentration Inside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "NFT Concentration Inside Over 10 Years, Stochastic Vs Deterministic",col="gold",ylim = c(3.3*10^-10,8.5*10^-10))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('gold','goldenrod','goldenrod1','goldenrod2','goldenrod3','gold1','gold2','khaki','khaki1','lightgoldenrod','lightgoldenrod1','lightgoldenrod2','palegoldenrod','navajowhite2','yellow','yellow1','yellow2','yellow3','darkgoldenrod1','darkgoldenrod2','darkgoldenrod3','khaki2','khaki3','gold3','gold4','darkgoldenrod')
  style=c(1:25)
  lines(out$tnow,out$Fival,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$Fi,lty=26)
legend('topright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('gold','goldenrod','goldenrod1','goldenrod2','goldenrod3','gold1','gold2','khaki','khaki1','lightgoldenrod','lightgoldenrod1','lightgoldenrod2','palegoldenrod','navajowhite2','yellow','yellow1','yellow2','yellow3','darkgoldenrod1','darkgoldenrod2','darkgoldenrod3','khaki2','khaki3','gold3','gold4','black')
       ,lty = c(1:26),ncol = 5 ,title = "Run #",bty = 'n')


##Stoch & Det Fo Graph
plot(out$tnow,out$Foval,type = 'l',lty=1,ylab = 'NFT Concentration Outside (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "NFT Concentration Outside Over 10 Years, Stochastic Vs Deterministic",col="darkorange")
#meanFo=c(mean(out$Foval))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  #meanFo[i]=c(mean(out$Foval))
  cols=c('darkorange','darkorange1','darkorange2','darkorange3','chocolate','chocolate1','chocolate2','orange','orange1','orange2','orange3','sienna1','sienna2','sienna3','sandybrown','salmon1','salmon2','coral','tan1','tan2','chocolate3','chocolate4','darkorange4','peru','lightsalmon2','lightsalmon3')
  style=c(1:25)
  lines(out$tnow,out$Foval,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$Fo,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('darkorange','darkorange1','darkorange2','darkorange3','chocolate','chocolate1','chocolate2','orange','orange1','orange2','orange3','sienna1','sienna2','sienna3','sandybrown','salmon1','salmon2','coral','tan1','tan2','chocolate3','chocolate4','darkorange4','peru','lightsalmon2','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #")

##Stoch & Det M1 Graph

plot(out$tnow,out$M1val,type = 'l',lty=1,ylab = 'Type 1 Microglia Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Type 1 Microglia Concentration Over 10 Years, Stochastic Vs Deterministic",col="darkmagenta",ylim = c(.02,.05))
#meanM1=c(mean(out$M1val))
for (i in 2:25) {
  out=ALStochModelV5(3650)
  #meanM1[i]=c(mean(out$M1val))
  cols=c('darkmagenta','darkviolet','darkslateblue','darkorchid','darkorchid1','darkorchid2','darkorchid3','darkorchid4','mediumslateblue','mediumpurple','mediumpurple1','mediumpurple2','mediumpurple3','mediumpurple4','mediumorchid','mediumorchid1','mediumorchid2','mediumorchid3','mediumorchid4','magenta4','orchid','orchid3','orchid4','plum4','purple','purple3')
  style=c(1:25)
  lines(out$tnow,out$M1val,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$M1,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('darkmagenta','darkviolet','darkslateblue','darkorchid','darkorchid1','darkorchid2','darkorchid3','darkorchid4','mediumslateblue','mediumpurple','mediumpurple1','mediumpurple2','mediumpurple3','mediumpurple4','mediumorchid','mediumorchid1','mediumorchid2','mediumorchid3','mediumorchid4','magenta4','orchid','orchid3','orchid4','plum4','purple','black')
,lty = c(1:26),ncol = 2 ,title = "Run #")

##Stoch & Det M2 Graph

plot(out$tnow,out$M2val,type = 'l',lty=1,ylab = 'Type 2 Microglia Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Type 2 Microglia Concentration Over 10 Years, Stochastic Vs Deterministic",col="aquamarine")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('aquamarine','aquamarine1','aquamarine2','aquamarine3','aquamarine4','cyan2','darkturquoise','darkslategray2','darkslategray1','lightblue','lightblue3','lightblue4','lightcyan3','lightcyan4','paleturquoise2','paleturquoise3','paleturquoise4','lightskyblue3','lightskyblue4','lightsteelblue3','lightsteelblue4','lightsteelblue','steelblue','slategray2','slategray3','slategray4')
  style=c(1:25)
  lines(out$tnow,out$M2val,lty=style[i],col=cols[i],lwd=1.5)
}
lines(output$time,output$M2,lty=26)
legend('topright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('aquamarine','aquamarine1','aquamarine2','aquamarine3','aquamarine4','cyan2','darkturquoise','darkslategray2','darkslategray1','lightblue','lightblue3','lightblue4','lightcyan3','lightcyan4','paleturquoise2','paleturquoise3','paleturquoise4','lightskyblue3','lightskyblue4','lightsteelblue3','lightsteelblue4','lightsteelblue','steelblue','slategray2','slategray3','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #")

##Stoch & Det Astrocyte

plot(out$tnow,out$Aval,type = 'l',lty=1,ylab = 'Astrocyte Concentration (g/ml)', xlab = 'Time (days)',lwd=1.5,main = "Astrocyte Concentration Over 10 Years, Stochastic Vs Deterministic",col="burlywood4")
for (i in 2:25) {
  out=ALStochModelV5(3650)
  cols=c('burlywood4','burlywood3','burlywood2','lightgoldenrod4','lightsalmon4','navajowhite2','navajowhite3','navajowhite4','khaki4','orange4','chocolate4','darkorange4','sienna','sienna4','sandybrown','salmon4','saddlebrown','lemonchiffon4','rosybrown4','rosybrown','tan1','tan2','tan3','tan4','peru','goldenrod4')
  style=c(1:25)
  lines(out$tnow,out$Aval,lty=style[i],col=cols[i],lwd=1.5)
  #lines(output$time,output$A)
}
lines(output$time,output$A,lty=26)
legend('bottomright',c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th','12th','13th','14th','15th','16th','17th','18th','19th','20th','21st','22nd','23rd','24th','25th','Det'),col =c('burlywood4','burlywood3','burlywood2','lightgoldenrod4','lightsalmon4','navajowhite2','navajowhite3','navajowhite4','khaki4','orange4','chocolate4','darkorange4','sienna','sienna4','sandybrown','salmon4','saddlebrown','lemonchiffon4','rosybrown4','rosybrown','tan1','tan2','tan3','tan4','peru','black')
       ,lty = c(1:26),ncol = 2 ,title = "Run #")

