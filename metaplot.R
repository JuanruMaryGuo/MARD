library(latex2exp)
library(xtable)
library(readxl)
library(plotrix)
# Create data:
rm(final1)
rm(final2)
rm(final3)
sd_na<-function(x){sd(x,na.rm = TRUE)}
for (i in c("10.310","10.330","10.610","10.630")){
  t_normal1 = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[2:11]
  t_normal1lower = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[12:21]
  t_normal1upper = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[22:31]
  #print(summary(t_beta1))
  t_normal2 = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "result",col_names = TRUE)
  t_normal2upper = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultupper",col_names = TRUE)
  t_normal2lower = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultlower",col_names = TRUE)
  #print(summary(t_beta2))
  t_normal = cbind(t_normal1,t_normal2)
  t_normallower = cbind(t_normal1lower,t_normal2lower)
  t_normalupper = cbind(t_normal1upper,t_normal2upper)
  t_nomal = round(cbind(colSums(is.na(t_normal)),(colMeans(t_normal,na.rm = TRUE)+0.2),std.error(t_normal,na.rm = TRUE)**2, colMeans((t_normal+0.2)**2,na.rm = TRUE),
                        colMeans((t_normallower<=-0.2)*(t_normalupper>=-0.2),na.rm = TRUE)),10)
  
  
  final1 = tryCatch(cbind(final1,t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100), error=function(e){cbind(t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100)})
}




##############################################################################################
# -0.02

sd_na<-function(x){sd(x,na.rm = TRUE)}
# results for -0.02
for (i in c("20.310","20.330","20.610","20.630")){
  t_normal1 = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[2:29]
  t_normal1lower = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[30:57]
  t_normal1upper = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[58:85]
  
  
  t_normal2 = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "result",col_names = TRUE)
  t_normal2upper = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultupper",col_names = TRUE)
  t_normal2lower = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultlower",col_names = TRUE)
  t_normal = cbind(t_normal1,t_normal2)
  t_normallower = cbind(t_normal1lower,t_normal2lower)
  t_normalupper = cbind(t_normal1upper,t_normal2upper)
  t_nomal = round(cbind(colSums(is.na(t_normal)),(colMeans(t_normal,na.rm = TRUE)+0.02),std.error(t_normal,na.rm = TRUE)**2, colMeans((t_normal+0.02)**2,na.rm = TRUE),
                        colMeans((t_normallower<=-0.02)*(t_normalupper>=-0.02),na.rm = TRUE)),10)
  

  final2 = tryCatch(cbind(final2,t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100), error=function(e){cbind(t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100)})
}


# -0.02

sd_na<-function(x){sd(x,na.rm = TRUE)}
# results for -0.02
for (i in c("30.310","30.330","30.610","30.630")){
  t_normal1 = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[2:29]
  t_normal1lower = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[30:57]
  t_normal1upper = read.csv(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\N",i,".csv",sep = ""))[58:85]
  
  
  t_normal2 = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "result",col_names = TRUE)
  t_normal2upper = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultupper",col_names = TRUE)
  t_normal2lower = read_excel(paste("C:\\Users\\mary\\Desktop\\metacodes\\SIMULATE_NORMAL\\",i,".xlsx",sep = ""),sheet = "resultlower",col_names = TRUE)
  t_normal = cbind(t_normal1,t_normal2)
  t_normallower = cbind(t_normal1lower,t_normal2lower)
  t_normalupper = cbind(t_normal1upper,t_normal2upper)
  t_nomal = round(cbind(colSums(is.na(t_normal)),(colMeans(t_normal,na.rm = TRUE)+0.02),std.error(t_normal,na.rm = TRUE)**2, colMeans((t_normal+0.02)**2,na.rm = TRUE),
                        colMeans((t_normallower<=-0.02)*(t_normalupper>=-0.02),na.rm = TRUE)),10)
  
  
  final3 = tryCatch(cbind(final3,t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100), error=function(e){cbind(t_nomal[,2]*1000,t_nomal[,4]*1000,t_nomal[,5]*100)})
}



width = 4600
height = 2200
h = 180

pdf("normal.pdf", width = width/h, height = height/h)
par(mar=c(4.1, 5.1, 4.1, 2.1)*1.2)
layout(matrix(c(1,2,2,3,3,4,5,5,6,6,7,8,8,9,9), nrow = 3, ncol = 5, byrow = TRUE))


a11=c("I1","I2","II1","II2","II3","II4","II5","II6","II7","II8")
a12=c("I1a","I1b","I1c","I2",
     "II1a","II2a","II3a","II4a","II5a","II6a","II7a","II8a",
     "II1b","II2b","II3b","II4b","II5b","II6b","II7b","II8b",
     "II1c","II2c","II3c","II4c","II5c","II6c","II7c","II8c")
a2=c("III1","III2","IV")

###############################################################################################

b1=final1[,1]
b2=final1[,4]
b3=final1[,7]
b4=final1[,10]

# Make a basic graph
plot( 1:13,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("Bias($10^{-3}$)") , col=rgb(0.2,0.4,0.1,0.7), 
      lwd=3 , pch=17 , ylim=c(-5,5), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:10, labels=a11)
axis(1, at=11:13, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 1", "Setting 2","Setting 3","Setting 4"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

mtext("(A)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)

b1=final2[,1]
b2=final2[,4]
b3=final2[,7]
b4=final2[,10]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("Bias($10^{-3}$)") , col=rgb(0.2,0.4,0.1,0.7), 
      lwd=3 , pch=17 , ylim=c(-10,10), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 5", "Setting 6","Setting 7","Setting 8"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))

mtext("(B)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)



b1=final3[,1]
b2=final3[,4]
b3=final3[,7]
b4=final3[,10]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("Bias($10^{-3}$)") , col=rgb(0.2,0.4,0.1,0.7), 
      lwd=3 , pch=17 , ylim=c(-15,15), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 9", "Setting 10","Setting 11","Setting 12"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))

mtext("(C)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)
#####################################################################################################################

b1=final1[,2]
b2=final1[,5]
b3=final1[,8]
b4=final1[,11]

# Make a basic graph
plot( 1:13,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("MSE($10^{-3}$)") , 
      col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17 , ylim=c(-2,2), cex.lab=1.8  )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:10, labels=a11)
axis(1, at=11:13, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 1", "Setting 2","Setting 3","Setting 4"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

mtext("(D)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)

b1=final2[,2]
b2=final2[,5]
b3=final2[,8]
b4=final2[,11]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("MSE($10^{-3}$)") , 
      col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17 , ylim=c(-0.2,0.2), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 5", "Setting 6","Setting 7","Setting 8"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))
mtext("(E)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)


b1=final3[,2]
b2=final3[,5]
b3=final3[,8]
b4=final3[,11]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab=TeX("MSE($10^{-3}$)"), 
      col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=17 , ylim=c(-10,10), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 9", "Setting 10","Setting 11","Setting 12"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))
mtext("(F)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)

########################################################################
b1=final1[,3]
b2=final1[,6]
b3=final1[,9]
b4=final1[,12]

# Make a basic graph
plot( 1:13,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab="CP(%)" , col=rgb(0.2,0.4,0.1,0.7) , 
      lwd=3 , pch=17 , ylim=c(0,100), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:10, labels=a11)
axis(1, at=11:13, labels=a2,font = 2)
# Add a legend
legend("bottomleft", 
       legend = c("Setting 1", "Setting 2","Setting 3","Setting 4"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

mtext("(G)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)

b1=final2[,3]
b2=final2[,6]
b3=final2[,9]
b4=final2[,12]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab="CP(%)" , col=rgb(0.2,0.4,0.1,0.7) , 
      lwd=3 , pch=17 , ylim=c(0,100), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)

# Add a legend
legend("bottomleft", 
       legend = c("Setting 5", "Setting 6","Setting 7","Setting 8"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))

mtext("(H)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)


b1=final3[,3]
b2=final3[,6]
b3=final3[,9]
b4=final3[,12]

# Make a basic graph
plot( 1:31,b1 ,xaxt = "n", type="b" , bty="l" , xlab="Different Models" , ylab="CP(%)" , col=rgb(0.2,0.4,0.1,0.7) , 
      lwd=3 , pch=17 , ylim=c(0,100), cex.lab=1.8 )
lines(b2 , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=16 , type="b" )
lines(b3 , col=rgb(0.1,0.4,0.8,0.7) , lwd=3 , pch=18 , type="b" )
lines(b4 , col=rgb(0.8,0.1,0.1,0.7) , lwd=3 , pch=15 , type="b" )
axis(1, at=1:28, labels=a12)
axis(1, at=29:31, labels=a2,font = 2)

# Add a legend
legend("bottomleft", 
       legend = c("Setting 9", "Setting 10","Setting 11","Setting 12"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7),
               rgb(0.1,0.4,0.8,0.7),
               rgb(0.8,0.1,0.1,0.7)), 
       pch = c(17,16,18,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.2, 0.1))

mtext("(I)",side=3,line=1, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.6)



dev.off()




