if(!("foreign" %in% installed.packages()[,1])){install.packages("foreign")}
if(!("foreign" %in% installed.packages()[,1])){install.packages("CCA")}
library(CCA)
library(foreign)
setwd("C:/Users/MACIEK/Desktop/analizaWW/Canonical Anlysis STATA oraz R")
zadowo<-as.data.frame(read.dta("zadowolenie.dta",convert.factors = FALSE))
nazwy01<-c("rodzina" ,"przyjaciele" ,"zdrowie" ,"sukces")
nazwy02<-c("dep_wyglad", "dep_zapal", "dep_zdrowie", "dep_sen", "dep_meczenie","plec", "wiek2011" ,"kontakty", "aktywnosc", "edukacja", "zaufanie", "zaleznosc", "stanciv", "dochod")
sap1<-sapply(c(4,6,13,14),function(x) zadowo[,x]<<-as.factor(zadowo[,x]))
zadowo<-as.data.frame(zadowo)

##ROZBICIE ZMIENNYCH FACTOROWYCH

for(i in c("plec","zaleznosc","stanciv","zaufanie")){
  sap2<-sapply(2:length(levels(zadowo[,i])),function(x) zadowo[paste0(i,x)]<<-as.numeric(zadowo[,i]==levels(zadowo[,i])[x]))
}

#WYROZNIENIE MACIERZY X ORAZ Y

zadowo<-zadowo[,!names(zadowo) %in% c("plec","zaleznosc","stanciv","zaufanie")]
nazwy1<-c("rodzina" ,"przyjaciele" ,"zdrowie" ,"sukces")
nazwy2<-c("dep_wyglad", "dep_zapal", "dep_zdrowie", "dep_sen", "dep_meczenie","plec2", "wiek2011" ,"kontakty", "aktywnosc", "edukacja", "zaufanie2", "zaleznosc2", "stanciv2", "stanciv3", "stanciv4", "stanciv5", "dochod")
matY<-zadowo[,nazwy1]
matX<-zadowo[,nazwy2]

matcor(matX,matY)

#ANALIZA KANONICZNA

cc1<-cc(matX,matY)

#WYSTANDARYZOWANE WYNIKI

wyniki<-list()
wyniki[[1]] <- diag(sqrt(diag(cov(matY)))) %*% cc1$ycoef
rownames(wyniki[[1]])<-nazwy1
wyniki[[2]] <- diag(sqrt(diag(cov(matX)))) %*% cc1$xcoef
rownames(wyniki[[2]])<-nazwy2
wyniki

#FUNKCJA  PIERWSZA WILKS

WILKSL<-function(matX,matY,cc1){
  ev <- (1 - cc1$cor^2)
  n <- dim(matX)[1]
  p <- length(matX)
  q <- length(matY)
  k <- min(p, q)
  m <- n - 3/2 - (p + q)/2
  w <- rev(cumprod(rev(ev)))
  # initialize
  d1 <- d2 <- f <- vector("numeric", k)
  
  for (i in 1:k) {
    s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    si <- 1/s
    d1[i] <- p * q
    d2[i] <- m * s - p * q/2 + 1
    r <- (1 - w[i]^si)/w[i]^si
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1
    q <- q - 1
  }
  pv <- pf(f, d1, d2, lower.tail = FALSE)
  dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv)
  return(dmat)
  ## source: http://www.ats.ucla.edu/stat/r/dae/canonical.htm
}

WILKSL(matX,matY,cc1)

#FUNKCJA DRUDA WILKS

WILKS2<-function(cc1){
  dmat2<-matrix(0,nrow=ncol(matY),ncol=2)
  sapply(1:ncol(matY),function(i) dmat2[i,1]<<-prod((1-cc1$cor^2)[i:(ncol(matY))]))
  return(dmat2)
}

WILKS2(cc1)

#FUNKCJA REDUNDANCJA 

REDUNT<-function(matX,matY,cc1){
  eigenmatY<-cc1$cor
  vector1<-vector(,length(matY))
  sapply(1:ncol(matY),function(i) vector1[i]<<-(eigenmatY[i])^2)
  names1<-c("opposite variance","own variance")
  names2<-c("prop stdvar v","prop stdvar u")
  matim<<-list()
  for(i in (1:ncol(matY))){
    a<-sum((cc1$scores$corr.Y.xscores[,i])^2)/ncol(matY)
    b<-a/vector1[i]
    c<-sum((cc1$scores$corr.X.yscores[,i])^2)/ncol(matX)
    d<-c/vector1[i]
    assign(paste0("amat",i),matrix(c(c,d,a,b),byrow=TRUE,nrow=2,ncol=2,dimnames=list(names2,names1)),inherit=TRUE)
    matim[[i]]<<-get(paste0("amat",i))
  }
  return(matim)
}

REDUNT(matX,matY,cc1)

