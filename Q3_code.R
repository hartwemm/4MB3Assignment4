library(ggplot2)
library(viridis)
# Gillespie Algorithm
SI.Gillespie <- function (beta, N , I0, tmax) {
  t=0
  I=I0
  iter=tmax
  time <- c(t)
  prev <- c(I)
  while (t <= tmax && I <= N) {
    a0=beta*(N-I)*I
    t <- t + 1/a0*log(1/(1-runif(1)))
    I <- I+1
    time <- c(time,t)
    prev <- c(prev,I)
  }
  df <- data.frame(Time=time,Prevalence=prev)
  return(df)
}
df <- SI.Gillespie(1,1000,1,tmax=1)
ggplot(df,aes(x=Time,y=Prevalence))+
  geom_point(col="dodgerblue")+
  ggtitle("Gillespie Approximation to SI Model")


## Part b  
SI.exact <- function(beta,N,I0,tmax){
  t=seq(0,tmax,length.out=201)
  I=I0*N/(I0+(N-I0)*exp(-beta*N*t))
  df <- data.frame(Time=t,Prevalence=I)
  return(df)
}

N=c(32,100,1000,10000)
B=1
I0=1
tmax=1
set.seed(2018)
for (j in 1:4) {
  df <- SI.Gillespie(B,N[j],I0,tmax)
  Tmax <- df[(nrow(df)-1),"Time"]
  p <- ggplot(df,aes(x=Time,y=Prevalence)) + 
    geom_line(color=viridis(30)[1],size=0.5) +
    ggtitle(paste("Sol'ns of SI Model with Pop'n",N[j]))+
    theme(plot.title = element_text(hjust = 0.5,size=14))
  for (i in c(2:30)) {
    df2 <- SI.Gillespie(B,N[j],I0,tmax)
    Tmax <- max(Tmax,df2[(nrow(df2)-1),"Time"])
    p <- p + geom_line(data=df2,color=viridis(30)[i],size=1)
  }
  df <- SI.exact(B,N[j],I0,Tmax)
  p <- p + geom_line(data=df,col="red",size=2) +
    coord_cartesian(xlim=c(0,Tmax))
  if (j==1) {
    p1 <- p
  } else if (j==2) {
    p2 <- p
  } else if (j==3){
    p3 <- p
  } else {
    p4 <- p
  }
}
grid.arrange(p1,p2,p3,p4,nrow=2)

  