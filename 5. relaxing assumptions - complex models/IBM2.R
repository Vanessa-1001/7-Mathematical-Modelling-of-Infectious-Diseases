initS <- 99
initI <- 1
beta <- 0.5
gamma<-1
tau<-1
no.of.timesteps <- 40
pop <- c(rep(0,initS),rep(1,initI))
sim.table <- as.data.frame(matrix(NA, no.of.timesteps, 7))
colnames(sim.table) <- c('WeekNo','S','E','I','R', 'lambda', 'incidence')

incidence <- 0
for(t in 1:no.of.timesteps)
{
  #Simulate what happens at each timestep
  sim.table$WeekNo[t] <- t
  lambda <- beta*sum(pop==1| pop==2)/length(pop)
  print(lambda)
  random.no <- runif(length(pop))
  # incidence <- 0
  for(i in 1:length(pop))
  {
    #Simulate what happens to each individual
    if(pop[i]==0)
    {
      if(random.no[i]<lambda)
      {
        pop[i] <- 1
        incidence <- incidence + 1
      }
    }
    else if(pop[i]==1 | pop[i]==2)
    {
      pop[i] <- pop[i] + 1
    }
    
  }
  sim.table$S[t]  <- sum(pop==0)
  sim.table$E[t] <- sum(pop==1)
  sim.table$I[t] <- sum(pop==2)
  sim.table$R[t]  <- sum(pop==3)
  sim.table$lambda[t] <- lambda
  sim.table$incidence[t] <- incidence
  
}


plot(sim.table$WeekNo,sim.table$S,type='l',ylim=c(0,length(pop)),col='black',xlab="Weeks",ylab=("# People"))
lines(sim.table$WeekNo,sim.table$E+sim.table$I,col='red')
lines(sim.table$WeekNo,sim.table$R,col='blue')
lines(sim.table$WeekNo,sim.table$incidence, col = "orange2")
legend("topright", legend = c("S","E+I","R","incidence"), lwd = 1, col=c("black","red","blue","orange2"))


# plot(sim.table$lambda, type = "l")

