##******************************************************************##
# Agent Based Simulation of the Galam Majority Model in a simple way
# obs:  agents are chosen sequentially
# Galam, Serge. "Sociophysics: a review of Galam models." 
# International Journal of Modern Physics C 19.03 (2008): 409-440.
##******************************************************************##

# Main rule: Iif there is a local majority (2 × 1) in favor of
# one of the two possible opinions, the individual with minority opinion
# will follow the local majority.

nao     = 550  # number of initial A
nbo     = 450  # number of initial B

tmax    = 15   # number of time steps 
nsims   = 5    # number of independent simulations

N       = nao+nbo

ids     = 1:N

out.A.num = matrix(0,nrow=tmax, ncol=nsims)
out.B.num = matrix(0,nrow=tmax, ncol=nsims)

# Loop over diferent configurations
for (simu in 1:nsims) {
  
  # initialization
  A.num  = c(nao)
  B.num  = c(nbo)
    
  opi  = c( rep(1, A.num), rep(-1, B.num) )
  opi2 = opi
  
  # Dynamics 
  for (t in 2:tmax) {
    
    # Visit each agent
    for (node in 1:N ) 
    { 
            
      if( opi[ node ] == 1 ) # if node A has opinion +1
      {
        B.den = B.num[t-1]/N
        Cab    = B.den*B.den 
        if( runif(1) <= Cab ) 
        {
          opi2[ node ] = -1
        }
      }
      else
        if( opi[ node ] == -1 ) 
        {
          A.den = A.num[t-1]/N
          Cba = A.den*A.den
          if( runif(1) <= Cba ) 
          {
            opi2[ node ] = +1  
          }      
        }
    } # end of one MC step
    
    opi = opi2              # simultaneous-parallel updating
    
    # Opi of all nodes
    A.num[t] = sum(opi  == +1)
    B.num[t] = sum(opi  == -1)    
  }
  
  # save the results for each configuration(simulation, simu)
  out.A.num[, simu] = A.num
  out.B.num[, simu] = B.num
  cat(simu, " ")
}


###***************** Plot  I ********************#####
A.mean = apply(out.A.num, 1, mean)/N
plot(1:tmax, A.mean,xlim=c(1,tmax), ylim=c(0,1),pch=16,
     ylab=expression(p[A]), xlab="t", type="p", col="blue" )

######**********   interactive map *******************##
p = c( nao/N )

for( t in 2:tmax){
  p[t] = p[t-1]*p[t-1]*p[t-1] + 3*p[t-1]*p[t-1]*(1-p[t-1])
}


###***************** Plot II ********************#####
mtext( sprintf("%s%0.2f", " Initial proportion of agents with opinion A=",p[1] ), 
       side=3, cex=0.8, line=0)

lines(p, col="red", type="l",pch=16)

legend( "bottomleft", 0.95 , cex=0.7, 
        lty = c(1,0), 
        pch = c(NA, 16), 
        col = c( "red","blue"), 
        legend=c( expression(p[t+1]==p[t]^3+3*p[t]*(1-p[t])) ,expression( "Monte Carlo Simulation" ) ) )

