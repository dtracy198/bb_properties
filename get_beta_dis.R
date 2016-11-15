library(deSolve)
library(rootSolve)
setwd("/Users/sherriexie/Dropbox/Levy_Lab/bb_properties")

# Final figures will show how total cost due to disclosure and prevalence 
# varies over a range of baseline prevalence values (0.02-0.10)
# A separate figure will be made for low, medium, and high rates of migration 
# x partial and maximal disclosure effects 

migration <- 1/(10*365)
# high migration rate = 1/(2*365)
# moderate = 1/(5*365)
# low = 1/(10*365)

maxc <- 1
# maximal disclosure effect: c = 1
# partial disclosure effect: c = 0.5

# =============================================================================
# Set parameter ranges (i.e. boundaries of the uniform distributions)
# =============================================================================

delta.range <- c(2/365, 26/365)  # 0.005479452 - 0.07123288
gamma.range <- c(1/365, 12/365)  # 0.002739726 - 0.03287671
i.range <- c(0, 0.5)
k.range <- c(0, 0.5)
b.range <- c(1,3)

# =============================================================================
# Functions
# =============================================================================

# ----------------------------------------------------------------------------
# SetODEs function:
# - Contains the ODEs for the model that will be fed to the ode solver in the
#   RunODEs function below.
# - Note dtrt.dt tracks the flow of Ir -> Sr & Er and Iv -> Sv & Ev, which will
#   be useful later when we will want to determine the number of treatments in 
#   our cost calculations.
# - dtov.dt tracks the flow of Xv -> Xr (where X = S, E, I), i.e. all instances
#   where a vacant property is rented out again. The amount of turnover will
#   also be used in our cost calculations.

SetODEs<-function(t,y,p){
  Sr <- y[1]
  Er <- y[2]
  Ir <- y[3]
  Sv <- y[4]
  Ev <- y[5]
  Iv <- y[6]
  Sv2 <- y[7]
  Ev2 <- y[8]
  trt <- y[9]
  tov <- y[10]
  with(as.list(p),{
    dSr.dt <- -beta*Sr*Ir/N + (1-i)*gamma*Ir + (1-c)*(1-k*b*Ir/(Sr+Er+b*Ir))*n*Sv2 +(1-k*b*Ir/(Sr+Er+b*Ir))*n*Sv - m*Sr
    dEr.dt <- -delta*Er + i*gamma*Ir + beta*Sr*Ir/N + (k*b*Ir/(Sr+Er+b*Ir))*n*Sv + (1-c)*(k*b*Ir/(Sr+Er+b*Ir))*n*Sv2 + (1-c)*n*Ev2 + n*Ev - m*Er 
    dIr.dt <- -gamma*Ir + delta*Er - b*m*Ir + (1-c)*n*Iv 
    dSv.dt <- m*Sr - n*Sv + c*(1/D)*Sv2
    dEv.dt <- m*Er - n*Ev + c*(1/D)*Ev2
    dIv.dt <- b*m*Ir - gamma*Iv - (1-c)*n*Iv
    dSv2.dt <- (1-i)*gamma*Iv - (1-c)*n*Sv2 - c*(1/D)*Sv2
    dEv2.dt <- i*gamma*Iv - (1-c)*n*Ev2 - c*(1/D)*Ev2
    dtrt.dt <- gamma*Ir + gamma*Iv
    dtov.dt <- n*Sv + (1-c)*n*Sv2 + n*Ev + (1-c)*n*Ev2 + (1-c)*n*Iv
    return(list(c(dSr.dt, dEr.dt, dIr.dt, dSv.dt, dEv.dt, dIv.dt, dSv2.dt, dEv2.dt, dtrt.dt, dtov.dt)))
  })
}

# ----------------------------------------------------------------------------
# GetParameters function: 
# - Draws parameters from uniform distributions
# - Returns a vector of randomly drawn parameter values to feed to GetBeta 
#   function below

GetParameters <- function(){
  # This function draws parameter values from uniform distributions, with
  # boundaries set above
  delta <- runif(1, delta.range[1], delta.range[2])
  gamma <- runif(1, gamma.range[1], gamma.range[2])
  i <- runif(1, i.range[1], i.range[2])
  k <- runif(1, k.range[1], k.range[2])
  b <- runif(1, b.range[1], b.range[2])
  return(c(delta,gamma,i,k,b))
}

# ----------------------------------------------------------------------------
# RunODEs function: 
# - Function uses an ode solver (function "ode" from the deSolve library) to
#   solve the system of differential equations outlined in SetODEs
# - The output of interest is the bed bug prevalence which is equal to the
#   number in Ir and Iv at equilibrium
# - This function returns "prevalence - PPP" so that the GetBeta function can
#   later solve for the beta value associated with a particular prevalence 

RunODEs <- function(beta.p, p.set, PPP){
  # Set initial conditions and time interval
  Sr0 <- 900*(1-2*PPP)
  Er0 <- 900*PPP
  Ir0 <- 900*PPP
  Sv0 <- 100*.8  # Vacancy rate (% of homes vacant) = 9.8% in Philly
  Ev0 <- 100*.05
  Iv0 <- 100*.05
  Sv20 <- 100*.05
  Ev20 <- 100*.05
  trt0 <- 0
  tov0 <- 0
  y0 <- c(Sr0, Er0, Ir0, Sv0, Ev0, Iv0, Sv20, Ev20, trt0, tov0)
  # Set t = 3650 days (10 years) to give the system enough time to reach 
  # equilibrium
  t <- seq(from=0, to=3650, by=1)

  
  # Set parameter values
  beta <- beta.p
  delta <- p.set[1] 
  gamma <- p.set[2]
  i <- p.set[3]
  k <- p.set[4]
  b <- p.set[5]
  c <- 0
  D <- 365  # length of disclosure (in days)
  m <- migration  # migration rate (rate at which properties are vacated)
  n <- 4/365  # reoccupancy rate (rate at which vacant properties get rented out 
  # again)
  N <- Sr0+Er0+Ir0+Sv0+Ev0+Iv0+Sv20+Ev20
  p <- list(beta=beta, delta=delta, gamma=gamma, b=b, c=c, D=D, i=i, k=k, m=m, n=n, N=N)
  
  out <- ode(y=y0, times=t, func=SetODEs, parms=p)
  
  # prevalence is the sum of Ir + Iv at the last time point (i.e. at equilibrium)
  # divided by N=1000
  prevalence <- (out[,4][length(t)]+out[,7][length(t)])/1000
  return(prevalence - PPP)
}

# ----------------------------------------------------------------------------
# GetBeta function: 
# - Input for GetBeta are randomly drawn parameter values (which will be 
#   obtained from the GetParameters function) and a desired prevalence
# - GetBeta uses a solver function (uniroot from the RootSolve library) to 
#   solve for the beta value that gives you the desired equilibrium prevalence 
#   given a particular set of other parameter values 

GetBeta <- function(p.set, PPP){
  fun <- function(x) RunODEs(x, p.set, PPP)
  return(uniroot(fun, c(0.00001, 0.15))$root)
}

# ----------------------------------------------------------------------------
# GetCost function 
# - This function takes a set of beta and other parameter values and uses
#   the ode solver to determine the number in Sr, Er, etc. at equilibrium
# - Then sums over all the treatment, turnover and vacancy over the last year
#   which can then be used to calculate total costs/unit/yr
# - The total cost is reported as the difference in the cost at a given 
#   disclosure effect subtracted by the "baseline" cost with all else
#   being the same but no disclosure effect

GetCost <- function(beta.p,p.set,PPP){
  c.range <- c(0,maxc) 
  cost.vec <- as.numeric()
  for (ii in 1:length(c.range)){
    # Set initial conditions and time interval
    Sr0 <- 900*(1-2*PPP)
    Er0 <- 900*PPP
    Ir0 <- 900*PPP
    Sv0 <- 100*.8  # Vacancy rate (% of homes vacant) = 9.8% in Philly
    Ev0 <- 100*.05
    Iv0 <- 100*.05
    Sv20 <- 100*.05
    Ev20 <- 100*.05
    trt0 <- 0
    tov0 <- 0
    y0 <- c(Sr0, Er0, Ir0, Sv0, Ev0, Iv0, Sv20, Ev20, trt0, tov0)
    t <- seq(from=0, to=3650, by=1)
    
    # Set parameter values
    beta <- beta.p
    delta <- p.set[1]
    gamma <- p.set[2]
    i <- p.set[3]
    k <- p.set[4]
    b <- p.set[5]
    c <- c.range[ii]
    D <- 365  # length of disclosure (in days)
    m <- migration  # migration rate (rate at which properties are vacated)
    n <- 4/365  # reoccupancy rate (rate at which vacant properties get rented out 
    # again)
    N <- Sr0+Er0+Ir0+Sv0+Ev0+Iv0+Sv20+Ev20
    p <- list(beta=beta, delta=delta, gamma=gamma, b=b, c=c, D=D, i=i, k=k, m=m, n=n, N=N)
    
    # To calculate the costs, we will sum over the number of treatments, number
    # of turnover events, and total vacancy over the last year
    out <- ode(y=y0, times=t, func=SetODEs, parms=p)
    start <- length(t)-365
    end <- length(t)
    trt <- out[,10][end]-out[,10][start] # total treatments in the last year
    tov <- out[,11][end]-out[,11][start] # total turnover in the last year
    # trt and tov units = number/year/1000 units
    # The cost per trt and per tov is $1000 so we don't need to multiply by 
    # anything to convert units to $/year/unit
    
    # below gives total vacancy (unit = units*days/year) among the 1000 units in 
    # the last year
    vac <- sum(out[,5][start:end]) + sum(out[,6][start:end]) + 
      sum(out[,7][start:end]) + sum(out[,8][start:end]) + sum(out[,9][start:end])
    
    # average monthly rent is $920 (i.e. $920/30days) and divide by 1000 to get
    # per unit (rather than per 1000 units)
    vac <- vac/(30*1000)*920 # convert units from units*days/year to $/year/unit
    
    cost.vec[ii] <- sum(trt + tov + vac)
    
    # beta was solved for the case where disclosure has no effect (i.e. c=0)
    # we are interested to see how prevalence changes with disclosure (i.e.
    # at c=1/maximal disclsoure effect and c=0.5/partial disclosure effect)
    if (c!= 0){
      prev <- (out[t[length(t)],4] + out[t[length(t)],7])/1000
    }
  }
  # total cost is the difference between total cost at c=maxc (1 or 0.5)
  # amd total cost at c=0
  tot.cost <- cost.vec[2] - cost.vec[1]
  return(c(tot.cost, prev))
}

# ----------------------------------------------------------------------------
# SampleMatrix function 
# - This function calls on all the other functions to calculate the total cost
#   at a chosen "baseline" prevalence value.
# - The steps are as follows:
#   1. GetParameters() randomly draws parameter values from their uniform 
#      distributions
#   2. GetBeta() solves for beta given the set of parameter values drawn
#      and the bed bug prevalence we are interested in 
#   3. GetCost() calculates the total cost of disclosure for the given
#      parameter/beta values and prevalence level
# - Output is a data frame of drawn parameter values, total cost, and
#   prevalence value in the face of disclosure (which will be compared
#   to the baseline prevalence)
# - Each row of the data frame represents a different draw

SampleMatrix <- function(samples, PPP){
  delt <- gam <- iv <- kv <- bv <- bet <- tot.cost <- prev <- as.numeric()
  for (ii in 1:samples){
    parameter.set <- GetParameters()
    beta.set <- GetBeta(parameter.set, PPP)
    delt[ii] <- parameter.set[1]
    gam[ii] <- parameter.set[2]
    iv[ii] <- parameter.set[3]
    kv[ii] <- parameter.set[4]
    bv[ii] <- parameter.set[5]
    bet[ii] <- beta.set
    cost.out <- GetCost(beta.set, parameter.set, PPP)
    tot.cost[ii] <- cost.out[1]
    prev[ii] <- cost.out[2]
  }
  return(data.frame(cbind(delt,gam,iv,kv,bv,bet,tot.cost,prev)))
}

# Run functions over a range of bed bug prevalence levels: 2%, 4%, 6%, 8% and 10%
data02 <- SampleMatrix(100,.02)
data04 <- SampleMatrix(100,.04)
data06 <- SampleMatrix(100,.06)
data08 <- SampleMatrix(100,.08)
data10 <- SampleMatrix(100,.10)

summary(data02$tot.cost)
summary(data04$tot.cost)
summary(data06$tot.cost)
summary(data08$tot.cost)
summary(data10$tot.cost)

# Save all data for a given migration rate and a different disclosure effect (c)
# as its own list.
# Each item of the list contains a separate data frame for each bed bug prevalence.
data<-list(NA)
data[[1]] <- data02
data[[2]] <- data04
data[[3]] <- data06
data[[4]] <- data08
data[[5]] <- data10

saveRDS(data, "cost_mig10yr_maxc.rds")


  