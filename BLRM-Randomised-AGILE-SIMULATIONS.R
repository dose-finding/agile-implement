# Defining the MCMC model
library("rjags")
model1.string <-"
model {
for (i in 1:m){
logit(p[i]) <- alpha0[1] + alpha1[1] * sdose[i]
s[i] ~ dbin(p[i], n[i])		
}

theta[1:2] ~ dmnorm(priorMean[1:2],
priorPrec[1:2,
1:2])
## extract actual coefficients
alpha0<- theta[1]
alpha1 <- exp(theta[2])
}
"
model1.spec<-textConnection(model1.string)

# Defining function to find standarditised doses for the given skeleton ptox and parameters alpha
find.x <- function(ptox, alpha ) {
  alpha<-matrix(alpha,ncol=2)
  x <- (qlogis(ptox)-(alpha[,1]))/alpha[,2]
  return( x )
}

# Defining doses (1 stands for SoC and 2:5 stand for experimental doses)
doses<-c(1,2,3,4,5)
# SoC
D<-doses[1]


# Defining Scenarios in Table 2
true<-c(0.10,0.30,0.45,0.60,0.70) # Sc 1
# true<-c(0.10,0.15,0.30,0.45,0.60) # Sc 2
# true<-c(0.10,0.12,0.15,0.30,0.45) # Sc 3
# true<-c(0.10,0.11,0.12,0.15,0.30) # Sc 4

# Number of MCMC Samples used to approximate the posterior distribution
iter<-10000
# Number of Simulations used to produce OC
nsims<-2000

# Cohort size for experimental group
cohort<-4
# Cohort size for SoC/Control group
cohort.control<-2
# Total number of patients
N<-30
# Starting dose
firstdose<-2

# Target increase in the toxicity (over the control)
target.increase<-0.20
# Half-width of the tolerance interval around target.increase
delta<-0.05

# Prior Probability of AE at the SoC
p0.control<-0.10
# Overdosing Threshold
overdose<-0.25

# Calibrated prior parameters
var1<-1.10
var2<-0.30       
slope<-(-0.05) 
spacing<-0.075


#Defining Skeleton and standartised dose levels corresponding to this skeleton
p.tox0<-c(p0.control,p0.control + spacing* seq(1,length(doses)-1)) # finding the skeleton
priorMean<-c(log(p0.control/(1-p0.control)),slope)  
priorVar<-matrix(c(var1,0.0,0.0,var2),2,2)  
priorPrec<-solve(priorVar)
alpha.prior.plug<-c(priorMean[1],exp(priorMean[2]+diag(priorVar)[2]/2))
sdose<-find.x(p.tox0,alpha=alpha.prior.plug) # standartised dose levels

# Defining matrices to store the results
ss<-mat.or.vec(nsims,1)
selection<-mat.or.vec(nsims,length(doses))
p<-mat.or.vec(iter,length(doses))
  
# Running Simulations        
          for (z in 1:nsims){
            nextdose<-firstdose
            counter<-0
            stop<-0
            n<-rep(0,length(doses))
            s<-rep(0,length(doses))
            
            
            while(sum(n)<N){
              
              n[1]<-n[1]+cohort.control
              n[nextdose]<-n[nextdose]+cohort
              
              #Assigning the patients and evaluating DLTs
              s[1]<-s[1]+sum(rbinom(cohort.control,1,true[1]))
              s[nextdose]<-s[nextdose]+sum(rbinom(cohort,1,true[nextdose]))
              
              #Fitting the Bayesian model
              model1.spec<-textConnection(model1.string)
              mydata <- list(n=n,s=s,m=length(doses),sdose=sdose,priorMean=priorMean,priorPrec=priorPrec)
              jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
              update(jags, iter,progress.bar="none")
              tt<-jags.samples(jags,c('alpha0','alpha1'),iter,progress.bar="none")
              
              # Extracting vectors of posterior samples of the model parameters
              a0<-tt$alpha0[1,,]
              a1<-tt$alpha1[1,,]
              
              #Fitting the model with these parameters
              for (j in 1:length(doses)){
                logit <- a0 + a1 * sdose[j]
                p[,j]<-exp(logit)/(1+exp(logit))
              }
              
              # Finding the probability of being in the target interval and overdosing probability
              prob.next<-mat.or.vec(length(doses),1)
              for (j in 2:length(doses)){
                y<-p[,j]-p[,1]
                prob.next[j]<-mean(y <=(target.increase+delta) & (y>=target.increase-delta))
                if(mean(y>=(target.increase+2*delta))>overdose){
                  prob.next[j]<-0
                }
              }
              
              
              # If all unsafe - stop the trial, otherwise assign to the max Prob of Target dose (subject to no skipping constraint)
              if(all(prob.next==0)){
                stop<-1
                break()
              }else{
                nextdose<-min(nextdose+1,which.max(prob.next))
              }
              
            }
            
            # Storing results of the simulation
            if(stop==0){
              selection[z,nextdose]<-1
              ss[z]<-sum(n)
            }else{
              counter<-counter+1
              ss[z]<-sum(n)
            }
            cat(z,"\n")
          }

# Proportion of Each Dose Selection
          colMeans(selection)
# Mean Sample size
          mean(ss)

