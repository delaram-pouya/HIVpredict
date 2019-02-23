library("ctbn")
rm(list=ls(all=TRUE))

setwd("C:\\Users\\Delaram\\Desktop\\HIV_ctbn2")
#setwd('/mnt/c/Users/Delaram/Desktop/HIV_ctbn2')
# variables
vars  <- read.table( "HIV_Vars.txt", sep="\t",   header=TRUE)
# static part
bn.str <- read.table("HIV_BN_Struct.txt" , sep=",", header=TRUE)
bn.cpts  <- ReadProbStruct("HIV_CPTs3.txt")
# dynamic part
dyn.str  <- read.table( "HIV_Dyn_Struct.txt", sep=",", header=TRUE)
dyn.cims <- ReadProbStruct("HIV_CIMs3.txt")

head(bn.str)
head(bn.cpts)
head(dyn.str)
head(dyn.cims)

# a6. create the CTBN object
xpCtbn <- NewCtbn(vars,bn.str,bn.cpts,dyn.str,dyn.cims)
#xpCtbn <- NewCtbn(vars,bn.str,NULL,dyn.str,NULL)

# check inputs
GetCtbnVars  (xpCtbn)
GetBnStruct  (xpCtbn)
GetBnCPTs    (xpCtbn)
GetDynStruct (xpCtbn)
GetDynIntMats(xpCtbn)
# save ctbn
# SaveRCtbn (xpCtbn,"C:\\Users\\Delaram\\Desktop\\HIV_CTBN.rctbn")
### sample from CTBN ###
par(mfrow = c(3,4))
# sample full trajectories
samplesFull <- SampleFullTrjs(xpCtbn, num=100)
sampleFull  <- samplesFull[[1]]
head(sampleFull)
lapply(seq(1,length(sampleFull),by = 3), function(i) plot(sampleFull[1:nrow(sampleFull)-1,c(1,i)],type="s",main="Full Trajectory"))


# sample partial trajectories
samplesPart <- SamplePartialTrjs(xpCtbn, num=10, rem.frac=0.1)
samplePart  <- samplesPart[[1]]
head(samplePart)
lapply(seq(1,length(samplePart),by = 3), function(i) plot(sampleFull[1:nrow(samplePart)-1,c(1,i)],type="s",main="Partial Trajectory"))

#######################################################
### parameters learning from fully observed data ###
samples  <- SampleFullTrjs(xpCtbn,num=100)
#new ctbn and parameters learning
xpCtbn2  <- CloneCtbn(xpCtbn)
LearnCtbnParams(xpCtbn2,samples)

# get parameters
cimsOut  <- GetDynIntMats(xpCtbn)
cimsOut2 <- GetDynIntMats(xpCtbn2)
cptsOut  <- GetBnCPTs(xpCtbn)
cptsOut2 <- GetBnCPTs(xpCtbn2)

# delete new ctbn
xpCtbn2  <- DeleteCtbn(xpCtbn2)
garbage  <- gc()

# compare true vs new parameters
print("bn: d.NRTI.1$`VL=1`")
print( cptsOut$d.NRTI.1$`VL=1`)
print( cptsOut2$d.NRTI.1$`VL=1`)
### calculated sum squared error
print(sqrt(sum((cptsOut$d.NRTI.1$`VL=1`-cptsOut2$d.NRTI.1$`VL=1`)^2)))

print("dynamics: d.NRTI.1$`VL=1,g.NRTI.1=0=1`")
print( cimsOut$d.NRTI.1$`VL=1,g.NRTI.1=1`)
print(cimsOut2$d.NRTI.1$`VL=1,g.NRTI.1=1`)
### calculated sum squared error
print(sqrt(sum((cimsOut$d.NRTI.1$`VL=1,g.NRTI.1=2`-cimsOut2$d.NRTI.1$`VL=1,g.NRTI.1=2`)^2)))

################################################
### parameters learning from partially observed data ###

# choose inference engine (default: exact inference)
inf <- "importance";
# sample partial trajectories
samples <- SamplePartialTrjs(xpCtbn,num=10,rem.frac=0.01)

# new ctbn and parameters learning
xpCtbn3 <- CloneCtbn(xpCtbn)
LearnCtbnParams(xpCtbn3,samples,inf.type=inf,num.samples=100)


# ptm <- proc.time()
# available inference engines and relative parameters
# switch(inf, 
#        importance={ LearnCtbnParams(xpCtbn3,samples,inf.type=inf,num.samples=100)},
#        gibbs=     { LearnCtbnParams(xpCtbn3,samples,inf.type=inf,num.samples=100,burn.iters=100)},
#        gibbsaux=  { LearnCtbnParams(xpCtbn3,samples,inf.type=inf,num.samples=100,burn.iters=100)},
#        meanfield= { LearnCtbnParams(xpCtbn3,samples,inf.type=inf,eps=0.001)},
#        { LearnCtbnParams(xpCtbn3,samples) })
# proc.time()-ptm

#  get parameters
cimsOut  <- GetDynIntMats(xpCtbn)
cimsOut3 <- GetDynIntMats(xpCtbn3)
cptsOut  <- GetBnCPTs(xpCtbn)
cptsOut3 <- GetBnCPTs(xpCtbn3)

print("delete new ctbn")
xpCtbn3  <- DeleteCtbn(xpCtbn3)
garbage  <- gc()

cptsOut$drug1$`VL=1`
print("compare true vs new parameters")
print("bn: drug1$`VL=1`")
print( cptsOut$drug1$`VL=1`)
print(cptsOut3$drug1$`VL=1`)
print(sqrt(sum((cptsOut$drug1$`VL=1`-cptsOut3$drug1$`VL=1`)^2)))

print("dynamics: drug1$`VL=1,gss1=2`")
print(cimsOut$drug1$`VL=1,gss1=2`)
print(cimsOut3$drug1$`VL=1,gss1=2`)
print(sqrt(sum((cimsOut$drug1$`VL=1,gss1=2`-cimsOut3$drug1$`VL=1,gss1=2`)^2)))

# delete ctbn
#xpCtbn  <- DeleteCtbn(xpCtbn)
#garbage <- gc()


################################################

### c. query CTBN ###

samplesFull  <- SampleFullTrjs(xpCtbn,num=100)
sampleFull = samplesFull[[1]]
# simulate partial data > prepare partial trajectory to query
trj  <- sampleFull
var  <- c("CD4")
val  <- 2
# put random noise
rem  <- round(runif(10, 1,length(trj[,var])-2))
trj[rem,var]=-1

# c.2: expected time and number of transitions
full <- QueryCtbnStats(xpCtbn,sampleFull,var)
part <- QueryCtbnStats(xpCtbn,trj,       var)
print("Time that CD4 stays in 0")
print(full[[var]][1,1])
print("Expected time that CD4 stays in 0")
print(part[[var]][1,1])
print("Transitions of CD4 from 0 to 2")
print(full[[var]][1,val+1])
print("Expected transitions of CD4 from 0 to 2")
print(part[[var]][1,val+1])

# c.3: expected time that the process stays in some states
insts <- data.frame("CD4"=c(0, 1, 2))
print(insts)
times <- QueryCtbnTime(xpCtbn,trj,insts)
print("Expected time that CD4=2")
print(times[val+1,1])

# c.4: filtering
time  <- seq(0,10,0.01)
vals  <- rep(val,length(time))
insts <- data.frame(Time=time,CD4=val)
x     <- QueryCtbnFilter(xpCtbn,trj,insts)
tlt   <- paste (paste(var,val,sep="="),sep=", ")
plot(x,main=tlt,xlab="Time",ylab="Probability",type="l")

# c.5: smoothing
x     <- QueryCtbnSmooth(xpCtbn,trj,insts)
plot(x,main=tlt,xlab="Time",ylab="Probability",type="l")


# delete ctbn
xpCtbn  <- DeleteCtbn(xpCtbn)
garbage <- gc()








Names = c("VL", "CD4", 
          "d.NRTI.1", "d.NRTI.2", "d.NRTI.3", "d.NRTI.4",
          "d.NNRTI.1", "d.NNRTI.2", "d.NNRTI.3", "d.NNRTI.4",
          "d.PI.1", "d.PI.2", "d.PI.3", "d.PI.4",
          "d.IntI.1", "d.IntI.2", "d.IntI.3", "d.IntI.4",
          
          "g.NRTI.1", "g.NRTI.2", "g.NRTI.3", "g.NRTI.4",
          "g.NNRTI.1", "g.NNRTI.2", "g.NNRTI.3", "g.NNRTI.4",
          "g.PI.1", "g.PI.2", "g.PI.3", "g.PI.4",
          "g.IntI.1", "g.IntI.2", "g.IntI.3", "g.IntI.4")

