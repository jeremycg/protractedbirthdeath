# this is the meat of the functions, actually runs the sims. takes a
# list of 5 parameters and a time
pbdsim2 <- function(pars, totaltime = 15) {
  overload=0
  good <- list()  #list of good species
  incipient <- list()  #incipients
  good[[1]] <- c(1, totaltime, totaltime, -1, 0)  #sim starts with one good species
  incipient[[1]] <- c(2, totaltime, -1, -1, 1)  #and one incipient species
  # an individual looks like: taxaid, birth time, time at 'good'
  # transition,time at death, parent if one of these hasn't happened yet,
  # -1 is used
  taxaid <- 3  #the taxa label of the next species
  deadgood <- list()  #all dead good individuals
  deadincipient <- list()  #all dead incipient
  t <- 0  #we start at t=0
  while (t <= totaltime) {
    # until t is bigger than time (should never hit this, break out below)
    numgood <- length(good)  #we want to know how many good species
    numincipient <- length(incipient)  #and incipient species
    if(numgood+numincipient>100000){
      overload<-1
      break
    }
    if (numgood == 0 && numincipient == 0) {
      break
    }  #if everything is extinct, we break
    probs <- c(pars[1] * numgood, pars[2] * numincipient, pars[3] *
                 numincipient, pars[4] * numgood, pars[5] * numincipient)
    # probs are the probability of each event happening so its prob good
    # speciation times number of good species etc
    denom <- sum(probs)  #we sum them for total rate per time point
    probs <- probs/denom  #and then correct to 1 for the weighting
    t <- t - log(runif(1))/denom  #this increases time by doob gillespie
    # it's equivalent to t=t+rexp(1,denom)
    if (t >= totaltime) {
      break
    }  #if we just went over time, break
    event <- sample(1:5, 1, prob = probs)  #5 things that can happen, weighted by probs
    if (event == 1) {
      # new incipient from good
      take <- sample(1:numgood, 1)  #choose a parent
      incipient[[numincipient + 1]] <- c(taxaid, totaltime - t, -1,
                                         -1, good[[take]][1])  #add incipient with time parent etc
      taxaid <- taxaid + 1  #increase taxaid for parental tracking
    } else if (event == 2) {
      # new good from incipient
      take <- sample(1:numincipient, 1)  #choose one
      good[[numgood + 1]] <- incipient[[take]]  #add it to good
      good[[numgood + 1]][3] <- totaltime - t  #change its time to be good rate to t
      incipient[[take]] <- NULL  #delete it from incipient
    } else if (event == 3) {
      # new from incipient
      take <- sample(1:numincipient, 1)  #choose one
      incipient[[numincipient + 1]] <- c(taxaid, totaltime - t, -1,
                                         -1, incipient[[take]][1])  #add it in
      taxaid <- taxaid + 1
    } else if (event == 4) {
      # dead good
      totaldeadgood <- length(deadgood)  #gets length
      take <- sample(1:numgood, 1)  #chooses one to die
      deadgood[[totaldeadgood + 1]] <- good[[take]]  #adds it
      deadgood[[totaldeadgood + 1]][4] <- totaltime - t  #adds death time
      good[[take]] <- NULL  #removes it
    } else if (event == 5) {
      # dead incipient
      totaldeadincipient <- length(deadincipient)  #gets length
      take <- sample(1:numincipient, 1)  #takes it
      deadincipient[[totaldeadincipient + 1]] <- incipient[[take]]  #adds it
      deadincipient[[totaldeadincipient + 1]][4] <- totaltime - t  #adds death time
      incipient[[take]] <- NULL  #removes it
    }
  }
  output<-c(good, incipient, deadgood, deadincipient)
  if(overload==1){
    attributes(output)<-list(overloaded="1")
  } else {
    attributes(output)<-list(overloaded="0")
  }
  return(output)  #returns a list with all taxa as above
}

# a function to repeat the above a ton of times, and store the output
# in a data frame
repsim2 <- function(pars, n, time = 15) {
    x <- c()  #holding for output
    i <- 1  #repeat number
    while (i <= n) {
        y <- c()  #temp holding
        for (j in pbdsim2(pars, time)) {
            y <- rbind(y, j)
        }  #run the sim
        y <- cbind(y, i)  #add on the run number column
        x <- rbind(x, y)  #add it all to output
        i <- i + 1  #increase repeats
    }
    x <- as.data.frame(x, row.names = F)  #converts to dataframe
    names(x) <- c("taxalabel", "timeatbirth", "speciationcomplete", "timeofdeath",
        "parent", "run")  #names outputs
    return(x)  #output
}

# a function to make a newick tree text from a given output is not of
# class phy - need to do read.tree(text=treemaker(x)) for that
treemaker <- function(x) {
    for (i in 1:length(x[, 1])) {
        # for each taxa
        if (x[i, ]$timeofdeath == -1)
            {
                x[i, ]$timeofdeath <- 0
            }  #changes -1s to 0s
    }
    x$label <- x$taxalabel  #makes a new column to hold the tree strings
    while (length(x[, 1]) >= 2) {
        # while there are more than two taxa
        birthtime <- min(x$timeatbirth)  #find the youngest taxa
        if (length(x[, 1]) == 2) {
            # if there are two
            child1 <- which(x$parent == 1)  #it breaks as we started with two taxa at t=end
            parent1 <- which(x$parent == 0)  #so manually set the parent/child
            # as the rest is done by stepwise functions we should never have more
            # than one child at the same timepoint unless we hit machine precision,
            # which is very very unlikely
        } else {
            child1 <- which(x$timeatbirth == birthtime)  #the child is the one with the birthtime=now
            parent1 <- which(x$taxalabel == x[child1, ]$parent)  #the parent is the parent of that one
        }
        x[parent1, ]$label <- paste("(", x[parent1, ]$label, ":", birthtime -
            x[parent1, ]$timeofdeath, ",", x[child1, ]$label, ":", birthtime -
            x[child1, ]$timeofdeath, ")", sep = "")
        # this changes the label of the parent, to newick format should look
        # like (parent:distance,child:distance) but, newick can nest, so the
        # parent or child part can be a previous label so we can recursively
        # make a tree
        # ((parent1:distance,child1:distance):distance,(parent:distance,child:distance))
        # etc etc.
        x[parent1, ]$timeofdeath <- birthtime  #collapses the node, so the parent node now dies at the split
        x <- x[-child1, ]  #removes the child
    }
    return(paste(x$label, ";", sep = ""))  #returns the string, with the required tailing ;
}

# so far: how close is each species sister taxa? as a plyrable function
# looks at all extant taxa finds closest relative branch time - single,
# not both run on a single run from repsim2
sisterlengths <- function(working) {
    phy <- treemaker(working)  #make a tree
    phy1 <- read.tree(text = phy)  #make it of class phy
    if (sum(which(working$timeofdeath == -1)) <= 1) {
        # if we had no taxa at the end
        return("no extant species")
    }
    phy2 <- drop.tip(phy1, as.character(working[which(working$timeofdeath !=
        -1), ]$taxalabel))  #this drops all dead taxa
    return(phy2$edge.length[phy2$edge[, 2] <= Ntip(phy2)])  #returns the distance to node for each taxa
    # this works as R stores edges in order - we only want the edges which
    # lead to species, not nodes so we want the ones that are less than the
    # number of tips
}
# output is a list - it depends on number of species, so can't be a
# dataframe

# next, a function for 'time slicing' we want to go at say 10 million
# y.a then today and see what percentage of species still exists input
# is single run, then time 1, and time 2
countpersistance <- function(df, t1 = 10, t2 = 5) {
    for (i in 1:length(df$timeofdeath)) {
        # for each taxa correct time of death to now
        if (df$timeofdeath[i] == -1) {
            df$timeofdeath[i] <- 0
        }
    }
    z1 <- df[which(df$timeatbirth >= t1 & df$timeofdeath <= t1), ]  #z1 is alive at t1
    z2 <- z1[which(z1$timeofdeath <= t2), ]  #z2 finds those of z1 that are alive at t2
    z3 <- df[which(df$timeatbirth >= t2 & df$timeofdeath <= t2), ]  #z2 is those alive at t2
    output <- (c(length(z1[, 1]), length(z2[, 1]), length(z3[, 1])))  #binds them together
    names(output) <- c("alivet1", "aliveboth", "alivet2")  #names them
    return(output)
}

# ok, so now let's try given it is a good species, how long did it
# take? input is a single run
timegivengood <- function(df) {
    output <- c()  #clears output to hold
    for (i in 1:length(df$speciationcomplete)) {
        # for each one if it's a good species, but not our starting one
        if (df$speciationcomplete[i] != 15 && df$speciationcomplete[i] !=
            -1) {
            output <- c(output, df$timeatbirth[i] - df$speciationcomplete[i])  #how long did it take?
        }
    }
    return(output)
}

# ok, now so ratio of good vs incipient. Let's just return the number
# of each so it's usable in other stuff
numgoodincip <- function(df) {
    liveincip <- 0
    livegood <- 0
    deadincip <- 0
    deadgood <- 0
    for (i in 1:length(df$timeatbirth)) {
        # for each one if it didnt die if it completed speciation
        if (df$timeofdeath[i] == -1) {
            if (df$speciationcomplete[i] != -1) {
                livegood <- livegood + 1  #its a good species
            } else if (df$speciationcomplete[i] == -1) {
                # if it didnt complete speciation
                liveincip <- liveincip + 1  #its an incipient species
            }
        } else if (df$timeofdeath[i] != -1) {
            # if it died and it completed speciation
            if (df$speciationcomplete[i] != -1) {
                deadgood <- deadgood + 1  #its a dead good species
            } else if (df$speciationcomplete[i] == -1) {
                # if it didnt complete speciation
                deadincip <- deadincip + 1  #its a dead incipient
            }
        }
    }
    return(c(liveincip, livegood, deadincip, deadgood))  #returns it as a list
}

pbdsim=function(pars,totaltime){
  good=list()
  incipient=list()
  good[[1]]=c(1,totaltime,totaltime,-1)
  incipient[[1]]=c(2,totaltime,-1,-1)
  taxaid=3
  deadgood=list()
  deadincipient=list()
  t=0
  while(t<=totaltime){
    numgood=length(good)
    numincipient=length(incipient)
    if(numgood==0 && numincipient==0){break}
    probs=c(pars[1]*numgood,pars[2]*numincipient,pars[3]*numincipient,pars[4]*numgood,pars[5]*numincipient)
    denom=sum(probs)
    probs=probs/denom
    t=t-log(runif(1))/denom
    if(t>=totaltime){break}
    event=sample(1:5,1,prob=probs)
    if(event==1){
      #new incipient from good
      incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1)
      taxaid=taxaid+1
    } else if (event==2){
      #new good from incipient
      take=sample(1:numincipient,1)
      good[[numgood+1]]=incipient[[take]]
      good[[numgood+1]][3]=totaltime-t
      incipient[[take]]=NULL
    } else if (event==3){
      #new from incipient
      incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1)
      taxaid=taxaid+1
    } else if (event==4){
      #dead good
      totaldeadgood=length(deadgood)
      take=sample(1:numgood,1)
      deadgood[[totaldeadgood+1]]=good[[take]]
      deadgood[[totaldeadgood+1]][4]=totaltime-t
      good[[take]]=NULL
    } else if (event==5){
      #dead incipient
      totaldeadincipient=length(deadincipient)
      take=sample(1:numincipient,1)
      deadincipient[[totaldeadincipient+1]]=incipient[[take]]
      deadincipient[[totaldeadincipient+1]][4]=totaltime-t
      incipient[[take]]=NULL
    }
  }
  return(c(good,incipient,deadgood,deadincipient))
}


repsim <- function(pars,n,time=15){
  x=c()
    i=1
    while(i <= n){
    y=c()
    for(j in pbdsim(pars,time)){y=rbind(y,j)}
    y=cbind(y,i)
        x=rbind(x,y)
        i=i+1
    }
  x=as.data.frame(x,row.names=F)
  names(x)=c("taxalabel","timeatbirth","speciationcomplete","timeofdeath","run")
    return(x)
}

sumfunct<-function(x,time){
  z=data.frame(0:time,rep(0,time+1),rep(0,time+1),rep(0,time+1),rep(0,time+1),rep(0,time+1))
  names(z)=c("time", "livingspecies", "livingincipient", "extinctspecies", "extinctincipient", "alltaxa")
  for(i in 1:length(x[,1])){
    zz=floor(x[i,])
    if(zz$speciationcomplete==-1){
      if(zz$timeofdeath==-1){
        z$livingincipient[(zz$timeatbirth+1):1]=z$livingincipient[(zz$timeatbirth+1):1]+1}
      else if(zz$timeofdeath!=-1){
        if(zz$timeatbirth!=zz$timeofdeath){
          z$livingincipient[(zz$timeatbirth+1):(zz$timeofdeath+2)]=z$livingincipient[(zz$timeatbirth+1):(zz$timeofdeath+2)]+1
          z$extinctincipient[(zz$timeofdeath+1):1]=z$extinctincipient[(zz$timeofdeath+1):1]+1}
        else {
          z$extinctincipient[(zz$timeofdeath+1):1]=z$extinctincipient[(zz$timeofdeath+1):1]+1}
      }
    }
    else if(zz$speciationcomplete!=-1){
      if(zz$timeofdeath==-1){
        if(zz$timeatbirth!=zz$speciationcomplete){
          z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]=z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]+1
          z$livingspecies[(zz$speciationcomplete+1):1]=z$livingspecies[(zz$speciationcomplete+1):1]+1}
        else if(zz$timeatbirth==zz$speciationcomplete){
          z$livingspecies[(zz$speciationcomplete+1):1]=z$livingspecies[(zz$speciationcomplete+1):1]+1}
      }
      else if(zz$timeofdeath!=-1){
        if(zz$timeatbirth!=zz$speciationcomplete && zz$speciationcomplete!=zz$timeofdeath){
          z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]=z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]+1
          z$livingspecies[(zz$speciationcomplete+1):(zz$timeofdeath+2)]=z$livingspecies[(zz$speciationcomplete+1):(zz$timeofdeath+2)]+1
          z$extinctspecies[(zz$timeofdeath+1):1]=z$extinctspecies[(zz$timeofdeath+1):1]+1}
        else if(zz$timeatbirth!=zz$speciationcomplete && zz$speciationcomplete==zz$timeofdeath){
          z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]=z$livingincipient[(zz$timeatbirth+1):(zz$speciationcomplete+2)]+1
          z$extinctspecies[(zz$timeofdeath+1):1]=z$extinctspecies[(zz$timeofdeath+1):1]+1}
        else if(zz$timeatbirth==zz$speciationcomplete && zz$speciationcomplete!=zz$timeofdeath){
          z$livingspecies[(zz$speciationcomplete+1):(zz$timeofdeath+2)]=z$livingspecies[(zz$speciationcomplete+1):(zz$timeofdeath+2)]+1
          z$extinctspecies[(zz$timeofdeath+1):1]=z$extinctspecies[(zz$timeofdeath+1):1]+1}
        else if(zz$timeatbirth==zz$speciationcomplete && zz$speciationcomplete==zz$timeofdeath){
          z$extinctspecies[(zz$timeofdeath+1):1]=z$extinctspecies[(zz$timeofdeath+1):1]+1}
      }
    }

  }
  z$alltaxa=z$livingspecies+z$livingincipient
  return(z)
}
#loop it over all the repeats
#gives a data frame with all the values from the output mapped onto each run
loopfunct<-function(x,n,time=15){
  holder=c()
  for(i in 1:n){
    tempout=sumfunct(x[which(x$run==i),],time)
    tempout$run=i
    holder=rbind(holder,tempout)
    }
  return(holder)
  }

#ddply solution - much slower? weird
#loopfunct<-function(x,time=15){
#	ddply(x,.(run),function(df){sumfunct(df,time)})
#	}

#get means and sd!
dplyframe<-function(x){
  z=ddply(x,.(time),function(df){
  c(mean(df$livingspecies),mean(df$livingincipient),mean(df$extinctspecies),mean(df$extinctincipient),mean(df$alltaxa),
  sd(df$livingspecies),sd(df$livingincipient),sd(df$extinctspecies),sd(df$extinctincipient),sd(df$alltaxa))
  })
  names(z)=c("time","meanlivingsp","meanlivingin","meanextsp","meanextin","meantaxa",
    "sdlivingsp","sdlivingin","sdextsp","sdextin","sdtaxa")
  return(z)
  }

#put it all together
summaryrepsim<-function(pars,n,time){
  dplyframe(loopfunct(repsim(pars,n,time),n,time))
  }

subsetdata<-function(data,var1,val1,var2,val2,var3,val3){
  data=as.data.frame(data)
  data2=data[data[,var1]==val1&data[,var2]==val2&data[,var3]==val3,]
  return(data2)
}

titleplot<-function(x){
  if(x=="a"){return ("good speciation rate")}
  else if (x=="b"){return("speciation completion rate")}
  else if (x=="d"){return("incipient speciation rate")}
  else if (x=="e"){return("good species extinction rate")}
  else if (x=="f"){return("incipient extinction rate")}
}

plotcontour<-function(data,variable,xx,yy,var1,val1,var2,val2,var3,val3,logged=F,numbins){
  holding<<-subsetdata(data,var1,val1,var2,val2,var3,val3)
  if(logged==T){holding[,variable]=log(holding[,variable])}
  return(ggplot(holding,aes(holding[,xx],holding[,yy],z=holding[,variable]),environment=environment()) + stat_contour(aes(colour = ..level..),environment=environment(),bins=numbins)+ xlim(0, 1)+ylim(0,1)+xlab(titleplot(xx))+ylab(titleplot(yy)))
}

plotsim<-function(x){
  ggplot()+
    geom_line(data = x, aes(x = -time, y = meanlivingsp,colour=brewer.pal(5, "Accent")[1]))+
    geom_line(data = x, aes(x = -time, y = meanlivingin,colour=brewer.pal(5, "Accent")[2]))+
    geom_line(data = x, aes(x = -time, y = meanextsp,colour=brewer.pal(5, "Accent")[3]))+
    geom_line(data = x, aes(x = -time, y = meanextin,colour=brewer.pal(5, "Accent")[4]))+
    geom_line(data = x, aes(x = -time, y = meantaxa,colour=brewer.pal(5, "Accent")[5]))+
    geom_ribbon(data = x,aes(x = -time,ymax=meanlivingsp+sdlivingsp, ymin=meanlivingsp-sdlivingsp, fill=brewer.pal(5, "Accent")[1]),alpha = 0.2)+
    geom_ribbon(data = x,aes(x = -time,ymax=meanlivingin+sdlivingin, ymin=meanlivingin-sdlivingin, fill=brewer.pal(5, "Accent")[2]),alpha = 0.2)+
    geom_ribbon(data = x,aes(x = -time,ymax=meanextsp+sdextsp, ymin=meanextsp-sdextsp, fill=brewer.pal(5, "Accent")[3]),alpha = 0.2)+
    geom_ribbon(data = x,aes(x = -time,ymax=meanextin+sdextin, ymin=meanextin-sdextin, fill=brewer.pal(5, "Accent")[4]),alpha = 0.2)+
    geom_ribbon(data = x,aes(x = -time,ymax=meantaxa+sdtaxa, ymin=meantaxa-sdtaxa, fill=brewer.pal(5, "Accent")[5]),alpha = 0.2)+
    xlab("Time (Mya)")+
    ylab("Number of taxa")+
    scale_fill_discrete(name="",
                        labels=c("Total taxa", "Living species", "Living Incipient","Extinct species","Extinct incipient"))+
    scale_colour_hue(name="",
                     labels=c("Total taxa", "Living species", "Living Incipient","Extinct species","Extinct incipient"))
}

tau<-function(l2,l3,m2){
  D=sqrt((l2+l3)^2+2*((l2-l3)*m2)+m2^2)
  return((2/(D-l2+l3-m2))*log(2/(1+((l2-l3+m2)/D))))
}



tauloop<-function(x,variable){
  holding=c()
  if(variable=="incipext"){
    for(i in seq(from=0.05, to=1,by=0.01)){
      for(j in seq(from=0.05, to=1,by=0.01)){
        holding=rbind(holding,c(i,j,tau(i,j,x)))
      }
    }
  }	else if(variable=="incipsp"){
    for(i in seq(from=0.05, to=1,by=0.01)){
      for(j in seq(from=0.05, to=1,by=0.01)){
        holding=rbind(holding,c(i,j,tau(i,x,j)))
      }
    }
  }	else if(variable=="speccomp"){
    for(i in seq(from=0.05, to=1,by=0.01)){
      for(j in seq(from=0.05, to=1,by=0.01)){
        holding=rbind(holding,c(i,j,tau(x,i,j)))
      }
    }
  }
  holding=as.data.frame(holding)
  names(holding)=c("x","y","z")
  holding=na.omit(holding)
  holding=holding[which(holding$z!=Inf),]
  return(holding)
}

plottau<-function(holding,vars,max){
  if(vars=="incipext"){
    d="speccomp"
    e="incipsp"}
  else if (vars=="incipsp"){
    d="speccomp"
    e="incipext"}
  else if (vars=="speccomp"){
    d="incipsp"
    e="incipext"}
  v <- ggplot(holding, aes(x, y, z = z))
  print(v+geom_tile(aes(fill = z)) + stat_contour(bins=20)+xlab(d)+ylab(e)+ scale_fill_gradient(limits=c(0, max)))
}