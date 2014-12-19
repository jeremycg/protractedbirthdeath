# this is the meat of the functions, actually runs the sims. takes a
# list of 5 parameters and a time
pbdsim2 <- function(pars, totaltime = 15) {
  overload=0
  good <- list()  #list of good species
  incipient <- list()  #incipients
  good[[1]] <- c(1, totaltime, totaltime, -1, 0,0)  #sim starts with one good species
  incipient[[1]] <- c(2, totaltime, -1, -1, 1,1)  #and one incipient species
  # an individual looks like: taxaid, birth time, time at 'good'
  # transition,time at death, parent, effective parent if one of these hasn't happened yet,
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
    if(sum(probs)==0){break}
    denom <- sum(probs)  #we sum them for total rate per time point
    probs <- probs/denom  #and then correct to 1 for the weighting
    t=t+rexp(1,denom)  #this increases time by doob gillespie
    # it's equivalent to t<-t-log(runif(1))/denom, but a little faster
    if (t >= totaltime) {
      break
    }  #if we just went over time, break
    event <- sample(1:5, 1, prob = probs)  #5 things that can happen, weighted by probs
    if (event == 1) {
      # new incipient from good
      take <- sample(1:numgood, 1)  #choose a parent
      incipient[[numincipient + 1]] <- c(taxaid, totaltime - t, -1,
                                         -1, rep(good[[take]][1],2))  #add incipient with time parent etc
      taxaid <- taxaid + 1  #increase taxaid for parental tracking
    } else if (event == 2) {
      # new good from incipient
      take <- sample(1:numincipient, 1)  #choose one
      good[[numgood + 1]] <- incipient[[take]]  #add it to good
      good[[numgood + 1]][3] <- totaltime - t  #change its time to be good to t
      incipient[[take]] <- NULL  #delete it from incipient
    } else if (event == 3) {
      # new from incipient
      take <- sample(1:numincipient, 1)  #choose one
      incipient[[numincipient + 1]] <- c(taxaid, totaltime - t, -1,
                                         -1, incipient[[take]][1],incipient[[take]][6])  #add it in
      taxaid <- taxaid + 1
    } else if (event == 4) {
      # dead good
      totaldeadgood <- length(deadgood)  #gets length
      take <- sample(1:numgood, 1)  #chooses one to die
      deadgood[[totaldeadgood + 1]] <- good[[take]]  #adds it
      deadgood[[totaldeadgood + 1]][4] <- totaltime - t  #adds death time
      if(numincipient>0){
        testparent<-good[[take]][6]
        offspringofdead<-matrix(unlist(incipient), nrow=numincipient, byrow=T)[,6]==testparent
        if(sum(offspringofdead)!=0){
          ordered<-sample(which(offspringofdead))
          chosen<-ordered[1]
          neweffectiveparent<-incipient[[chosen]][1]
          for(i in 1:length(ordered)-1){
            incipient[[ordered[i+1]]][6]<-neweffectiveparent
          }
          #z %>% map_if(function(df){df[6]==testparent},function(df2){df2[6]<-neweffectiveparent;df2})
          #potential loop replacement
          #but much slower for now, see how the lowliner package speeds up over time
          good[[numgood + 1]] <- incipient[[chosen]]
          good[[numgood + 1]][3] <- totaltime - t
          incipient[[chosen]] <- NULL
        }
      }
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
  output<-c()
  x <- rlply(n,pbdsim2(pars,time))
  xx<-lapply(x,combinelists)
  lengths<-lapply(xx,nrow)
  zzz<-rbind_all(xx)
  zzz[,7]<-unlist(mapply(rep,1:n,lengths,SIMPLIFY = F))
  names(zzz)<- c("taxalabel", "timeatbirth", "speciationcomplete", "timeofdeath",
                           "parent", "effective parent","run")
  return(zzz)
}

combinelists<-function(x){
  data.frame(matrix(unlist(x), nrow=length(x), byrow=T))
}

# a function to make a newick tree text from a given output is not of
# class phy - need to do read.tree(text=treemaker(x)) for that
treemaker <- function(x) {
    x<-as.data.frame(x)
    x$timeofdeath[x$timeofdeath==-1]<-0
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
  df[df$timeofdeath==-1,]$timeofdeath<-0
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
  out=df$timeatbirth - df$speciationcomplete
  return(out[df$speciationcomplete != 15 & df$speciationcomplete !=-1])
}

# ok, now so ratio of good vs incipient. Let's just return the number
# of each so it's usable in other stuff
numgoodincip <- function(df) {
    liveincip <- sum(df$timeofdeath==-1&df$speciationcomplete==-1)
    livegood <- sum(df$timeofdeath==-1&df$speciationcomplete!=-1)
    deadincip <- sum(df$timeofdeath!=-1&df$speciationcomplete==-1)
    deadgood <- sum(df$timeofdeath!=-1&df$speciationcomplete!=-1)
    return(c(liveincip, livegood, deadincip, deadgood))  #returns it as a list
}

#deprecated
#use pddsim2
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

#this is deprecated
#slow, and leaves stuff out
#use repsim2
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

sumfunctpart1<-function(x,time=15){
  output<-c(rep(0,time+1))
  output[(time+1-x[[2]]):(time+2)]<-1
  if(x[[3]]>=x[[4]]){
    output[(time+1-x[[3]]):(time+2)]<-2
    output[(time+1-x[[4]]):(time+2)]<-3
  } else{
    output[(time+1-x[[4]]):(time+2)]<-4
  }
  return(output[0:time+1])
}

sumfunctpart2<-function(downcol){
  incipient<-sum(downcol==1)
  good<-sum(downcol==2)
  deadgood<-sum(downcol==3)
  deadincipient<-sum(downcol==4)
  return(c(good,incipient,deadgood,deadincipient))
}


sumfunct<-function(df,time=15){
  df[2:4]<-floor(df[2:4])
  out<-as.data.frame(t(apply(apply(df,1,sumfunctpart1,time=time),1,sumfunctpart2)))
  names(out)<-c("livingspecies","livingincipient","extinctspecies","extinctincipient")
  out$time=seq(from=15,to=0)
  out$alltaxa<-out[,1]+out[,2]
  return(out)
}

loopfunct<-function(x,time=15){
 ddply(x,.(run),function(df){sumfunct(df,time)})
}

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
  dplyframe(loopfunct(repsim2(pars,n,time),time))
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



taufunct<-function(var2,x,y,fixed){
  if(var2=="incipext"){
    return(tau(x,y,fixed))
  } else if(var2=="incipsp"){
    return(tau(x,fixed,y))
  } else if(var2=="speccomp"){
    return(tau(fixed,x,y))
  }
}

tauloop<-function(x,var){
  a<-seq(from=0.05,to=1,by=0.01)
  b<-seq(from=0.05,to=1,by=0.01)
  holding<-melt(outer(a,b,FUN="taufunct",var2=var,fixed=x))
  holding[,1]<-0.01*holding[,1]+0.04
  holding[,2]<-0.01*holding[,2]+0.04
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
