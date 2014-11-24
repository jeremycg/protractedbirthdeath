#load a bunch of stuff
libraries=c("DDD","PBD","plyr","dplyr","ggplot2","RColorBrewer","compiler","shiny","data.table")
for(i in libraries){if(i %in% rownames(installed.packages()) == FALSE) {install.packages(i)}}
lapply(libraries, require, character.only=T)

pbdsim2=function(pars,totaltime=15){
	good=list()
	incipient=list()
	good[[1]]=c(1,totaltime,totaltime,-1,0)
	incipient[[1]]=c(2,totaltime,-1,-1,1)
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
			take=sample(1:numgood,1)
			incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1,good[[take]][1])
			taxaid=taxaid+1
		} else if (event==2){
			#new good from incipient
			take=sample(1:numincipient,1)
			good[[numgood+1]]=incipient[[take]]
			good[[numgood+1]][3]=totaltime-t
			incipient[[take]]=NULL
		} else if (event==3){
			#new from incipient
			take=sample(1:numincipient,1)
			incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1,incipient[[take]][1])
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


repsim2 <- function(pars,n,time=15){
	x=c()
    i=1
    while(i <= n){
		y=c()
		for(j in pbdsim2(pars,time)){y=rbind(y,j)}
		y=cbind(y,i)
        x=rbind(x,y)
        i=i+1
    }
	x=as.data.frame(x,row.names=F)
	names(x)=c("taxalabel","timeatbirth","speciationcomplete","timeofdeath","parent","run")
    return(x)
}

treemaker <- function(x){
  for(i in 1:length(x[,1])){
      if(x[i,]$timeofdeath==-1){x[i,]$timeofdeath=0}
  }
  x$label=x$taxalabel
  while(length(x[,1])>=2){
    birthtime=min(x$timeatbirth)
    if(length(x[,1])==2){
      child1=which(x$parent==1)
      parent1=which(x$parent==0)
    }
    else{
      child1=which(x$timeatbirth==birthtime)
      parent1=which(x$taxalabel==x[child1,]$parent)
    }
    x[parent1,]$label<-paste("(",x[parent1,]$label,":",birthtime-x[parent1,]$timeofdeath,",",x[child1,]$label,":",
                        birthtime-x[child1,]$timeofdeath,")",sep="")
    x[parent1,]$timeofdeath=birthtime
    x=x[-child1,]
  }
  return(paste(x$label,";",sep=""))
}

#so far:
#how close is each species sister taxa?
#as a plyrable function

sisterlengths<-function(working){
	phy<-treemaker(working)
	phy1<-read.tree(text=phy)
	if(sum(which(working$timeofdeath==-1))<=1){
		return("no extant species")
		}
	phy2<-drop.tip(phy1,as.character(working[which(working$timeofdeath!=-1),]$taxalabel))
	return(phy2$edge.length[phy2$edge[,2] <= Ntip(phy2)])
	}

dlply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun<-sisterlengths)

#fine,,, a little cumbersome in the output, but whatever

#next, a function for "time slicing"
#we want to go at say 10 million y.a
#then today and see what percentage of species still exists
countpersistance<-function(df,t1=10,t2=5){
	for(i in 1:length(df$timeofdeath)){
		if(df$timeofdeath[i]==-1){
			df$timeofdeath[i]<-0
		}
	}
	z1<-df[which(df$timeatbirth>=t1&df$timeofdeath<=t1),]
	z2<-z1[which(z1$timeofdeath<=t2),]
	z3<-df[which(df$timeatbirth>=t2&df$timeofdeath<=t2),]
	output<-(c(length(z1[,1]),length(z2[,1]),length(z3[,1])))
	names(output)<-c("alivet1","aliveboth","alivet2")
	return(output)
	}
#outputs num alive at t1, number of those alive at t1 that made it to t2, and number alive at t2
ddply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun=countpersistance)

#ok, so now let's try given it is a good species, how long did it take?
#should be easy!
timegivengood<-function(df){
	output<-c()
	for(i in 1:length(df$speciationcomplete)){
		if(df$speciationcomplete[i]!=15&&df$speciationcomplete[i]!=-1){
			output<-c(output,df$timeatbirth[i]-df$speciationcomplete[i])
		}
	}
	return(output)
}
#and loop it
dlply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun<-timegivengood)

#ok, now so ratio of good vs incipient. Let's just return the number of each so it's usabl;e in other stuff
numgoodincip<-function(df){
	liveincip <- 0
	livegood  <- 0
	deadincip <- 0
	deadgood  <- 0
	for(i in 1:length(df$timeatbirth)){
		if(df$timeofdeath[i] == -1){
			if(df$speciationcomplete[i] != -1){
				livegood <- livegood + 1
			} else if (df$speciationcomplete[i] == -1){
				liveincip <- liveincip + 1
			}
		} else if(df$timeofdeath[i] != -1){
			if(df$speciationcomplete[i] != -1){
				deadgood <- deadgood + 1
			} else if (df$speciationcomplete[i] == -1){
				deadincip <- deadincip + 1
			}
		}
	}
	return(c(liveincip,livegood,deadincip,deadgood))
}
				
#last is a pertubation one, can just run normal sim