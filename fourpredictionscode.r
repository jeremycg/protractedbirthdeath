#load a bunch of stuff
libraries=c("DDD","PBD","plyr","dplyr","ggplot2","RColorBrewer","compiler","shiny","data.table")
for(i in libraries){if(i %in% rownames(installed.packages()) == FALSE) {install.packages(i)}}
lapply(libraries, require, character.only=T)
#terrible way to load packages...


#this is the meat of the functions, actually runs the sims.
pbdsim2=function(pars,totaltime=15){#takes a list of 5 parameters and a time
	good=list()#list of good species
	incipient=list()#incipients
	good[[1]]=c(1,totaltime,totaltime,-1,0)#sim starts with one good species
	incipient[[1]]=c(2,totaltime,-1,-1,1)#and one incipient species
	#an individual looks like taxaid, birth time, time at "good" transition,time at death, parent
	#if one of these hasn't happened yet, -1 is used
	taxaid=3#the taxa label of the next species
	deadgood=list()#all dead good individuals
	deadincipient=list()#all dead incipient
	t=0#we start at t=0
	while(t<=totaltime){#until t is bigger than time (should never hit this, break out below)
		numgood=length(good)#we want to know how many good species
		numincipient=length(incipient)#and incipient species
		if(numgood==0 && numincipient==0){break}#if everything is extinct, we break
		probs=c(pars[1]*numgood,pars[2]*numincipient,pars[3]*numincipient,pars[4]*numgood,pars[5]*numincipient)
		#probs are the probability of each event happening
		#so its prob good speciation times number of good species etc
		denom=sum(probs)#we gotta sum them for total rate per time point
		probs=probs/denom#and then correct to 1 for the weighting
		t=t-log(runif(1))/denom#this increases time by doob gillespie
		#its equivalent to t=t+rexp(1,denom)
		if(t>=totaltime){break}#if we just went over time, break
		event=sample(1:5,1,prob=probs)#we have 5 things that can happen, weighted by probs
		if(event==1){
			#new incipient from good
			take=sample(1:numgood,1)#choose a parent
			incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1,good[[take]][1])#add incipient with time parent etc 
			taxaid=taxaid+1#increase taxaid for parental tracking
		} else if (event==2){
			#new good from incipient
			take=sample(1:numincipient,1)#choose one
			good[[numgood+1]]=incipient[[take]]#add it to good
			good[[numgood+1]][3]=totaltime-t#change its time to be good rate to t
			incipient[[take]]=NULL#delete it from incipient
		} else if (event==3){
			#new from incipient
			take=sample(1:numincipient,1)#choose one
			incipient[[numincipient+1]]=c(taxaid,totaltime-t,-1,-1,incipient[[take]][1])#add it in
			taxaid=taxaid+1
		} else if (event==4){
			#dead good
			totaldeadgood=length(deadgood)#gets length
			take=sample(1:numgood,1)#chooses one to die
			deadgood[[totaldeadgood+1]]=good[[take]]#adds it
			deadgood[[totaldeadgood+1]][4]=totaltime-t#adds death time
			good[[take]]=NULL#removes it
		} else if (event==5){
			#dead incipient
			totaldeadincipient=length(deadincipient)#gets length
			take=sample(1:numincipient,1)#takes it
			deadincipient[[totaldeadincipient+1]]=incipient[[take]]#adds it
			deadincipient[[totaldeadincipient+1]][4]=totaltime-t#adds death time
			incipient[[take]]=NULL#removes it
		}
	}
	return(c(good,incipient,deadgood,deadincipient))#returns a list with all taxa as above
}

#a function to repeat the above a ton of times, and store the output in a data frame
repsim2 <- function(pars,n,time=15){
	x=c()#holding for output
    i=1#repeat number
    while(i <= n){
		y=c()#temp holding
		for(j in pbdsim2(pars,time)){y=rbind(y,j)}#run the sim
		y=cbind(y,i)#add on the run number column
        x=rbind(x,y)#add it all to output
        i=i+1#increase repeats
    }
	x=as.data.frame(x,row.names=F)#converts to dataframe
	names(x)=c("taxalabel","timeatbirth","speciationcomplete","timeofdeath","parent","run")#names outputs
    return(x)#output
}

#a function to make a newick tree text from a given output
#is not of class phy - need to do read.tree(text=treemaker(x)) for that
treemaker <- function(x){
  for(i in 1:length(x[,1])){#for each taxa
      if(x[i,]$timeofdeath==-1){x[i,]$timeofdeath=0}#changes -1s to 0s
  }
  x$label=x$taxalabel#makes a new column to hold the tree strings
  while(length(x[,1])>=2){#while there are more than two taxa
    birthtime=min(x$timeatbirth)#find the yongest taxa
    if(length(x[,1])==2){#if there are two
      child1=which(x$parent==1)#it breaks as we started with two taxa at t=end
      parent1=which(x$parent==0)#so manually set the parent/child
	  #as the rest is done by stepwise functions we should never have more than one child at the same timepoint
	  #unless we hit machine precision, which is very very unlikely
    }
    else{
      child1=which(x$timeatbirth==birthtime)#the child is the one with the birthtime=now
      parent1=which(x$taxalabel==x[child1,]$parent)#the parent is the parent of that one
    }
    x[parent1,]$label<-paste("(",x[parent1,]$label,":",birthtime-x[parent1,]$timeofdeath,",",x[child1,]$label,":",
                        birthtime-x[child1,]$timeofdeath,")",sep="")
	#this changes the label of the parent, to newick format
	#should look like (parent:distance,child:distance)
	#but, newick can nest, so the parent or child part can be a previous label
	#so we can recursively make a tree
    x[parent1,]$timeofdeath=birthtime#collapses the node, so the parent node now dies at teh split
    x=x[-child1,]#removes the child
  }
  return(paste(x$label,";",sep=""))#returns the string, with the required tailing ;
}

#so far:
#how close is each species sister taxa?
#as a plyrable function
#looks at all extant taxa
#finds closest relative branch time - single, not both
sisterlengths<-function(working){#run on a single run from repsim2
	phy<-treemaker(working)#make a tree
	phy1<-read.tree(text=phy)#make it of class phy
	if(sum(which(working$timeofdeath==-1))<=1){#if we had no taxa at the end
		return("no extant species")
		}
	phy2<-drop.tip(phy1,as.character(working[which(working$timeofdeath!=-1),]$taxalabel))#this drops all dead taxa
	return(phy2$edge.length[phy2$edge[,2] <= Ntip(phy2)])#returns the distance to node for each taxa
	#this works as R stores edges in order - we only want the edges which lead to species, not nodes
	#so we want the ones that are less than the number of tips
	}
#output is a list - it depends on number of species, so can't be a dataframe

dlply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun=sisterlengths)
#needs dlply die to list. migth be worht making a dataframe containing lists, but lists are good for now

#next, a function for "time slicing"
#we want to go at say 10 million y.a
#then today and see what percentage of species still exists
#input is single run, then time 1, and time 2
countpersistance<-function(df,t1=10,t2=5){
	for(i in 1:length(df$timeofdeath)){#for each taxa
		if(df$timeofdeath[i]==-1){#correct time of death to now
			df$timeofdeath[i]<-0
		}
	}
	z1<-df[which(df$timeatbirth>=t1&df$timeofdeath<=t1),]#z1 is alive at t1
	z2<-z1[which(z1$timeofdeath<=t2),]#z2 finds those of z1 that are alive at t2
	z3<-df[which(df$timeatbirth>=t2&df$timeofdeath<=t2),]#z2 is thos alive at t2
	output<-(c(length(z1[,1]),length(z2[,1]),length(z3[,1])))#binds them together
	names(output)<-c("alivet1","aliveboth","alivet2")#names them
	return(output)
	}
#outputs num alive at t1, number of those alive at t1 that made it to t2, and number alive at t2
ddply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun=countpersistance)
#dataframe out, as it will always give 3 numbers

#ok, so now let's try given it is a good species, how long did it take?
#inout is a single run
timegivengood<-function(df){
	output<-c()#clears output to hold
	for(i in 1:length(df$speciationcomplete)){#for each one
		if(df$speciationcomplete[i]!=15&&df$speciationcomplete[i]!=-1){#if it's a good species, but not our starting one
			output<-c(output,df$timeatbirth[i]-df$speciationcomplete[i])#how long did it take?
		}
	}
	return(output)
}
#and loop it
dlply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun=timegivengood)
#again, its a list as it is one number peroutput species

#ok, now so ratio of good vs incipient. Let's just return the number of each so it's usabl;e in other stuff
numgoodincip<-function(df){
	liveincip <- 0
	livegood  <- 0
	deadincip <- 0
	deadgood  <- 0
	for(i in 1:length(df$timeatbirth)){#for each one
		if(df$timeofdeath[i] == -1){#if it didnt die
			if(df$speciationcomplete[i] != -1){#if it completed speciation
				livegood <- livegood + 1#its a good species
			} else if (df$speciationcomplete[i] == -1){#if it didnt complete speciation
				liveincip <- liveincip + 1#its an incipient species
			}
		} else if(df$timeofdeath[i] != -1){#if it died
			if(df$speciationcomplete[i] != -1){#and it completed speciation
				deadgood <- deadgood + 1#its a dead good species
			} else if (df$speciationcomplete[i] == -1){#if it didnt complete speciation
				deadincip <- deadincip + 1#its a dead incipient
			}
		}
	}
	return(c(liveincip,livegood,deadincip,deadgood))#returns it as a list
}
#loopit
ddply(repsim2(c(0.1,0.1,0.1,0.1,0.1),20),.(run),.fun=numgoodincip)
#again, it can be a dataframe as it always returns 4 values.
#last is a pertubation one, can just run normal sim