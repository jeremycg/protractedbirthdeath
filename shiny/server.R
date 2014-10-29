libraries=c("DDD","PBD","plyr","dplyr","ggplot2","RColorBrewer","compiler","shiny","data.table")
for(i in libraries){if(i %in% rownames(installed.packages()) == FALSE) {install.packages(i)}}
lapply(libraries, require, character.only=T)
setwd(choose.dir())	
		
	

		
x=fread("compiled.csv")
y=x[which(x$time==0),]
y$cvlivingsp=exp(y$sdlivingsp/(y$meanlivingsp+0.00001))
y$cvlivingin=exp(y$sdlivingin/(y$meanlivingin+0.00001))
y$cvextsp=exp(y$sdextsp/(y$meanextsp+0.00001))
y$cvextin=exp(y$sdextin/(y$meanextin+0.00001))
y$cvtaxa=exp(y$sdtaxa/(y$meantaxa+0.00001))
y1=as.matrix(y)
 
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
 
shinyServer(function(input, output) {
 
 
	output$Plot1 <- renderPlot({
	data2<-x[which(x$a==input$goodrate&x$b==input$compprate&x$d==input$inciprate&x$e==input$extgood&x$f==input$extincip),]
	print(plotsim(data2))
  })
   output$Plot2 <- renderPlot({
	data2<-x[which(x$a==input$goodrate2&x$b==input$compprate2&x$d==input$inciprate2&x$e==input$extgood2&x$f==input$extincip2),]
	print(plotsim(data2))
  })
 
	output$Plot3 <- renderPlot({
	data2<-summaryrepsim(c(input$goodrate,input$compprate,input$inciprate,input$extgood,input$extincip),15,15)
	print(plotsim(data2))
	})
	output$Plot4 <- renderPlot({
	data2<-summaryrepsim(c(input$goodrate2,input$compprate2,input$inciprate2,input$extgood2,input$extincip2),15,15)
	print(plotsim(data2))
	})
	output$Plot5 <- renderPlot({
	boxplot(log(y1[,which(names(y)==input$variable)]+1)~y$a,main="speciation rate of good species:")
	})
	output$Plot6 <- renderPlot({
	boxplot(log(y1[,which(names(y)==input$variable)]+1)~y$b,main="speciation completion rate:")
	})
	output$Plot7 <- renderPlot({
	boxplot(log(y1[,which(names(y)==input$variable)]+1)~y$d,main="speciation rate of incipient species:")
	})
	output$Plot8 <- renderPlot({
	boxplot(log(y1[,which(names(y)==input$variable)]+1)~y$e,main="extinction rate of good species:")
	})
	output$Plot9 <- renderPlot({
	boxplot(log(y1[,which(names(y)==input$variable)]+1)~y$f,main="extinction rate of incipient species:")
	})
	output$Plot10<- renderPlot({
	print(plotcontour(y,input$variable1,input$plotx,input$ploty,input$fixed1,input$fixed1rate,input$fixed2,input$fixed2rate,input$fixed3,input$fixed3rate,input$lograte1,input$binrate1))
	})
	output$Plot11<- renderPlot({
	print(plotcontour(y,input$variable2,input$plotx2,input$ploty2,input$fixed12,input$fixed1rate2,input$fixed22,input$fixed2rate2,input$fixed32,input$fixed3rate2,input$lograte2,input$binrate2))
	})
	output$Plot12<- renderPlot({
	x1<-repsim2(c(input$goodrate,input$compprate,input$inciprate,input$extgood,input$extincip),1,15)
	x2<-treemaker(x1)
	plot(read.tree(text=x2))
	})
	output$Plot13<- renderPlot({
	x1<-repsim2(c(input$goodrate2,input$compprate2,input$inciprate2,input$extgood2,input$extincip2),1,15)
	x2<-treemaker(x1)
	plot(read.tree(text=x2))
	})
	output$Plot14<-renderPlot({
		plottau(tauloop(input$taurate,input$tauvariable),input$tauvariable,input$maxcolour)
	})
})