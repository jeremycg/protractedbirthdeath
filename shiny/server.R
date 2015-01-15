libraries=c("DDD","PBD","plyr","dplyr","ggplot2","RColorBrewer","compiler","shiny","data.table")
for(i in libraries){if(i %in% rownames(installed.packages()) == FALSE) {install.packages(i)}}
lapply(libraries, require, character.only=T)
library(protractedbirthdeath)
setwd(choose.dir())

x=fread("compiled.csv")
y=x[which(x$time==0),]
y$cvlivingsp=exp(y$sdlivingsp/(y$meanlivingsp+0.00001))
y$cvlivingin=exp(y$sdlivingin/(y$meanlivingin+0.00001))
y$cvextsp=exp(y$sdextsp/(y$meanextsp+0.00001))
y$cvextin=exp(y$sdextin/(y$meanextin+0.00001))
y$cvtaxa=exp(y$sdtaxa/(y$meantaxa+0.00001))
y1=as.matrix(y)

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
	data2<-summaryrepsim2(c(input$goodrate2,input$compprate2,input$inciprate2,input$extgood2,input$extincip2),15,15)
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
	x2<-treemaker2(x1)
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
