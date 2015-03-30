#' Simulate speciation by a protracted birth death process.
#' @param pars a vector of 5 parameters: c(speciation rate of good species,
#'   speciation rate of incipient species, completion rate, death rate of
#'   good species and death rate of incipient species).
#' @param totaltime time to run simulation, defaults to 15
#' @return a list with each taxa as its own entry, with birth and death times
#'   has an atrribute overloaded=1 if too many species were formed
#' @seealso \code{\link{repsim2}} which wraps this function and does repeats
#' @export
#' @examples
#'\dontrun{
#' pbdsim2(c(0.1,0.1,0.1,0.1,0.1),15)
#' pbdsim2(c(0.2,0.2,0.2,0.1,0.1),15)
#'}
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

#' Repeat \code{\link{pbdsim2}} multiple times and output in a dataframe.
#' @param pars a vector of 5 parameters: c(speciation rate of good species,
#'   speciation rate of incipient species, completion rate, death rate of
#'   good species and death rate of incipient species).
#' @param n a number of times to repeat the simulation
#' @param time time to run simulation, defaults to 15
#' @return a dataframe with each taxa as its own entry, with birth and death times
#' @seealso \code{\link{pbdsim2}} which is the single repeat version,
#'   \code{\link{summaryrepsim}} which takes averages of repeats
#' @export
#' @examples
#'\dontrun{
#' repsim2(c(0.1,0.1,0.1,0.1,0.1),15,15)
#' repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15)
#' summaryrepsim(c(0.1,0.1,0.1,0.1,0.1),15,15)
#'}
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

#' An internal function to quickly combine lists
#' @param x a list
#' @return a dataframe containing sensibly formatted data for use in \code{\link{repsim2}}
combinelists<-function(x){
  data.frame(matrix(unlist(x), nrow=length(x), byrow=T))
}

#' Takes a single inputted run and makes a (text) tree.
#' deprecated - much slower than treemaker2
#' @param x a single run from repsim2
#' @return a string which can be read into ape using read.tree
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @seealso \code{\link{treemaker2}} a much faster implementation
#' @export
#' @examples
#'\dontrun{
#' x<-treemaker(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#' phy<-read.tree(text=x)
#' plot(phy)
#'}
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

#' Takes a single inputted run and makes a (text) tree.
#' much faster than treemaker
#' @param z a single run from repsim2
#' @return a string which can be read into ape using read.tree
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
#' @examples
#'\dontrun{
#' x<-treemaker2(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#' phy<-read.tree(text=x)
#' plot(phy)
#'}
treemaker2<-function(z){
  z$timeofdeath[z$timeofdeath==-1]<-0
  z$label<-z$taxalabel
  while(nrow(z)>1){
    z<-removeyoungest(z)
  }
  paste(z$label, ";", sep = "")
}
#' Remove youngest taxa from tree
#' the recurssive function used to prune trees for tree construction in
#' \code{\link{treemaker2}}. Works by taking the youngest offspring from each parent,
#' making sure it isnt a parent. If it isnt, it is trimmed and removed from the list
#' and its parent is relabelled. Probably my hardest thought out function.
#' Significantly faster than simply removing the youngest taxa as this prunes all at once.
#' @param z a single run from repsim2
#' @return dataframe with the youngest offsprings all removed
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @examples
#'\dontrun{
#' x<-removeyoungest(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#' }
removeyoungest<-function(z){
  zdeadtest<-z%>%group_by(parent) %>% filter(timeatbirth==min(timeatbirth))
  zdead<-zdeadtest[!(zdeadtest$taxalabel %in% z$parent),]
  z<-z[!(z$taxalabel %in% zdead$taxalabel),]
  tmppar<-z$taxalabel %in% zdead$parent
  parents<-z[tmppar,]
  z<-z[!tmppar,]
  zdead<-zdead[with(zdead, order(parent)), ]
  parents<-parents[with(parents, order(taxalabel)), ]
  parents$label <- paste("(", parents$label, ":", zdead$timeatbirth -
  parents$timeofdeath, ",", zdead$label, ":", zdead$timeatbirth -
  zdead$timeofdeath, ")", sep = "")
  parents$timeofdeath<-zdead$timeatbirth
  z<-rbind.fill(z,parents)
  z
}

#' Finds branchlength to closest extant sister taxa.
#' @param working a single run from repsim2
#' @return a list of branch lengths
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
#' @examples
#'\dontrun{
#' sisterlengths(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#'}
sisterlengths <- function(working) {
    phy <- treemaker2(working)  #make a tree
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


#' Finds persistance of taxa over time.
#' gives a number of taxa alive at t1 and t2 or both
#' @param df single run from repsim2
#' @param t1 first time point
#' @param t2 second time point
#' @return a vector of 3 values - those alive at t1, both and t2
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
#' @examples
#'\dontrun{
#' countpersistance(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15),10,5)
#'}
countpersistance <- function(df, t1 = 10, t2 = 5) {
  df<-as.data.frame(df)
  df$timeofdeath[df$timeofdeath==-1]<-0
  z1 <- df[which(df$timeatbirth >= t1 & df$timeofdeath <= t1), ]  #z1 is alive at t1
  z2 <- z1[which(z1$timeofdeath <= t2), ]  #z2 finds those of z1 that are alive at t2
  z3 <- df[which(df$timeatbirth >= t2 & df$timeofdeath <= t2), ]  #z2 is those alive at t2
  output <- (c(length(z1[, 1]), length(z2[, 1]), length(z3[, 1])))  #binds them together
  names(output) <- c("alivet1", "aliveboth", "alivet2")  #names them
  return(output)
}

#' Finds time to be a good species.
#' @param df a single run from repsim2
#' @return a vector of times - equal to number of good species
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
#' @examples
#'\dontrun{
#' timegivengood(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#'}
timegivengood <- function(df) {
  out=df$timeatbirth - df$speciationcomplete
  return(out[df$speciationcomplete != 15 & df$speciationcomplete !=-1])
}

#' Finds numbers of each possible outcome at end of simulation.
#' @param df a single run from repsim2
#' @return a vector of 4 named numbers - dead and live good and incipient taxa
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
#' @examples
#'\dontrun{
#' numgoodincip(repsim2(c(0.2,0.2,0.2,0.1,0.1),1,15))
#'}
numgoodincip <- function(df) {
    liveincip <- sum(df$timeofdeath==-1&df$speciationcomplete==-1)
    livegood <- sum(df$timeofdeath==-1&df$speciationcomplete!=-1)
    deadincip <- sum(df$timeofdeath!=-1&df$speciationcomplete==-1)
    deadgood <- sum(df$timeofdeath!=-1&df$speciationcomplete!=-1)
    return(c(liveincip, livegood, deadincip, deadgood))  #returns it as a list
}
#' sumfunctpart1
#' an internal function to determine living species at each time point
#' run across (apply,1) on a repsim result to give each lines time of death etc
#' @param x a single line of a run from repsim2
#' @param time the time the simulation was run for
#' @return a vector with 0,1,2,3,4 at each integer t depending on taxa state
#' summed to give summaries
#' @seealso \code{\link{summaryrepsim}} which links everything together
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
#' sumfunctpart2
#' an internal function to determine living species at each time point
#' run across (apply,1) on a sumfunctpart1 result to give each lines time of death etc
#' @param downcol a single line of a run from sumfunct1
#' @return a matrix with taxa state for each n
#' @seealso \code{\link{summaryrepsim}} which links everything together
sumfunctpart2<-function(downcol){
  incipient<-sum(downcol==1)
  good<-sum(downcol==2)
  deadgood<-sum(downcol==3)
  deadincipient<-sum(downcol==4)
  return(c(good,incipient,deadgood,deadincipient))
}
#' sumfunct
#' an internal function to determine living species at each time point
#' @param df a single run of repsim2
#' @param time the time repsim2 was run for
#' @return a dataframe with number of species in each class at each integer time point
#' @seealso \code{\link{summaryrepsim}} which loops this over multiple runs for means and sds
sumfunct<-function(df,time=15){
  df[2:4]<-floor(df[2:4])
  out<-as.data.frame(t(apply(apply(df,1,sumfunctpart1,time=time),1,sumfunctpart2)))
  names(out)<-c("livingspecies","livingincipient","extinctspecies","extinctincipient")
  out$time=seq(from=15,to=0)
  out$alltaxa<-out[,1]+out[,2]
  return(out)
}
#' loopfunct
#' an internal function to determine living species at each time point
#' @param x a single run of repsim2
#' @param time the time repsim2 was run for
#' @return a dataframe with number of species in each class at each integer time point run over all repeats
#' @seealso \code{\link{summaryrepsim}} which loops this over multiple runs for means and sds
loopfunct<-function(x,time=15){
 ddply(x,.(run),function(df){sumfunct(df,time)})
}

#' dplyframe
#' an internal function to determine living species at each time point
#' @param x a single run of loopfunct
#' @return a dataframe with mean and sd number of species in each class at each integer time point run over all repeats
#' @seealso \code{\link{summaryrepsim}} which loops this over multiple runs for means and sds
dplyframe<-function(x){
  z=ddply(x,.(time),function(df){
  c(mean(df$livingspecies),mean(df$livingincipient),mean(df$extinctspecies),mean(df$extinctincipient),mean(df$alltaxa),
  sd(df$livingspecies),sd(df$livingincipient),sd(df$extinctspecies),sd(df$extinctincipient),sd(df$alltaxa))
  })
  names(z)=c("time","meanlivingsp","meanlivingin","meanextsp","meanextin","meantaxa",
    "sdlivingsp","sdlivingin","sdextsp","sdextin","sdtaxa")
  return(z)
}

#' summaryrepsim
#' A summary of repeats with means and sds.
#' @param pars a vector of 5 parameters: c(speciation rate of good species,
#'   speciation rate of incipient species, completion rate, death rate of
#'   good species and death rate of incipient species).
#' @param n a number of times to repeat the simulation
#' @param time time to run simulation, defaults to 15
#' @return a data frame with means and sds of numbers of taxa at each timepoint
#' @seealso \code{\link{repsim2}} which produces the inputs \code{\link{plotsim}}
#'   which plots this functions output
#' @export
#' @examples
#'\dontrun{
#' summaryrepsim(c(0.2,0.2,0.2,0.1,0.1),15,15)
#'}
summaryrepsim<-function(pars,n,time){
  dplyframe(loopfunct(repsim2(pars,n,time),time))
}
#' subsetdata
#' subsets data based on parameters
#' @param data a loaded frame of multiple results see example data for formatting
#' @param var1 the name of variable 1
#' @param val1 the selected value of variable 1
#' @param var2 the name of variable 2
#' @param val2 the selected value of variable 2
#' @param var3 the name of variable 3
#' @param val3 the selected value of variable 3
subsetdata<-function(data,var1,val1,var2,val2,var3,val3){
  data=as.data.frame(data)
  data2=data[data[,var1]==val1&data[,var2]==val2&data[,var3]==val3,]
  return(data2)
}

#' titles a plot
#' @param x a variable taken from the example data
titleplot<-function(x){
  if(x=="a"){return ("good speciation rate")}
  else if (x=="b"){return("speciation completion rate")}
  else if (x=="d"){return("incipient speciation rate")}
  else if (x=="e"){return("good species extinction rate")}
  else if (x=="f"){return("incipient extinction rate")}
}

#' A contourplot of number of species across differing parameters
#' chooses a variable - mean sd of alltaxa, good taxa etc
#' and plots it for fixed values of three chosen variables
#' @param data a loaded frame of multiple results see example data for formatting
#' @param variable a value to plot - usually mean
#' @param xx the x axis variable
#' @param yy the y axis variable
#' @param var1 the name of variable 1
#' @param val1 the selected value of variable 1
#' @param var2 the name of variable 2
#' @param val2 the selected value of variable 2
#' @param var3 the name of variable 3
#' @param val3 the selected value of variable 3
#' @param logged whether to plot logged data, defaults to F
#' @param numbins the number of bins to use in plotting contours
#' @return a ggplot with the two non named variables and countours for the chosen variable
#' @seealso \code{\link{repsim2}} which produces the inputs
#' @export
plotcontour<-function(data,variable,xx,yy,var1,val1,var2,val2,var3,val3,logged=F,numbins){
  holding<<-subsetdata(data,var1,val1,var2,val2,var3,val3)
  if(logged==T){holding[,variable]=log(holding[,variable])}
  return(ggplot(holding,aes(holding[,xx],holding[,yy],z=holding[,variable]),environment=environment()) + stat_contour(aes(colour = ..level..),environment=environment(),bins=numbins)+ xlim(0, 1)+ylim(0,1)+xlab(titleplot(xx))+ylab(titleplot(yy)))
}

#' A plot of repeats with means and sds.
#' @param x an output from summaryrepsim
#' @return a ggplot
#' @seealso \code{\link{repsim2}} which produces the inputs \code{\link{plotsim}}
#'   which plots this functions output
#' @export
#' @examples
#'\dontrun{
#' plotsim(summaryrepsim(c(0.2,0.2,0.2,0.1,0.1),15,15))
#'}
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

#' calculates tau - time to complete speciation
#' @param l2 speciation completion rate
#' @param l3 incipient speciation rate
#' @param m2 incipient species death rate
#' @return a single tau value
tau<-function(l2,l3,m2){
  D=sqrt((l2+l3)^2+2*((l2-l3)*m2)+m2^2)
  return((2/(D-l2+l3-m2))*log(2/(1+((l2-l3+m2)/D))))
}


#' calculates tau with fixed values
#' allows choosing of variables to plot
#' @param var2 variable to fix
#' @param x first other variable
#' @param y second other variable
#' @param fixed the fixed value to give var2
#' @return a single tau value
taufunct<-function(var2,x,y,fixed){
  if(var2=="incipext"){
    return(tau(x,y,fixed))
  } else if(var2=="incipsp"){
    return(tau(x,fixed,y))
  } else if(var2=="speccomp"){
    return(tau(fixed,x,y))
  }
}

#' A matrix of tau values.
#' @param var value to fix
#' @param x the value to give it
#' @return a matrix of tau, read to plot
#' @seealso \code{\link{plottau}} which plots the outputs
#' @export
#' @examples
#'\dontrun{
#' tauloop(0.5,"incipext")
#'}
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

#' A plot of tau (time to speciation)
#' @param holding an output of tauloop
#' @param vars variable to fix
#' @param max the value to clip tau on the plot at
#' @return a countour ggplot of tau
#' @seealso \code{\link{tauloop}} which produces the inputs
#'   which this function plots
#' @export
#' @examples
#'\dontrun{
#' z<-tauloop(0.5,"incipext")
#' plottau(z,"incipext",5)
#'}
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

#' A shiny app to view the simulations
#' @return a shiny app
#' @export
#' @examples
#'\dontrun{
#' pbdshiny()
#'}
pbdshiny <- function() {
  shinyApp(
    ui = pageWithSidebar(
      headerPanel("Protracted birth death"),
      sidebarPanel(sliderInput("goodrate", "speciation rate of good species:",
          min = 0.00,max = 1.0, value = 0.4,step=0.01,ticks=T),
        sliderInput("compprate","speciation completion rate:",
          min = 0.00,max = 1.0,value = 1.0,step=0.01,ticks=T),
        sliderInput("inciprate","speciation rate of incipient species:",
          min = 0.0,max = 1.0,value = 0.4,step=0.01,ticks=T),
        sliderInput("extgood","extinction rate of good species:",
          min = 0.0,max = 1.0,value = 1.0,step=0.01,ticks=T),
        sliderInput("extincip","extinction rate of incipient species:",
          min = 0.0,max = 1.0,value = 1.0,step=0.01,ticks=T),
        sliderInput("goodrate2","2nd speciation rate of good species:",
          min = 0.0,max = 1.0,value = 0.4,step=0.01,ticks=T),
        sliderInput("compprate2","2nd speciation completion rate:",
          min = 0.0,max = 1.0,value = 1.0,step=0.01,ticks=T),
        sliderInput("inciprate2","2nd speciation rate of incipient species:",
          min = 0.0,max = 1.0,value = 0.4,step=0.01,ticks=T),
        sliderInput("extgood2","2nd extinction rate of good species:",
          min = 0.0,max = 1.0,value = 1.0,step=0.01,ticks=T),
        sliderInput("extincip2","2nd extinction rate of incipient species:",
          min = 0.0,max = 1.0,value = 1.0,step=0.01,ticks=T)
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("On the fly", plotOutput("Plot3"),plotOutput("Plot4")),
          tabPanel("phylograms",plotOutput("Plot12"),plotOutput("Plot13")),
          tabPanel("tau",
            selectInput("tauvariable", "Choose a variable to fix:",
              choices = c("speccomp","incipsp","incipext")),
            sliderInput("taurate","value",min = 0.0,max = 1.0,
              value = 0.4,step=0.01,ticks=T),
            sliderInput("maxcolour","max colour",min=0.0,max=100.0,
              value=8,step=0.01,ticks=T),
            plotOutput("Plot14")
          )
        )
      )
    ),
    server = function(input, output) {
      output$Plot3 <- renderPlot({
        data2<-summaryrepsim(c(input$goodrate,input$compprate,input$inciprate,
          input$extgood,input$extincip),15,15)
        print(plotsim(data2))
        })
      output$Plot4 <- renderPlot({
        data2<-summaryrepsim(c(input$goodrate2,input$compprate2,
          input$inciprate2,input$extgood2,input$extincip2),15,15)
        print(plotsim(data2))
        })
      output$Plot12<- renderPlot({
        x1<-repsim2(c(input$goodrate,input$compprate,input$inciprate,
          input$extgood,input$extincip),1,15)
        x2<-treemaker2(x1)
        plot(read.tree(text=x2))
        })
      output$Plot13<- renderPlot({
        x1<-repsim2(c(input$goodrate2,input$compprate2,input$inciprate2,
          input$extgood2,input$extincip2),1,15)
        x2<-treemaker(x1)
        plot(read.tree(text=x2))
        })
      output$Plot14<-renderPlot({
        plottau(tauloop(input$taurate,input$tauvariable),input$tauvariable,
          input$maxcolour)
        })
    }
  )
}
