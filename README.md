#protractedbirthdeath
====================

##modelling speciation with incipient species

This shiny app uses a modification of the PBD package by Etienne and Rosindell (2014) allowing for 0s as parameters, and extinction. It also treats "orphaned" incipient species in a more intelligent way.

The main core of the project is now a package. To install, run:

```R
# install.packages("devtools")
devtools::install_github("jeremycg/protractedbirthdeath/protractedbirthdeath")
```

##Examples

The main part of the package lies in the pbdsim2 function.
Output is a list, parse it (and run multiple repeats) using repsim2

```R
library(protractedbirthdeath)
run<-pbdsim2(c(0.7,0.7,0.7,0.1,0.1),15)
runs<-repsim2(c(0.7,0.7,0.7,0.1,0.1),5,15)
```
pars[1] is good speciation rate
pars[2] is speciation completion rate
pars[3] is incipient speciation rate
pars[4] is good extinction rate
pars[5] is incipient extinction rate

the following numbers are repeats and time (defaulted to 15 million years)
pbdsim2 has a limit on numbers of species at 100000. this can be altered in the code.

summaryrepsim will take means and sem of repeats.
Plotsim will make a plot of these summaries.

```R
summary1<-summaryrepsim(c(0.7,0.7,0.7,0.1,0.1),5,15)
plotsim(summary1)
```

##Shiny

The package began life as a shiny package, and this is still extant.

To run, simply download the exampledata folder, install and load shiny (install.packages(shiny), library(shiny))
and run runGitHub("protractedbirthdeath","jeremycg",subdir="shiny")
The r function choose.dir() is currently windows only - you can manually modify the source to your directory.

Output is in multiple tabs - first tab reads from the example data and plots previous data,
second does the same but on the fly - some parameter choices might take a long time.
Third tab gives boxplots of number of species at t=0 for differing single parameters.
Fourth draws two trees based on the inputs from the sliders
Fifth gives contour plots of number of species, based on three fixed variables, and two to plot. This won't work with example data
Sixth is tau - mean time to speciation again based on contour plots.
