protractedbirthdeath
====================

modelling speciation with incipient species

This shiny app uses a modification of the PBD package by Etienne and Rosindell (2014) allowing for 0s as parameters, and extinction.
To run, simply download the exampledata folder, install and load shiny (install.packages(shiny), library(shiny))
and run runGitHub("protractedbirthdeath","jeremycg",subdir="shiny")
The r function choose.dir() is currently windows only - you can manually modify the source to your directory.

Output is in multiple tabs - first tab reads from the example data and plots previous data,
second does the same but on the fly - some parameter choices might take a long time.
Third tab gives boxplots of number of species at t=0 for differing single parameters.
Fourth draws two trees based on the inputs from the sliders
Fifth gives contour plots of number of species, based on three fixed variables, and two to plot. This won't work with example data
Sixth is tau - mean time to speciation again based on contour plots.
