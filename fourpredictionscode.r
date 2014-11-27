# load a bunch of stuff
libraries <- c("DDD", "PBD", "plyr", "dplyr", "ggplot2", "RColorBrewer",
    "compiler", "shiny", "data.table")
for (i in libraries) {
    if (i %in% rownames(installed.packages()) == FALSE) {
        install.packages(i)
    }
}
lapply(libraries, require, character.only = T)
# terrible way to load packages...
library(protractedbirthdeath)

#sisterspecies lengths
dlply(repsim2(c(0.1, 0.1, 0.1, 0.1, 0.1), 20), .(run), .fun = sisterlengths)
#persistance through time
ddply(repsim2(c(0.1, 0.1, 0.1, 0.1, 0.1), 20), .(run), .fun = countpersistance)
#timegivengood
dlply(repsim2(c(0.1, 0.1, 0.1, 0.1, 0.1), 20), .(run), .fun = timegivengood)
#numberof species in different classes
ddply(repsim2(c(0.1, 0.1, 0.1, 0.1, 0.1), 20), .(run), .fun = numgoodincip)
