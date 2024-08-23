library(reshape2)

## loading data
onema_community <- read.csv("inst/onema/community_data.csv")
env             <- read.csv("inst/onema/environment.csv")
protocol        <- read.csv("inst/onema/fishing_protocol.csv") # link between opcods & stations

onema_community$species <- factor(onema_community$species,
                                  levels = unique(onema_community$species))
onema_community$opcod   <- factor(onema_community$opcod,
                                  levels = unique(onema_community$opcod))

# filter data for stations with environmental info
protocol        <- protocol[protocol$station %in% env$station,]
onema_community <- onema_community[onema_community$opcod %in% protocol$opcod,] # does not contain NA
species_order <- unique(onema_community$species)
opcod_order   <- unique(onema_community$opcod)
stations        <- protocol[match(opcod_order, protocol$opcod), , drop=FALSE]$station
idx             <- match(stations, env$station)
stations_env    <- env[idx, ,drop = FALSE] # does not contain NA

# select and format data to run normalblockr on it
X    <- stations_env[ , 2:ncol(stations_env)]
Y    <- apply(as.matrix(dcast(onema_community, species ~ opcod, value.var = "biomass", fill = 0))[,-1], 1, as.numeric)
Y    <- log(1 + Y)
Y500 <- Y[order(rowSums(Y == 0))[1:500],]
X500 <- as.matrix(X[order(rowSums(Y == 0))[1:500],])
rownames(X500) <- NULL
colnames(X500) <- NULL
Y500 <- Y500[, - which(colSums(Y500 == 0) >= 400)]
n = nrow(Y500)
Xones <- matrix(rep(1, n), nrow = n)

## ZI inflated normal with diagonal covariance
zi_normal <- normalblockr::normal_zi$new(Y500, Xones)
zi_normal$optimize()

myModel <- normalblockr::NB_unknown_zi$new(Y500, Xones, 2:10)
myModel$optimize()
