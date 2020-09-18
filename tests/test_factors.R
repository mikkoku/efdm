library(efdm)

statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
statespace$split <- statespace$ds == "sp" & statespace$region == "n"
pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3), region=c("n", "s"))
pairdata$split <- pairdata$ds == "sp" & pairdata$region == "n"
actprob <- state0 <- statespace
actprob$test <- 1
state0$area <- c(1,1,1,1, 0,0,0,0, 0,0,0,0)

act <- define_activity("test", c("vol"))
act <- build_statespace(act, statespace, factors=c("region"), by=c("split", "ds"))
act1 <- estimatetransprobs(act, pairdata, "nochange")
r1 <- runEFDM(state0, actprob, list(act1), 10)

statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3)
statespace$split <- statespace$ds == "sp" & statespace$region == "n"
pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3), region=c("n", "s"))
pairdata$split <- pairdata$ds == "sp" & pairdata$region == "n"
actprob <- state0 <- statespace
actprob$test <- 1
state0$area <- c(1,1,1,1, 0,0,0,0, 0,0,0,0)

act <- define_activity("test", c("vol"))
act <- build_statespace(act, statespace, factors=c("region"), by=c("split", "ds"))
act1 <- estimatetransprobs(act, pairdata, "nochange")
r2 <- runEFDM(state0, actprob, list(act1), 10)
r1$ds <- factor(r1$ds, levels=c("sp", "pi"))
r1$region <- factor(r1$region, levels=c("n", "s"))
sortdf <- function(x) x[do.call(order, x),]
stopifnot(all(sortdf(r1)==sortdf(r2)))

