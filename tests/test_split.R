library(efdm)
statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3), region=c("n", "s"))
actprob <- state0 <- statespace
actprob$test <- 1
state0$area <- c(1,1,1,1, 0,0,0,0, 0,0,0,0)
orderresult <- function(res) {
  res <- res[do.call(order, res[c("time", "region", "ds", "vol")]),c("time", "region", "ds", "vol", "activity", "area")]
  row.names(res) <- NULL
  res
}

act <- define_activity("test", c("vol"))
act <- build_statespace(act, subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"))
act1 <- estimatetransprobs(act, NULL, "nochange")
act <- define_activity("test", c("vol"))
act <- build_statespace(act, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"))
act2 <- estimatetransprobs(act, NULL, "nochange")
orderresult(runEFDM(state0, actprob, list(act1, act2), 1))

act <- define_activity("test", c("vol"))
act <- build_statespace(act, subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"))
act1 <- estimatetransprobs(act, subset(pairdata, ds=="sp"&region=="n"), "nochange")
act <- define_activity("test", c("vol"))
act <- build_statespace(act, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"))
act2 <- estimatetransprobs(act, subset(pairdata, !(ds=="sp"&region=="n")), "nochange")
orderresult(runEFDM(state0, actprob, list(act1, act2), 10))
