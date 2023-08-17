data.table::setDTthreads(1)
library(testthat)
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

act1 <- define_activity("test", c("vol"))
transprobs(act1) <- estimatetransprobs("vol", NULL, subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"), prior="nochange")
act2 <- define_activity("test", c("vol"))
transprobs(act2) <- estimatetransprobs("vol", NULL, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"), prior="nochange")
orderresult(runEFDM(state0, actprob, list(act1, act2), 1))

act1 <- define_activity("test", c("vol"))
transprobs(act1) <- estimatetransprobs("vol", subset(pairdata, ds=="sp"&region=="n"), subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"), prior="nochange")
act2 <- define_activity("test", c("vol"))
transprobs(act2) <- estimatetransprobs("vol", subset(pairdata, !(ds=="sp"&region=="n")), subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"), prior="nochange")
orderresult(runEFDM(state0, actprob, list(act1, act2), 10))


statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
act <- define_activity("test", c("vol"))
statespace1 <- subset(statespace, !(ds=="sp"&region=="n"))
act <- efdm:::build_complex_statespace(act, statespace1, statespace1, factors=c("ds", "region"))
stopifnot(nrow(act$statespace[[1]])==3)

statespace <- expand.grid(a=1:2, b=c("f", "g"), c=c(4,9), vol=1:3, stringsAsFactors=FALSE)
ii <- sample(nrow(statespace), sample(nrow(statespace), 1))
act <- define_activity("test", c("vol"))
act <- efdm:::build_complex_statespace(act, statespace[ii,], statespace[ii,], factors=c("a", "b"), by="c")
expect_identical(sum(sapply(act$statespace, function(x) nrow(x$statespace0))), length(ii))
