library(efdm)
statespace <- expand.grid(a=1:2, b=1:2, vol=1:5)
pairdata <- data.frame(a=c(1,1,2,2), b=c(1,2,1,2), vol0=c(1,1,1,1), vol1=c(2,3,4,5))
state0 <- statespace
actprob <- statespace
actprob$test <- 1
state0$area <- 0
state0$area[1] <- 1

act1 <- define_activity("test", c("vol"))
transprobs(act1) <- estimatetransprobs("vol", pairdata, statespace, by=c("a", "b"), prior="nochange")
transprobs(act1)
runEFDM(state0, actprob, list(act1), 1)
act2 <- define_activity("test", c("vol"))
transprobs(act2) <- estimatetransprobs("vol", pairdata, statespace, factors="a", by="b", prior="nochange")
transprobs(act2)
runEFDM(state0, actprob, list(act2), 1)
act3 <- define_activity("test", c("vol"))
transprobs(act3) <- estimatetransprobs("vol", pairdata, statespace, factors="b", by="a", prior="nochange")
transprobs(act3)
runEFDM(state0, actprob, list(act3), 1)
act4 <- define_activity("test", c("vol"))
transprobs(act4) <- estimatetransprobs("vol", pairdata, statespace, factors=c("a", "b"), prior="nochange")
transprobs(act4)
runEFDM(state0, actprob, list(act4), 1)
act5 <- define_activity("test", c("vol"))
transprobs(act5) <- estimatetransprobs("vol", pairdata, statespace, factors=c("b", "a"), prior="nochange")
transprobs(act5)
runEFDM(state0, actprob, list(act5), 1)

act6 <- define_activity("test", c("vol"))
transprobs(act6) <- estimatetransprobs("vol", pairdata[c("vol0", "vol1")], data.frame(vol=1:5), prior="nochange")
transprobs(act6)
runEFDM(state0, actprob, list(act6), 1)

try(estimatetransprobs("vol", pairdata[c("vol0", "vol1")], data.frame(vol=1:4), prior="nochange"))
