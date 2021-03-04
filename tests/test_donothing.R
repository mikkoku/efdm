library(efdm)
library(testthat)

statespace <- expand.grid(a=1:2, b=1:2, vol=1:5)
state0 <- statespace
actprob <- statespace
actprob$test <- 1
state0$area <- 0
state0$area[1] <- 1


act1 <- define_activity("test", character())
expect_error(transprobs(act1) <- data.frame(vol0 = 1:4, vol1=2:5, prob=1))
state1 <- runEFDM(state0, actprob, list(act1), 5)
expect_equal(state1$vol, rep(1, 6))

act2 <- define_activity("test2", "vol")
transprobs(act2) <- data.frame(vol0=1:5, vol1=c(2:5,5), prob=1)
actprob$test <- 0.1
actprob$test2 <- 0.9
expect_known_hash(runEFDM(state0, actprob, list(act1, act2), 5), "1b23f7f782")

