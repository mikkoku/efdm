library(efdm)
library(testthat)
statespace <- expand.grid(a=1:2, b=1:2, vol=1:5)
pairdata <- data.frame(a=c(1,1,2,2), b=c(1,2,1,2), vol0=c(1,1,1,1), vol1=c(2,3,4,5))
state0 <- statespace
actprob <- statespace
actprob$test <- 1
state0$area <- 0
state0$area[1] <- 1

act <- define_activity("test", c("vol"))
act <- build_statespace(act, statespace, by=c("a", "b"))
act1 <- estimatetransprobs(act, pairdata, "nochange")
A <- transprobs(act1)
res0 <- runEFDM(state0, actprob, list(act1), 1)

transprobs(act1) <- A
res1 <- runEFDM(state0, actprob, list(act1), 1)
expect_identical(res0, res1)

act2 <- define_activity("test", c("vol"))
transprobs(act2) <- A
res2 <- runEFDM(state0, actprob, list(act2), 1)
expect_identical(res0, res2)

A$asdf <- 1
transprobs(act1) <- A
expect_error(runEFDM(state0, actprob, list(act1), 1, check=FALSE),
             "Variable 'asdf' required by activity 'test' not present in state.")
expect_error(runEFDM(state0, actprob, list(act1), 1),
             "Activity 'test' has variable asdf not in actprob.")


act3 <- define_activity("test", c("vol"))
transprobs(act3) <- data.frame(vol0 = 1:5, vol1=c(2:5, 5), prob=1)
res3 <- runEFDM(state0, actprob, list(act3), 3)
# Test something?
