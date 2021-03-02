library(efdm)
library(testthat)
statespace <- expand.grid(a=1:2, b=1:2, vol=1:5)
state0 <- statespace
actprob <- statespace
actprob$test <- 1
state0$area <- 0
state0$area[1] <- 1


act3 <- define_activity("test", c("vol"))
transprobs(act3) <- data.frame(vol0 = 1:4, vol1=2:5, prob=1)
expect_error(runEFDM(state0, actprob, list(act3), 5), "States with no activity")
expect_error(runEFDM(state0, data.frame(test=1), list(act3), 5), "No state variables in actprob")



transprobs(act3) <- data.frame(vol0 = 1:5, vol1=c(2:6), prob=1)
expect_error(runEFDM(state0, actprob, list(act3), 3), "Activity transitions to state not in actprob")


expect_error(runEFDM(state0, actprob, list(act3, act3), 3), "States with multiple activities")


act4 <- define_activity("test2", c("vol"))
transprobs(act4) <- data.frame(vol0 = 1:5, vol1=c(2:5,5), prob=1)
expect_error(runEFDM(state0, actprob, list(act4), 5), "Activity probabilities for test2 not given.")

act5 <- define_activity("test", c("vol"))
transprobs(act5) <- data.frame(vol0 = 1:6, vol1=2, prob=1)
expect_error(runEFDM(state0, actprob, list(act5), 5), "Activity transitions from state not in actprob")


actprob$test[3] <- NA
expect_error(runEFDM(state0, actprob, list(act3), 5), "Actprob should not have NAs.")

actprob$test2 <- actprob$test <- 0.5
expect_error(runEFDM(state0, actprob, list(act3), 5), "Activity probabilities for test2 given, but activities not.")

