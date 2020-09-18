library(efdm)
# Test that prior is not dependent on the order of statespace
statespace <- expand.grid(a=1:2, b=1:2, vol=1:15, age=1:35)
statespace <- statespace[sample(nrow(statespace), nrow(statespace)),]

act <- define_activity("test", c("vol", "age"))
act <- build_statespace(act, statespace, by=c("a", "b"))
act1 <- estimatetransprobs(act, NULL, prior_ff())
A <- extract_transitions(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))

act <- define_activity("test", c("vol", "age"))
act <- build_statespace(act, statespace, factors=c("a", "b"))
act1 <- estimatetransprobs(act, NULL, prior_ff())
A <- extract_transitions(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))

act <- define_activity("test", c("vol", "age"))
act <- build_statespace(act, statespace, factors=c("b"), by="a")
act1 <- estimatetransprobs(act, NULL, prior_ff())
A <- extract_transitions(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))
