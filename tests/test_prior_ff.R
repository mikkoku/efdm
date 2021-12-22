library(efdm)
# Test that prior is not dependent on the order of statespace
statespace <- expand.grid(a=1:2, b=1:2, vol=1:15, age=1:35)
statespace <- statespace[sample(nrow(statespace), nrow(statespace)),]

act <- define_activity("test", c("vol", "age"))
act1 <- estimatetransprobs(c("vol", "age"), NULL, statespace, by=c("a", "b"), prior=prior_ff())
A <- transprobs(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))

act <- define_activity("test", c("vol", "age"))
act1 <- estimatetransprobs(c("vol", "age"), NULL, statespace, factors=c("a", "b"), prior=prior_ff())
A <- transprobs(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))

act <- define_activity("test", c("vol", "age"))
act1 <- estimatetransprobs(c("vol", "age"), NULL, statespace, factors=c("b"), by="a", prior=prior_ff())
A <- transprobs(act1)
stopifnot(all(A$vol1==1 & A$age1 == 1))
