
data <- readRDS("costbyprev_mig10yr.rds")


d1 <- data[[1]]
d2 <- data[[2]]
d3 <- data[[3]]
d4 <- data[[4]]
d5 <- data[[5]]

par(mfrow=c(2,3),oma=c(0,0,2,0))
plot(d1$bet, d1$tot.cost, main = "2% baseline prev")
plot(d2$bet, d2$tot.cost, main = "4% baseline prev")
plot(d3$bet, d3$tot.cost, main = "6% baseline prev")
plot(d4$bet, d4$tot.cost, main = "8% baseline prev")
plot(d5$bet, d5$tot.cost, main = "10% baseline prev")
title("beta vs cost, low tov + partial dis", outer=TRUE, cex=1.5)

par(mfrow=c(2,3))
plot(d1$bet, d1$prev, main = "2% baseline prev")
plot(d2$bet, d2$prev, main = "4% baseline prev")
plot(d3$bet, d3$prev, main = "6% baseline prev")
plot(d4$bet, d4$prev, main = "8% baseline prev")
plot(d5$bet, d5$prev, main = "10% baseline prev")
title("beta vs new prev, low tov + partial dis", outer=TRUE)
