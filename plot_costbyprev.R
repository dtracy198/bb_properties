# Summarizing results from costbyprev.rds
library(ggplot2)

setwd("/Users/sherriexie/Dropbox/Levy_Lab/bb_properties")
data <- readRDS("costbyprev_mig10yr_test.rds")

d1 <- data[[1]]
d2 <- data[[2]]
d3 <- data[[3]]
d4 <- data[[4]]
d5 <- data[[5]]

prev <- cost <- redprev <- numeric()
prev <- c(rep(2, dim(data[[1]])[1]), rep(4, dim(data[[2]])[1]), 
          rep(6, dim(data[[3]])[1]), rep(8, dim(data[[4]])[1]),
          rep(10, dim(data[[5]])[1]))
# remember to multiply costs by -1 for _mig2yr.rds & _mig5yr
cost <- c(data[[1]]$tot.cost, data[[2]]$tot.cost, data[[3]]$tot.cost,
          data[[4]]$tot.cost, data[[5]]$tot.cost)
costbyprev <- data.frame(prev=prev, cost=cost)

sstat <- data.frame(matrix(rep(0,7*5),nrow=5,ncol=7))
names(sstat) <- c("prev", "min", "q1", "med", "mean", "q3", "max")
sstat$prev <- seq(2,10,2)
sstat[1,2:7] <- summary(costbyprev$cost[which(costbyprev$prev==2)])
sstat[2,2:7] <- summary(costbyprev$cost[which(costbyprev$prev==4)])
sstat[3,2:7] <- summary(costbyprev$cost[which(costbyprev$prev==6)])
sstat[4,2:7] <- summary(costbyprev$cost[which(costbyprev$prev==8)])
sstat[5,2:7] <- summary(costbyprev$cost[which(costbyprev$prev==10)])

redprev <- c(data[[1]]$prev*100, data[[2]]$prev*100, data[[3]]$prev*100,
             data[[4]]$prev*100, data[[5]]$prev*100)
prevbyprev <- data.frame(prev=prev, redprev=redprev)

sprev <- data.frame(matrix(rep(0,7*5),nrow=5,ncol=7))
names(sprev) <- c("prev", "min", "q1", "med", "mean", "q3", "max")
sprev$prev <- seq(2,10,2)
sprev[1,2:7] <- summary(prevbyprev$redprev[which(prevbyprev$prev==2)])
sprev[2,2:7] <- summary(prevbyprev$redprev[which(prevbyprev$prev==4)])
sprev[3,2:7] <- summary(prevbyprev$redprev[which(prevbyprev$prev==6)])
sprev[4,2:7] <- summary(prevbyprev$redprev[which(prevbyprev$prev==8)])
sprev[5,2:7] <- summary(prevbyprev$redprev[which(prevbyprev$prev==10)])

ggplot(sstat, aes(x=prev, y=med)) + 
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.1) +
  geom_line() +
  geom_point() +
  xlab("Prevalence (%)") +
  ylab("Total Cost ($/unit/year)") +
  theme(axis.text.x=element_text(margin = margin(b=10), size=12),
        axis.text.y=element_text(size=12, margin=margin(l=10)),
        axis.title=element_text(face="bold", size =12))

ggplot(sprev, aes(x=prev, y=med)) + 
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.1) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(0,10,2), name = "Baseline Prevalence (%)",
                     limits = c(1.5,10)) +
  scale_y_continuous(breaks=seq(0,10,2), name = "Prevalence with Disclosure (%)",
                     limits = c(0,10)) +
  theme(axis.text.x=element_text(margin = margin(b=10), size=12),
        axis.text.y=element_text(size=12, margin=margin(l=10)),
        axis.title=element_text(face="bold", size =12)) +
  geom_segment(aes(x = 2, y = 2, xend = 10, yend = 10), linetype=2, colour="red")

