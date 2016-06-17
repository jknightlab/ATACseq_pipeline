# "temp" should contain three columns -- parameters, FScore and fraction_on_target
# 
# parameters	F0.5	fraction_on_target
# l_10.t_0.001	0.0341355849	0.9275687549
# l_10.t_0.01	0.0378076052	0.9225974445
# l_10.t_0.1	0.097440015	0.9145559285
# l_10.t_0.99	0.28203335	0.8233712572
# l_100.t_0.001	0.1913228994	0.9052032701
# l_100.t_0.01	0.2444097117	0.8908444699


data <- read.table("temp", sep='\t', header=TRUE)

data1 <- read.table("temp1", sep='\t', header=TRUE)


plot (data$fraction_on_target, pch=19, cex=2, col="blue", ylim=c(0,1))
lines (data$fraction_on_target, col="blue")
points (data$F0.5, pch=19, cex=2)
lines (data$F0.5)

# legend(28.1, 0.92, c("Fraction on target", "F-Score"), lty=c(1,1), lwd=c(3,3), col=c("blue","black"), cex=2)


plot (data1$fraction_on_target, pch=19, cex=2, col="blue", ylim=c(0,1))
lines (data1$fraction_on_target, col="blue")
points (data1$F0.5, pch=19, cex=2)
lines (data1$F0.5)

# legend(49, 0.92, c("Fraction on target", "F-Score"), lty=c(1,1), lwd=c(3,3), col=c("blue","black"), cex=2)

