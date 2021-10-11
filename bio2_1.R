#
# bioinfomatices II 
# Day2
#
# 2021/10/04 Shiho AOKI

# install packages.
# install.packages("tidyverse");
# install.packages("ggfortify");

library(tidyverse);
library(ggfortify);

# load Indometh dataset
Data <- Indometh;
summary(Data);

# create histogram about conc
v = conlnames(Data);
hist(as.matrix(Data[, 3], main=v[3]));

# create plot
plot(Data)

# ################# Analysis###########
library(psych);
psych::pairs.panels(Data);
conAndTime <- lm(conc~time, Data);
autoplot(conAndTime);

plot(prcomp(Data[,-1]));
P <- prcomp(Data[,-1]);
plot(P$x[, 1:2], col=Data[, -1]);
