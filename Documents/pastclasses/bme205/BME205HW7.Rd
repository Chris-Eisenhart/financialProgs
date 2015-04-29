mydata <- read.csv("~/Desktop/output")
keysAndValues = subset(mydata, select = c("key","value"))
plot(keysAndValues)
title(main="ORF Occurences")