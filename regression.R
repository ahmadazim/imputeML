data = read.csv("/Users/ahmadazim/Documents/Research/imputeML/gene_21053_100.csv",header = TRUE)
data[1:6,1:6]
data[,1] = NULL

summary(lm(y21053 ~ ., data = data))

data2 = read.csv("/Users/ahmadazim/Documents/Research/imputeML/gene_30933_100.csv",header = TRUE)
data2[1:6, 1:6]
data2[,1] = NULL

summary(lm(y30933 ~ ., data = data2))

data[,102] = 1
data2[,102] = 2

allData = rbind(data, data2)
