install.packages("mlbench")
library("mlbench")

install.packages("caret")
library("caret")

install.packages("xlsx")
library(xlsx)

install.packages("purrr")
library(purrr)


data <- read.csv('/AAC-CT-PAAA Birlesik.csv',sep=";")


x=data[,2:414]
y=data$label
#for train function vectors should be represent as factor. 
y <- as.factor(y)

control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(x, y, "lvq", trControl = trainControl(method = "cv"))
# estimate variable importance
importance <- varImp(model, scale=FALSE)

#'importance' is a list. 'print(importance)' function lists the most significant 20 features
features=map_df(importance, ~as.data.frame(t(.)))

write.xlsx(features, '/features.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)


print(importance)

