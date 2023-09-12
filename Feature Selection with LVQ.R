install.packages("mlbench")
library("mlbench")

install.packages("caret")
library("caret")

install.packages("xlsx")
library(xlsx)

install.packages("purrr")
library(purrr)


#AAC, CT ve PAAC birleþtirince xlsx formatý için exceed memory hatasý verdi o yüzden csv ye çevirdim. csv formatýný virgülle deðil
#noktalý virgül ile ayrýlmýþ olarak dönüþtürdü.
data <- read.csv('C:/Users/zeyne/Desktop/Covid-PPI Proje/Makale1/AAC-CT-PAAC/AAC-CT-PAAA Birlesik.csv',sep=";")


x=data[,2:414]
y=data$label
#train fonksiyonu için bu vektörün factor olarak çevrilmesini istedi yöntem
y <- as.factor(y)

control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(x, y, "lvq", trControl = trainControl(method = "cv"))
# estimate variable importance
importance <- varImp(model, scale=FALSE)

#importance bir list verisi. onu data frame e dönüþtürüyor. print(importance) most significant 20 özelliði listeliyor. 
#Daha fazlasýný listeletmek için data frame e çevirip excel formatýnda yazdýrdým. 
features=map_df(importance, ~as.data.frame(t(.)))

write.xlsx(features, 'C:/Users/zeyne/Desktop/features.xlsx', sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)


print(importance)

