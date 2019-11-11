library(tm)
library(SnowballC)
library(e1071)
Data<-read.csv("~/MyData.csv",stringsAsFactors = FALSE)
str(Data)
Data<-Data[,3:4]
Data$therme_key <- factor(Data$therme_key)
str(Data)
table<-table(Data$therme_key)


DataCorpus  <- Corpus(VectorSource(Data$therme_value))

#on prefere garder les num la ponctuation et tout ce qui caract?rise les donn?s ayant des classes avec G# 
#le reste des mots qui n'ont pas de tags ont ?t? deja trait? dans process.R 

dtm <- DocumentTermMatrix(DataCorpus)

dtm_train <- dtm[1:134053,]
dtm_test <- dtm[134054:223423,]

train_labels <- Data[1:134053,]$therme_key
test_labels <- Data[134054:223423,]$therme_key

frequent_terms <- findFreqTerms(dtm_train,5)

dtm_freq_train <- dtm_train[,frequent_terms]
dtm_freq_test <- dtm_test[,frequent_terms]

convert_counts <- function(x){
  x <- ifelse(x > 0,"Yes","No")
}

train <- apply(dtm_freq_train,MARGIN = 2,
                   convert_counts)
test <- apply(dtm_freq_test,MARGIN = 2,
                  convert_counts)

words_classifier <- naiveBayes(dtm_freq_train,train_labels)
words_test_pred <- predict(words_classifier,dtm_freq_test)



library(gmodels)
CrossTable(words_test_pred,test_labels,
           prop.chisq = FALSE,
           prop.t = FALSE,
           dnn = c("predicted","actual"))
