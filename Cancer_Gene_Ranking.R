library(Biobase)
library(GEOquery)
library(igraph)
library(e1071)
library(caret)
library(ggplot2)

df_ggplt <- function(df){
  dt <- data.frame(GENE=character(),EXPR=double(),Condition=factor())
  for(i in 1:nrow(df)){
    exprs <- unlist(df[i,])
    names(exprs) <- NULL
    dtt <- data.frame(GENE=rep(rownames(df)[i]),EXPR=exprs,Condition=as.factor(c(rep('Cancer',60),rep('Normal',60))))
    dt <- rbind(dt,dtt)
  }
  ggplot(data = dt, aes(x=GENE, y=EXPR)) + geom_boxplot(aes(fill=Condition)) + xlab("Genes") +
    ylab("Expression")
}

#gds <- getGEO(filename = "GDS3837.soft.gz")
#Meta(gds)$sample_count
#tbl_gds <- Table(gds)[,]
tbs <- read.csv("GDS3837.csv",header = TRUE)
size <- dim(tbs)
tbs <- cbind(tbs,AVGE = apply(tbs[3:size[2]],1,mean))
tbs <- tbs[tbs$AVGE>=6,]
sd_flag <- apply(tbs[3:size[2]],1,sd) >= 0.5
tbs <- tbs[sd_flag,]
tbs <- tbs[order(tbs$IDENTIFIER,-tbs$AVGE),]
tbs <- tbs[!duplicated(tbs$IDENTIFIER),]
GEP <- tbs[,3:size[2]]
rownames(GEP) <- tbs$IDENTIFIER
pin <- read.delim("HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt",sep = "\t",header = FALSE)

indx <- c()
for(i in 1:dim(pin)[1]){
  print(i)
  symbols <- unlist(pin[i,c(1,4)])
  if((sum(rownames(GEP)==symbols[1]) + sum(rownames(GEP)==symbols[2])) == 2){
    indx <- c(indx,i)
  }
}

pEP <- pin[indx,c(1,4)]
colnames(pEP) <- c("Source_A","Source_B")

common_genes <- intersect(rownames(GEP),union(pEP$Source_A,pEP$Source_B))
cgep <- GEP[common_genes,]

write.csv(pEP,"RPIN.csv")
write.csv(cgep,"RGEP.csv")

system('python Rank.py')

rgenes <- read.csv("RankedGenes.csv",header = TRUE)
rgenes <- rgenes[order(-rgenes$Rank),]

cgep <- cgep[rgenes$Gene,]
df_ggplt(cgep[1:10,])

#target <- c(rep(0,40),rep(1,13))
target <- c(rep(1,60),rep(0,60))

set.seed(20)
train <- c(sample(1:60,45),sample(61:120,45))
test <- -train

train_data <- cgep[,train]
train_data <- t(train_data)
train_data <- cbind(train_data,class = target[train])
test_data <- cgep[,test]
test_data <- t(test_data)


train_data_reduced <- as.matrix(train_data[,c(1:10,4278)])
test_data_reduced <- test_data[,1:10]

grid <- expand.grid(C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_model <- train(as.factor(class) ~ .,data = train_data_reduced, method = "svmLinear",
                 trControl=trctrl,tuneGrid = grid_radial,tuneLength = 10)
plot(svm_model)

X <- 1:50
Y <- c()
F1 <- c()
for(i in 1:50){
  if(i == 1){
    train_data_reduced <- as.matrix(train_data[,c(i,4278)],ncol=2)
    test_data_reduced <- as.matrix(test_data[,i],ncol=1)
    colnames(test_data_reduced) <- colnames(train_data_reduced)[1]
  }
  else{
    train_data_reduced <- as.matrix(train_data[,c(1:i,4278)])
    test_data_reduced <- test_data[,1:i]
  }
  svm_model <- svm(as.factor(class) ~ .,data = train_data_reduced, method = "svmLinear",cost = 0.5)

  test_pred <- predict(svm_model, newdata = test_data_reduced)
  cm <- confusionMatrix(test_pred,as.factor(target[test]))
  Y <- c(Y,cm$overall['Accuracy'])
  F1 <- c(F1,cm$byClass['F1'])
}

names(Y) <- NULL
names(F1) <- NULL

ggplot(data.frame(Features = 1:50,Accr = Y, F1S = F1), aes(Features)) + 
  geom_line(aes(y = Y, colour = "Accuracy")) + 
  geom_line(aes(y = F1, colour = "F1-Score")) + ylab("Performance")

