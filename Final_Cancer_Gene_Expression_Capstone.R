############################################################################
############################################################################
###                                                                      ###
###                  SECTION 0: LOAD REQUIRED LIBRARIES                  ###
###                                                                      ###
############################################################################
############################################################################

if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(nnet)) install.packages("nnet")
if(!require(knitr)) install.packages("knitr")
if(!require(rpart)) install.packages("rpart")
if(!require(rpart.plot)) install.packages("rpart.plot")

library(tidyverse)
library(caret)
library(data.table)
library(lubridate)
library(bannerCommenter)
library(nnet)
library(knitr)
library(rpart)
library(rpart.plot)

############################################################################
############################################################################
###                                                                      ###
###                        SECTION 1: CREATE DATA                        ###
###                                                                      ###
############################################################################
############################################################################

# Download data from github
ytrain <- scan("https://raw.github.com/skylershapiro/Cancer_Gene_Expression/master/14cancer.ytrain.txt")
xtrain <- read.table("https://raw.github.com/skylershapiro/Cancer_Gene_Expression/master/14cancer.xtrain.txt", header=FALSE, sep="")

ytest <- scan("https://raw.github.com/skylershapiro/Cancer_Gene_Expression/master/14cancer.ytest.txt")
xtest <- read.table("https://raw.github.com/skylershapiro/Cancer_Gene_Expression/master/14cancer.xtest.txt", header=FALSE, sep="")

# Create cancer gene expression dataset from pre-made test and train sets to create my own random partition
cancer_gene_expression <- as.data.frame(t(cbind(xtrain, xtest)))
rownames(cancer_gene_expression) <- 1:nrow(cancer_gene_expression)
cancer_type_nums <- c(ytrain, ytest)

# Convert cancer names from number labels to text
cancer_names <- c("breast", "prostate", "lung", "colorectal", "lymphoma", "bladder", "melanoma", "uterus", "leukemia", "renal", "pancreas", "ovary", "meso", "cns")
cancer_types <- factor(cancer_type_nums, levels = c(1:14), labels = cancer_names)

###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 2: DATA EXPLORATION                     ###
###                                                                     ###
###########################################################################
###########################################################################

# Check dimensions of train set
dim(cancer_gene_expression)

# Check for missing values
sum(is.na(cancer_gene_expression))

# Heatmap 
#Takes long time to run!!!!
heatmap3(as.matrix(cancer_gene_expression[,1:50]))

# Table of cancer type occurences
sort(summary(cancer_types),decreasing = TRUE) %>% kable()




############################################################################
############################################################################
###                                                                      ###
###                SECTION 3: PCA FOR DIMENSION REDUCTION                ###
###                                                                      ###
############################################################################
############################################################################

# Calculate principal components
pc <- prcomp(cancer_gene_expression)
X <- pc$x

# Create partition 
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
test_index <- createDataPartition(y = cancer_types, times = 1, p = 0.2, list = FALSE)

# Partition matrix of principal components of gene expression data
pc_test <- X[test_index,]
pc_train <- X[-test_index,]

# Partition cancer type labels
ctype_test <- cancer_types[test_index]
ctype_train <- cancer_types[-test_index]

##################################################################
##                     PCA Data Exploration                     ##
##################################################################

# Table of distribution of cancer types in test and train sets
sort(summary(ctype_test), decreasing = TRUE) %>% kable()
sort(summary(ctype_train), decreasing = TRUE) %>% kable()

# Plot of first 2 principal components
data.frame(pc_1 = X[,1], pc_2 = X[,2], cancer = cancer_types) %>%
  ggplot(aes(pc_1, pc_2, col = cancer)) +
  geom_point() 

# Cumulative percent of variation explained 
var_explained <- cumsum(pc$sdev^2/sum(pc$sdev^2))
plot(var_explained, ylim = c(0,1),main = "Cumulative Proportion of Variation Explained",xlab = "PC Number", ylab = "Variation Explained")
abline(h=var_explained[37],col='red',v=37)

# We will use 30 principal components because it is the lowest number of components that explains 90% of the variability
var_explained[37]

# Table of cumulative percent of variation explained (firgure this out)
var_explained_labels <- paste("First",paste(1:length(var_explained), "Principal Components"))
var_pcs <- data.frame(number_pcs = var_explained_labels, variance_explained = var_explained, stringsAsFactors = FALSE)
var_pcs$number_pcs[1] <- "First Principle Component"

# Use top ... amount of pcs for models
# pc1 pc2 faceted plot by cancer type
data.frame(pc_1 = pc$x[,1], pc_2 = pc$x[,2], cancer = cancer_types) %>%
  ggplot(aes(pc_1, pc_2,col=cancer)) +
  geom_point() +
  stat_ellipse() +
  facet_grid(cols=vars(cancer)) +
  theme(axis.text.x = element_text(angle = 40, size = 5), strip.text = element_text(size=7,face = "bold")) 

###########################################################################
###########################################################################
###                                                                     ###
###                    SECTION 4: MODEL CONSTRUCTION                    ###
###                                                                     ###
###########################################################################
###########################################################################

#################################################################
##            Data Wrangling for Model Construction            ##
#################################################################

# Wrangle data with 37 principal components for modeling, attach cancer type labels to training set
newdat_train <- cbind(pc_train[,1:37], ctype_train)
temp <- as.data.frame(newdat_train)
newdat_train <- temp
newdat_train$ctype_train <- ctype_train
# Wrangle data with 37 principal components for modeling, attach cancer type labels to testing set
newdat_test <- cbind(pc_test[,1:37], ctype_test)
temp <- as.data.frame(newdat_test)
newdat_test <- temp
newdat_test$ctype_test <- ctype_test

# Wrangle data using only 15 principal components for training set
fift_newdat_train <- cbind(pc_train[,1:15], ctype_train)
temp <- as.data.frame(fift_newdat_train)
fift_newdat_train <- temp
fift_newdat_train$ctype_train <- ctype_train
# Wrangle data using only 15 components for testing set
fift_newdat_test <- cbind(pc_test[,1:15], ctype_test)
temp <- as.data.frame(fift_newdat_test)
fift_newdat_test <- temp
fift_newdat_test$ctype_test <- ctype_test

##################################################################
##             Model 1: K - Nearest Neighbors (KNN)             ##
##################################################################
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
# Build KNN model
control <- trainControl(method = "cv", number = 10, p = .9)
train_knn <- train(ctype_train ~ ., method = "knn", 
                      data = newdat_train,
                      tuneGrid = data.frame(k = seq(1, 50, 2)),
                      trControl = control)

# Plot k values
ggplot(train_knn,highlight = TRUE) + 
  ggtitle( "K-values versus Model Accuracy")

# Choose optimal k
k = train_knn$bestTune
k

# Accuracy on train set using optimal k-value
knn_train_acc <- max(train_knn$results$Accuracy)

# Predict on test set
y_hat_knn <- predict(train_knn, newdat_test, type = "raw")

# Check final accuracy
knn_test_acc <- confusionMatrix(data = y_hat_knn, reference = ctype_test)$overall["Accuracy"]

# Add to table of train and test accuracy
tab_knn <- data.frame(Method="KNN", train_accuracy=knn_train_acc, test_accuracy=knn_test_acc) 
tab_knn %>% kable()

##****************************************************************
##               KNN with 15 principal components               **
##****************************************************************
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
# Build model
control <- trainControl(method = "cv", number = 10, p = .9)
train_knn <- train(ctype_train ~ ., method = "knn", 
                   data = fift_newdat_train,
                   tuneGrid = data.frame(k = seq(1, 50, 2)),
                   trControl = control)

# Predict on test set
y_hat_knn <- predict(train_knn, fift_newdat_test, type = "raw")

# Check final accuracy
fift_knn_test_acc <- confusionMatrix(data = y_hat_knn, reference = ctype_test)$overall["Accuracy"]

# Add results to table
model_results <- data_frame(Method = "1: K - Nearest Neighbors (KNN)", 
                            Accuracy_37_principal_components = knn_test_acc,
                            Accuracy_15_principal_components = fift_knn_test_acc)



#################################################################
##                 Model 2: Random Forest (RF)                 ##
#################################################################

# Use basic rf model to generate decision tree plot
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
fit <- rpart(ctype_train ~ ., data = newdat_train)
rpart.plot(fit,fallen.leaves = FALSE)

# Build rf model using Rborist package
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3)
train_rf <- train(ctype_train ~ .,
                    method = "Rborist",
                    tuneGrid = data.frame(predFixed = 2, minNode = c(3, 50)),
                    data = newdat_train,
                    trControl=control)

# Accuracy on train set
train_rf_acc <- max(train_rf$results$Accuracy)
train_rf_acc
# Model Accuracy on test set
rf_acc <- as.numeric(confusionMatrix(predict(train_rf, newdat_test), ctype_test)$overall["Accuracy"])
rf_acc

# Table comparing train and test accuracy
tab_rf <- data.frame(Method="RF", train_accuracy=train_rf_acc, test_accuracy=rf_acc) 
tab_rf %>% kable()
##****************************************************************
##               rf with 15 principal components                **
##****************************************************************

# Build rf model using Rborist package
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
control <- trainControl(method='repeatedcv', number=10, repeats=3)
fift_train_rf <- train(ctype_train ~ .,
                  method = "Rborist",
                  tuneGrid = data.frame(predFixed = 2, minNode = c(3, 50)),
                  data = fift_newdat_train,
                  trControl=control)

# Model Accuracy
fift_rf_acc <- as.numeric(confusionMatrix(predict(fift_train_rf, fift_newdat_test), ctype_test)$overall["Accuracy"])

# Add results to table
model_results <- bind_rows(model_results,
                        data_frame(Method = "2: Random Forest (RF)",
                        Accuracy_37_principal_components = rf_acc,
                        Accuracy_15_principal_components=fift_rf_acc))

#################################################################
##           Model 3: Multinomial Logistic Regression (MLR)    ##
#################################################################

# Create reference
newdat_train$ctype_train <- relevel(newdat_train$ctype_train, ref = "breast")

# Build model with newly made reference
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
train_multi_log <- multinom(ctype_train ~ ., data = newdat_train)

# Summary of model
summary(train_multi_log)

# Table of predicted probabilities for first 6 people in training set
head(round(fitted(train_multi_log), 2)) %>% kable(format = "simple")

# Predicting the values for train dataset
newdat_train$ctypepredicted <- predict(train_multi_log, newdata = newdat_train, "class")

# Building classification table
train_tab <- table(newdat_train$ctype_train, newdat_train$ctypepredicted)

# Calculating accuracy - sum of diagonal elements divided by total obs
train_tab_acc <- (round((sum(diag(train_tab))/sum(train_tab))*100,2))/100

# Predicting the class for test dataset
newdat_test$ctypepredicted <- predict(train_multi_log, newdata = newdat_test, "class")

# Building classification table
tab <- table(newdat_test$ctype_test, newdat_test$ctypepredicted)
tab %>% kable(caption = "Classification table of actual (rows) versus predicted (columns) cancer types")

# Compute final accuracy
mlr_acc <- (round((sum(diag(tab))/sum(tab))*100,2))/100

# Table comparing train and test accuracy
tab_MLR <- data.frame(Method="MLR", train_accuracy=train_tab_acc, test_accuracy=mlr_acc)
tab_MLR %>% kable()


##****************************************************************
##               MLR with 15 principal components               **
##****************************************************************

# Create reference
fift_newdat_train$ctype_train <- relevel(fift_newdat_train$ctype_train, ref = "breast")

# Build model with newly made reference
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
fift_train_multi_log <- multinom(ctype_train ~ ., data = fift_newdat_train)

# Summary of model
summary(fift_train_multi_log)

# Table of predicted probabilities for first 6 people in training set
head(round(fitted(fift_train_multi_log), 2)) %>% kable()

# Predicting the values for train dataset
fift_newdat_train$ctypepredicted <- predict(fift_train_multi_log, newdata = fift_newdat_train, "class")

# Predicting the class for test dataset
fift_newdat_test$ctypepredicted <- predict(fift_train_multi_log, newdata = fift_newdat_test, "class")

# Building classification table
fift_tab <- table(fift_newdat_test$ctype_test, fift_newdat_test$ctypepredicted)

# Compute final accuracy
fift_mlr_acc <- (round((sum(diag(fift_tab))/sum(fift_tab))*100,2))/100

# Add to results table
model_results <- bind_rows(model_results,
                           data_frame(Method = "3: Multinomial Logistic Regression (MLR)",
                                      Accuracy_37_principal_components=mlr_acc,
                                      Accuracy_15_principal_components=fift_mlr_acc))
