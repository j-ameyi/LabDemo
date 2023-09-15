rm(list= ls())

rnorm(10)

error.1


edit something

# For the multiple linear regression slides, 
prostate<-read.delim("prostate.txt", row.names=1)
prostate$gleason<-as.factor(ifelse(prostate$gleason>=8, 8, prostate$gleason))
prostate.test <- prostate[!prostate$train,]
prostate <- prostate[prostate$train,]

# For the logistic regression slides:
heart<-read.csv("heartdisease.csv", row.names=1)


#####
# LAB 10

install.packages(c("glmnet", "pROC", "bootstrap"))


# 1 -3
library(Biobase)
load("Module10.RData")

dim(kidney)
kidney

# 4: Penalized mulitple linear regression
library(glmnet)
lasso.kdpi <- glmnet( x = t(exprs(kidney)), y = pData(kidney)$KDPI, 
                      family = gaussian, alpha = 1)
# number of lambda's
length(lasso.kdpi$lambda)


# 5. Identifying the "best" model
# 10-fold cross validation
library(bootstrap)

# Next, in order to use the crossval function we must first define our 
# predictive model as a function:
lasso.fit <- function(x, y, lambda) {glmnet(x = x, y = y, lambda = lambda, 
                                            family = gaussian, alpha = 1)}

# The crossval function also requires that we define a function the takes the
# model fit and outputs a prediction
lasso.predict <- function(fit, x) { predict(fit, newx = x, type = "response") }

# A "list" object to store our cross-validation results since predict() gives a list
lasso.cv <- list()

# calculate prediction mean squared error for each lambda
# set seed cos we want each lambda evaluated on the same random partitions of the training set
for (s in 1:length(lasso.kdpi$lambda)) {
     set.seed(4)
     lasso.cv[[s]]<-crossval(x = t(exprs(kidney)), y = pData(kidney)$KDPI, 
                              lambda = lasso.kdpi$lambda[s], 
                                        lasso.fit, lasso.predict, ngroup=10)
}

lasso.cv[[1]]
lasso.cv[[1]]$leave.out


# 6
names(lasso.cv[[1]])

# These below store the indexes (row numbers in the training dataset) of the observations 
# assigned to that fold.
lasso.cv[[1]]$groups
lasso.cv[[2]]$groups
lasso.cv[[3]]$groups
lasso.cv[[100]]$groups
all.equal(lasso.cv[[1]]$groups, lasso.cv[[2]]$groups)

# 7
# The predictions for each value of lambda (hence each LASSO penalized model) 
length(lasso.cv[[100]]$cv.fit)


# 8: Prediction mean squared error 
# First, set up an object to store the prediction mean squared error.
mse <- numeric()

# for each value of lambda.  sum(y - yhat)^2 /N
for (i in 1:length(lasso.cv)) {
      mse[i]<-sum((pData(kidney)$KDPI - lasso.cv[[i]]$cv.fit)^2)/190
}

which.min(mse)

# 9
# Plot the cross-validated prediction mean squared error against lambda.
plot(lasso.kdpi$lambda, mse)
lasso.kdpi$lambda[which.min(mse)]
abline(v=lasso.kdpi$lambda[which.min(mse)], lty=3, col="red")


# 10
# Now that we know what model is optimal, let's examine the coefficients. 
#First, extract the coefficient estimates for each value of lambda
coefficients <- coef(lasso.kdpi)

# This could equivalently be achieved using 
coefficients2 <- predict(lasso.kdpi, type="coefficients")
dim(coefficients2)

# Our coefficients object includes estimated parameters for our variables (genes)
# in rows and one column for each value of lambda. To extract our optimal model,
final.model <- coefficients[, which.min(mse)]

# How many non-zero coefficient estimates are in this final model?
sum(final.model != 0)


# 11-12: Penalized logistic regression

sum(pData(kidney)$DGF == 1)

# 13
# Elastic net penalized logistic regression
# x must have subjects in rows, genes in columns
# Perhaps we want to encourage model parsimony so we set alpha = 0.80 which
# places more weight on the LASSO penalty
enet.dgf <- glmnet( x = t(exprs(kidney)), y = pData(kidney)$DGF, 
                    family = binomial, alpha = 0.80)


# How many elastic net models (that is, how many lambda values) were fit?
length(enet.dgf$lambda)


# 14
# Identifying the "best" model among all that were fit. We do not have a test set
# Let's use V-fold cross-validation to help us identify the optimal lambda.
# Specifically, we will use 10-fold cross-validation.

# Next, in order to use the crossval function we must first define our 
# predictive model as a function
enet.fit <- function(x, y, lambda) {glmnet(x = x, y = y, lambda = lambda, 
                                           family = binomial, alpha = 0.8)}

# The crossval function also requires that we define a function the takes the
# model fit and outputs a prediction. 
enet.predict <- function(fit, x) { predict(fit, newx = x, type = "response") }

# predict on glmnet gives a list so our object for storing the cross-validation results is
enet.cv <- list()

# For a binary response, we could choose to minimize prediction error or maximize
# the area under the ROC curve (AUC).
# Let's maximize the AUC for our model selection criterion
for (s in 1:length(enet.dgf$lambda)) {
  set.seed(16)
  enet.cv[[s]]<-crossval(x = t(exprs(kidney)), y = pData(kidney)$DGF, 
                         lambda = enet.dgf$lambda[s],  enet.fit, 
                         enet.predict, ngroup=10)
}

enet.cv[[1]]
enet.cv[[1]]$leave.out


# 15
# We need to calculate the AUC for each of these models using the results from 
# our 10-fold cross-validation procedure. First, set up an object to store the AUC.
AUC <- numeric()

# Next, we need to calculate the AUC for each value of lambda. To do that, load the pROC R package
library(pROC)

# We then wrap the necessary calculation in a for loop.
for (i in 1:length(enet.cv)) {
      AUC[i]<-roc(pData(kidney)$DGF ~ enet.cv[[i]]$cv.fit)$auc
}

which.max(AUC)


# 16
# Plot the cross-validated AUC error against lambda.
plot(enet.dgf$lambda, AUC)
enet.dgf$lambda[which.max(AUC)]
abline(v=enet.dgf$lambda[which.max(AUC)], col="red")


# 17
# Now that we know what model is optimal, let's examine the coefficients. 
# First, extract the coefficient estimates for each value of lambda
coefficients <- coef(enet.dgf)

# To extract our optimal model,
final.model <- coefficients[,which.max(AUC)]
sum(final.model != 0)

# 18
# We may be interested in the predictions associated with our final model.
# We can extract these from our original model using
predictions <- predict(enet.dgf, newx = t(exprs(kidney)), 
                       type = "response")[,which.max(AUC)]

# Examining histogram below, we recognize that the probabilities are returned
hist(predictions)

# Suppose we want a binary classification from this model. 
# Use the coords function in the pROC package and inspect the output from this call:
thresh <- coords(roc(pData(kidney)$DGF ~ predictions))
thresh
thresh[thresh$sensitivity >= 0.90, ]


# 19-20
# Suppose instead we want to simultaneously maximize Sensitivity and Specificity 
# so we use Youden's index as our threshold.
Youden <- coords(roc(pData(kidney)$DGF ~ predictions), x = "best")
Youden

predicted.class <- ifelse(predictions<= Youden$threshold, 0, 1)

# How many subjects with DGF were classified using this rule as having DGF?
table <- table(predicted.class, pData(kidney)$DGF)
table

# 21
sensitivity <- table[2,2] / colSums(table)[2]
sensitivity

specificity <- table[1,1] / colSums(table)[1]
specificity










