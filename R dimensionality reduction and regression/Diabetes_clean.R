# Title: 
# Integrated modelling of labile and glycated hemoglobin with glucose 
# for enhanced diabetes detection and short-term monitoring

# Author: David G. Aragones
# Email: david.aragones@uclm.es

# Corresponding author: Gabriel F. Calvo
# Email: gabriel.fernandez@uclm.es

# Description: This code contains the logistic models and the UMAP 
# developed for dimensionality reduction



# Loading data ##############################################################

# Read data from general participants
data.general <- openxlsx::read.xlsx("General.xlsx")

# Read data from pediatric participants
data.pediatric <- openxlsx::read.xlsx("Pediatric.xlsx")

# Join general and pediatric data side by side
X <- mapply(c, data.general, data.pediatric)

# Rescale values in [0, 1] for each column
Xrs <- sapply(1:ncol(X), 
              function(j) scales::rescale(X[, j]))

# Extract the target variable (first column of data.general)
Y <- data.general[, 1]

# Convert "No" / "Yes" to 0 / 1 in the target variable
Y <- ifelse(Y == "No", 0, 1)



# Logistic models ###########################################################

# Proportion for training data (p), the rest is for testing (1-p)
training.p <- 0.8

# Randomly sample indices for training data
sample.indices <- sample(1:nrow(data.general), training.p * nrow(data.general))

# Create training and testing datasets based on the sampled indices
train.data <- data.general[sample.indices, ]
test.data <- data.general[-sample.indices, ]

# Fit a logistic regression model using training data
logmodel <- glm(Y[sample.indices] ~ ., data = train.data, family = "binomial")

# Predict probabilities on the test data
predicted.probabilities <- predict(logmodel, newdata = test.data, 
                                   type = "response")

# Convert predicted probabilities to binary classes (0 / 1)
predicted.classes <- as.factor(ifelse(predicted.probabilities > 0.5, 1, 0))

# Evaluate model performance using confusion matrix
confusion.matrix <- caret::confusionMatrix(data = predicted.classes, 
                                           reference = as.factor(Y[-sample.indices]))

# Calculate ROC curve
roc_obj <- pROC::roc(Y[-sample.indices], predicted.probabilities)

# Calculate AUC score
auc_score <- pROC::auc(roc_obj)

# Set up 10-fold cross-validation
cv_indices <- caret::createFolds(Y[sample.indices], k = 10)

# Perform cross-validation
for (fold in seq_along(cv_indices)) {
  
  # Extract indices for the current fold
  cv_train_indices <- unlist(cv_indices[-fold])
  cv_test_indices <- cv_indices[[fold]]
  
  # Create training and testing datasets based on the current fold
  cv_train.data <- data.general[cv_train_indices, ]
  cv_test.data <- data.general[cv_test_indices, ]
  
  # Fit a logistic regression model using training data
  logmodel <- glm(Y[cv_train_indices] ~ ., data = cv_train.data, family = "binomial")
  
  # Predict probabilities on the test data
  predicted.probabilities <- predict(logmodel, newdata = cv_test.data, type = "response")
  
  # Convert predicted probabilities to binary classes (0 / 1)
  predicted.classes <- as.factor(ifelse(predicted.probabilities > 0.5, 1, 0))
  
  # Evaluate model performance using confusion matrix
  confusion.matrix <- caret::confusionMatrix(data = predicted.classes, 
                                      reference = as.factor(Y[cv_test_indices]))
  
  # Store confusion matrix for the current fold
  cv_confusion_matrices[[fold]] <- confusion.matrix
  
  # Calculate ROC curve for the current fold
  roc_obj <- pROC::roc(Y[cv_test_indices], predicted.probabilities)
  
  # Store ROC object for the current fold
  cv_roc_objects[[fold]] <- roc_obj
  
  # Calculate AUC score for the current fold
  auc_score <- pROC::auc(roc_obj)
  
  # Store AUC score for the current fold
  cv_auc_scores[fold] <- auc_score
}



# Dimensionality reduction ##################################################

# Perform PCA on the original data
pca_result <- prcomp(data.general)

# UMAP model using the rescaled data
umap.model <- uwot::umap(Xrs, n_neighbors = 200)








