library(arsenal)

setwd("eave_simulated_data/")

B <- read.csv("simulations/demographics.csv", header = TRUE)
B <- B[, c(1, 3, 4, 2)]
B$SCSIMD5 <- as.factor(B$SCSIMD5)
B$NumComorbidities <- as.factor(B$NumComorbidities)
levels(B$Sex) <- c("Female", "Male")
levels(B$AgeGroup) <- gsub("years", "", levels(B$AgeGroup))
levels(B$SCSIMD5) <- c("1 (Most deprived)", 2:4, "5 (Least deprived)")
levels(B$NumComorbidities) <- c(0:3, "4+")
tabledem <- tableby( ~ Sex + SCSIMD5 + NumComorbidities + AgeGroup, data = B)
labels(tabledem) <- c(
  Overall   = "Total", 
  Sex   = "Sex", 
  AgeGroup   = "Age Group", 
  SCSIMD5   = "SCSIMD5", 
  NumComorbidities   = "Number of comorbidities")
summary(tabledem, text = "latex")

