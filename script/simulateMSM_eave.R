library(parallel)
library(foreach)
library(iterators)
library(doParallel)
set.seed(1)
expit <- function(x) return( exp(x) / (1+exp(x)) )

setwd("/Users/victor/eave_simulated_data/")
# setwd("/home/victor/trial_emulation_project/seaman_paper/")

multimorbidity_table <- read.csv("data/multimorbidity.csv", stringsAsFactors = TRUE)
multimorbidity_table$SCSIMD5 <- as.factor(multimorbidity_table$SCSIMD5)
multimorbidity_table$NumComorbidities <- as.factor(multimorbidity_table$NumComorbidities)
multimorbidity_table$p <- multimorbidity_table$n/sum(multimorbidity_table$n)
multimorbidity_table <- multimorbidity_table[!is.na(multimorbidity_table$SCSIMD5), ]

# multimorbidity_table$sex <- ifelse(multimorbidity_table$sex == "F", 0, 1)

### multimorbidity_table$AgeGroup <- ifelse(multimorbidity_table$AgeGroup < "60", "18-59", "60+")
# multimorbidity_table$AgeGroup <- ifelse(multimorbidity_table$AgeGroup < "60", 0, 1)

### multimorbidity_table$SCSIMD5 <- ifelse(multimorbidity_table$SCSIMD5 < 3, "1-2", "3+")
# multimorbidity_table$SCSIMD5 <- ifelse(multimorbidity_table$SCSIMD5 < 3, 0, 1)

### multimorbidity_table$Number_count <- ifelse(multimorbidity_table$Number_count < 2, "0-1", "2+")
# multimorbidity_table$Number_count <- ifelse(multimorbidity_table$Number_count < 2, 0, 1)

# multimorbidity_table <- as.matrix(multimorbidity_table)

### multimorbidity_table$sex <- as.factor(multimorbidity_table$sex)
### multimorbidity_table$AgeGroup <- as.factor(multimorbidity_table$AgeGroup)
### multimorbidity_table$SCSIMD5<- as.factor(multimorbidity_table$SCSIMD5)
### multimorbidity_table$Number_count <- as.factor(multimorbidity_table$Number_count)

# browser()
# Parameters
N <- 100000 # Sample size
K <- 3 # Number of visits

# Initialise variables
i <- 1 # Index for individuals in the database i=1:N
B <- matrix(0, nrow = N, ncol = ncol(multimorbidity_table)-2)
colnames(B) <- colnames(multimorbidity_table[1:(ncol(multimorbidity_table)-2)])
A <- matrix(0, nrow = N, ncol = K)
Ys <- matrix(0, nrow = N, ncol = K)
  
# Parameters for treatment allocation
alphas <- c(1, -0.05, 0.1, -0.05, 0.25)  # The first alpha is for the baseline risk. The others for the four confounders

# Parameters for the marginal structural model
betas <- c(-2, -0.5, -0.5, -0.5) # The first alpha is for the baseline risk. The second and third are first and second vaccine dose

# Parameters for the confounding mechanism
gammas <- rep(0.1, 4)

# correlation parameter of Gaussian copula
rhos <- seq(0, 0.4, 0.1)


# Step 1. Simulate from p(multimorbidity_table). No covariates in the MSM for our purposes. Set k = 0.
B <- multimorbidity_table[
  sample(nrow(multimorbidity_table), size = N, replace = TRUE, prob = multimorbidity_table[, "n"]), 1:(ncol(multimorbidity_table)-2)
]
rownames(B) <- NULL
# Save demographics data
write.csv(B, "simulations/demographics.csv", quote = FALSE, row.names = FALSE)
B.numeric <- matrix(0, nrow = nrow(B), ncol = ncol(B), dimnames = list(NULL, colnames(B)))
for (j in 1:ncol(B)) {
  B.numeric[, j] <- as.numeric(B[, j]) - 1
}
# cdf for H
# hCDF <- as.list(cumsum(table(sort(as.vector(B.numeric %*% gammas))))/nrow(B.numeric))
uniqueHValues <- as.vector(B.numeric %*% gammas)
foreach (proc = 1:5) %dopar% {
  rho = rhos[proc]
  for (i in 1:N) {
    k <- 0
    # Step 2. Simulate from p(L_k | -). No time-varying confounders for our purposes.
    if (i %% 100 == 0) cat(proc, i, "\n")
    # browser()
    while (k < K) {
      # Step 3. Simulate from p(A_k | -). Treatment variable
      # logitp <- -1 + c(B[i, ] %*% rep(0.2, 4))
      prob.tmt <- expit(sum(alphas * c(1, as.vector(B.numeric[i, ]))) + 0.5 * sum(A[i, ]))
      A[i, k+1] <- runif(1) < prob.tmt
      
      # Step 4. Calculate the H_k.
      # Define h_0^{a_0} = -0.1 * sum(B)
      h <- c(as.vector(B.numeric[i, ]) %*% gammas)
      
      # Step 5. Calculate U_H, Z_H.
      u.h.min <- length(which(uniqueHValues <= h-1))/N
      u.h.max <- length(which(uniqueHValues <= h))/N
      u.h <- runif(1, u.h.min, u.h.max)
      z.h <- qnorm(u.h)
      
      # Step 6. Calculate Z_Y, U_Y.
      z.y <- rnorm(1, rho*z.h, sqrt(1-rho^2))
      u.y <- pnorm(z.y)
      
      # Step 7. Compute Y_{k+1}.
      g.score <- expit(as.vector(c(1, A[i, ]) %*% betas))
      Ys[i, k+1] <- u.y >= g.score
  
      # Step 8
      if (Ys[i, k+1] == 0) {
        break
      } else {
        k <- k + 1
      }
    }
  }
  colnames(A) <- paste0("A", 0:(K-1))
  colnames(Ys) <- paste0("Y", 1:K)
  write.csv(A, paste0("simulations/", rho, "/treatment.csv"), quote = FALSE, row.names = FALSE)
  write.csv(Ys, paste0("simulations/", rho, "/outcome.csv"), quote = FALSE, row.names = FALSE)
}
simulatedData <- cbind(B.numeric, A, Ys)
colnames(simulatedData) <- c(colnames(B), paste0("A", 0:(K-1)), paste0("Y", 1:K))
simulatedData <- as.data.frame(simulatedData)
head(simulatedData)

e0 <- fitted(
  glm(
    A0 ~ Sex + AgeGroup + SCSIMD5 + NumComorbidities,
    family = binomial,
    data = simulatedData
  )
)

wNum <- fitted(
  glm(
    A0 ~ 1,
    family = binomial,
    data = simulatedData
  )
)

wDen <- fitted(
  glm(
    A0 ~ Sex + AgeGroup + SCSIMD5 + NumComorbidities,
    family = binomial,
    data = simulatedData
  )
)
simulatedData$weights <- (1/e0) * simulatedData$A0 + (1/(1-e0)) * (1-simulatedData$A0)
simulatedData$weightsNum <- wNum * simulatedData$A0 + (1-wNum) * (1-simulatedData$A0)
simulatedData$weightsDen <- wDen * simulatedData$A0 + (1-wDen) * (1-simulatedData$A0)
modelfit <- glm((Y1 == 0) ~ as.factor(A0), family = binomial, data = simulatedData, weights = weightsNum/weightsDen)
summary(modelfit)

modelfit <- glm((Y1 == 0) ~ as.factor(A0), family = binomial, data = simulatedData, weights = weights)
summary(modelfit)
