pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  parallel
)

# modified from David Benkeser file
# simulate data from linear/logistic models
simulation_function <- function(iteration){
  
  n <- 1e3
  W <- runif(n)
  A <- rbinom(n, 1, plogis(W))
  logit_Y <- rbinom(n, 1, plogis(A + W + A*W))
  continuous_Y <- A + W + A*W +rnorm(n)
  linear_Y <- as.numeric(continuous_Y > median(continuous_Y))
  data <- data.frame(monte_carlo = iteration, W, A, logit_Y, linear_Y)
  # these are used throughout to compute g-comp
  data_Ais1 <- data; data_Ais1$A <- 1
  data_Ais0 <- data; data_Ais0$A <- 0
  
  analytic_function <- function(outcome_model = model1){
    pred_Ais1 <- predict(
      outcome_model, newdata = data_Ais1, type = "response"
    )
    pred_Ais0 <- predict(
      outcome_model, newdata = data_Ais0, type = "response"
    )
    
    # E[Y(1)] estimate (for ATE)
    EY1_hat <- mean(pred_Ais1)
    # E[Y(0)] estimate (for ATT)
    EY0_hat <- mean(pred_Ais0)
    # E[Y(1) | A = 1] estimate (for ATT)
    EY1_Ais1_hat <- mean(pred_Ais1[A == 1])
    # E[Y(0) | A = 1] estimate (for ATT)
    EY0_Ais1_hat <- mean(pred_Ais0[A == 1])
    
    # contrasts 1
    ATE_hat <- EY1_hat - EY0_hat
    ATT_hat <- EY1_Ais1_hat - EY0_Ais1_hat
    
    # contrasts 2
    ATE_logOR_hat <- log((EY1_hat/(1 - EY1_hat))/(EY0_hat/(1 - EY0_hat)))
    ATT_logOR_hat <- log((EY1_Ais1_hat/(1 - EY1_Ais1_hat))/(EY0_Ais1_hat/(1 - EY0_Ais1_hat)))
    
    res <- data.frame(comparisonRD = ATE_hat - ATT_hat,
                      comparisonOR = ATE_logOR_hat - ATT_logOR_hat)
    
    return(res)
  }
  
  # case 1, linear DGM, single linear regression, no interactions
  model1 <- glm(
    linear_Y ~ A + W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 2, linear DGM, single linear regression, interactions
  model2 <- glm(
    linear_Y ~ A*W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 3, logit DGM, single linear regression, no interactions
  model3 <- glm(
    logit_Y ~ A+W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 4, logit DGM, single linear regression, interactions
  model4 <- glm(
    logit_Y ~ A*W, 
    data = data,
    family = gaussian(link = "identity")
  )
  
  # case 5, linear DGM, single logit regression, no interactions
  model5 <- glm(
    linear_Y ~ A + W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 6, linear DGM, single logit regression, interactions
  model6 <- glm(
    linear_Y ~ A*W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 7, logit DGM, single logit regression, no interactions
  model7 <- glm(
    logit_Y ~ A+W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 8, logit DGM, single logit regression, interactions
  model8 <- glm(
    logit_Y ~ A*W, 
    data = data,
    family = binomial(link = "logit")
  )
  
  model_list <- list(model1, model2, model3, model4, model5, model6, model7, model8)
  
  simulation_results <- lapply(model_list, function (x) analytic_function(outcome_model = x))
  
  sim_res <- do.call(rbind, simulation_results)
  
  sim_res$scenario <- 1:8
  
  return(cbind(iteration, sim_res))

}


## TODO: ADD RISK RATIO
set.seed(123)
sim_res <- mclapply(1:1e4, 
                    function(x) simulation_function(iteration = x),
                    mc.cores = detectCores() - 2)

res <- do.call(rbind, sim_res)

head(res, 10)

results <- res %>% 
  group_by(scenario) %>% 
  summarize(RD_diff = mean(comparisonRD),
            OR_diff = mean(comparisonOR))

round(results, 4)

aug_res <- data.frame(
  DGM = rep(c("Linear", "Linear", "Logit", "Logit"),2),
  AnalyticModel = c(rep("Linear", 4),rep("Logit", 4)),
  InteractionInModel = rep(c("No", "Yes"), 4)
)

reasons <- data.frame(
  reasonRD = c(),
  reasonOR = c()
)

cbind(aug_res, round(results, 4))
