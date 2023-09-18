
pval_khi2 <- function(var1,var2){
  t <- table(var1,var2)
  return(chisq.test(t)$p.value)
}


pval_compute_student <- function(data_list,quanti_names =c("age","vp","PF","lactat"),quali_names = c("sdra","centre","immunodep")){
  t(sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[,k],data$treatment_allocation)),sapply(quanti_names,function(k) t.test(data[,k][data$treatment_allocation == 1],data[,k][data$treatment_allocation == 0])$p.value)))) 
}

pval_compute_wilcox <- function(data_list,quanti_names =c("age","vp","PF","lactat"),quali_names = c("sdra","centre","immunodep")){
  t(sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[,k],data$treatment_allocation)),sapply(quanti_names,function(k) wilcox.test(data[,k][data$treatment_allocation == 1],data[,k][data$treatment_allocation == 0])$p.value)))) 
}

pval_compute_subgroup_wilcox <- function(data_list,strat_name = "sdra",quanti_names =c("age","vp","PF","lactat"),quali_names = c("centre","immunodep")){
  ARDS0 = sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[data[,strat_name] == 0,k],data[data[,strat_name] == 0,"treatment_allocation"])),sapply(quanti_names,function(k) wilcox.test(data[data$treatment_allocation == 0 & data[,strat_name] == 0,k],data[data$treatment_allocation == 1 & data[,strat_name] == 0,k])$p.value))) 
  ARDS1 = sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[data[,strat_name] == 1,k],data[data[,strat_name] == 1,"treatment_allocation"])),sapply(quanti_names,function(k) wilcox.test(data[data$treatment_allocation == 0 & data[,strat_name] == 1,k],data[data$treatment_allocation == 1 & data[,strat_name] == 1,k])$p.value)))
  rownames(ARDS0) = c(quali_names,quanti_names)
  rownames(ARDS1) = c(quali_names,quanti_names)
  return(list(t(ARDS0),t(ARDS1)))
}

pval_compute_subgroup_student <- function(data_list,strat_name = "sdra",quanti_names =c("age","vp","PF","lactat"),quali_names = c("centre","immunodep")){
  ARDS0 = sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[data[,strat_name] == 0,k],data[data[,strat_name] == 0,"treatment_allocation"])),sapply(quanti_names,function(k) t.test(data[data$treatment_allocation == 0 & data[,strat_name] == 0,k],data[data$treatment_allocation == 1 & data[,strat_name] == 0,k])$p.value))) 
  ARDS1 = sapply(data_list, function(data) c(sapply(quali_names,function(k) pval_khi2(data[data[,strat_name] == 1,k],data[data[,strat_name] == 1,"treatment_allocation"])),sapply(quanti_names,function(k) t.test(data[data$treatment_allocation == 0 & data[,strat_name] == 1,k],data[data$treatment_allocation == 1 & data[,strat_name] == 1,k])$p.value)))
  rownames(ARDS0) = c(quali_names,quanti_names)
  rownames(ARDS1) = c(quali_names,quanti_names)
  return(list(t(ARDS0),t(ARDS1)))
}


MSE <- function(data_list){
  n = length(data_list)
  real_v = numeric(n)
  estimated_v = numeric(n)
  for(k in 1:n){
    data_list[[k]]$treatment_allocation = factor(data_list[[k]]$treatment_allocation, levels = c("0","1"))
    data_list[[k]]$outcome = factor(data_list[[k]]$outcome, levels = c("0","1"))
    data_list[[k]]$centre = factor(data_list[[k]]$centre, levels = c("5","1","3","4","7","8","9")) 
    model = glm(data_list[[k]][,c("centre",
                                  "age",
                                  "sdra",
                                  "PF",
                                  "immunodep",
                                  "vp",
                                  "lactat",
                                  "treatment_allocation",
                                  "outcome")],
                formula = outcome ~ . -1, family = "binomial")
    coefficients = coef(model)
    estimated_treatment = coefficients["treatment_allocation1"]
    estimated_v[k] = estimated_treatment
    real_v[k] = data_list[[k]]$beta_treatment[1] 
  }
  return(list("RMSE" = sqrt((1/n) * sum((real_v - estimated_v)**2)),
              "real_beta" = real_v,
              "estimated_beta" = estimated_v))
}

power <- function(nA,nB,nA_1,nB_1,OR = pA*(1-pB)/pB/(1-pA),alpha =0.05){
  pA = nA_1/nA
  pB = nB_1/nB
  kappa=nA/nB
  alpha=alpha
  z=log(OR)*sqrt(nB)/sqrt(1/(kappa*pA*(1-pA))+1/(pB*(1-pB)))
  return(pnorm(z-qnorm(1-alpha/2))+pnorm(-z-qnorm(1-alpha/2)))
}

unadjusted_power <- function(data_list,alpha =0.05){
  n = nrow(data_list[[1]])
  res = numeric(length(data_list))
  for(index in 1:length(data_list)){
    data = data_list[[index]]
    nA = length(data[which(data$treatment_allocation == 1),1])
    nB = length(data[which(data$treatment_allocation == 0),1])
    nA_1 = length(data[which(data$outcome == 1 & data$treatment_allocation == 1),1])
    nB_1 = length(data[which(data$outcome == 1 & data$treatment_allocation == 0),1])
    res[index] = power(nA,nB,nA_1,nB_1,alpha = alpha)
  }
  list_res = list("mean" = mean(res,na.rm = TRUE), "sd" = sd(res))
  return(list_res)
}

adjusted_power <- function(data_list,RMSE_list,alpha =0.05){
  n = nrow(data_list[[1]])
  res = numeric(length(data_list))
  for(index in 1:length(data_list)){
    data = data_list[[index]]
    adjusted_OR = exp(RMSE_list$estimated_beta[[index]])
    nA = length(data[which(data$treatment_allocation == 1),1])
    nB = length(data[which(data$treatment_allocation == 0),1])
    nA_1 = length(data[which(data$outcome == 1 & data$treatment_allocation == 1),1])
    nB_1 = length(data[which(data$outcome == 1 & data$treatment_allocation == 0),1])
    res[index] = power(nA,nB,nA_1,nB_1,alpha = alpha,OR = adjusted_OR)
  }
  list_res = list("mean" = mean(res,na.rm = TRUE), "sd" = sd(res))
  return(list_res)
}

significant_betas <- function(data_list,alpha = 0.05){
  n = length(data_list)
  n_beta_signif = 0
  n_beta_sdra0_signif = 0
  n_beta_sdra1_signif = 0
  n_close_testing_signif = 0
  for(k in 1:n){
    data_list[[k]]$treatment_allocation = factor(data_list[[k]]$treatment_allocation, levels = c("0","1"))
    data_list[[k]]$outcome = factor(data_list[[k]]$outcome, levels = c("0","1"))
    data_list[[k]]$centre = factor(data_list[[k]]$centre, levels = c("5","1","3","4","7","8","9")) 
    
    model_tot = glm(data_list[[k]][,c("centre",
                                  "age",
                                  "sdra",
                                  "PF",
                                  "immunodep",
                                  "vp",
                                  "lactat",
                                  "treatment_allocation",
                                  "outcome")],
                formula = outcome ~ . -1, family = "binomial")
    
    model_sdra0 = glm(data_list[[k]][data_list[[k]]$sdra == 0,c("centre",
                                        "age",
                                        "PF",
                                        "immunodep",
                                        "vp",
                                        "lactat",
                                        "treatment_allocation",
                                        "outcome")],
                      formula = outcome ~ . -1, family = "binomial")
    
    model_sdra1 = glm(data_list[[k]][data_list[[k]]$sdra == 1,c("centre",
                                                                "age",
                                                                "PF",
                                                                "immunodep",
                                                                "vp",
                                                                "lactat",
                                                                "treatment_allocation",
                                                                "outcome")],
                      formula = outcome ~ . -1, family = "binomial")
    pval_tot = summary(model_tot)$coefficients["treatment_allocation1", "Pr(>|z|)"]
    pval_sdra0 = summary(model_sdra0)$coefficients["treatment_allocation1", "Pr(>|z|)"]
    pval_sdra1 = summary(model_sdra1)$coefficients["treatment_allocation1", "Pr(>|z|)"]
    if(pval_tot <= alpha){
      n_beta_signif = n_beta_signif + 1
    }
    if(pval_sdra0 <= alpha){
      n_beta_sdra0_signif = n_beta_sdra0_signif + 1
    }
    if(pval_sdra1 <= alpha){
      n_beta_sdra1_signif = n_beta_sdra1_signif + 1
    }
    if(pval_sdra0 <= alpha & pval_sdra1 <= alpha){
      n_close_testing_signif = n_close_testing_signif + 1
    }
  }
  return(list("beta_signif" = round(n_beta_signif * 100/n,1),
              "beta_sdra0_signif" = round(n_beta_sdra0_signif * 100 / n,1),
              "beta_sdra1_signif" = round(n_beta_sdra1_signif * 100 / n,1),
              "close_testing_signif" = round(n_close_testing_signif * 100 / n,1)))
}

