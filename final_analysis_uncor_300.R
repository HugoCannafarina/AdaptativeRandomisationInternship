rm(list=ls())
source("desc_functions.R")
library(ggplot2)

#######################################
### Correlated datasets of size 300 ###
#######################################

load(file = "res_uncor_300/res_blocks_uncor_300.rda")
load(file = "res_uncor_300/res_MSB_uncor_300.rda")
load(file = "res_uncor_300/res_pocock_2_uncor_300.rda")
load(file = "res_uncor_300/res_pocock_4_uncor_300.rda")


########### RMSE ##########

mse_blocks <- MSE(res_blocks)
mse_MSB <- MSE(res_MSB)
mse_p2 <- MSE(res_pocock_2)
mse_p4 <- MSE(res_pocock_4)

round(mse_MSB$RMSE,2)
round(mse_blocks$RMSE,2)
round(mse_p2$RMSE,2)
round(mse_p4$RMSE,2)

round(mse_MSB$RMSE**2,2)
round(mse_blocks$RMSE**2,2)
round(mse_p2$RMSE**2,2)
round(mse_p4$RMSE**2,2)



ggplot(data = data.frame(mse_MSB), aes(x = estimated_beta, y = real_beta)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),  # Taille des nombres sur les axes (augmentée à 14 ici)
    axis.title = element_blank()  # Supprimer les titres des axes
  )

ggplot(data = data.frame(mse_blocks), aes(x = estimated_beta, y = real_beta)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),  # Taille des nombres sur les axes (augmentée à 14 ici)
    axis.title = element_blank()  # Supprimer les titres des axes
  )

ggplot(data = data.frame(mse_p2), aes(x = estimated_beta, y = real_beta)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),  # Taille des nombres sur les axes (augmentée à 14 ici)
    axis.title = element_blank()  # Supprimer les titres des axes
  )

ggplot(data = data.frame(mse_p4), aes(x = estimated_beta, y = real_beta)) +
  geom_point() +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),  # Taille des nombres sur les axes (augmentée à 14 ici)
    axis.title = element_blank()  # Supprimer les titres des axes
  )

res_RMSE = c("MSB" = mse_MSB$RMSE, "Stratified-blocks" = mse_blocks$RMSE, "Pocock (2 classes)" = mse_p2$RMSE, "Pocock (4 classes)" = mse_p4$RMSE )



########### Imbalance ,stratified student ##########

pval_sg_s_MSB <- pval_compute_subgroup_student(res_MSB)
pval_sg_s_blocks <- pval_compute_subgroup_student(res_blocks)
pval_sg_s_p2 <- pval_compute_subgroup_student(res_pocock_2)
pval_sg_s_p4 <- pval_compute_subgroup_student(res_pocock_4)


round(sapply(1:ncol(pval_sg_s_MSB[[1]]),function(k) length(which(pval_sg_s_MSB[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_s_MSB[[2]]),function(k) length(which(pval_sg_s_MSB[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_s_blocks[[1]]),function(k) length(which(pval_sg_s_blocks[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_s_blocks[[2]]),function(k) length(which(pval_sg_s_blocks[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_s_p2[[1]]),function(k) length(which(pval_sg_s_p2[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_s_p2[[2]]),function(k) length(which(pval_sg_s_p2[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_s_p4[[1]]),function(k) length(which(pval_sg_s_p4[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_s_p4[[2]]),function(k) length(which(pval_sg_s_p4[[2]][,k] <= 0.05)))*100/5000,2)


########### Imbalance ,non-stratified student ##########

pval_s_MSB <- pval_compute_student(res_MSB)
pval_s_blocks <- pval_compute_student(res_blocks)
pval_s_p2 <- pval_compute_student(res_pocock_2)
pval_s_p4 <- pval_compute_student(res_pocock_4)

round(sapply(1:ncol(pval_s_MSB),function(k) length(which(pval_s_MSB[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_s_blocks),function(k) length(which(pval_s_blocks[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_s_p2),function(k) length(which(pval_s_p2[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_s_p4),function(k) length(which(pval_s_p4[,k] <= 0.05)))*100/5000,2)

########### Imbalance ,stratified wilcoxon ##########

pval_sg_w_MSB <- pval_compute_subgroup_wilcox(res_MSB)
pval_sg_w_blocks <- pval_compute_subgroup_wilcox(res_blocks)
pval_sg_w_p2 <- pval_compute_subgroup_wilcox(res_pocock_2)
pval_sg_w_p4 <- pval_compute_subgroup_wilcox(res_pocock_4)


round(sapply(1:ncol(pval_sg_w_MSB[[1]]),function(k) length(which(pval_sg_w_MSB[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_w_MSB[[2]]),function(k) length(which(pval_sg_w_MSB[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_w_blocks[[1]]),function(k) length(which(pval_sg_w_blocks[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_w_blocks[[2]]),function(k) length(which(pval_sg_w_blocks[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_w_p2[[1]]),function(k) length(which(pval_sg_w_p2[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_w_p2[[2]]),function(k) length(which(pval_sg_w_p2[[2]][,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_sg_w_p4[[1]]),function(k) length(which(pval_sg_w_p4[[1]][,k] <= 0.05)))*100/5000,2)
round(sapply(1:ncol(pval_sg_w_p4[[2]]),function(k) length(which(pval_sg_w_p4[[2]][,k] <= 0.05)))*100/5000,2)

########### Imbalance ,non-stratified wilcoxon ##########

pval_w_MSB <- pval_compute_wilcox(res_MSB)
pval_w_blocks <- pval_compute_wilcox(res_blocks)
pval_w_p2 <- pval_compute_wilcox(res_pocock_2)
pval_w_p4 <- pval_compute_wilcox(res_pocock_4)

round(sapply(1:ncol(pval_w_MSB),function(k) length(which(pval_w_MSB[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_w_blocks),function(k) length(which(pval_w_blocks[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_w_p2),function(k) length(which(pval_w_p2[,k] <= 0.05)))*100/5000,2)

round(sapply(1:ncol(pval_w_p4),function(k) length(which(pval_w_p4[,k] <= 0.05)))*100/5000,2)

########### Unadjusted power ##########


round(as.numeric(unadjusted_power(res_MSB)),1)
round(as.numeric(unadjusted_power(res_blocks)),1)
round(as.numeric(unadjusted_power(res_pocock_2)),1)
round(as.numeric(unadjusted_power(res_pocock_4)),1)

########### Adjusted power ##########

round(as.numeric(adjusted_power(res_MSB,mse_MSB)),1)
round(as.numeric(adjusted_power(res_blocks,mse_blocks)),1)
round(as.numeric(adjusted_power(res_pocock_2,mse_p2)),1)
round(as.numeric(adjusted_power(res_pocock_4,mse_p4)),1)

########### Significant betas and adjusted power ##########

res_betas_msb <- significant_betas(res_MSB)
res_betas_blocks <- significant_betas(res_blocks)
res_betas_p2 <- significant_betas(res_pocock_2)
res_betas_p4 <- significant_betas(res_pocock_4)

round(as.numeric(res_betas_msb),1)
round(as.numeric(res_betas_blocks),1)
round(as.numeric(res_betas_p2),1)
round(as.numeric(res_betas_p4),1)

