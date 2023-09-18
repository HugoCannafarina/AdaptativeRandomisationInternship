rm(list=ls())
library(dplyr)
library(ggplot2)
library(corrplot)
library(carat)
library(parallel)
source("sim_functions.R")
load(file = "source_data/uncor_820.rda")


RNGkind("L'Ecuyer-CMRG")
set.seed(92)
mc.reset.stream()

param_MSB = list("xi" = 0.7,
                 "center_name" = "centre",
                 "strat_name" = "sdra")

param_blocks = list("block_size" = 4,
                    "strat_variables" = c("sdra","centre"))

param_pocock_2 = list("xi" = 0.7,
                      "nb_class" = 2,
                      "strat_name" = "sdra")

param_pocock_4 = list("xi" = 0.7,
                      "nb_class" = 4,
                      "strat_name" = "sdra")

nb_coeur <- detectCores() - 1


res_MSB <- mclapply(X = uncor_820,
                    FUN = function(data) MSB_strat(data = data,
                                                   xi = param_MSB$xi,
                                                   center_name = param_MSB$center_name,
                                                   strat_name = param_MSB$strat_name),
                    mc.set.seed = TRUE,
                    mc.cores = nb_coeur)

res_blocks <- lapply(X = uncor_820,
                     FUN = function(data) rando_blocks(data = data,
                                                       block_size = param_blocks$block_size,
                                                       strat_variables = param_blocks$strat_variables))



res_pocock_2 <- lapply(X = uncor_820,
                       FUN = function(data) rando_pocock_strat(data = data,
                                                               xi = param_pocock_2$xi,
                                                               nb_class = param_pocock_2$nb_class,
                                                               strat_name = param_pocock_2$strat_name))


res_pocock_4 <- lapply(X = uncor_820,
                       FUN = function(data) rando_pocock_strat(data = data,
                                                               xi = param_pocock_4$xi,
                                                               nb_class = param_pocock_4$nb_class,
                                                               strat_name = param_pocock_4$strat_name))

load(file = "beta.rda")
res_MSB <- lapply(res_MSB,function(data) sim_outcome(data,coeff))
res_blocks  <- lapply(res_blocks ,function(data) sim_outcome(data,coeff))
res_pocock_2  <- lapply(res_pocock_2 ,function(data) sim_outcome(data,coeff))
res_pocock_4 <- lapply(res_pocock_4 ,function(data) sim_outcome(data,coeff))

save(res_MSB, file = "res_uncor_820/res_MSB_uncor_820.rda")
save(res_blocks, file = "res_uncor_820/res_blocks_uncor_820.rda")
save(res_pocock_2, file = "res_uncor_820/res_pocock_2_uncor_820.rda")
save(res_pocock_4, file = "res_uncor_820/res_pocock_4_uncor_820.rda")