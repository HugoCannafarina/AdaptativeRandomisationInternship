rm(list=ls())
library(openxlsx)
library(dplyr)
library(MASS)
library(ggplot2)
library(corrplot)
library(pracma)
library(carat)
library(VGAM)
source("cov_functions.R")
set.seed(20)


##################################################################
### Sensitivity analysis to choose which correlations to model ###
##################################################################

sc1 <- read.xlsx("SC1.xlsx")

sc1 <- sc1 %>%
  mutate(centre = as.factor(centre),
         immunodep = as.factor(immunodep),
         sdra = as.factor(sdra))

sc1_subset = sc1[which(sc1$lactat > 2),c("idpat","centre","age","sdra","PF","immunodep","vp","lactat")]


cor_matrix <- table_cor(sc1_subset[,c("age","sdra","PF","immunodep","vp","lactat")])


corrplot(cor_matrix,tl.col = "black",addgrid.col = TRUE)

#Creation of a plot to display the p-values of the rank correlation test


pval_rank <- table_pval_rank(sc1_subset[,-c(1,2)])
dist_pval_rank <- sort(as.vector(pval_rank)[!duplicated(as.vector(pval_rank))])
dist_pval_rank <- dist_pval_rank[dist_pval_rank!=0]
id_rank = lapply(dist_pval_rank, function(n) which(pval_rank == n,arr.ind=TRUE)[1:2])
name_rank  = sapply(id_rank,function(index) paste(rownames(pval_rank)[index[1]],colnames(pval_rank)[index[2]], sep = "/"))
dist_pval_rank <- data.frame(index = 1:length(dist_pval_rank), value = dist_pval_rank)
dist_pval_rank$name = name_rank

ggplot(dist_pval_rank, aes(x = index, y = value)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.3, linetype = "dashed") +
  geom_text(aes(label = name), hjust = -0.1, size = 8, angle = 90) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16)
  ) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, 1.5))


###########################################################################
### Research of the appropriate correlation matrix for the NORTA method ###
###########################################################################

var_list = list(c("PF","sdra"), #Specify the correlation that you want to model (ordinal - ordinal not yet implemented)
                c("sdra","age"),
                c("vp","PF"),
                c("immunodep","sdra"),
                c("PF","age"),
                c("lactat","vp"),
                c("immunodep","age"),
                c("immunodep","PF"))

fct_list = list("age" = "quantile", #Specify the quantile functions that you want to use for each variable
                "sdra" = "qbinom", 
                "vp" = "quantile",
                "immunodep" = "qbinom",
                "lactat" = "quantile",
                "PF" = "quantile")

mean_sdra = mean(as.numeric(as.character(sc1_subset$sdra)),na.rm = TRUE) #Proportion of patients with ARDS
mean_immunodep = mean(as.numeric(as.character(sc1_subset$immunodep)),na.rm = TRUE) #Proportion of immunodeficient patients
mean_CDF_sdra = integrate(esp_sdra, lower = -Inf, upper = Inf ,abs.tol = 10**-15)$value #mean of the cumulative distribution function of a Bernoulli with mean equal mean_sdra and after NORTA transformation of parameters
var_CDF_sdra = integrate(esp_squared_sdra, lower = -Inf, upper = Inf,abs.tol = 10**-15)$value - mean_CDF_sdra**2 #variance of the cumulative distribution function of a Bernoulli with mean equal mean_sdra and after NORTA transformation of parameters
mean_CDF_immunodep = integrate(esp_immunodep, lower = -Inf, upper = Inf ,abs.tol = 10**-15)$value #mean of the cumulative distribution function of a Bernoulli with mean equal mean_immunodep and after NORTA transformation of parameters
var_CDF_immunodep = integrate(esp_squared_immunodep, lower = -Inf, upper = Inf,abs.tol = 10**-15)$value - mean_CDF_immunodep**2 #variance of the cumulative distribution function of a Bernoulli with mean equal mean_immunodep and after NORTA transformation of parameters
mean_CDF_quanti = integrate(esp_quanti, lower = -Inf, upper = Inf ,abs.tol = 10**-15)$value #variance of the CDF of a gaussian with gaussian distribution 
var_CDF_quanti = integrate(esp_squared_quanti, lower = -Inf, upper = Inf,abs.tol = 10**-15 )$value - mean_CDF_quanti**2 #variance of the CDF of a gaussian with gaussian distribution 

param_list = list( #prob, mean_cdf and sd_cdf needed for discrete variables, none for quantitative ones.
  "sdra" = list("prob" = mean_sdra, "mean_cdf" = mean_CDF_sdra, "sd_cdf" = sqrt(var_CDF_sdra)),
  "immunodep" = list("prob" = mean_immunodep, "mean_cdf" = mean_CDF_immunodep, "sd_cdf" =  sqrt(var_CDF_immunodep)),
  "quanti" = list("mean_cdf" = mean_CDF_quanti, "sd_cdf" = sqrt(var_CDF_quanti))
)

sigma <- find_cor(data = sc1_subset,
                  cor_matrix = cor_matrix, 
                  var_list = var_list, 
                  param_list = param_list ) #Optimal correlation matrix for the Norta method

cor_820 <- simul_data_cor(data = sc1_subset,
                          sigma = sigma,
                          var_center = "centre",
                          fct_list = fct_list,
                          n_rep = 820, 
                          n_sim = 5000) #820 patients, correlated data

uncor_820 <- simul_data(data=sc1_subset,
                        n_rep = 820,
                        n_sim = 5000) #820 patients, uncorrelated data

cor_300 <- simul_data_cor(data = sc1_subset,
                          sigma = sigma,
                          var_center = "centre",
                          fct_list = fct_list,
                          n_rep = 300, 
                          n_sim = 5000)#300 patients, correlated data

uncor_300 <- simul_data(data = sc1_subset,
                          n_rep = 300, 
                          n_sim = 5000)#300 patients, uncorrelated data

save(cor_820, file = "cor_820.rda")
save(uncor_820, file = "uncor_820.rda")
save(cor_300, file = "cor_300.rda")
save(uncor_300, file = "uncor_300.rda")



















