
R <- function(var1,var2){
  ecdf1 = ecdf(as.numeric(as.character(var1)))
  ecdf2 = ecdf(as.numeric(as.character(var2)))
  cor = cor(ecdf1(as.numeric(as.character(var1))),ecdf2(as.numeric(as.character(var2))), use = "complete.obs",method = "pearson")
  return(cor)
}


table_cor <- function(data){
  m = ncol(data)
  table = matrix(rep(NA,m**2),ncol= m,nrow=m)
  for(j in 1:m){
    for(k in 1:m){
      table[j,k] = R(data[,j],data[,k])
    }
  }
  colnames(table) = colnames(data)
  rownames(table) = colnames(data)
  return(table)
}

table_pval_rank <- function(data){
  m = ncol(data)
  table = matrix(rep(NA,m**2),ncol= m,nrow=m)
  for(j in 1:m){
    for(k in 1:m){
      var1 = data[,j]
      var2 = data[,k]
      ecdf1 = ecdf(as.numeric(as.character(var1)))
      ecdf2 = ecdf(as.numeric(as.character(var2)))
      table[j,k] = cor.test(ecdf1(as.numeric(as.character(var1))),ecdf2(as.numeric(as.character(var2))), use = "complete.obs",method = "pearson")$p.value
    }
  }
  colnames(table) = colnames(data)
  rownames(table) = colnames(data)
  return(table)
}

table_pval <- function(data){
  m = ncol(data)
  table = matrix(rep(NA,m**2),ncol= m,nrow=m)
  for(j in 1:m){
    if(is.numeric(data[,j])){
      for(k in 1:m){
        if(k == j){
          pval = 0
        }
        else if(is.numeric(data[,k])){
          pval = pval_cor(data[,k],data[,j])
        }
        else{
          pval = pval_anova(data[,j],data[,k])
        }
        table[j,k] = pval
      }}
    else{
      for(k in 1:m){
        if(k == j){
          pval = 0
        }
        else if(is.numeric(data[,k])){
          pval = pval_anova(data[,k],data[,j])
        }
        else{
          pval = pval_khi2(data[,k],data[,j])
        }
        table[j,k] = pval
      }
      
    }
    
  }
  table_num = table
  colnames(table) = colnames(data)
  rownames(table) = colnames(data)
  return(list(format(round(table,3), scientific = FALSE),table_num))
}


table_pval_np <- function(data){
  m = ncol(data)
  table = matrix(rep(NA,m**2),ncol= m,nrow=m)
  for(j in 1:m){
    if(is.numeric(data[,j])){
      for(k in 1:m){
        if(k == j){
          pval = 0
        }
        else if(is.numeric(data[,k])){
          pval = pval_wilcoxon(data[,k],data[,j])
        }
        else{
          pval = pval_kruskal(data[,j],data[,k])
        }
        table[j,k] = pval
      }}
    else{
      for(k in 1:m){
        if(k == j){
          pval = 0
        }
        else if(is.numeric(data[,k])){
          pval = pval_kruskal(data[,k],data[,j])
        }
        else{
          pval = pval_khi2(data[,k],data[,j])
        }
        table[j,k] = pval
      }
      
    }
    
  }
  table_num = table
  colnames(table) = colnames(data)
  rownames(table) = colnames(data)
  return(list(format(round(table,3), scientific = FALSE),table_num))
}






cor_point = function(cor,fct1,param1,fct2,param2,n_iter = 1,n_rep = 20000){
  sigma <- matrix(c(1,cor,cor,1), nrow=2, ncol=2)
  mu <- c(0,0)
  res_cor = numeric(n_iter)
  for(rep in 1:n_iter){
    rawvars <- mvrnorm(n=n_rep, mu=mu, Sigma=sigma)
    param1$p = pnorm(rawvars[,1])
    param2$p = pnorm(rawvars[,2])
    sample1 <- do.call(fct1,param1)
    sample2 <- do.call(fct2,param2)
    res_cor[rep] <- R(sample1,sample2)
  }
  return(mean(res_cor))
}

dist_cor <- function(cor,wanted_cor,fct1,param1,fct2,param2,n_iter = 3000,n_rep = 300){
  res = cor_point(cor,fct1,param1,fct2,param2,n_iter,n_rep)
  return(abs(res-wanted_cor))
}


find_cor_sim <- function(var_list,fct_list,param_list,cor_matrix,n_iter = 3500,n_rep = 300){
   res_matrix = matrix(0,ncol = ncol(cor_matrix),nrow = nrow(cor_matrix))
   colnames(res_matrix) = colnames(cor_matrix)
   rownames(res_matrix) = rownames(cor_matrix)
   for(i in 1:length(var_list)){
     varname1 = var_list[[i]][1]
     varname2 = var_list[[i]][2]
     wanted_cor = cor_matrix[which(rownames(cor_matrix) == varname1),which(colnames(cor_matrix) == varname2)]
     if(wanted_cor < 0){
       inf = -1
       sup = 0
     }
     else{
       inf = 0
       sup = 1
     }
     param1 = param_list[[varname1]]
     param2 = param_list[[varname2]]
     fct1 = fct_list[[varname1]]
     fct2 = fct_list[[varname2]]
     best_cor = optimise(dist_cor,wanted_cor = wanted_cor,interval = c(inf,sup),fct1 = fct1,param1 = param1,fct2 = fct2,param2 = param2,n_iter = n_iter,n_rep = n_rep)$minimum
     res_matrix[which(rownames(res_matrix) == varname1),which(colnames(res_matrix) == varname2)] = best_cor
     res_matrix[which(rownames(res_matrix) == varname2),which(colnames(res_matrix) == varname1)] = best_cor
    }
    res_matrix = res_matrix + diag(1,nrow=nrow(sigma),ncol=ncol(sigma))
    return(res_matrix)
 }



g_bern <- function(rho, p){
  integral_0 = integrate(f_integral,lower = 0, upper = 1,abs.tol = 10**-15, rho = rho, l = 0, p = p)$value
  integral_1 = integrate(f_integral,lower = 0, upper = 1,abs.tol = 10**-15, rho = rho, l = 1, p = p)$value
  res = pbinom(0,size = 1, prob = p)*integral_0 + pbinom(1,size = 1, prob = p)*integral_1
  return(res)
}

f_integral <- function(u,rho,l,p){
  z_l0 = qnorm(pbinom(l,size = 1, prob = p))
  z_l1 = qnorm(pbinom(l-1,size = 1, prob = p))
  res = u * (pnorm((z_l0 - rho * qnorm(u)) / sqrt(1 - rho**2)) - pnorm((z_l1 - rho * qnorm(u)) / sqrt(1-rho**2)))
  return(res)
}


g_discrete <- function(rho,p1,p2){
  z_1_m1 = qnorm(pbinom(-1,size = 1, prob = p1))
  z_1_0 = qnorm(pbinom(0,size = 1, prob = p1))
  z_2_m1 = qnorm(pbinom(-1,size = 1, prob = p2))
  z_2_0 = qnorm(pbinom(0,size = 1, prob = p2))
  res <- (1-p1)*(1-p2)*pbinorm(-z_1_m1,-z_2_m1,cov12 = rho) + 
    p1*p2*pbinorm(-z_1_0,-z_2_0,cov12 = rho) +
    (1-p1)*p2*pbinorm(-z_1_m1,-z_2_0,cov12 = rho) +
    p1*(1-p2)*pbinorm(-z_1_0,-z_2_m1,cov12 = rho)
  return(res)
}


esp_sdra <- function(z){
  pbinom(qbinom(p= pnorm(z),size = 1,prob = 0.5048544) ,size = 1, prob = 0.5048544)* dnorm(z)
}

esp_squared_sdra <- function(z){
  ((pbinom(qbinom(p= pnorm(z),size = 1,prob = 0.5048544) ,size = 1, prob = 0.5048544))**2) * dnorm(z)
}

esp_immunodep <- function(z){
  pbinom(qbinom(p= pnorm(z),size = 1,prob = 0.2716049) ,size = 1, prob = 0.2716049)* dnorm(z)
}

esp_squared_immunodep <- function(z){
  ((pbinom(qbinom(p= pnorm(z),size = 1,prob = 0.2716049) ,size = 1, prob = 0.2716049))**2) * dnorm(z)
}

esp_quanti <- function(z){
  pnorm(z) * dnorm(z)
}

esp_squared_quanti <- function(z){
  pnorm(z)**2 * dnorm(z)
}

objective_f <- function(rho,wanted_cor,prob,mean_cdf_1,mean_cdf_2,sd_cdf_1,sd_cdf_2){
  return(g_bern(rho,prob) - mean_cdf_1*mean_cdf_2 - wanted_cor*sd_cdf_1*sd_cdf_2)
}

objective_f_discrete <- function(rho,wanted_cor,prob1,prob2,mean_cdf_1,mean_cdf_2,sd_cdf_1,sd_cdf_2){
  return(g_discrete(rho = rho,p1 = prob1, p2 = prob2) - mean_cdf_1*mean_cdf_2 - wanted_cor*sd_cdf_1*sd_cdf_2)
}



find_cor <- function(data,cor_matrix,var_list,param_list){
  res_matrix = matrix(0,ncol = ncol(cor_matrix),nrow = nrow(cor_matrix))
  colnames(res_matrix) = colnames(cor_matrix)
  rownames(res_matrix) = rownames(cor_matrix)
  for(i in 1:length(var_list)){
    varname1 = var_list[[i]][1]
    varname2 = var_list[[i]][2]
    wanted_cor = cor_matrix[which(rownames(cor_matrix) == varname1),which(colnames(cor_matrix) == varname2)]
    if(wanted_cor < 0){
      inf = -1
      sup = 0
    }
    else{
      inf = 0
      sup = 1
    }
    if(is.factor(data[,varname1]) & is.numeric(data[,varname2])){
      res_cor <- brentDekker(objective_f,
                             a = inf,
                             b= sup,
                             tol = 1e-6,
                             wanted_cor = wanted_cor,
                             prob = param_list[[varname1]]$prob,
                             mean_cdf_1 = param_list[[varname1]]$mean_cdf,
                             sd_cdf_1 = param_list[[varname1]]$sd_cdf,
                             mean_cdf_2 = param_list[["quanti"]]$mean_cdf,
                             sd_cdf_2 = param_list[["quanti"]]$sd_cdf)$root
    }
    if(is.factor(data[,varname2]) & is.numeric(data[,varname1])){
      res_cor <- brentDekker(objective_f,
                             a = inf,
                             b= sup,
                             tol = 1e-6,
                             wanted_cor = wanted_cor,
                             prob = param_list[[varname2]]$prob,
                             mean_cdf_1 = param_list[[varname2]]$mean_cdf,
                             sd_cdf_1 = param_list[[varname2]]$sd_cdf,
                             mean_cdf_2 = param_list[["quanti"]]$mean_cdf,
                             sd_cdf_2 = param_list[["quanti"]]$sd_cdf)$root
      
    }
    if(is.factor(data[,varname1]) & is.factor(data[,varname2])){
      res_cor <- brentDekker(objective_f_discrete,
                             a = inf,
                             b= sup,
                             tol = 1e-6,
                             wanted_cor = wanted_cor,
                             prob1 = param_list[[varname1]]$prob,
                             prob2 = param_list[[varname2]]$prob,
                             mean_cdf_1 = param_list[[varname1]]$mean_cdf,
                             sd_cdf_1 = param_list[[varname1]]$sd_cdf,
                             mean_cdf_2 = param_list[[varname2]]$mean_cdf,
                             sd_cdf_2 = param_list[[varname2]]$sd_cdf)$root
    }
    if(is.numeric(data[,varname1]) & is.numeric(data[,varname2])){
      res_cor = 2*sin((wanted_cor * pi)/6)
    }
    res_matrix[which(rownames(res_matrix) == varname1),which(colnames(res_matrix) == varname2)] = res_cor
    res_matrix[which(rownames(res_matrix) == varname2),which(colnames(res_matrix) == varname1)] = res_cor
  }
  res_matrix = res_matrix + diag(1,nrow=nrow(res_matrix),ncol=ncol(res_matrix))
  return(res_matrix)
}


noise_cov <- function(sigma,var_list,eps){
  for(tuple in var_list){
    varname1 = tuple[[1]]
    varname2 = tuple[[2]]
    sigma[varname1,varname2] = sigma[varname1,varname2] + rnorm(n = 1, sd = eps)
    sigma[varname2,varname1] = sigma[varname1,varname2]
  }
  return(sigma)
}


simul_data_cor <- function(data,sigma,var_center,fct_list,n_sim,n_rep ){
  prop_center = prop.table(table(data[,var_center]))
  result = vector("list",n_sim)
  mu <- numeric(ncol(sigma))
  factor_index <- which(sapply(data, is.factor))
  num_index <- which(sapply(data, is.numeric))
  decimals = sapply(data,n_decimals)
  i = 1
  while(i <= n_sim){
    sigma_noise <- noise_cov(sigma,var_list,eps = 0.025) #Noise is added to the objective correlation matrix because f behaves linearly far from the bounds
    if(is_positive_definite(sigma_noise)){
      data_cor = data.frame(matrix(ncol=ncol(data),nrow = n_rep))
      colnames(data_cor) = colnames(data)
      data_cor[,1] = 1:n_rep
      data_cor[,2] = sample(as.numeric(names(prop_center)),
                            n_rep,
                            prob = prop_center,
                            replace = TRUE)
      rawvars <- mvrnorm(n=n_rep, mu=mu, Sigma=sigma_noise)
      colnames(rawvars) = colnames(sigma)
      for(varname in colnames(sigma)){
        fct_name = fct_list[[varname]]
        if(fct_name == "quantile"){
          param = list("x" = data[,varname]
                       ,"na.rm" = TRUE
                       , "p" = pnorm(rawvars[,varname] ))
          data_cor[,varname] = do.call(quantile,param)
        }
        if(fct_name == "qbinom"){
          prob = mean(as.numeric(as.character(data[,varname])),na.rm = TRUE)
          param = list("size" = 1,
                       "prob" = prob, 
                       "p" =  pnorm(rawvars[,varname] ))
          data_cor[,varname] = do.call(qbinom,param)
        }
      }
      n_pval = n_wilcox_pvals(data_cor[,3:8],data[,colnames(sigma)],threshold = 10**-4)
      if(n_pval < 1){
        for(index in factor_index){
          data_cor[,index] <-  as.factor(data_cor[,index])
        }
        #data_cor[,3:ncol(data_cor)] <- noise(data_cor[,3:ncol(data_cor)])
        data_cor[,num_index] <- sapply(num_index,function(index) round(data_cor[,index],decimals[index]))
        result[[i]] <- data_cor
        i = i+1
      }
    }
  }
  return(result)
}

simul_data <- function(data,n_sim,n_rep ){
  result = vector("list",n_sim)
  factor_index <- which(sapply(data, is.factor))
  num_index <- which(sapply(data, is.numeric))
  decimals = sapply(data,n_decimals)
  i = 1
  while(i <= n_sim){
    data_sim = data.frame(matrix(ncol=ncol(data),nrow = n_rep))
    colnames(data_sim) = colnames(data)
    data_sim[,1] = 1:n_rep
    for(varname in colnames(data)[2:ncol(data)]){
      data_sim[,varname] = sample(as.numeric(as.character(data[!is.na(data[,varname]),varname])),n_rep,replace = TRUE)
    }
    for(index in factor_index){
      data_sim[,index] <-  as.factor(data_sim[,index])
    }
    #data_sim[,3:ncol(data_sim)] <- noise(data_sim[,3:ncol(data_sim)])
    n_pval = n_wilcox_pvals(data_sim[,3:8],data[,3:8],threshold = 10**-4)
    if(n_pval < 1){
      data_sim[,num_index] <- sapply(num_index,function(index) round(data_sim[,index],decimals[index]))
      result[[i]] <- data_sim
      i = i+1
    }
  }
  return(result)
}

n_wilcox_pvals <- function(data1,data2,threshold){
  n_pvals = 0
  for(varname in colnames(data1)){
    pval <- wilcox.test(as.numeric(as.character(data1[,varname])),as.numeric(as.character(data2[,varname])))$p.value
    if(pval < threshold){
      n_pvals = n_pvals + 1
    }
  }
  return(n_pvals)
}


# Function to count the number of decimals in a numeric variable
n_decimals <- function(var) {
  # Remove NA values from the variable
  var_wo_NA = var[!is.na(var)]
  
  # Extract the first non-NA value
  number = var_wo_NA[1]
  
  # Convert the number to a character string
  number_string <- as.character(number)
  
  # Split the number string into integer and decimal parts
  parts <- strsplit(number_string, ".", fixed = TRUE)[[1]]
  
  # Check if a decimal part exists
  if (length(parts) > 1) {
    decimal_part <- parts[2]
    
    # Calculate the number of decimal places
    number_of_decimals <- nchar(decimal_part)
    
    return(number_of_decimals)
  } else {
    return(0)
  }
}




name_format <- function(name_vect){
  name_vect <- data.frame("name"= name_vect)
  name_vect <- name_vect %>%
    mutate(name = case_when(
      name == "age" ~ "Age",
      name == "sdra" ~ "SDRA",
      name =="PF" ~ "Pression partielle O2",
      name == "immunodep" ~ "Immunodéficience",
      name == "vp" ~ "Dose vasopresseurs",
      name == "lactat" ~ "Lactate"
    ))
  return(as.vector(name_vect$name))
}






is_positive_definite <- function(matrix_input) {
  # Check if the input is a square matrix
  if (!is.matrix(matrix_input) || nrow(matrix_input) != ncol(matrix_input)) {
    return(FALSE)
  }
  
  # Check if the matrix is symmetric
  if (!identical(matrix_input, t(matrix_input))) {
    return(FALSE)
  }
  
  # Calculate the eigenvalues of the matrix
  eigenvalues <- eigen(matrix_input)$values
  
  # Check if all eigenvalues are strictly positive
  if (all(eigenvalues > 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



