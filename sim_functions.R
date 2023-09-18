MSB <- function(data,xi = 0.7,center_name){
  n = nrow(data)
  k = ncol(data) - 1 #number of covariates
  center_id = which(colnames(data) == center_name)
  vote_register = matrix(NA,nrow = n, ncol = k) #One column for each covariate
  treatment_allocation = rep(NA,n)
  treatment_allocation[1] = rbinom(1, size = 1, prob = 0.5) #The first patient in assigned randomly
  for(i in 2:20){ #Random allocation rule for the first 20 patients
    N1 = length(treatment_allocation[which(treatment_allocation == 1)])
    p = (20/2 -N1)/(20-(i-1))
    treatment_allocation[i] = rbinom(1, size = 1, prob = p)
  }
  for(i in 21:n){
    vote_register[i,] = vote_f(data[1:i,],treatment_allocation[1:i],center_id) #Vote function is applied sequentially to each row after the coinflips
    treatment_allocation[i] = vote_result(vote_register[i,],xi) #Vote is casted using the results of the vote function
  }
  data_result = cbind(data,vote_register,treatment_allocation)
  colnames(data_result)[(k+2):(ncol(data_result) - 1)] <- paste("vote_",colnames(data[2:(k+1)]))
  return(data_result)
}

# MSB function ------------------------------------------------------------


MSB_strat <- function(data,xi,center_name,strat_name){
  strat_id = which(colnames(data) == strat_name)
  categories = levels(data[,strat_id])
  data_result = data.frame()
  for(category in categories){
    data_result = rbind(data_result,MSB(data[which(data[,strat_id] == category),-strat_id],xi,center_name)) #A dataframe subset is created using each category of the stratified covariate, then MSB is applied on these dataframes
  }
  data_result = data_result[order(data_result[,1]),] #The id column is ordered
  data_result = cbind(data_result[,1],data[,strat_id],data_result[,2:ncol(data_result)]) #The second row of the result dataframe is set as the stratified covariate
  colnames(data_result)[1:2] = colnames(data[c(1,strat_id)]) 
  return(data_result)
}




vote_f <- function(data,treatment_allocation,center_id){ #Computes the votes for each covariate of the last individual in the dataframe
  k = ncol(data) -1 #Number of covariates
  vote_register_i = numeric(k) #Vector of votes of the ith row
  for(j in 1:k){
    data_subset = data.frame(data[,j+1],treatment_allocation) #Only the corresponding covariate and the treatment allocation column are selected
    if(is.numeric(data[,j+1])){ #Quantitative variable
      vote_register_i[j] = vote_f_quanti(data_subset) 
    }
    else if(j+1 == center_id){#Clinical center variable
      vote_register_i[j] = vote_f_center(data_subset)
    }
    else{ #Qualitative variable
      vote_register_i[j] = vote_f_quali(data_subset)
    }
  }
  return(vote_register_i)
}

statistics_quanti_f <- function(x){
  return(c(n = length(x),mean = mean(x),var = var(x))) #All the statstics needed to compute the Welch test in the two groups for quantitative variables
}


vote_f_quanti <- function(data){
  n = nrow(data)
  statistics = aggregate(data[1:n-1, 1], by = list(data[1:n-1, 2]), statistics_quanti_f) #Data aggregated by the second variable, without taking the last row, which corresponds to the current individual
  n0 = statistics[[2]][1] #sample size for treatment = 0
  n1 = statistics[[2]][2] #sample size for treatment = 1
  mean0 = statistics[[2]][3] #mean for treatment = 0
  mean1 = statistics[[2]][4] #mean for treatment = 1
  var0 = statistics[[2]][5] #variance for treatment = 0
  var1 = statistics[[2]][6] #variance for treatment = 1
  test_statistic = (mean0 - mean1)/sqrt((var0/n0) + (var1/n1)) #Welch test statistic
  df = ((var0/n0) + (var1/n1))**2 /((var0**2)/(n0**2 * (n0-1)) + (var1**2)/(n1**2 * (n1-1))) #Degrees of freedom of the law of the statistic under null hypothesis (student law)
  limit_statistic = qt(1-0.3/2,df) #statistic corresponding to a pvalue of 0.03 with this data
  if((test_statistic < -limit_statistic & data[n,1] > mean1) | (test_statistic > limit_statistic & data[n,1] < mean1)){
    result = 0 #Vote for treatment 0
  }
  else if((test_statistic < -limit_statistic & data[n,1] < mean0) | (test_statistic > limit_statistic & data[n,1] > mean0)){
    result = 1 #Vote for treatment 1
  }
  else{
    result = 2
  }
  return(result) 
  
}



vote_f_quali <- function(data){
  n = nrow(data)
  category = data[n,1]
  if(category %in% levels(droplevels(data[1:n-1,1])) & length(levels(droplevels(data[,1]))) > 1){#If this is the first time this category appears, the vote is set to neutral
    data[,1] = droplevels(data[,1]) #drop factor levels not in the current data
    t = table(data[1:n-1,1],data[1:n-1,2]) 
    khi2 = chisq.test(t) 
    test_statistic = khi2$statistic
    pvalue = khi2$p.value
    expected = khi2$expected[rownames(khi2$expected) == category,] #expected values for category
    observed = khi2$observed[rownames(khi2$observed) == category,] #observed values for category
    if(pvalue < 0.3 & expected[1] > observed[1] ){
      result = 0 #Vote for treatment 0 
    }
    else if(pvalue < 0.3 & expected[2] > observed[2]){
      result = 1 #Vote for treatment1
    }
    else{
      result = 2 #Neutral
    }
  }
  else{
    result = 2
  }
  return(result)
}

vote_f_center <- function(data){
  n = nrow(data)
  category = data[n,1]
  if(category %in% levels(droplevels(data[1:n-1,1]))){  #drop factor levels not in the current data and if there are other observed values (otherwise the khi2 test is impossible)
    #data[,1] = droplevels(data[,1])
    nc0 = length(data[which(data[1:n-1,1]== category  & data[1:n-1,2] == 0 ),1])
    nc = length(data[which(data[1:n-1,1]== category),1])
    n0 = length(data[which(data[1:n-1,2] == 0 ),1])
    n1 = length(data[which(data[1:n-1,2] == 1 ),1])
    ntot = n0 + n1
    if(nc>= 20){
      test_statistic = ((nc0/nc) - 1/2)/sqrt(1/2*(n1/ntot)*(1/nc))
      pvalue = 2 * (1 - pnorm(abs(test_statistic))) 
    }
    else{
      if((nc0/nc) < 1/2 ){
        #pvalue = 2 * sum(sapply(0:nc0, function(i) choose(nc,i) * ((n0/ntot)**i) * (n1/ntot)**(nc0-i)))
        #pvalue = sum(sapply(0:nc0, function(i) choose(nc,i) * ((n0/ntot)**i) * (n1/ntot)**(nc-i)))
        pvalue <- binom.test(nc0,nc,1/2,alternative = "two.sided" )$p.value
        #print(c((nc0/nc),(n0/ntot),pvalue,nc))
      }
      else{
        #pvalue = 2 * sum(sapply(nc0:nc, function(i) choose(nc,i) * ((n0/ntot)**i) * (n1/ntot)**(nc0-i)))
        #pvalue = sum(sapply(nc0:nc, function(i) choose(nc,i) * ((n0/ntot)**i) * (n1/ntot)**(nc-i)))
        pvalue <- binom.test(nc0,nc,1/2,alternative = "two.sided" )$p.value
        #print(c((nc0/nc),(n0/ntot),pvalue,nc))
      }
    }
    if(pvalue < 0.3 & (nc0/nc) < 1/2){
      result = 0 #Vote for treatment 0 
    }
    else if(pvalue < 0.3 & (nc0/nc) > 1/2){
      result = 1 #Vote for treatment1
    }
    else{
      result = 2 #Neutral
    }
  }
  else{
    result = 2 #Neutral if the category does not appear yet
  }
  return(result)
}

vote_result <- function(vote,xi){
  n0 = length(vote[which(vote == 0)])
  n1 = length(vote[which(vote == 1)])
  if(n1 > n0){
    result = rbinom(1, size = 1, prob = xi)
  }
  else if(n0 > n1){
    result = 1 - rbinom(1, size = 1, prob = xi)
  }
  else{
    result = rbinom(1, size = 1, prob = 0.5)
  }
  return(result)
}



rando_blocks <- function(data,block_size,strat_variables = c("sdra","centre")){
  data_strat = data.frame(sapply(strat_variables, function(name) data[,name]))
  res <- StrPBR(data_strat, bsize = block_size)
  treatment_allocation = ifelse(res$assignments == "B",1,0)
  data$treatment_allocation = treatment_allocation
  return(data)
}

rando_pocock_strat <- function(data,xi,nb_class,strat_name){
  num_index <- which(sapply(data[,-1], is.numeric)) + 1
  data$treatment_allocation = NA
  data_copy = data
  data_copy[,num_index] <- sapply(num_index, function(index) as.numeric(cut_number(data_copy[,index],nb_class )))
  res_0 <- PocSimMIN(data_copy[which(data_copy[,strat_name] == "0"),-1], p = xi)$assignments
  res_1 <- PocSimMIN(data_copy[which(data_copy[,strat_name] == "1"),-1], p = xi)$assignments
  data[which(data[,strat_name] == "0"),]$treatment_allocation = ifelse(res_0 == "B",1,0)
  data[which(data[,strat_name] == "1"),]$treatment_allocation = ifelse(res_1 == "B",1,0)
  return(data)
}


sim_outcome <- function(data,beta_v,mean_effect = -0.09461064 ,sd_effect = 0.5752965){
  data$beta_treatment = NA
  beta_treatment = -abs(rnorm(1, mean = mean_effect, sd = sd_effect))
  data$beta_treatment[1] = beta_treatment
  data$outcome = rbinom(nrow(data), size = 1, prob = inv_logit(beta_v["centre1"]*I(data$centre ==1) +  
                                                                 beta_v["centre3"]*I(data$centre ==3) +
                                                                 beta_v["centre5"]*I(data$centre ==5) +
                                                                 beta_v["centre4"]*I(data$centre ==4) +
                                                                 beta_v["centre7"]*I(data$centre ==6) +
                                                                 beta_v["centre8"]*I(data$centre ==8) +  
                                                                 beta_v["centre9"]*I(data$centre ==9) +
                                                                 beta_v["age"] * data$age +
                                                                 beta_v["immunodep1"] * I(data$immunodep == 1) +
                                                                 beta_v["vp"] * data$vp +
                                                                 beta_v["lactat"] * data$lactat +
                                                                 beta_v["sdra1"] * I(data$sdra == 1) +
                                                                 beta_v["PF"] * data$PF +
                                                                 data$beta_treatment[1] * I(data$treatment_allocation == 1)))
  return(data)
}

inv_logit <- function(x){
  return(exp(x) / (1+exp(x)))
}






