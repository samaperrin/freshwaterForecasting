

# AUC generation
if (validate_full_model == TRUE) {
  whole_model_output <- readRDS(file=paste0("./Data/whole_model_output_kfold_pooled.RDS"))} else{
    whole_model_output <- readRDS(file=paste0("./Data/whole_model_output_pseudoAbs_kfold_pooled.RDS"))
  }
# Test 1 species first

k_part <- 5

auc_full <- matrix(NA,5,5, dimnames = list(species_list))
variation_explained <- matrix(NA,5,5, dimnames = list(species_list))
variation_explained_adjust <- matrix(NA,5,5, dimnames = list(species_list))
sensitivity_table <- matrix(NA,5,5, dimnames = list(species_list))
specificity_table <- matrix(NA,5,5, dimnames = list(species_list))

for (k in 1:k_part) {
  
  
  params_test <- whole_model_output$model_ouput[[k]]
  
  beta_test <- params_test$beta
  alpha_test <- params_test$alpha
  draws_test <- params_test$draws
  
  beta_all <- matrix(apply(as.matrix(calculate(beta_test, values = draws_test)),2, median),
                     nrow=nrow(beta_test),ncol=ncol(beta_test))
  alpha_all <- matrix(apply(as.matrix(calculate(alpha_test, values = draws_test)),2, median),
                      nrow=nrow(alpha_test),ncol=ncol(alpha_test))
  
  for (s in 1:length(species_list)) {
  
      data_test <- whole_model_output$split_data[[species_list[s]]][[k]]$validation_data
    
    
    # Test to run new data
    new_data <- data_test[,c(5:7,9,11:13)]
    new_data$interaction1 = data_test$HFP * data_test$nearby_weighted
    new_data$interaction2 = data_test$upstream_presence * data_test$nearby_weighted
    new_data$interaction3 = data_test$upstream_presence * data_test$downstream_presence
    
    
    attempt <- as.matrix(new_data) %*% as.matrix(beta_all[,s])
    eta <- greta::sweep(attempt, 2, alpha_all[s], "+")
    expeta <- exp(eta)
    probabilities <- expeta/(1+expeta)
    
    actuals <- data_test$introduced
    
    threshold_sen <- mean(actuals)
    
    probabilities_for_actuals <- probabilities[actuals == 1]
    probabilities_for_actual_absences <- probabilities[actuals == 0]
    true_presences <- which(probabilities_for_actuals > threshold_sen)
    true_absences <- which(probabilities_for_actual_absences < threshold_sen)

    
    number_presences <- sum(actuals)
    number_absences <- length(actuals)-sum(actuals)
    
    roc_all <- roc(actuals, probabilities[,1])
    
    variation_explained[s,k] <- Dsquared(obs = actuals, pred = probabilities[,1], family = 'binomial')
    variation_explained_adjust[s,k] <- Dsquared(obs = actuals, pred = probabilities[,1],
                                                adjust=TRUE,family = 'binomial',npar=10)
    

    auc_full[s,k] <- auc(roc_all)
    sensitivity_table[s,k] <- length(true_presences)/length(probabilities_for_actuals)
    specificity_table[s,k] <- length(true_absences)/length(probabilities_for_actual_absences)
    
  }

  }

variation_explained_adjust <- as.data.frame(variation_explained_adjust)

variation_explained_adjust$mean <- apply(variation_explained_adjust,1,mean)
variation_explained_adjust$SD <- apply(variation_explained_adjust,1,sd)

auc_full <- as.data.frame(auc_full)

auc_full$mean <- apply(auc_full,1,mean)
auc_full$SD <- apply(auc_full,1,sd)

sensitivity_table <- as.data.frame(sensitivity_table)

sensitivity_table$mean <- apply(sensitivity_table,1,mean)
sensitivity_table$SD <- apply(sensitivity_table,1,sd)

specificity_table <- as.data.frame(specificity_table)

specificity_table$mean <- apply(specificity_table,1,mean)
specificity_table$SD <- apply(specificity_table,1,sd)

model_fit_stats <- list(AUC = auc_full, Dsquared = variation_explained_adjust, sensitivity = sensitivity_table,
                        specificity = specificity_table)
  
  
  
  
