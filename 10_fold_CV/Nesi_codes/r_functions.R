### function to make data splits based on site
split_func  <- function(fold_input = fold_names[1],data_input)
{ 
  train_indices <- which(data_input$fold != fold_input )
  test_indices <- which(data_input$fold == fold_input)
  
  indices <-
    list(analysis   = train_indices, 
         assessment = test_indices)
  split <- make_splits(indices, data_input)
  return(split)}
## scale train and test seperately and apply the mean and sd of the baseline  to the followup

scale_train <- function(baseline_data,followup_data){
  ### select features
  baseline_data_tibble <- baseline_data %>% select_if(is.numeric)
  followup_data_tibble <- followup_data %>% select_if(is.numeric)
  
  mean_baseline_data <- baseline_data_tibble %>% summarise(across(everything(), mean))
  sd_baseline_data <- baseline_data_tibble %>% summarise(across(everything(), sd))
  ### find out which of the features has 0 sd
  if(length(which(sd_baseline_data==0)) >0){
    zero_sd_vec <- which(sd_baseline_data==0)
    ### only scale the variable with none zero sd
    scaled_baseline_data <- baseline_data_tibble[-zero_sd_vec] %>% scale()
    scaled_baseline_data_all <- bind_cols(scaled_baseline_data,baseline_data_tibble[zero_sd_vec])
    
    
    # Standardize the test sample
    scaled_followup_data <- followup_data_tibble[-zero_sd_vec] %>%
      rowwise() %>%
      mutate((across(everything())-mean_baseline_data[-zero_sd_vec])/sd_baseline_data[-zero_sd_vec]) %>%
      ungroup()
    scaled_followup_data_all <- bind_cols(scaled_followup_data,followup_data_tibble[zero_sd_vec])
  }
  else{
    {
      ### only scale the variable with none zero sd
      scaled_baseline_data <- baseline_data_tibble %>% scale()
      # Standardize the test sample
      scaled_followup_data <- followup_data_tibble %>%
        rowwise() %>%
        mutate((across(everything())-mean_baseline_data)/sd_baseline_data) %>%
        ungroup()
      scaled_followup_data_all <- bind_cols(scaled_followup_data,followup_data_tibble)
    }
  }
  output_baseline_data <- baseline_data %>% dplyr::select(where(is.character))%>%
    bind_cols(scaled_baseline_data)
  output_followup_data <- followup_data %>% dplyr::select(where(is.character))%>%
    bind_cols(scaled_followup_data)
  return(list(output_baseline_data = output_baseline_data,
              output_followup_data=output_followup_data))
}

scale_seperate <- function(baseline_data,followup_data){
  ### select features
  baseline_data_tibble <- baseline_data %>% select_if(is.numeric)
  followup_data_tibble <- followup_data %>% select_if(is.numeric)
  
  scaled_baseline_data <- baseline_data_tibble %>% scale()
  
  # Standardize the test sample
  scaled_followup_data <- followup_data_tibble %>% scale()
  
  output_baseline_data <- baseline_data %>% dplyr::select(where(is.character))%>%
    bind_cols(scaled_baseline_data)
  output_followup_data <- followup_data %>% dplyr::select(where(is.character))%>%
    bind_cols(scaled_followup_data)
  return(list(output_baseline_data = output_baseline_data,
              output_followup_data=output_followup_data))
}

## compute and processing gractor
gfactor_cross_sites_seperate <- function(split_input){
  ### select the variables
  
  train_baseline <- training(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  test_baseline <- testing(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  train_followup <- training(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  test_followup <- testing(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na() 
  
  scaled_train_list <- scale_train(baseline_data = train_baseline,
                                   followup_data = train_followup)
  scaled_test_list <- scale_train(baseline_data = test_baseline,
                                  followup_data = test_followup)
  
  
  NeuroCog2ndOrder.Fit <- lavaan::cfa(model = NeuroCog2ndOrder, 
                                      data = scaled_train_list[["output_baseline_data"]],estimator="MLR")
  
  second_order_train_baseline <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                    newdata = scaled_train_list[["output_baseline_data"]])
  
  second_order_test_baseline <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                   newdata = scaled_test_list[["output_baseline_data"]])
  
  
  second_order_train_followup <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                    newdata = scaled_train_list[["output_followup_data"]])
  
  second_order_test_followup <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                   newdata = scaled_test_list[["output_followup_data"]])
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
  
  second_order_pred <- list(output_train_baseline = second_order_train_baseline, 
                            output_test_baseline = second_order_test_baseline,
                            output_train_followup=second_order_train_followup,
                            output_test_followup=second_order_test_followup)
  
  ###use this gfactor as response for all the other combat models
  ### response do not needs to be combat
  gfactor_data_list <- map2(.x =scaled_data_list,.y = second_order_pred, ~mutate(.x,gfactor = .y[,4])%>%drop_na())
  
  ### get the gfactor for the final output
  gfactor_list <- gfactor_data_list %>% map(.,~dplyr::select(.,all_of(c("SRC_SUBJECT_ID","gfactor"))))
  
  
  scaled_train_gfactor <- scale_train(baseline_data = gfactor_list[["output_train_baseline"]],
                                      followup_data = gfactor_list[["output_train_followup"]])
  scaled_test_gfactor <- scale_train(baseline_data = gfactor_list[["output_test_baseline"]],
                                     followup_data = gfactor_list[["output_test_followup"]])
  
  
  scaled_gfactor_list <- list(output_train_baseline = scaled_train_gfactor[["output_baseline_data"]], 
                              output_test_baseline = scaled_test_gfactor[["output_baseline_data"]],
                              output_train_followup=scaled_train_gfactor[["output_followup_data"]],
                              output_test_followup=scaled_test_gfactor[["output_followup_data"]])
  
  return(scaled_gfactor_list)
}


### scale gfactor completely between every data set

gfactor_cross_sites_individual <- function(split_input){
  ### select the variables
  
  train_baseline <- training(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  test_baseline <- testing(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  train_followup <- training(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  
  test_followup <- testing(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na() 
  
  scaled_train_list <- scale_seperate(baseline_data = train_baseline,
                                   followup_data = train_followup)
  scaled_test_list <- scale_seperate(baseline_data = test_baseline,
                                  followup_data = test_followup)
  
  
  NeuroCog2ndOrder.Fit <- lavaan::cfa(model = NeuroCog2ndOrder, 
                                      data = scaled_train_list[["output_baseline_data"]],estimator="MLR")
  
  second_order_train_baseline <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                    newdata = scaled_train_list[["output_baseline_data"]])
  
  second_order_test_baseline <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                   newdata = scaled_test_list[["output_baseline_data"]])
  
  
  second_order_train_followup <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                    newdata = scaled_train_list[["output_followup_data"]])
  
  second_order_test_followup <- lavaan::lavPredict(NeuroCog2ndOrder.Fit, 
                                                   newdata = scaled_test_list[["output_followup_data"]])
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
  
  second_order_pred <- list(output_train_baseline = second_order_train_baseline, 
                            output_test_baseline = second_order_test_baseline,
                            output_train_followup=second_order_train_followup,
                            output_test_followup=second_order_test_followup)
  
  ###use this gfactor as response for all the other combat models
  ### response do not needs to be combat
  gfactor_data_list <- map2(.x =scaled_data_list,.y = second_order_pred, ~mutate(.x,gfactor = .y[,4])%>%drop_na())
  
  ### get the gfactor for the final output
  gfactor_list <- gfactor_data_list %>% map(.,~dplyr::select(.,all_of(c("SRC_SUBJECT_ID","gfactor"))))
  
  
  scaled_train_gfactor <- scale_seperate(baseline_data = gfactor_list[["output_train_baseline"]],
                                      followup_data = gfactor_list[["output_train_followup"]])
  scaled_test_gfactor <- scale_seperate(baseline_data = gfactor_list[["output_test_baseline"]],
                                     followup_data = gfactor_list[["output_test_followup"]])
  
  
  scaled_gfactor_list <- list(output_train_baseline = scaled_train_gfactor[["output_baseline_data"]], 
                              output_test_baseline = scaled_test_gfactor[["output_baseline_data"]],
                              output_train_followup=scaled_train_gfactor[["output_followup_data"]],
                              output_test_followup=scaled_test_gfactor[["output_followup_data"]])
  
  return(scaled_gfactor_list)
}




IQR_remove_vec <- function(data_split,x){
  outlier_tibble <- data_split%>%
    dplyr::select(all_of(x))%>%
    mutate_at(vars(all_of(x)), ~ ifelse(
      .x > quantile(.x, na.rm = TRUE)[4] + 3 * IQR(.x, na.rm = TRUE) |
        .x < quantile(.x, na.rm = TRUE)[2] - 3 * IQR(.x, na.rm = TRUE), 
      TRUE, FALSE))
  data_split <- data_split %>% mutate(outlier_indicator = rowMeans(outlier_tibble))
  data_split_outlier <- filter(data_split,outlier_indicator > 0.05) ##threshold to remove outliers
  data_split_select <- data_split %>% filter(! SRC_SUBJECT_ID %in% data_split_outlier$SRC_SUBJECT_ID)%>%
    dplyr::select(-outlier_indicator)
  return(data_split_select)
}


IQR_remove_num <- function(data_split,x){
  numeric_features <- data_split %>% select_if(is.numeric)%>% colnames()

  outlier_tibble <- data_split%>%
    dplyr::select(all_of(numeric_features))%>%
    mutate_if(is.numeric, ~ ifelse(
      .x > quantile(.x, na.rm = TRUE)[4] + 3 * IQR(.x, na.rm = TRUE) |
        .x < quantile(.x, na.rm = TRUE)[2] - 3 * IQR(.x, na.rm = TRUE), 
      TRUE, FALSE))
  data_split <- data_split %>% mutate(outlier_indicator = rowMeans(outlier_tibble))
  data_split_outlier <- filter(data_split,outlier_indicator > 0.05) ##threshold to remove outliers
  data_split_select <- data_split %>% filter(! SRC_SUBJECT_ID %in% data_split_outlier$SRC_SUBJECT_ID)%>%
    dplyr::select(-outlier_indicator)
  return(data_split_select)
}

batch_adjustOne_vec <- function(data_fold, x){
  ols_data <- dplyr::select(data_fold,all_of(x))
  ols_matrix <- as.matrix(t(ols_data))
  dimnames(ols_matrix)<- NULL
  data_fold$SITE_ID_L <- as.factor(data_fold$SITE_ID_L)##have to be changed into factor before combat
  data_fold$SITE_ID_L <- droplevels(data_fold$SITE_ID_L)##drop the empty levels of a factor
  ols_matrix_com <- ComBat(dat = ols_matrix,batch = data_fold$SITE_ID_L)
  ols_data_com <- data.frame(t(ols_matrix_com))
  #%>%
  # scale()### scale the observations after combat
  colnames(ols_data_com) <- colnames(ols_data)
  data_resp <- data_fold%>%dplyr::select(-starts_with(x))
  data_output <- bind_cols(data_resp,ols_data_com)
  
  return(data_output)
}


batch_adjustOne_ref_vec <- function(data_fold, x,test_data_fold){
  ### create a new fold variable and combine train and test together
  data_fold <- data_fold %>% mutate(fold = "train")
  test_data_fold <- test_data_fold %>% mutate(fold = "test")
  data_fold_all <- bind_rows(data_fold,test_data_fold)
  
  ols_data <- dplyr::select(data_fold_all,all_of(x))
  ols_matrix <- as.matrix(t(ols_data))
  dimnames(ols_matrix)<- NULL
  data_fold_all$fold <- as.factor(data_fold_all$fold)##have to be changed into factor before combat
  data_fold_all$fold <- droplevels(data_fold_all$fold)##drop the empty levels of a factor
  
  ols_matrix_com <- ComBat(dat = ols_matrix,batch = data_fold_all$fold, ref.batch = "train")
  ols_data_com <- data.frame(t(ols_matrix_com))
  #%>%
  # scale()### scale the observations after combat
  colnames(ols_data_com) <- colnames(ols_data)
  data_resp <- data_fold_all%>%dplyr::select(-starts_with(x))
  data_output <- bind_cols(data_resp,ols_data_com)
  
  data_output_test <- filter(data_output, fold == "test")%>%
    dplyr::select(-fold)
  
  return(data_output_test)
}



### only processing the features, gfactor, the respone, is computed in the other function
data_processing_cross_sites_seperate <- function(split_input,features_input= features){
  ### select the variables
  
  train_baseline <- training(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_vec(x = features_input)
  
  test_baseline <- testing(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_vec(x = features_input)
  
  train_followup <- training(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_vec(x = features_input)
  
  test_followup <- testing(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_vec(x = features_input)
  
  scaled_train_list <- scale_train(baseline_data = train_baseline,
                                   followup_data = train_followup)
  scaled_test_list <- scale_train(baseline_data = test_baseline,
                                  followup_data = test_followup)
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
  
  
  ### scale and combat features
  
  
  ### scale and combat
  ## only select features and subject ID
  
  scaled_data_list_select <-scaled_data_list%>% map(.,~                                                     
                                                      dplyr::select(.,all_of(c("SRC_SUBJECT_ID","SITE_ID_L")),all_of(features_input))) 
  
  train_baseline_combat_train <- batch_adjustOne_vec(data_fold = scaled_data_list_select[["output_train_baseline"]],
                                                     x = features_input) %>% 
                                 mutate(EVENTNAME ="baseline_year_1_arm_1")%>% 
    mutate(fold = "train")
  
  
  train_followup_combat_train <- batch_adjustOne_vec(data_fold = scaled_data_list_select[["output_train_followup"]],
                                                     x = features_input) %>% 
        mutate(EVENTNAME="2_year_follow_up_y_arm_1")%>% 
    mutate(fold = "train")
  
  
  test_baseline_combat_test <- batch_adjustOne_vec(data_fold = scaled_data_list_select[["output_test_baseline"]],
                                                     x = features_input)%>% 
    mutate(EVENTNAME ="baseline_year_1_arm_1")%>% 
    mutate(fold = "test")
  
  
  test_followup_combat_test <- batch_adjustOne_vec(data_fold = scaled_data_list_select[["output_test_followup"]],
                                                     x = features_input) %>% 
    mutate(EVENTNAME="2_year_follow_up_y_arm_1")%>% 
    mutate(fold = "test")

  
  output_data <- bind_rows(train_baseline_combat_train,train_followup_combat_train,test_baseline_combat_test,test_followup_combat_test)
  return(output_data)
}



### only processing the features, gfactor, the respone, is computed in the other function
data_processing_cross_sites_seperate_no_combat <- function(split_input,features_input=features){
  ### select the variables
  
  train_baseline <- training(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_num(x = features_input)
  
  test_baseline <- testing(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_num(x = features_input)
  
  train_followup <- training(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_num(x = features_input)
  
  test_followup <- testing(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()%>%
    IQR_remove_num(x = features_input)
  
  scaled_train_list <- scale_train(baseline_data = train_baseline,
                                   followup_data = train_followup)
  scaled_test_list <- scale_train(baseline_data = test_baseline,
                                  followup_data = test_followup)
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
  
  
  output_data_train_baseline <- scaled_data_list[["output_train_baseline"]]%>% 
    mutate(fold = "train")
  
  output_data_test_baseline <- scaled_data_list[["output_test_baseline"]]%>% 
    mutate(fold = "test")
  
  output_data_train_followup <- scaled_data_list[["output_train_followup"]]%>% 
    mutate(fold = "train")
  
  output_data_test_followup <- scaled_data_list[["output_test_followup"]]%>% 
    mutate(fold = "test")
  
  output_data <- bind_rows(output_data_train_baseline,output_data_test_baseline,output_data_train_followup,output_data_test_followup)
  return(output_data)
}
### used to process social econumic status and psychopathology models

### only processing the features, gfactor, the respone, is computed in the other function
data_processing_cross_sites_seperate_no_combat_no_iqr <- function(split_input,features_input=features){
  ### select the variables
  
  train_baseline <- training(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  #%>%
  #  IQR_remove_num(x = features_input)
  
  test_baseline <- testing(split_input)%>%
    filter(EVENTNAME=="baseline_year_1_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  #%>%
  #  IQR_remove_num(x = features_input)
  
  train_followup <- training(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  #%>%
  #  IQR_remove_num(x = features_input)
  
  test_followup <- testing(split_input)%>%
    filter(EVENTNAME=="2_year_follow_up_y_arm_1") %>%
    dplyr::select(all_of(features_input),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
  #%>%
  #  IQR_remove_num(x = features_input)
  
  scaled_train_list <- scale_train(baseline_data = train_baseline,
                                   followup_data = train_followup)
  scaled_test_list <- scale_train(baseline_data = test_baseline,
                                  followup_data = test_followup)
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
  
  
  output_data_train_baseline <- scaled_data_list[["output_train_baseline"]]%>% 
    mutate(fold = "train")
  
  output_data_test_baseline <- scaled_data_list[["output_test_baseline"]]%>% 
    mutate(fold = "test")
  
  output_data_train_followup <- scaled_data_list[["output_train_followup"]]%>% 
    mutate(fold = "train")
  
  output_data_test_followup <- scaled_data_list[["output_test_followup"]]%>% 
    mutate(fold = "test")
  
  output_data <- bind_rows(output_data_train_baseline,output_data_test_baseline,output_data_train_followup,output_data_test_followup)
  return(output_data)
}


recipe_prep <- function(train_input=train_gfactor_scan_enet,features_input=features){
  norm_recipe <- recipe( as.formula("gfactor~."), data = train_input) %>%
    update_role(all_of(features_input), new_role = "predictor")%>%
    update_role("gfactor", new_role = "outcome" )%>%
    step_dummy(all_nominal()) 
    preped_recipe <- norm_recipe%>%
    prep(training = train_input, retain = TRUE)
  return(list(norm_recipe= norm_recipe,preped_recipe=preped_recipe))
}

recipe_prep_scale <- function(train_input=train_gfactor_scan_enet,features_input=features){
  norm_recipe <- recipe( as.formula("gfactor~."), data = train_input) %>%
    update_role(all_of(features_input), new_role = "predictor")%>%
    update_role("gfactor", new_role = "outcome" )%>%
    step_normalize(all_numeric()) 
    preped_recipe <- norm_recipe %>%
    prep(training = train_input, retain = TRUE)
  return(list(norm_recipe=norm_recipe,preped_recipe=preped_recipe))
}
### fit the elastic net model 
### fit the elastic net model 
enet_tuning <- function(recipe_input){
  set.seed(123) 
  train_input <- recipe_input[["preped_recipe"]] %>% bake(new_data=NULL)
  tuning_cv_folds <- train_input  %>%
    vfold_cv(v = 10)
  
  ## mtry is the number of predictors to sample at each split
  ## min_n (the number of observations needed to keep splitting nodes)
  model_tune <-linear_reg(penalty =tune(),  
                          mixture = tune()) %>%
    set_mode("regression") %>%
    set_engine("glmnet")
  
  tune_wf <- workflow() %>%
    add_recipe(recipe_input[["norm_recipe"]]) %>% 
    add_model(model_tune)
  
  ## automate generate grid for hyperparameters
  
  model_grid <- 
    model_tune %>% 
    extract_parameter_set_dials(tune_wf)%>% 
    ## update the values of the parameters
    update(penalty =penalty(range = c(-10,1), trans = log10_trans()))%>%
    update(mixture =mixture(range = c(0,1)))%>%
    grid_regular(levels = c(200,11))
  
  tune_ctrl <- control_grid(save_pred = TRUE, verbose = TRUE
                            ,parallel_over = "everything"
  )
  
  
  #start <- Sys.time()
  tune_res <- tune_grid(
    tune_wf,
    resamples = tuning_cv_folds,
    metrics = metric_set(rmse),
    grid = model_grid,
    control= tune_ctrl
  )
  
  best_tune <- select_best(tune_res, 
                           metric = "rmse")
  
  best_tuned_param <- show_best(tune_res, 
                                metric="rmse")
  
  enet_final_wf <- tune_wf %>% finalize_workflow(best_tune)
  return(list(enet_wf_final = enet_final_wf, 
              best_enet_model = best_tune,
              best_enet_forest_param = best_tuned_param))
}



### predict the enet model
model_final_fit <- function(recipe_input,
                            wf_input,
                            test_data){
  train_input <- recipe_input[["preped_recipe"]] %>% 
    bake(new_data=NULL) %>% drop_na()
  
  ##baked recipe scale the test data with the mean and sd in the training data
  test_input <-  bake( recipe_input[["preped_recipe"]],
                       new_data=test_data) %>% drop_na()
  
  model_final_fit <- 
    wf_input%>%
    parsnip::extract_spec_parsnip()%>%
    parsnip::fit(data = train_input, formula= as.formula("gfactor~."))
  
  model_predict <- predict(model_final_fit, 
                           new_data = test_input %>% 
                             drop_na() ) %>%
    rename(model_predict = .pred) %>% 
    bind_cols(test_input%>% drop_na())  
  
  ##processing output
  
  output_list <- vector("list",length=2)
  names(output_list) <- c(paste0("model","_final_fit"),
                          paste0("model","_predict"))
  
  output_list[[paste0("model","_final_fit")]] <- model_final_fit
  output_list[[paste0("model","_predict")]] <- model_predict
  return(output_list)
}
### random forest tuning

random_forest_tuning <- function(recipe_input){

  train_input <- recipe_input[["preped_recipe"]] %>% bake(new_data=NULL)
  tuning_cv_folds <- train_input  %>%
    vfold_cv(v = 10)

  
  ## mtry is the number of predictors to sample at each split
  ## min_n (the number of observations needed to keep splitting nodes)
  model_tune <-rand_forest(mtry = tune(),
                           trees = 500,
                           min_n = tune()) %>%
    set_mode("regression") %>%
    set_engine("ranger")
  
  tune_wf <- workflow() %>%
    add_recipe(recipe_input) %>%
    add_model(model_tune)
  
  
  ## automate generate grid for hyperparameters
  
  model_grid <- 
    model_tune %>% 
    extract_parameter_set_dials(tune_wf)%>% 
    ## update the values of the parameters
    update(min_n =min_n(range = c(2,2000)))%>%
    update(mtry =mtry(range = c(1, 90)))%>%
    grid_latin_hypercube(size = 3000)
  
  
  tune_ctrl <- control_grid(save_pred = TRUE, verbose = TRUE
                            ,parallel_over = "everything"
  )
  #library(doFuture)
  #registerDoFuture()
  #plan(multisession(workers = 30))
  start <- Sys.time()
  
  #library(doParallel)
  #registerDoParallel(strtoi(Sys.getenv('SLURM_CPUS_PER_TASK')))    
  
  #doParallel::registerDoParallel(cores = 20)
  
  tune_res <- tune_grid(
    tune_wf,
    resamples = tuning_cv_folds,
    metrics = metric_set(rmse),
    grid = model_grid,
    control= tune_ctrl
  )
  end <- Sys.time()
  
  ##check how much of the grid it went through
  #tune_res %>% 
  #  collect_metrics()%>%
  #    print()
  
  
  
  best_tune <- select_best(tune_res, 
                           metric = "rmse")
  
  best_tuned_param <- show_best(tune_res, 
                                metric="rmse")
  
  rf_final_wf <- tune_wf %>% finalize_workflow(best_tune)
  
  return(list(rf_wf_final = rf_final_wf, 
              best_rf_model = best_tune,
              best_rf_param = best_tuned_param))
}

metric_compute_site <- function(data_input,site_input){
  cor_model <- cor(data_input$model_predict,
                   data_input$gfactor,
                   use = "pairwise.complete.obs")
  
  tradrsq_model <- yardstick::rsq_trad(data=data_input, 
                                       truth=.data$gfactor, 
                                       estimate=.data$model_predict)
  
  mae_model <- yardstick::mae(data=data_input, 
                              truth=.data$gfactor, 
                              estimate=.data$model_predict)
  
  rmse_model <- yardstick::rmse(data=data_input, 
                                truth=.data$gfactor, 
                                estimate=.data$model_predict)
  return(tibble(correlation=cor_model,  tradrsq= tradrsq_model$.estimate ,MAE= mae_model$.estimate, RMSE=rmse_model$.estimate,
                site = unique(site_input$SITE_ID_L)))
} 


metric_compute <- function(data_input){
  cor_model <- cor(data_input$model_predict,
                   data_input$target,
                   use = "pairwise.complete.obs")
  
  tradrsq_model <- yardstick::rsq_trad(data=data_input, 
                                       truth=.data$target, 
                                       estimate=.data$model_predict)
  
  mae_model <- yardstick::mae(data=data_input, 
                              truth=.data$target, 
                              estimate=.data$model_predict)
  
  rmse_model <- yardstick::rmse(data=data_input, 
                                truth=.data$target, 
                                estimate=.data$model_predict)
  return(tibble(correlation=cor_model,  tradrsq= tradrsq_model$.estimate ,MAE= mae_model$.estimate, RMSE=rmse_model$.estimate))
} 


average_metric <- function(metric_list,pred_names){
  metric_average <- metric_list %>% dplyr::select(-site)%>% colMeans()
  metric_sd <- metric_list %>% dplyr::select(-site)%>% as.matrix()%>% matrixStats::colSds()
  
  title_name <- plotting_names$plotting_name[which(plotting_names$Original_name==pred_names)]
  output_tibble <-tibble(correlation= metric_average[1],
                         cor_sd = metric_sd[1],
                         tradrsq= metric_average[2],
                         rsq_sd = metric_sd[2],
                         MAE = metric_average[3],
                         mae_sd = metric_sd[3],
                         RMSE =metric_average[4],
                         rmse_sd = metric_sd[4],
                         modality=title_name)
}



average_metric_one_mod <- function(metric_list){
  metric_average <- metric_list %>% dplyr::select(-"site")%>% colMeans()
  metric_sd <- metric_list %>% dplyr::select(-"site")%>% as.matrix()%>% matrixStats::colSds()
    output_tibble <-tibble(correlation= metric_average[1],
                         cor_sd = metric_sd[1],
                         tradrsq= metric_average[2],
                         rsq_sd = metric_sd[2],
                         MAE = metric_average[3],
                         mae_sd = metric_sd[3],
                         RMSE =metric_average[4],
                         rmse_sd = metric_sd[4])
}

### partial least square regression models

pls_tune <- function(recipe_input,feature_input=features){
  set.seed(123) 
  train_input <- recipe_input[["preped_recipe"]] %>% bake(new_data=NULL)
  tuning_cv_folds <- train_input  %>%
    vfold_cv(v = 10)
  
  pls_model <- parsnip::pls(num_comp = tune()) %>% 
    set_mode("regression") %>% 
    set_engine("mixOmics")
  
  pls_workflow <- workflow() %>% 
    add_recipe(recipe_input[["norm_recipe"]]) %>% 
    add_model(pls_model)
  
  # create grid
  pls_grid <- expand.grid(num_comp = seq (from = 1, to = length(feature_input), by = 1))
  
  tune_ctrl <- control_grid(save_pred = TRUE, verbose = TRUE
                            ,parallel_over = "everything")
  
  
  #start <- Sys.time()
  tune_res <- tune_grid(
    pls_workflow,
    resamples = tuning_cv_folds,
    metrics = metric_set( rmse),
    grid = pls_grid,
    control= tune_ctrl
  )
  
  
  model_results <- tune_res %>% collect_metrics()
  
  
  # Now find the least complex model that has no more than a 0.1% loss of RMSE:
  best_tune <- select_by_pct_loss(tune_res, 
                           metric = "rmse",
                          limit = 0.1,
                            num_comp)
  
  best_tuned_param <- show_best(tune_res, 
                                metric="rmse")
  
  
  pls_final_wf <- pls_workflow %>% finalize_workflow(best_tune)
  
  return(list(pls_final_wf = pls_final_wf, 
              best_pls_model = best_tune,
              best_pls_param = best_tuned_param,
              pls_grid = tune_res))
}

### functions for feature importance in enet

rename_coef_site <- function(data_input,site_input){
  names(data_input) <- c("term",    paste0(site_input,"_estimate") )
  return(data_input)
}

coef_processing <- function(list_input){
  cleaned_list <- list_input %>% map(.,~filter(.x,term !="(Intercept)")%>%
                                       dplyr::select(-penalty))
  renamed_list <- map2(.x = cleaned_list,.y=site_char,~rename_coef_site(data_input =.x,site_input=.y))
  joined_tibble <- plyr::join_all(renamed_list,by = "term")
  coefs_tibbles <- joined_tibble %>% dplyr::select(ends_with("_estimate"))
  output_tibble <-joined_tibble %>% mutate(mean_estimate = rowMeans(coefs_tibbles))
  return(output_tibble)
}


### functions for scatterplot


scatter_plot_gfactor <- function(data_input,name_input){
  corr_metric <- cor(data_input$model_predict,
                     data_input$gfactor,
                     use = "pairwise.complete.obs")
  
  scatter_plot <-  ggplot(data_input,aes(x = scale(.data[["model_predict"]]) , 
                                         y = scale(.data[["gfactor"]]))) +
    geom_jitter(height = 0.1, width = 0.1, size = 1, col = 'grey60') +
    geom_smooth(method = 'lm', se = FALSE, col = 'black')  +
    labs(x = NULL,
         y = NULL,
         title = paste (name_input,'\nr = ',round(corr_metric, 3)))+
    theme(axis.text.x = element_text(size = 12),                      
          axis.text.y = element_text(size = 12),                     
          plot.title = element_text(size=16)) + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank())
  return(scatter_plot)
  
}

scatter_plot_gfactor_string <- function(data_input,name_input,corr_string_input){
  scatter_plot <-  ggplot(data_input,aes(x = scale(.data[["model_predict"]]) , 
                                         y = scale(.data[["gfactor"]]))) +
    geom_pointdensity(size = 1) +
    scale_color_viridis()+
    geom_smooth(method = 'lm', se = FALSE, col = 'black')  +
    labs(x = NULL,
         y = NULL,
         title = paste (name_input,'\nr = ',corr_string_input))+
    scale_x_continuous(limits=c(-5,5))+
    scale_y_continuous(limits=c(-5,5))+
    theme_classic() + 
    theme(axis.text.x = element_text(size = 35),                      
          axis.text.y = element_text(size = 35),                     
          plot.title = element_text(size=35)) + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "none")
  return(scatter_plot)
  
}

### functions to analyse dummy variables

### only processing the features, gfactor, the respone, is computed in the other function
recipe_prep_dummy_no_scale <- function(train_input=train_gfactor_scan_enet){
  norm_recipe <- recipe( as.formula("gfactor~."), data = train_input) %>%
    update_role(all_of(features), new_role = "predictor")%>%
    update_role("gfactor", new_role = "outcome" )%>% 
    #update_role("SRC_SUBJECT_ID", new_role = "id variable")%>%
    # Impute missing for categorical
    step_impute_mode(all_nominal(), -all_outcomes()) %>%
    # change all nominal into dummy variables
    step_dummy(all_nominal(), -all_outcomes()) %>%
    # Impute missing for numeric
    step_impute_knn(all_numeric(), -all_outcomes(),neighbors = 5) %>%
    # remove na from outcome
    step_naomit(all_outcomes(),all_predictors())
  #%>%
  # step_normalize(all_predictors(),means = 0, sds = 1)
  return(norm_recipe)
}

dummy_data_table_recipe_processing <- function(data_input){
  data_input_select <- data_input%>%dplyr::select(-all_of(subj_info))
  ### select the variables
  output_data  <-  recipe_prep_dummy_no_scale(train_input = data_input_select)%>%
    prep(data_input_select)%>%
    juice()
  
  
  output_data_with_info <- mutate(output_data,SRC_SUBJECT_ID = data_input$SRC_SUBJECT_ID,
                                  EVENTNAME = data_input$EVENTNAME,
                                  SITE_ID_L=data_input$SITE_ID_L) 
  return(output_data_with_info)
}

data_processing_cross_sites_seperate_dummy <- function(baseline_train,
                                                 baseline_test,
                                                 followup_train,
                                                 followup_test,
                                                 features_input= features){
 ### use recipe to deal with dummy variables
  processed_baseline_train <- dummy_data_table_recipe_processing(data_input = baseline_train)
  processed_baseline_test <- dummy_data_table_recipe_processing(data_input = baseline_test)
  processed_followup_train <- dummy_data_table_recipe_processing(data_input = followup_train)
  processed_followup_test <- dummy_data_table_recipe_processing(data_input = followup_test)
  
 dummy_feature_vec <-  processed_baseline_train %>% 
   dplyr::select(-all_of(c(subj_info)))%>%
   dplyr::select(-"gfactor")%>%
                 colnames()
   ## do not process gfactor in this function 
 ## gfactor is processed in another script
 ## use the same response variable across all the analyses
 
 processed_baseline_train_iqr <- processed_baseline_train%>%
   dplyr::select(all_of(dummy_feature_vec),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
 #%>%
#    IQR_remove_vec(x = dummy_feature_vec)
  
 processed_baseline_test_iqr <- processed_baseline_test %>%
   dplyr::select(all_of(dummy_feature_vec),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
 #%>%
#    IQR_remove_vec(x = dummy_feature_vec)
  
 processed_followup_train_iqr <- processed_followup_train %>%
   dplyr::select(all_of(dummy_feature_vec),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
 #%>%
#    IQR_remove_vec(x = dummy_feature_vec)
  
 processed_followup_test_iqr <- processed_followup_test %>%
   dplyr::select(all_of(dummy_feature_vec),all_of(subj_info))%>%
    distinct(SRC_SUBJECT_ID, .keep_all = TRUE)%>% 
    drop_na()
 #%>%
  #  IQR_remove_vec(x = dummy_feature_vec)
  
  scaled_train_list <- scale_train(baseline_data = processed_baseline_train_iqr,
                                   followup_data = processed_followup_train_iqr)
  scaled_test_list <- scale_train(baseline_data = processed_baseline_test_iqr,
                                  followup_data = processed_followup_test_iqr)
  
  
  scaled_data_list <- list(output_train_baseline = scaled_train_list[["output_baseline_data"]], 
                           output_test_baseline = scaled_test_list[["output_baseline_data"]],
                           output_train_followup=scaled_train_list[["output_followup_data"]],
                           output_test_followup=scaled_test_list[["output_followup_data"]])
 
   gfacor_baseline_train <- baseline_train %>% dplyr::select(all_of(c(subj_info)),"gfactor")
  gfacor_baseline_test <- baseline_test %>% dplyr::select(all_of(c(subj_info)),"gfactor")
  gfacor_followup_train <- followup_train %>% dplyr::select(all_of(c(subj_info)),"gfactor")
  gfacor_followup_test <- followup_test %>% dplyr::select(all_of(c(subj_info)),"gfactor")
  
  
  scaled_data_list[["output_train_baseline"]] <- full_join(scaled_data_list[["output_train_baseline"]],
                                                          gfacor_baseline_train, by= subj_info)%>% 
    drop_na()
  scaled_data_list[["output_test_baseline"]] <- full_join(scaled_data_list[["output_test_baseline"]],
                                                          gfacor_baseline_test, by= subj_info)%>% 
    drop_na()
  scaled_data_list[["output_train_followup"]] <- full_join(scaled_data_list[["output_train_followup"]],
                                                           gfacor_followup_train, by= subj_info)%>% 
    drop_na()
  scaled_data_list[["output_test_followup"]] <- full_join(scaled_data_list[["output_test_followup"]],
                                                          gfacor_followup_test, by= subj_info)%>% 
    drop_na()
  
  return(scaled_data_list)
}

### another workflow: use recipe to  process and scale the dummy variables then put into model directly
### The recipe to deal with the dummy variables.
recipe_prep_dummy <- function(train_input=train_gfactor_scan_enet){
  norm_recipe <- recipe( as.formula("gfactor~."), data = train_input) %>%
    update_role(all_of(features), new_role = "predictor")%>%
    update_role("gfactor", new_role = "outcome" )%>%
    # Impute missing for categorical
    step_impute_mode(all_nominal(), -all_outcomes()) %>%
    # change all nominal into dummy variables
    step_dummy(all_nominal(), -all_outcomes()) %>%
    # normalize numeric predictors and outcome
    step_normalize(all_predictors())%>%
    # Impute missing for numeric
    step_impute_knn(all_numeric(), -all_outcomes(),neighbors = 5) %>%
    # remove na from outcome
    step_naomit(all_outcomes(),all_predictors())
  #%>%
  # step_normalize(all_predictors(),means = 0, sds = 1)
  return(norm_recipe)
}

##see what the recipe is like
#trial_data <- gfactor_ses_baseline_train_select[[1]]
#trial_recipe <- recipe_prep_dummy(train_input = trial_data)
#trial_tibble <- trial_recipe%>%
#      prep(training = trial_data, retain = TRUE)
#trial_data_table <- trial_tibble%>% bake(new_data = NULL)

#enet tuning function changed due to the new recipe

### fit the elastic net model 
enet_tuning_dummy <- function(train_input,recipe_input){
  set.seed(123) 
  tuning_cv_folds <- train_input  %>%
    vfold_cv(v = 10)
  
  ## mtry is the number of predictors to sample at each split
  ## min_n (the number of observations needed to keep splitting nodes)
  model_tune <-linear_reg(penalty =tune(),  
                          mixture = tune()) %>%
    set_mode("regression") %>%
    set_engine("glmnet")
  
  tune_wf <- workflow() %>%
    add_recipe(recipe_input) %>%
    add_model(model_tune)
  
  ## automate generate grid for hyperparameters
  
  model_grid <- 
    model_tune %>% 
    extract_parameter_set_dials(tune_wf)%>% 
    ## update the values of the parameters
    update(penalty =penalty(range = c(-10,1), trans = log10_trans()))%>%
    update(mixture =mixture(range = c(0,1)))%>%
    grid_regular(levels = c(200,11))
  
  tune_ctrl <- control_grid(save_pred = TRUE, verbose = TRUE
                            ,parallel_over = "everything"
  )
  
  
  #start <- Sys.time()
  tune_res <- tune_grid(
    tune_wf,
    resamples = tuning_cv_folds,
    metrics = metric_set(rmse),
    grid = model_grid,
    control= tune_ctrl
  )
  
  best_tune <- select_best(tune_res, 
                           metric = "rmse")
  
  best_tuned_param <- show_best(tune_res, 
                                metric="rmse")
  
  enet_final_wf <- tune_wf %>% finalize_workflow(best_tune)
  return(list(enet_wf_final = enet_final_wf, 
              best_enet_model = best_tune,
              best_enet_forest_param = best_tuned_param))
}
### performance of baseline

model_final_fit_dummy <- function(recipe_input,
                                  wf_input,
                                  train_data,
                                  test_data){
  train_input <- prep(recipe_input, train_data)%>% juice() 
  
  ##baked recipe scale the test data with the mean and sd in the training data
  test_input <-  prep(recipe_input, test_data)%>% juice()
  
  model_final_fit <- 
    wf_input%>%
    parsnip::extract_spec_parsnip()%>%
    parsnip::fit(data = train_input, formula= as.formula("gfactor~."))
  
  model_predict <- predict(model_final_fit, 
                           new_data = test_input %>% 
                             drop_na() ) %>%
    rename(model_predict = .pred) %>% 
    bind_cols(test_input%>% drop_na())  
  
  ##processing output
  
  output_list <- vector("list",length=2)
  names(output_list) <- c(paste0("model","_final_fit"),
                          paste0("model","_predict"))
  
  output_list[[paste0("model","_final_fit")]] <- model_final_fit
  output_list[[paste0("model","_predict")]] <- model_predict
  return(output_list)
}


# The VennDiagram package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION TO DRAW VENN DIAGRAM WITH TWO SETS ###################################################
draw.pairwise.venn_invert <- function(
    area1,
    area2,
    cross.area,
    category = rep('', 2),
    euler.d = TRUE,
    scaled = TRUE,
    inverted = FALSE,
    ext.text = TRUE,
    ext.percent = rep(0.05, 3),
    lwd = rep(2, 2),
    lty = rep('solid', 2),
    col = rep('black', 2),
    fill = NULL,
    alpha = rep(0.5, 2),
    label.col = rep('black', 3),
    cex = rep(1, 3),
    fontface = rep('plain', 3),
    fontfamily = rep('serif', 3),
    cat.pos = c(-50, 50),
    cat.dist = rep(0.025, 2),
    cat.cex = rep(1, 2),
    cat.col = rep('black', 2),
    cat.fontface = rep('plain', 2),
    cat.fontfamily = rep('serif', 2),
    cat.just = rep(list(c(0.5, 0.5)), 2),
    cat.default.pos = 'outer',
    cat.prompts = FALSE,
    ext.pos = rep(0, 2),
    ext.dist = rep(0, 2),
    ext.line.lty = 'solid',
    ext.length = rep(0.95, 2),
    ext.line.lwd = 1,
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE,
    sep.dist = 0.05,
    offset = 0,
    cex.prop=NULL,
    print.mode = 'raw',
    sigdigs=3,
    ...
) {
  
  # area1 > area2 OR area1 < area2 plots the same Venn diagram.  Invert using the 'inverted' argument.
  # check parameter lengths and plausibility of Venn diagram
  if (length(category) == 1) { category <- rep(category, 2); }
  else if (length(category) != 2) { flog.error('Unexpected parameter length for "category"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "category"'); }
  
  if (length(ext.percent) == 1) { ext.percent <- rep(ext.percent, 3); }
  else if (length(ext.percent) != 3) { flog.error('Unexpected parameter length for "ext.percent"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "ext.percent"'); }
  
  if (length(ext.pos) == 1) { ext.pos <- rep(ext.pos, 2); }
  else if (length(ext.pos) != 2) { flog.error('Unexpected parameter length for "ext.pos"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "ext.pos"'); }
  
  if (length(ext.dist) == 1) { ext.dist <- rep(ext.dist, 2); }
  else if (length(ext.dist) != 2) { flog.error('Unexpected parameter length for "ext.dist"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "ext.dist"'); }
  
  if (length(ext.length) == 1) { ext.length <- rep(ext.length, 2); }
  else if (length(ext.length) != 2) { flog.error('Unexpected parameter length for "ext.length"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "ext.length"'); }
  
  if (length(lwd) == 1) { lwd <- rep(lwd, 2); }
  else if (length(lwd) != 2) { flog.error('Unexpected parameter length for "lwd"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lwd"'); }
  
  if (length(lty) == 1) { lty <- rep(lty, 2); }
  else if (length(lty) != 2) { flog.error('Unexpected parameter length for "lty"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lty"'); }
  
  if (length(col) == 1) { col <- rep(col, 2); }
  else if (length(col) != 2) { flog.error('Unexpected parameter length for "col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "col"'); }
  
  if (length(label.col) == 1) { label.col <- rep(label.col, 3); }
  else if (length(label.col) != 3) { flog.error('Unexpected parameter length for "label.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "label.col"'); }
  
  if (length(cex) == 1) { cex <- rep(cex, 3); }
  else if (length(cex) != 3) { flog.error('Unexpected parameter length for "cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cex"'); }
  
  if (length(fontface) == 1) { fontface <- rep(fontface, 3); }
  else if (length(fontface) != 3) { flog.error('Unexpected parameter length for "fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontface"'); }
  
  if (length(fontfamily) == 1) { fontfamily <- rep(fontfamily, 3); }
  else if (length(fontfamily) != 3) { flog.error('Unexpected parameter length for "fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontfamily"'); }
  
  if (length(fill) == 1) { fill <- rep(fill, 2); }
  else if (length(fill) != 2 & length(fill) != 0) { flog.error('Unexpected parameter length for "fill"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fill"'); }
  
  if (length(alpha) == 1) { alpha <- rep(alpha, 2); }
  else if (length(alpha) != 2 & length(alpha) != 0) { flog.error('Unexpected parameter length for "alpha"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "alpha"'); }
  
  if (length(ext.line.lwd) != 1) { flog.error('Unexpected parameter length for "ext.line.lwd"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "ext.line.lwd"'); }
  
  if (length(cat.pos) == 1) { cat.pos <- rep(cat.pos, 2); }
  else if (length(cat.pos) != 2) { flog.error('Unexpected parameter length for "cat.pos"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.pos"'); }
  
  if (length(cat.dist) == 1) { cat.dist <- rep(cat.dist, 2); }
  else if (length(cat.dist) != 2) { flog.error('Unexpected parameter length for "cat.dist"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.dist"'); }
  
  if (length(cat.col) == 1) { cat.col <- rep(cat.col, 2); }
  else if (length(cat.col) != 2) { flog.error('Unexpected parameter length for "cat.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.col"'); }
  
  if (length(cat.cex) == 1) { cat.cex <- rep(cat.cex, 2); }
  else if (length(cat.cex) != 2) { flog.error('Unexpected parameter length for "cat.cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.cex"'); }
  
  if (length(cat.fontface) == 1) { cat.fontface <- rep(cat.fontface, 2); }
  else if (length(cat.fontface) != 2) { flog.error('Unexpected parameter length for "cat.fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontface"'); }
  
  if (length(cat.fontfamily) == 1) { cat.fontfamily <- rep(cat.fontfamily, 2); }
  else if (length(cat.fontfamily) != 2) { flog.error('Unexpected parameter length for "cat.fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontfamily"'); }
  
  if (length(offset) != 1) { flog.error('Unexpected parameter length for "Offset". Try using "rotation.degree" to achieve non-vertical offsets',name='VennDiagramLogger')
    stop('Unexpected parameter length for "Offset". Try using "rotation.degree" to achieve non-vertical offsets'); }
  
  if (!(is.list(cat.just) && length(cat.just) == 2 && length(cat.just[[1]]) == 2 && length(cat.just[[2]]) == 2)) {
    flog.error('Unexpected parameter format for "cat.just"',name='VennDiagramLogger')
    stop('Unexpected parameter format for "cat.just"');
  }
  
  # check uninterpretable parameters
  if (!euler.d & scaled) {
    flog.error('Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.',name='VennDiagramLogger')
    stop('Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.');
  }
  if (offset > 1 | offset < 0) {
    flog.error('"Offset" must be between 0 and 1.  Try using "rotation.degree = 180" to achieve offsets in the opposite direction.',name='VennDiagramLogger')
    stop('"Offset" must be between 0 and 1.  Try using "rotation.degree = 180" to achieve offsets in the opposite direction.');
  }
  
  if (cross.area > area1 | cross.area > area2) { flog.error('Impossible: cross section area too large.',name='VennDiagramLogger')
    stop('Impossible: cross section area too large.'); }
  cat.pos <- cat.pos + rotation.degree;
  
  # check category label defaults
  if (((cat.default.pos != 'outer') & (cat.default.pos != 'text')) & cat.prompts) {
    # PHH: removed this check from the if, so that code works with expressions: & isTRUE(category != rep("', 2))
    flog.info('No default location recognized.  Automatically changing to "outer"',name='VennDiagramLogger');
    cat.default.pos <- 'outer';
  }
  if ((cat.default.pos == 'outer') & cat.prompts) {
    flog.info('Placing category labels at default outer locations.  Use "cat.pos" and "cat.dist" to modify location.',name='VennDiagramLogger');
    flog.info(paste('Current "cat.pos":', cat.pos[1], 'degrees,', cat.pos[2], 'degrees'),name='VennDiagramLogger');
    flog.info(paste('Current "cat.dist":', cat.dist[1], ',', cat.dist[2]),name='VennDiagramLogger');
  }
  if ((cat.default.pos == 'text') & cat.prompts) {
    flog.info('Placing category labels at default text locations.  Use "cat.pos" and "cat.dist" to modify location.',name='VennDiagramLogger');
    flog.info(paste('Current "cat.pos":', cat.pos[1], 'degrees,', cat.pos[2], 'degrees'),name='VennDiagramLogger');
    flog.info(paste('Current "cat.dist":', cat.dist[1], ',', cat.dist[2]),name='VennDiagramLogger');
  }
  
  max.circle.size = 0.2;
  
  # initialize logical variables to hold special conditions
  special.coincidental <- FALSE;
  special.inclusion <- FALSE;
  special.exclusion <- FALSE;
  list.switch <- FALSE;
  
  # initialize gList to hold all Grobs generated
  grob.list <- gList();
  
  if (!inverted) {
    tmp1 <- max(area1, area2);
    tmp2 <- min(area1, area2);
    if (tmp1 != area1) { list.switch <- TRUE; }
    area1 <- tmp1;
    area2 <- tmp2;
    r1 <- sqrt(area1 / pi);
    r2 <- sqrt(area2 / pi);
    if (r2 == 0) {r2 <- 0.5*r1 }
    shrink.factor <- max.circle.size / r1;
  }
  else {
    tmp1 <- min(area1, area2);
    tmp2 <- max(area1, area2);
    if (tmp1 != area1) { list.switch <- TRUE; }
    area1 <- tmp1;
    area2 <- tmp2;
    r1 <- sqrt(area1 / pi);
    r2 <- sqrt(area2 / pi);
    if (r1 == 0) {r1 <- 0.5*r2 }
    shrink.factor <- max.circle.size / r2;
  }
  
  # reverse the list if the order is backwards OR inverted is called (both just reverts to normal)
  if (xor(list.switch, inverted)) {
    category <- rev(category);
    lwd <- rev(lwd);
    lty <- rev(lty);
    col <- rev(col);
    fill <- rev(fill);
    alpha <- rev(alpha);
    label.col <- rev(label.col);
    cex <- rev(cex);
    fontface <- rev(fontface);
    fontfamily <- rev(fontfamily);
    cat.pos <- rev(cat.pos);
    cat.dist <- rev(cat.dist);
    cat.col <- rev(cat.col);
    cat.cex <- rev(cat.cex);
    cat.fontface <- rev(cat.fontface);
    cat.fontfamily <- rev(cat.fontfamily);
    cat.just <- rev(cat.just);
    ext.pos <- rev(ext.pos);
    # ext.dist <- rev(ext.dist); # ext.dist intentionally not swapped
    ext.length <- rev(ext.length);
  }
  
  # convert radii to Grid dimensions
  r1 <- r1 * shrink.factor;
  r2 <- r2 * shrink.factor;
  
  # check special conditions
  if (area1 == area2 & area2 == cross.area) { special.coincidental <- TRUE; }
  if (cross.area != 0 & (cross.area == area2 | cross.area == area1)) { special.inclusion <- TRUE; }
  if (0 == cross.area) { special.exclusion <- TRUE; }
  
  denom <- area1+area2-cross.area;
  
  wrapLab <- function(num){
    stri = '';
    if(print.mode[1] == 'percent'){
      stri <- paste(signif(num*100/denom,digits=sigdigs),'%',sep='');
      if(isTRUE(print.mode[2] == 'raw'))
      {
        stri <- paste(stri,'\n(',num,')',sep='');
      }
    }
    if(print.mode[1] == 'raw')
    {
      stri <- num;
      if(isTRUE(print.mode[2] == 'percent'))
      {
        stri <- paste(stri,'\n(',paste(signif(num*100/denom,digits=sigdigs),'%)',sep=''),sep='');
      }
    }
    return(stri);
  }
  
  #	flog.info(c(area1,area2,cross.area),name='VennDiagramLogger');
  
  #	altCross <- cross.area;
  #	altArea1 <- area1;
  #	altArea2 <- area2;
  
  #	#Do processing on the areas and the cross.area to turn them into the required numbers for printing
  #	if(print.mode[1] == 'percent')
  #	{
  #		denom <- area1+area2-cross.area;
  #		area1 <- area1*100/denom;
  #		area2 <- area2*100/denom;
  #		cross.area <- cross.area*100/denom;
  #	}
  #	else #print.mode[1] == 'raw'
  #	{
  #		denom <- area1+area2-cross.area;
  #		altArea1 <- area1*100/denom;
  #		altArea2 <- area2*100/denom;
  #		altCross <- cross.area*100/denom;
  #	}
  
  #	flog.info(c(area1,area2,cross.area),name='VennDiagramLogger');
  
  # plot scaled, generic pairwise Venn diagram with or without external texts
  # ALL OF THE BELOW SECTIONS HAVE A SIMILAR STRUCTURE TO THIS IF BRACKET
  # IF YOU ARE TRYING TO FIGURE OUT WHAT A CERTAIN SECTION DOES, REFER TO THE ANALOGOUS SECTION INSIDE THIS IF BRACKET
  if (scaled & !special.inclusion & !special.exclusion & !special.coincidental) {
    
    # calculate centres of circles
    d <- find.dist(area1, area2, cross.area, inverted = inverted);
    d <- d * shrink.factor;
    x.centre.1 <- (1 + r1 - r2 - d) / 2;
    x.centre.2 <- x.centre.1 + d;
    
    # draw both circles and their borders
    tmp <- VennDiagram::ellipse(
      x = x.centre.1,
      y = 0.5,
      a = ifelse(!inverted, r1, r2),
      b = ifelse(!inverted, r1, r2),
      gp = gpar(
        lty = 0,
        fill = fill[1],
        alpha = alpha[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = x.centre.2,
      y = 0.5,
      a = ifelse(inverted, r1, r2),
      b = ifelse(inverted, r1, r2),
      gp = gpar(
        lty = 0,
        fill = fill[2],
        alpha = alpha[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = x.centre.1,
      y = 0.5,
      a = ifelse(!inverted, r1, r2),
      b = ifelse(!inverted, r1, r2),
      gp = gpar(
        lwd = lwd[1],
        lty = lty[1],
        col = col[1],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = x.centre.2,
      y = 0.5,
      a = ifelse(inverted, r1, r2),
      b = ifelse(inverted, r1, r2),
      gp = gpar(
        lwd = lwd[2],
        lty = lty[2],
        col = col[2],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    
    ## rescaling area labels to be proportional to area
    if(length(cex.prop) > 0){
      
      if(length(cex.prop) != 1){ 
        flog.error('Value passed to cex.prop is not length 1',name='VennDiagramLogger')
        stop('Value passed to cex.prop is not length 1')
      }
      
      ## figure out what function to use
      func = cex.prop
      if (!is(cex.prop, 'function')) {
        if(cex.prop == 'lin'){
          func = function(x) x
        }
        else if(cex.prop == 'log10'){
          func = log10
        }
        else flog.error(paste0('Unknown value passed to cex.prop: ', cex.prop),name='VennDiagramLogger')
        stop(paste0('Unknown value passed to cex.prop: ', cex.prop))
      }
      
      ## rescale areas
      areas = c(area1 - cross.area, cross.area, area2 - cross.area)
      maxArea = max(areas)            
      for(i in 1:length(areas)){                
        cex[i] = cex[i] * func(areas[i]) / func(maxArea)
        if(cex[i] <= 0) stop(paste0('Error in rescaling of area labels: the label of area ',
                                    i, ' is less than or equal to zero'))
      }
    }
    
    
    # if labels are to be placed outside circles
    if (ext.text) {
      area.1.pos <- x.centre.1 + ifelse(!inverted, -r1 + ( (2 * r1 - (r1 + r2 - d)) / 2), -r2 + ( (2 * r2 - (r2 + r1 - d)) / 2));
      area.2.pos <- x.centre.2 + ifelse(!inverted, r2 - ( (2 * r2 - (r1 + r2 - d)) / 2), r1 - ( (2 * r1 - (r2 + r1 - d)) / 2));
      # distinct area1 is more than the given percentage (label stays inside circle)
      if ( (area1 - cross.area) / area1 > ext.percent[1] & (area1 - cross.area) / area2 > ext.percent[1]) {
        # draw label normally
        tmp <- textGrob(
          label = wrapLab(ifelse(!inverted, area1, area2) - cross.area),
          x = area.1.pos,
          y = 0.5,
          gp = gpar(
            col = label.col[1],
            cex = cex[1],
            fontface = fontface[1],
            fontfamily = fontfamily[1]
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      # percentage is small enough to move label outside circle
      else {
        label.pos <- find.cat.pos(area.1.pos, 0.5, ext.pos[1], ext.dist[1], r1);
        area.1.xpos <- label.pos$x;
        area.1.ypos <- label.pos$y
        # draw label outside
        tmp <- textGrob(
          label = wrapLab(ifelse(!inverted, area1, area2) - cross.area),
          x = area.1.xpos,
          y = area.1.ypos,
          gp = gpar(
            col = label.col[1],
            cex = cex[1],
            fontface = fontface[1],
            fontfamily = fontfamily[1]
          )
        );
        grob.list <- gList(grob.list, tmp);
        # draw line from circle to label
        tmp <- linesGrob(
          x = c(area.1.pos + ext.length[1] * (area.1.xpos - area.1.pos), area.1.pos),
          y = c(0.5 + ext.length[1] * (area.1.ypos - 0.5), 0.5),
          gp = gpar(
            col = label.col[1],
            lwd = ext.line.lwd,
            lty = ext.line.lty
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      
      # distinct area2 is more than the given percentage (label stays inside the circle)
      if ((area2 - cross.area) / area2 > ext.percent[2] & (area2 - cross.area) / area1 > ext.percent[2]) {
        # draw label normally
        tmp <- textGrob(
          label = wrapLab(ifelse(inverted, area1, area2) - cross.area),
          x = area.2.pos,
          y = 0.5,
          gp = gpar(
            col = label.col[3],
            cex = cex[3],
            fontface = fontface[3],
            fontfamily = fontfamily[3]
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      # percentage is small enough to move label outside circle
      else {
        label.pos <- find.cat.pos(area.2.pos, 0.5, ext.pos[2], ext.dist[2], r2);
        area.2.xpos <- label.pos$x;
        area.2.ypos <- label.pos$y;
        # draw label outside
        tmp <- textGrob(
          label = wrapLab(ifelse(inverted, area1, area2) - cross.area),
          x = area.2.xpos,
          y = area.2.ypos,
          gp = gpar(
            col = label.col[3],
            cex = cex[3],
            fontface = fontface[3],
            fontfamily = fontfamily[3]
          )
        );
        grob.list <- gList(grob.list, tmp);
        # draw line from circle to label
        tmp <- linesGrob(
          x = c(area.2.pos + ext.length[1] * (area.2.xpos - area.2.pos), area.2.pos),
          y = c(0.5 + ext.length[1] * (area.2.ypos - 0.5), 0.5),
          gp = gpar(
            col = label.col[3],
            lwd = ext.line.lwd,
            lty = ext.line.lty
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      
      # if intersect area is more than the given percentage (label stays inside area)
      if (cross.area / area2 > ext.percent[3] & cross.area / area1 > ext.percent[3]) {
        # draw label normally
        tmp <- textGrob(
          label = wrapLab(cross.area),
          x = x.centre.1 + (d - ifelse(!inverted, r2, r1)) + (r1 + r2 - d) / 2,
          y = 0.5,
          gp = gpar(
            col = label.col[2],
            cex = cex[2],
            fontface = fontface[2],
            fontfamily = fontfamily[2]
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      # percentage is small enough to move label outside area
      else {
        cross.area.pos <- x.centre.1 + (d - r2) + (r1 + r2 - d) / 2;
        cross.pos <- find.cat.pos(cross.area.pos, 0.5, ext.pos[1], ext.dist[1], r1 + r2);
        cross.area.xpos <- cross.pos$x;
        cross.area.ypos <- cross.pos$y
        # draw label outside
        tmp <- textGrob(
          label = wrapLab(cross.area),
          x = cross.area.xpos,
          y = cross.area.ypos,
          gp = gpar(
            col = label.col[2],
            cex = cex[2],
            fontface = fontface[2],
            fontfamily = fontfamily[2]
          )
        );
        grob.list <- gList(grob.list, tmp);
        # draw line from area to label
        tmp <- linesGrob(
          x = c(cross.area.pos + ext.length[2] * (cross.area.xpos - cross.area.pos), cross.area.pos),
          y = c(0.5 + ext.length[2] * (cross.area.ypos - 0.5), 0.5),
          gp = gpar(
            col = label.col[2],
            lwd = ext.line.lwd,
            lty = ext.line.lty
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
    }
    
    # if the labels are not to be extended, draw them in their usual locations
    else {
      area.1.pos <-  x.centre.1 + ifelse(!inverted, -r1 + ( (2 * r1 - (r1 + r2 - d)) / 2), -r2 + ( (2 * r2 - (r2 + r1 - d)) / 2));
      tmp <- textGrob(
        label = wrapLab(ifelse(!inverted, area1, area2) - cross.area),
        x = area.1.pos,
        y = 0.5,
        gp = gpar(
          col = label.col[1],
          cex = cex[1],
          fontface = fontface[1],
          fontfamily = fontfamily[1]
        )
      );
      grob.list <- gList(grob.list, tmp);
      area.2.pos <- x.centre.2 + ifelse(!inverted, r2 - ( (2 * r2 - (r1 + r2 - d)) / 2), r1 - ( (2 * r1 - (r2 + r1 - d)) / 2));
      tmp <- textGrob(
        label = wrapLab(ifelse(inverted, area1, area2) - cross.area),
        x = area.2.pos,
        y = 0.5,
        gp = gpar(
          col = label.col[3],
          cex = cex[3],
          fontface = fontface[3],
          fontfamily = fontfamily[3]
        )
      );
      grob.list <- gList(grob.list, tmp);
      tmp <- textGrob(
        label = wrapLab(cross.area),
        x = x.centre.1 + (d - ifelse(!inverted, r2, r1)) + (r1 + r2 - d) / 2,
        y = 0.5,
        gp = gpar(
          col = label.col[2],
          cex = cex[2],
          fontface = fontface[2],
          fontfamily = fontfamily[2]
        )
      );
      grob.list <- gList(grob.list, tmp);
    }
    
    # find the location of the category labels
    if ('outer' == cat.default.pos) {
      cat.pos.1 <- find.cat.pos(x.centre.1, 0.5, (ifelse(!inverted, cat.pos[1], cat.pos[2]) + ifelse(xor(list.switch, inverted), 180, 0)) %% 360, cat.dist[1], ifelse(!inverted, r1, r2));
      cat.pos.2 <- find.cat.pos(x.centre.2, 0.5, (ifelse(!inverted, cat.pos[2], cat.pos[1]) + ifelse(xor(list.switch, inverted), 180, 0)) %% 360, cat.dist[2], ifelse(!inverted, r2, r1));
    }
    else if ('text' == cat.default.pos) {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], cat.dist[1]);
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], cat.dist[2]);
    }
    else {
      flog.error('Invalid value for "cat.default.pos", should be either "outer" or "text"',name='VennDiagramLogger')
      stop('Invalid value for "cat.default.pos", should be either "outer or "text"');
    }
    
    # draw category labels
    tmp <- textGrob(
      label = category[1],
      x = cat.pos.1$x,
      y = cat.pos.1$y,
      just = cat.just[[1]],
      gp = gpar(
        col = cat.col[1],
        cex = cat.cex[1],
        fontface = cat.fontface[1],
        fontfamily = cat.fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = category[2],
      x = cat.pos.2$x,
      y = cat.pos.2$y,
      just = cat.just[[2]],
      gp = gpar(
        col = cat.col[2],
        cex = cat.cex[2],
        fontface = cat.fontface[2],
        fontfamily = cat.fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
  }
  
  # plot scaled Venn diagram when one set is completely included in (but not exactly coincidental with) the other set
  # with or without external texts
  if (euler.d & special.inclusion & !special.coincidental) {
    
    if (inverted) {
      tmp1 <- area1;
      tmp2 <- area2;
      area1 <- tmp2;
      area2 <- tmp1;
    }
    
    if (!scaled & !inverted) {
      r1 <- 0.4;
      r2 <- 0.2;
    }
    
    if (!scaled & inverted) {
      r1 <- 0.2;
      r2 <- 0.4;
    }
    
    # draw circles and their borders
    tmp <- VennDiagram::ellipse(
      x = 0.5,
      y = 0.5,
      a = r1,
      b = r1,
      gp = gpar(
        lty = 0,
        fill = fill[1],
        alpha = alpha[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.5 - offset * (r1 - r2),
      y = 0.5,
      a = r2,
      b = r2,
      gp = gpar(
        lty = 0,
        fill = fill[2],
        alpha = alpha[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.5,
      y = 0.5,
      a = r1,
      b = r1,
      gp = gpar(
        lwd = lwd[1],
        lty = lty[1],
        col = col[1],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.5 - offset * (r1 - r2),
      y = 0.5,
      a = r2,
      b = r2,
      gp = gpar(
        lwd = lwd[2],
        lty = lty[2],
        col = col[2],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    # draw area labels in appropriate locations
    area.2.pos <- 0.5 - offset * (r1 - r2);
    tmp <- textGrob(
      label = wrapLab(area2),
      x = area.2.pos,
      y = 0.5,
      gp = gpar(
        col = label.col[2],
        cex = cex[2],
        fontface = fontface[2],
        fontfamily = fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    if (!ext.text | !scaled) {
      area.1.pos <- (1 + r1 + r2 - offset * (r1 - r2)) / 2;
      tmp <- textGrob(
        label = wrapLab(area1 - area2),
        x = area.1.pos,
        y = 0.5,
        gp = gpar(
          col = label.col[1],
          cex = cex[1],
          fontface = fontface[1],
          fontfamily = fontfamily[1]
        )
      );
      grob.list <- gList(grob.list, tmp);
    }
    
    if (ext.text & scaled) {
      # draw labels and lines if text is to be extended from areas
      if (area2 / area1 > 0.5) {
        area.1.pos <- (1 + r1 + r2 - offset * (r1 - r2)) / 2;
        area.pos <- find.cat.pos(area.1.pos, 0.5, ext.pos[1], ext.dist[1], r1);
        area.1.xpos <- area.pos$x;
        area.1.ypos <- area.pos$y;
        tmp <- textGrob(
          label = wrapLab(area1 - area2),
          x = area.1.xpos,
          y = area.1.ypos,
          gp = gpar(
            col = label.col[1],
            cex = cex[1],
            fontface = fontface[1],
            fontfamily = fontfamily[1]
          )
        );
        grob.list <- gList(grob.list, tmp);
        tmp <- linesGrob(
          x = c(area.1.pos + ext.length * (area.1.xpos - area.1.pos), area.1.pos),
          y = c(0.5 + ext.length * (area.1.ypos - 0.5), 0.5),
          gp = gpar(
            col = label.col[1],
            lwd = ext.line.lwd,
            lty = ext.line.lty
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
      else {
        area.1.pos <- (1 + r1 + r2 - offset * (r1 - r2)) / 2;
        tmp <- textGrob(
          label = wrapLab(area1 - area2),
          x = area.1.pos,
          y = 0.5,
          gp = gpar(
            col = label.col[1],
            cex = cex[1],
            fontface = fontface[1],
            fontfamily = fontfamily[1]
          )
        );
        grob.list <- gList(grob.list, tmp);
      }
    }
    
    # find the correct position of categories given default position and areas
    if (cat.default.pos == 'outer') {
      cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos[1], cat.dist[1], r1);
      cat.pos.2 <- find.cat.pos(0.5 - offset * (r1 - r2), 0.5, cat.pos[2], cat.dist[2], r2);
    }
    else if (cat.default.pos == 'text') {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], cat.dist[1]);
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], cat.dist[2]);
    }
    else {
      flog.error('Invalid value for "cat.default.pos", should be either "outer" or "text"',name='VennDiagramLogger')
      stop('Invalid value for "cat.default.pos", should be either "outer" or "text"');
    }
    
    # add category labels
    tmp <- textGrob(
      label = category[1],
      x = cat.pos.1$x,
      y = cat.pos.1$y,
      just = cat.just[[1]],
      gp = gpar(
        col = cat.col[1],
        cex = cat.cex[1],
        fontface = cat.fontface[1],
        fontfamily = cat.fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = category[2],
      x = cat.pos.2$x,
      y = cat.pos.2$y,
      just = cat.just[[2]],
      gp = gpar(
        col = cat.col[2],
        cex = cat.cex[2],
        fontface = cat.fontface[2],
        fontfamily = cat.fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
  }
  
  # plot scaled Venn diagrams when the two sets are coincidental
  if (euler.d & special.coincidental) {
    
    # draw the one circle and its border
    tmp <- VennDiagram::ellipse(
      x = 0.5,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lty = 0,
        fill = fill[1],
        alpha = alpha[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.5,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lwd = lwd[1],
        lty = lty[1],
        col = col[1],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    # draw labels on the same circle
    area.1.pos <- 0.46;
    tmp <- textGrob(
      label = wrapLab(area1),
      x = area.1.pos,
      y = 0.5,
      gp = gpar(
        col = label.col[2],
        cex = cex[2],
        fontface = fontface[2],
        fontfamily = fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    area.2.pos <- 0.54;
    tmp <- textGrob(
      label = wrapLab(area2),
      x = area.2.pos,
      y = 0.5,
      gp = gpar(
        col = label.col[2],
        cex = cex[2],
        fontface = fontface[2],
        fontfamily = fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = '(Coincidental)',
      x = 0.5,
      y = 0.45,
      gp = gpar(
        col = label.col[2],
        cex = cex[2],
        fontface = fontface[2],
        fontfamily = fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    if (cat.default.pos == 'outer') {
      cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos[1], cat.dist[1], max.circle.size);
      cat.pos.2 <- find.cat.pos(0.5, 0.5, cat.pos[2], cat.dist[2], max.circle.size);
    }
    else if (cat.default.pos == 'text') {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], cat.dist[1]);
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], cat.dist[2]);
    }
    else {
      flog.error('Invalid value for "cat.default.pos", should be either "outer" or "text"',name='VennDiagramLogger')
      stop('Invalid value for "cat.default.pos", should be either "outer" or "text"');
    }
    
    tmp <- textGrob(
      label = category[1],
      x = cat.pos.1$x,
      y = cat.pos.1$y,
      just = cat.just[[1]],
      gp = gpar(
        col = cat.col[1],
        cex = cat.cex[1],
        fontface = cat.fontface[1],
        fontfamily = cat.fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = category[2],
      x = cat.pos.2$x,
      y = cat.pos.2$y,
      just = cat.just[[2]],
      gp = gpar(
        col = cat.col[2],
        cex = cat.cex[2],
        fontface = cat.fontface[2],
        fontfamily = cat.fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
  }
  
  # plot scaled Venn diagrams when the two sets are mutually exclusive
  if (euler.d & special.exclusion) {
    
    if (!scaled) {
      r1 <- 0.2;
      r2 <- 0.2;
    }
    
    # determine centres of exclusive circles and draw them
    x.centre.1 <- (1 - 2 * (r1 + r2)) / 2 + r1 - sep.dist / 2;
    tmp <- VennDiagram::ellipse(
      x = x.centre.1,
      y = 0.5,
      a = r1,
      b = r1,
      gp = gpar(
        lty = 0,
        fill = fill[1],
        alpha = alpha[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    x.centre.2 <- 1 - (1 - 2 * (r1 + r2)) / 2 - r2 + sep.dist / 2;
    tmp <- VennDiagram::ellipse(
      x = x.centre.2,
      y = 0.5,
      a = r2,
      b = r2,
      gp = gpar(
        lty = 0,
        fill = fill[2],
        alpha = alpha[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = x.centre.1,
      y = 0.5,
      a = r1,
      b = r1,
      gp = gpar(
        lwd = lwd[1],
        lty = lty[1],
        col = col[1],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = x.centre.2,
      y = 0.5,
      a = r2,
      b = r2,
      gp = gpar(
        lwd = lwd[2],
        lty = lty[2],
        col = col[2],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    # draw area and category labels
    area.1.pos <- x.centre.1;
    tmp <- textGrob(
      label = wrapLab(area1),
      x = area.1.pos,
      y = 0.5,
      gp = gpar(
        col = label.col[1],
        cex = cex[1],
        fontface = fontface[1],
        fontfamily = fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    area.2.pos <- x.centre.2;
    tmp <- textGrob(
      label = wrapLab(area2),
      x = area.2.pos,
      y = 0.5,
      gp = gpar(
        col = label.col[3],
        cex = cex[3],
        fontface = fontface[3],
        fontfamily = fontfamily[3]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    if (cat.default.pos == 'outer') {
      cat.pos.1 <- find.cat.pos(x.centre.1, 0.5, cat.pos[1], cat.dist[1], r1);
      cat.pos.2 <- find.cat.pos(x.centre.2, 0.5, cat.pos[2], cat.dist[2], r2);
    }
    else if (cat.default.pos == 'text') {
      cat.pos.1 <- find.cat.pos(area.1.pos, 0.5, cat.pos[1], cat.dist[1]);
      cat.pos.2 <- find.cat.pos(area.2.pos, 0.5, cat.pos[2], cat.dist[2]);
    }
    else {
      flog.error('Invalid value for "cat.default.pos", should be either "outer" or "text"',name='VennDiagramLogger')
      stop('Invalid value for "cat.default.pos", should be either "outer" or "text"');
    }
    
    tmp <- textGrob(
      label = category[1],
      x = cat.pos.1$x,
      y = cat.pos.1$y,
      just = cat.just[[1]],
      gp = gpar(
        col = cat.col[1],
        cex = cat.cex[1],
        fontface = cat.fontface[1],
        fontfamily = cat.fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = category[2],
      x = cat.pos.2$x,
      y = cat.pos.2$y,
      just = cat.just[[2]],
      gp = gpar(
        col = cat.col[2],
        cex = cat.cex[2],
        fontface = cat.fontface[2],
        fontfamily = cat.fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
  }
  
  # plot non-scaled Venn diagram
  if ((!scaled & !euler.d) | (!scaled & euler.d & !special.inclusion & !special.exclusion & !special.coincidental)) {
    
    tmp <- VennDiagram::ellipse(
      x = 0.4,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lty = 0,
        fill = fill[1],
        alpha = alpha[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.6,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lty = 0,
        fill = fill[2],
        alpha = alpha[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.4,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lwd = lwd[1],
        lty = lty[1],
        col = col[1],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- VennDiagram::ellipse(
      x = 0.6,
      y = 0.5,
      a = max.circle.size,
      b = max.circle.size,
      gp = gpar(
        lwd = lwd[2],
        lty = lty[2],
        col = col[2],
        fill = 'transparent'
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = wrapLab(area1 - cross.area),
      x = 0.3,
      y = 0.5,
      gp = gpar(
        col = label.col[1],
        cex = cex[1],
        fontface = fontface[1],
        fontfamily = fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = wrapLab(area2 - cross.area),
      x = 0.7,
      y = 0.5,
      gp = gpar(
        col = label.col[3],
        cex = cex[3],
        fontface = fontface[3],
        fontfamily = fontfamily[3]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = wrapLab(cross.area),
      x = 0.5,
      y = 0.5,
      gp = gpar(
        col = label.col[2],
        cex = cex[2],
        fontface = fontface[2],
        fontfamily = fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    if (cat.default.pos == 'outer') {
      cat.pos.1 <- find.cat.pos(0.4, 0.5, cat.pos[1], cat.dist[1], max.circle.size);
      cat.pos.2 <- find.cat.pos(0.6, 0.5, cat.pos[2], cat.dist[2], max.circle.size);
    }
    else if (cat.default.pos == 'text') {
      cat.pos.1 <- find.cat.pos(0.3, 0.5, cat.pos[1], cat.dist[1]);
      cat.pos.2 <- find.cat.pos(0.7, 0.5, cat.pos[2], cat.dist[2]);
    }
    else {
      flog.error('Invalid value for "cat.default.pos", should be either "outer" or "text"',name='VennDiagramLogger')
      stop('Invalid value for "cat.default.pos", should be either "outer" or "text"');
    }
    
    tmp <- textGrob(
      label = category[1],
      x = cat.pos.1$x,
      y = cat.pos.1$y,
      just = cat.just[[1]],
      gp = gpar(
        col = cat.col[1],
        cex = cat.cex[1],
        fontface = cat.fontface[1],
        fontfamily = cat.fontfamily[1]
      )
    );
    grob.list <- gList(grob.list, tmp);
    
    tmp <- textGrob(
      label = category[2],
      x = cat.pos.2$x,
      y = cat.pos.2$y,
      just = cat.just[[2]],
      gp = gpar(
        col = cat.col[2],
        cex = cat.cex[2],
        fontface = cat.fontface[2],
        fontfamily = cat.fontfamily[2]
      )
    );
    grob.list <- gList(grob.list, tmp);
  }
  
  # rorate Venn if necessary and add other adjustments to plot
  grob.list <- adjust.venn(rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
  # draw the plot before returning the grob if specified
  if (ind) { grid.draw(grob.list); }
  return(grob.list);
}



draw.quad.venn_round <- function(
    area1,
    area2,
    area3,
    area4,
    n12,
    n13,
    n14,
    n23,
    n24,
    n34,
    n123,
    n124,
    n134,
    n234,
    n1234,
    category = rep('', 4),
    lwd = rep(2, 4),
    lty = rep('solid', 4),
    col = rep('black', 4),
    fill = NULL,
    alpha = rep(0.5, 4),
    label.col = rep('black', 15),
    cex = rep(1, 15),
    fontface = rep('plain', 15),
    fontfamily = rep('serif', 15),
    cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11),
    cat.col = rep('black', 4),
    cat.cex = rep(1, 4),
    cat.fontface = rep('plain', 4),
    cat.fontfamily = rep('serif', 4),
    cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE,
    cex.prop=NULL,
    print.mode = 'raw',
    sigdigs=3,
    direct.area = FALSE,
    area.vector = 0,
    ...
) {
  
  #area1 > area2 > area3 > area4
  # check parameter lengths
  if (length(category) == 1) { cat <- rep(category, 4); }
  else if (length(category) != 4) { flog.error('Unexpected parameter length for "category"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "category"'); }
  
  if (length(lwd) == 1) { lwd <- rep(lwd, 4); }
  else if (length(lwd) != 4) { flog.error('Unexpected parameter length for "lwd"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lwd"'); }
  
  if (length(lty) == 1) { lty <- rep(lty, 4); }
  else if (length(lty) != 4) { flog.error('Unexpected parameter length for "lty"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lty"'); }
  
  if (length(col) == 1) { col <- rep(col, 4); }
  else if (length(col) != 4) { flog.error('Unexpected parameter length for "col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "col"'); }
  
  if (length(label.col) == 1) { label.col <- rep(label.col, 15); }
  else if (length(label.col) != 15) { flog.error('Unexpected parameter length for "label.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "label.col"'); }
  
  if (length(cex) == 1) { cex <- rep(cex, 15); }
  else if (length(cex) != 15) { flog.error('Unexpected parameter length for "cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cex"'); }
  
  if (length(fontface) == 1) { fontface <- rep(fontface, 15); }
  else if (length(fontface) != 15) { flog.error('Unexpected parameter length for "fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontface"'); }
  
  if (length(fontfamily) == 1) { fontfamily <- rep(fontfamily, 15); }
  else if (length(fontfamily) != 15) { flog.error('Unexpected parameter length for "fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontfamily"'); }
  
  if (length(fill) == 1) { fill <- rep(fill, 4); }
  else if (length(fill) != 4 & length(fill) != 0) { flog.error('Unexpected parameter length for "fill"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fill"'); }
  
  if (length(alpha) == 1) { alpha <- rep(alpha, 4); }
  else if (length(alpha) != 4 & length(alpha) != 0) { flog.error('Unexpected parameter length for "alpha"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "alpha"'); }
  
  if (length(cat.pos) == 1) { cat.pos <- rep(cat.pos, 4); }
  else if (length(cat.pos) != 4) { flog.error('Unexpected parameter length for "cat.pos"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.pos"'); }
  
  if (length(cat.dist) == 1) { cat.dist <- rep(cat.dist, 4); }
  else if (length(cat.dist) != 4) { flog.error('Unexpected parameter length for "cat.dist"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.dist"'); }
  
  if (length(cat.col) == 1) { cat.col <- rep(cat.col, 4); }
  else if (length(cat.col) != 4) { flog.error('Unexpected parameter length for "cat.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.col"'); }
  
  if (length(cat.cex) == 1) { cat.cex <- rep(cat.cex, 4); }
  else if (length(cat.cex) != 4) { flog.error('Unexpected parameter length for "cat.cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.cex"'); }
  
  if (length(cat.fontface) == 1) { cat.fontface <- rep(cat.fontface, 4); }
  else if (length(cat.fontface) != 4) { flog.error('Unexpected parameter length for "cat.fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontface"'); }
  
  if (length(cat.fontfamily) == 1) { cat.fontfamily <- rep(cat.fontfamily, 4); }
  else if (length(cat.fontfamily) != 4) { flog.error('Unexpected parameter length for "cat.fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontfamily"'); }
  
  if (!(is.list(cat.just) && length(cat.just) == 4 && length(cat.just[[1]]) == 2 && length(cat.just[[2]]) == 2 && length(cat.just[[3]]) == 2 && length(cat.just[[4]]) == 2)) { flog.error('Unexpected parameter format for "cat.just"',name='VennDiagramLogger')
    stop('Unexpected parameter format for "cat.just"'); }
  cat.pos <- cat.pos + rotation.degree;
  
  if(direct.area){
    areas <- area.vector;
    #create the variables and assign their values from the area vector
    for(i in 1:15)
    {
      assign(paste('a',i,sep=''),area.vector[i]);
    }
  }
  else {
    # generate partial areas from given arguments
    a6  <- n1234;
    a12 <- n123 - a6;
    a11 <- n124 - a6;
    a5  <- n134 - a6;
    a7  <- n234 - a6;
    a15 <- n12 - a6 - a11 - a12;
    a4  <- n13 - a6 - a5 - a12;
    a10 <- n14 - a6 - a5 - a11;
    a13 <- n23 - a6 - a7 - a12;
    a8  <- n24 - a6 - a7 - a11;
    a2  <- n34 - a6 - a5 - a7;
    a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15;
    a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15;
    a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13;
    a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11;
    
    # check plausibility and 0 partial areas
    areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15);
  }
  areas.error <- c(
    'a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13',
    'a2  <- n34 - a6 - a5 - a7',
    'a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11',
    'a4  <- n13 - a6 - a5 - a12',
    'a5  <- n134 - a6',
    'a6  <- n1234',
    'a7  <- n234 - a6',
    'a8  <- n24 - a6 - a7 - a11',
    'a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15',
    'a10 <- n14 - a6 - a5 - a11',
    'a11 <- n124 - a6',
    'a12 <- n123 - a6',
    'a15 <- n12 - a6 - a11 - a12',
    'a13 <- n23 - a6 - a7 - a12',
    'a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15'
  );
  for (i in 1:length(areas)) {
    if (areas[i] < 0) {
      flog.error(paste('Impossible:', areas.error[i], 'produces negative area'),name='VennDiagramLogger')
      stop(paste('Impossible:', areas.error[i], 'produces negative area'));
    }
  }
  
  ## rescaling area labels to be proportional to area
  if(length(cex.prop) > 0){
    
    if(length(cex.prop) != 1) {
      flog.error('Value passed to cex.prop is not length 1',name='VennDiagramLogger')
      stop('Value passed to cex.prop is not length 1')
    }
    
    ## figure out what function to use
    func = cex.prop
    if (!is(cex.prop, 'function')) {
      if(cex.prop == 'lin'){
        func = function(x) x
      }
      else if(cex.prop == 'log10'){
        func = log10
      }
      else flog.error(paste0('Unknown value passed to cex.prop: ', cex.prop),name='VennDiagramLogger')
      stop(paste0('Unknown value passed to cex.prop: ', cex.prop))
    }
    
    ## rescale areas
    maxArea = max(areas)            
    for(i in 1:length(areas)){                
      cex[i] = cex[i] * func(areas[i]) / func(maxArea)
      if(cex[i] <= 0) stop(paste0('Error in rescaling of area labels: the label of area ',
                                  i, ' is less than or equal to zero'))
    }
  }
  
  # initialize gList to hold all Grobs generated
  grob.list <- gList();
  
  # plot the ellipses of the Venn diagram
  ellipse.positions <- matrix(
    nrow = 4,
    ncol = 7
  );
  colnames(ellipse.positions) <- c('x', 'y', 'a', 'b', 'rotation', 'fill.mapping', 'line.mapping');
  
  ellipse.positions[1,] <- c(0.65, 0.47, 0.35, 0.20,  45, 2, 2);
  ellipse.positions[2,] <- c(0.35, 0.47, 0.35, 0.20, 135, 1, 1);
  ellipse.positions[3,] <- c(0.50, 0.57, 0.33, 0.15,  45, 4, 4);
  ellipse.positions[4,] <- c(0.50, 0.57, 0.35, 0.15, 135, 3, 3);
  
  # draw the ellipses themselves
  for (i in 1:4) {
    grob.list <- gList(
      grob.list,
      VennDiagram::ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = gpar(
          lty = 0,
          fill = fill[ellipse.positions[i,'fill.mapping']],
          alpha = alpha[ellipse.positions[i,'fill.mapping']]
        )
      )
    );
  }
  
  # draw the ellipse borders
  for (i in 1:4) {
    grob.list <- gList(
      grob.list,
      ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = gpar(
          lwd = lwd[ellipse.positions[i,'line.mapping']],
          lty = lty[ellipse.positions[i,'line.mapping']],
          col = col[ellipse.positions[i,'line.mapping']],
          fill = 'transparent'
        )
      )
    );
  }
  
  # create the labels
  label.matrix <- matrix(
    nrow = 15,
    ncol = 3
  );
  colnames(label.matrix) <- c('label', 'x', 'y');
  
  label.matrix[ 1,] <- c(round(a1,2),  0.350, 0.77);
  label.matrix[ 2,] <- c(round(a2,2),  0.500, 0.69);
  label.matrix[ 3,] <- c(round(a3,2),  0.650, 0.77);
  label.matrix[ 4,] <- c(round(a4,2),  0.310, 0.67);
  label.matrix[ 5,] <- c(round(a5,2),  0.400, 0.58);
  label.matrix[ 6,] <- c(round(a6,2),  0.500, 0.47);
  label.matrix[ 7,] <- c(round(a7,2),  0.600, 0.58);
  label.matrix[ 8,] <- c(round(a8,2),  0.690, 0.67);
  label.matrix[ 9,] <- c(round(a9,2),  0.180, 0.58);
  label.matrix[10,] <- c(round(a10,2), 0.320, 0.42);
  label.matrix[11,] <- c(round(a11,2), 0.425, 0.38);
  label.matrix[12,] <- c(round(a12,2), 0.575, 0.38);
  label.matrix[13,] <- c(round(a13,2), 0.680, 0.42);
  label.matrix[14,] <- c(round(a14,2), 0.820, 0.58);
  label.matrix[15,] <- c(round(a15,2), 0.500, 0.28);
  
  processedLabels <- rep('',length(label.matrix[,'label']));
  if(print.mode[1] == 'percent'){
    processedLabels <- paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),'%',sep='');
    if(isTRUE(print.mode[2] == 'raw'))
    {
      processedLabels <- paste(processedLabels,'\n(',label.matrix[,'label'],')',sep='');
    }
  }
  if(print.mode[1] == 'raw'){
    processedLabels <- label.matrix[,'label'];
    if(isTRUE(print.mode[2] == 'percent'))
    {
      processedLabels <- paste(processedLabels,'\n(',paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),'%)',sep=''),sep='');
    }
  }
  
  
  for (i in 1:nrow(label.matrix)) {
    grob.list <- gList(
      grob.list,
      textGrob(
        label = processedLabels[i],
        x = label.matrix[i,'x'],
        y = label.matrix[i,'y'],
        gp = gpar(
          col = label.col[i],
          cex = cex[i],
          fontface = fontface[i],
          fontfamily = fontfamily[i]
        )
      )
    );
  }
  
  
  # find the location and plot all the category names
  cat.pos.x <- c(0.18, 0.82, 0.35, 0.65);
  cat.pos.y <- c(0.58, 0.58, 0.77, 0.77);
  
  for (i in 1:4) {
    
    # work out location of the category label
    this.cat.pos <- find.cat.pos(
      x = cat.pos.x[i],
      y = cat.pos.y[i],
      pos = cat.pos[i],
      dist = cat.dist[i]
    );
    
    # then print it
    grob.list <- gList(
      grob.list,
      textGrob(
        label = category[i],
        x = this.cat.pos$x,
        y = this.cat.pos$y,
        just = cat.just[[i]],
        gp = gpar(
          col = cat.col[i],
          cex = cat.cex[i],
          fontface = cat.fontface[i],
          fontfamily = cat.fontfamily[i]
        )
      )
    );
  }
  
  # adjust grob.list to fit and return grob.list
  grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
  
  # draw diagram before returning gList is specified by user
  if (ind) { grid.draw(grob.list); }
  return(grob.list);
}
