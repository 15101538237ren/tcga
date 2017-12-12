#install.packages(c("hash","fitdistrplus","dplyr"))
library(parallel)
library(hash)
library(fitdistrplus)
library(dplyr)
cancer_name_list = c("BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA") #% c("COAD")#
threads_number = length(cancer_name_list)

generate_pvalue_table_pipeline = function(cancer_name)
{
  cancer_names = c("BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "THCA")
  setwd("/disk/tcga/") # ~/PycharmProjects/tcga
  base_dir <- getwd()
  gene_idx_fp = file.path(base_dir,"global_files","gene_idx.txt")
  gene_names = read.table(gene_idx_fp, header=FALSE, stringsAsFactors = FALSE) 
  
  sig_level = -10
  
  gene_id_idx = 1
  gene_names_idx = 2
  gene_num = length(gene_names[[1]])
  gene_list_len = length(gene_names[[gene_names_idx]])
  data_frame_len = gene_list_len
  
  methy_data_dir = file.path(base_dir,"data","intermediate_file","methy_intermidiate","merged_stage")
  output_data_dir = file.path(base_dir,"data","intermediate_file","methy_pvalue","merged_stage")
  
  pvalue_positive_name_end = "_pp_value.dat"
  pvalue_negtive_name_end = "_pn_value.dat"
  
  score_file_positive_name_end = "_p_score.dat"
  score_file_negtive_name_end = "_n_score.dat"
  
  invalid_pvalue = 10 #set the invalid pvalue output into the dat file
  extreme_pvalue = -20 #set the extreme pvalue to represent p-value = 0, log(p-value)=-inf
  
  df_idx = 5 # df_col_index_start_of_data
  
  if(!file.exists(output_data_dir))
  {
    dir.create(output_data_dir,recursive = T)
    print(sprintf("create %s successful!", output_data_dir))
  }
  
  stage_name_list = c("normal", "i")
  get_cancer_idx = function(cancer_name)
  {
    for(i in 1:length(cancer_names))
    {
      if(cancer_name == cancer_names[i])
      {
        return (i)
      }
    }
    return (-1)
  }
  
  #less: not less than, >
  #great: not greater than, <
  beta_hypothesis_test = function(x, shape1, shape2, alternative)
  {
    beta_test_res = ks.test(x, "pbeta", shape1, shape2, alternative=alternative)
    #print(beta_test_res)
    
    return (beta_test_res$p.value)
  }
  
  p_value_calc = function(x, shape1, shape2, alternative)
  {
    if(alternative=='greater')
    {
      return (pbeta(x, shape1, shape2))
    }
    else if(alternative == 'less')
    {
      return (1.0 - pbeta(x, shape1, shape2))
    }
  }
  
  
  get_info_list = function(file_name)
  {
    str_list = strsplit(file_name, "[.]")
    file_name_pre = str_list[[1]][1]
    info_list = strsplit(file_name_pre, "_")
    return (info_list[[1]])
  }
  
  get_stage_name = function(file_name)
  {
    info_list = get_info_list(file_name)
    stage_idx = get_info_idx("stage")
    stage_name = info_list[stage_idx]
    return (stage_name)
  }
  
  get_cancer_name = function(file_name)
  {
    info_list = get_info_list(file_name)
    cancer_idx = get_info_idx("cancer")
    cancer_name = info_list[cancer_idx]
    return (cancer_name)
  }
  
  checkValue = function(methy_vector, num)
  {
    if(sum(methy_vector == -1) == num)
      return (TRUE)
    return (FALSE)
  }
  maefun <- function(pred, obs)
  {
    return (mean(abs(pred - obs)))
  }
  
  print(sprintf("start %s", cancer_name))
  if(get_cancer_idx(cancer_name) == -1)
    next
  
  normal_file_name = sprintf("%s/%s/%s_normal_methy_dat.dat", methy_data_dir, cancer_name, cancer_name)
  i_th_file_name = sprintf("%s/%s/%s_i_methy_dat.dat", methy_data_dir, cancer_name, cancer_name)
  
  if(file.exists(normal_file_name) && file.exists(i_th_file_name))
  {
    normal_data_frame = read.table(normal_file_name, header=TRUE)
    i_th_data_frame = read.table(i_th_file_name, header=TRUE)
    
    cancer_df_mp_score = data.frame(GeneIds=gene_names[[gene_id_idx]], LogNum=rep(0,data_frame_len), 
                                    TotalNums=rep(0,data_frame_len), Score=rep(0,data_frame_len))
    for(i in 2:length(names(i_th_data_frame)))
    {
      sample_id = names(i_th_data_frame)[i]
      cancer_df_mp_score[sample_id] = rep(0, data_frame_len)
    }
    #The difference of cancer_df_mp_score and cancer_df_mn_score is Score.
    #cancer_df_mp_score's score is M+ score, cancer_df_mn_score's score is M- score.
    cancer_df_mn_score = cancer_df_mp_score  
    
    #Length of cancer_sample_dataframe should be equal to length of cancer_df_mp_score 
    print(nrow(cancer_df_mp_score))
    
    for(i in 1:nrow(cancer_df_mp_score))
    {
      gene_id_item = cancer_df_mp_score[i, "GeneIds"]
      normal_frame_idx = gene_id_item
      i_th_frame_idx = gene_id_item
      normal_methy_frame = normal_data_frame[normal_frame_idx, 2:ncol(normal_data_frame)]   #List
      i_th_methy_frame = i_th_data_frame[i_th_frame_idx, 2:ncol(i_th_data_frame)]   #List
      
      normal_sample_num = ncol(normal_data_frame) - 1
      i_th_sample_num = ncol(i_th_data_frame) - 1
      
      #Change List to vector (origin vector has -1 values, just NA value)
      origin_normal_methy = as.numeric(normal_methy_frame)
      origin_i_th_methy = as.numeric(i_th_methy_frame)
      
      if(sum(is.na(origin_normal_methy)) > 0 || sum(is.na(origin_i_th_methy)) > 0)
      {
        next
      }
      
      if(checkValue(origin_normal_methy, normal_sample_num) || checkValue(origin_i_th_methy, i_th_sample_num) 
         || checkValue(origin_normal_methy, normal_sample_num - 1))
      {
        cancer_df_mp_score[i, df_idx : ncol(cancer_df_mp_score)] = rep(invalid_pvalue, ncol(cancer_df_mp_score) - df_idx + 1)
        cancer_df_mp_score[i, "LogNum"] = 1
        cancer_df_mp_score[i, "TotalNums"] = -1
        cancer_df_mp_score[i, "Score"] = -1
        
        cancer_df_mn_score[i, df_idx:ncol(cancer_df_mp_score)] = rep(invalid_pvalue, ncol(cancer_df_mp_score) - df_idx + 1)
        cancer_df_mn_score[i, "LogNum"] = 1
        cancer_df_mn_score[i, "TotalNums"] = -1
        cancer_df_mn_score[i, "Score"] = -1
        next
      }
      
      #calc effective sample num (the num of samples whose beta-value > 0)
      #guarantee the effective sample num > 0
      effective_normal_num = sum(origin_normal_methy != -1)
      effective_i_th_num = sum(origin_i_th_methy != -1)
      
      #calc the effective beta value vector
      normal_idx_subset = which(origin_normal_methy != -1)
      normal_methy = origin_normal_methy[normal_idx_subset]
      
      i_th_idx_subset = which(origin_i_th_methy != -1)
      i_th_methy = origin_i_th_methy[i_th_idx_subset]
      
      result = tryCatch ({
        fit_res = fitdist(normal_methy, "beta")
      },error=function(e){
        cat("fitdist error: i=", i, "\n")
      })
      fit_estimate_res = fit_res$estimate
      fit_shape1 = fit_estimate_res["shape1"] #p
      fit_shape2 = fit_estimate_res["shape2"] #q
      #plot(fit_res)
      #print(fit_res)
      
      i_th_array = array(i_th_methy, c(length(i_th_methy), 1))
      #print(i_th_array)
      alternative_up = "less"
      alternative_down = "greater"
      
      # the ks.test result p-value list is equal to pbeta value list
      # when ks.test return p-value = NaN, the pbeta result = 1.00000
      
      p_val_list_up = apply(i_th_array, 1, p_value_calc, shape1=fit_shape1, shape2=fit_shape2, alternative=alternative_up)
      p_val_list_down = apply(i_th_array, 1, p_value_calc, shape1=fit_shape1, shape2=fit_shape2, alternative=alternative_down)
      
      p_val_list_up = log10(p_val_list_up)
      p_val_list_down = log10(p_val_list_down)
      
      p_val_list_up[p_val_list_up == -Inf] = extreme_pvalue
      p_val_list_down[p_val_list_down == -Inf] = extreme_pvalue
      
      sub_num1 = sum(p_val_list_up < sig_level)
      sub_num2 = sum(p_val_list_down < sig_level)
      
      final_p_val_list_up = rep(invalid_pvalue, i_th_sample_num)
      final_p_val_list_up[i_th_idx_subset] = p_val_list_up
      
      final_p_val_list_down = rep(invalid_pvalue, i_th_sample_num)
      final_p_val_list_down[i_th_idx_subset] = p_val_list_down
      
      cancer_df_mp_score[i, df_idx:ncol(cancer_df_mp_score)] = final_p_val_list_up
      cancer_df_mp_score[i, "LogNum"] = sub_num1
      cancer_df_mp_score[i, "TotalNums"] = effective_i_th_num
      cancer_df_mp_score[i, "Score"] = sub_num1 / effective_i_th_num
      
      cancer_df_mn_score[i, df_idx:ncol(cancer_df_mp_score)] = final_p_val_list_down
      cancer_df_mn_score[i, "LogNum"] = sub_num2
      cancer_df_mn_score[i, "TotalNums"] = effective_i_th_num
      cancer_df_mn_score[i, "Score"] = sub_num2 / effective_i_th_num
      
      if(i %% 100 == 0)
      {
        cat(i,"/", data_frame_len,":", i/data_frame_len," is finished!\n")
      }
    }
    
    cancer_dir_path = sprintf("%s/%s", output_data_dir, cancer_name)
    
    if(!file.exists(cancer_dir_path))
    {
      dir.create(cancer_dir_path,recursive = T)
    }
    
    #---------------------------------------------
    cancer_pvalue_positive_file_name =  sprintf("%s/%s%s", cancer_dir_path, cancer_name, pvalue_positive_name_end)
    cancer_df_pp_values = cancer_df_mp_score[c(1, df_idx:ncol(cancer_df_mp_score))]
    names(cancer_df_pp_values) = seq(0, i_th_sample_num)
    write.table(cancer_df_pp_values, file=cancer_pvalue_positive_file_name, sep = "\t", row.names=F, quote=F)
    
    cancer_pvalue_negtive_file_name =  sprintf("%s/%s%s", cancer_dir_path, cancer_name,pvalue_negtive_name_end)
    cancer_df_pn_values = cancer_df_mn_score[c(1, df_idx:ncol(cancer_df_mp_score))]
    names(cancer_df_pn_values) = seq(0, i_th_sample_num)
    write.table(cancer_df_pn_values, file=cancer_pvalue_negtive_file_name, sep = "\t", row.names=F, quote=F)
    
    #
    score_file_positive_file_name =  sprintf("%s/%s%s", cancer_dir_path, cancer_name,score_file_positive_name_end)
    cancer_df_pp_scores = cancer_df_mp_score[c(1:df_idx-1)]
    write.table(cancer_df_pp_scores, file=score_file_positive_file_name, sep = "\t", row.names=F, col.names=F, quote=F)
    
    score_file_negtive_file_name =  sprintf("%s/%s%s", cancer_dir_path, cancer_name,score_file_negtive_name_end)
    cancer_df_pn_scores = cancer_df_mn_score[c(1:df_idx-1)]
    write.table(cancer_df_pn_scores, file=score_file_negtive_file_name, sep = "\t", row.names=F, col.names=F, quote=F)
  }
}

cl <- makeCluster(threads_number)
results <- parLapply(cl, cancer_name_list, generate_pvalue_table_pipeline)
stopCluster(cl)
