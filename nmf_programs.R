##############################################################
################## LOADING LIBRARIES #########################
##############################################################
{
  library(RcppML)
  library(tidyverse)
  library(extraoperators)
  library(ggcorrplot)
  library(fgsea)
  library(data.table)
  library('org.Hs.eg.db')
  library("clusterProfiler")
  library(msigdbr)
  library(qlcMatrix)
  library(data.table)
  library(scater)
  library(RColorBrewer)
  library(pheatmap)
  library(ComplexHeatmap)
  library(ggridges)
  library(ggplot2)
  library(ggstatsplot)
  library(palmerpenguins)
  library(yarrr)
  library(readxl)
  library(ggrepel)
}

##############################################################
################## LOADING DATA ##############################
##############################################################
{
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-05_clb22_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-06_clb25_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-06_clb36neg+tot_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-06-14_clb24neg+tot_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-06_clb29neg+tot_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-05_clb31_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-04-05_clb33_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-06-08_clb38neg+tot_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-06-08_clb39_Epi_SCE.RData")
  load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-06-08_clb42neg+tot_Epi_SCE.RData")
  #load("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/data/2022-07-13_clb42neg+tot_HighQuality_Epi_SCE.RData")
  
  
  sce_22 <- sce22Epi
  sce_25 <- sce25Epi
  sce_36 <- sceTum36
  sce_24 <- sceTum24
  sce_29 <- sceTum29
  sce_31 <- sce31Epi
  sce_33 <- sce33Epi
  sce_38 <- sceTum38
  sce_39 <- sce39Epi
  sce_42 <- sceTum42
  #sce_hq <- sceTum42
  
  
}

##############################################################
################## FUNCTIONS #################################
##############################################################
{
  #' @title Compute the jaccard index 
  #' @name jaccard
  #' @description Computes the jaccard index between two named vectors
  #' @param a a named  vector of genes
  #' @param b a named vector of genes
  #' @return A value corresponding to the jaccard index, bounded between 0 and 1.
  
  jaccard <- function(a, b) {
    intersection = length(intersect(names(a), names(b))) # compute the number of shared genes between a and b 
    union = length(names(a)) + length(names(b)) - intersection #compute genes from a and b minus the intersection 
    return (intersection/union)  #jaccard index defined by intersection/union
  }
  
  
  #' @title Compute the number of NMF programs 
  #' @name compute_num_programs
  #' @description Computes the number of programs in a certain list
  #' @param NMFs list of NMF programs classified by patient and k 
  #' @return An integer value corresponding to the number of NMF programs contained in the list 
  
  compute_num_programs <- function(NMFs){
    num_programs <- 0 #initial counter
    for( i in 1:length(NMFs)){ 
      for(j in 1:length(NMFs[[i]])){
        for(s in 1: length(NMFs[[i]][[j]])){
          num_programs <- num_programs + 1 #update counter
        } #close s loop
      } #close j loop
    } #close i loop 
    return(num_programs)
  }
  
  
  #' @title Compute a set of NMF programs
  #' @name compute_programs
  #' @description Computes a set of NMF programs based on a Non-Negative Matrix Factorization (NMF)
  #' @param sce A SingleCellExperiment object 
  #' @param k_range numeric vector that defines the range for k 
  #' @return A list containing all the NMF programs computed for each different value of k 
  
  compute_programs <- function(sce,k_range) {
    #Accessing the count data
    sce_counts = SummarizedExperiment::assay(sce,"logcounts")  #Accessing the count matrix from the single cell object 
    sce_counts_matrix <<- as.matrix(sce_counts)  #Casting object into a matrix 
    
    sce_counts_matrix <- sce_counts_matrix[-str_which(rownames(sce_counts_matrix),"^MT-"),] # Deleting mitochondrial genes from matrix 
    sce_counts_matrix <- sce_counts_matrix[-str_which(rownames(sce_counts_matrix),"^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"),] # Deleting ribosomal genes from matrix 
    
    
    # colnames(sce_counts_matrix) <- make.unique(colnames(sce_counts_matrix))
    # sce_counts_matrix <- as.data.frame(sce_counts_matrix) 
    # sce_counts_matrix <- sce_counts_matrix %>% mutate(Egenes = rowSums(.)/nrow(sce_counts_matrix)) %>% arrange(desc(Egenes)) %>% slice(1:7000) %>% select(-Egenes)
    # sce_counts_matrix <- normalizeCounts(sce_counts_matrix)
  
    list_h <- list()
    
    NMF_Programs <- list() #list that will store the programs for each k 
    for(i in k_range){
      #Apply NMF to a count matrix for a given K
      # if(ncol(sce_counts_matrix)<100){
      temp_nmf <- nmf(sce_counts_matrix,k=i,tol = 1e-5,seed=10)
      # }
      # else{
      #   set.seed(100)
      #   sce_counts_matrix <- sce_counts_matrix[,sample(ncol(sce_counts_matrix),size=100,replace=FALSE)]
      #   temp_nmf <- nmf(sce_counts_matrix,k=i,tol = 1e-5,seed=10) #seed set to 10 to be able to reproduce results.
      #}
      #Acceding the NMF Programs
      w <- temp_nmf$w
      h <- temp_nmf$h
      rownames(w) <- rownames(sce_counts_matrix) 
      # #Summarize each NMF Programs by the top 50 genes
      w <- lapply(seq_len(ncol(w)), function(i) w[,i])
      w <- lapply(w, function(x) head(sort(x, decreasing = TRUE), n = 50))
      names(w) <- paste0("cl",seq(from = 1, to =i)) #naming each of the k clusters
      NMF_Programs[[paste0("k",i)]] <- w  #naming each set of programs by their k 
      list_h[[paste0("k",i)]]<- h 
    }
    
    h <<- list_h
    return(NMF_Programs)
  }

  #' @title Computes a histogram of all jaccard indexes associated to each NMF Program. 
  #' @name diagnostic
  #' @description Computes a histogram of all jaccard indexes (one jaccard index for each NMF Program)
  #' @param jaccard_index_list A list with all jaccard indexes associated to all NMF Programs 
  
  diagnostic <- function(jaccard_index_list){
    listajacc <<- jaccard_index_list
    index <- numeric()
    for(i in 1:length(jaccard_index_list)){
      for(j in 1:length(jaccard_index_list[[i]])){
        for(k in 1:length(jaccard_index_list[[i]][[j]])){
          index <- c(index, jaccard_index_list[[i]][[j]][[k]])

          
        }
      }
    }
    h <- hist(index,main="Histogram of jaccard indexes",xlab="jaccard index", ylab="Frequency")
    print(h$counts)
  }

  
  #' @title Compute a robust within the sample set of NMF Programs 
  #' @name compute_robust_programs1
  #' @description Compares each program from each k to all the programs from the rest of k and keeps the ones that have a jaccard index greater than a given threshold.  
  #' @param NMF_Programs a list containing all the NMF Programs computed for all different k and for all samples. 
  #' @param k_range numeric vector that defines the range for k 
  #' @return A reduced list of all NMF programs robust within the sample for all samples. 
  
  
  #Defining a set of robust within the sample NMF programs:
  compute_robust_programs1 <- function(NMF_Programs,k_range,threshold){
    jaccard_index_list <- set_values_to_zero(NMF_Programs)  #copy data structure of the NMF_Programs list to associate a jaccard index to each NMF Program
    initial_nb_programs <- compute_num_programs(NMF_Programs)  #initial number of NMF Programs
    robust_NMF_Programs <- list() #list that will store the programs that are robust within the sample  
    for(p in 1:length(NMF_Programs)){    # p browses all different patients
      for(i in 1:(length(k_range)-1)){   # i browses all different k
        for(j in (i+1):length(k_range)){  # once a k chosen with i, j browses the rest of k. We do not compare programs within the same k 
          l1 = length(NMF_Programs[[p]][[i]])  # l1 keeps the number of clusters in i 
          l2 = length(NMF_Programs[[p]][[j]])  # l2 keeps the number of clusters in j 
          for(m in 1:l1){                 # m browses all different clusters of i 
            for(s in 1:l2){               # s browses all different clusters of j
              
              #Condition to keep an NMF Program within the sample: Jaccard index > 0.7
              #jaccard_index <- length(intersect(names(NMF_Programs[[p]][[i]][[m]]),names(NMF_Programs[[p]][[j]][[s]])))
              jaccard_index <- jaccard(NMF_Programs[[p]][[i]][[m]],NMF_Programs[[p]][[j]][[s]])
              jaccard_index_list[[p]][[i]][[m]] <- max(jaccard_index_list[[p]][[i]][[m]],jaccard_index)  #Update jaccard index value of program into jaccard index list. 
              jaccard_index_list[[p]][[j]][[s]] <- max(jaccard_index_list[[p]][[j]][[s]],jaccard_index)  #Update jaccard index value of program into jaccard index list. 
              if( jaccard_index > threshold){
                robust_NMF_Programs[[names(NMF_Programs)[[p]]]][[names(NMF_Programs[[p]])[i]]][[names(NMF_Programs[[p]][[i]])[m]]] <- NMF_Programs[[p]][[i]][[m]] #updating the list of robust programs, using names to keep the right k and cluster names
                robust_NMF_Programs[[names(NMF_Programs)[[p]]]][[names(NMF_Programs[[p]])[j]]][[names(NMF_Programs[[p]][[j]])[s]]] <- NMF_Programs[[p]][[j]][[s]] #updating the list of robust programs, using names to keep the right k and cluster names
              } # close if statement
            } #close s loop
          } #close m loop
        } #close j loop
      } #close i loop 
    } #close p loop
    diagnostic(jaccard_index_list)   #diagnostic function that computes the histogram of jaccard indexes   
    final_nb_programs <- compute_num_programs(robust_NMF_Programs)  # final number of NMF Programs 
    print(paste0("The initial number of NMF programs was ", initial_nb_programs, " , there are ", final_nb_programs, " NMF programs left, " , initial_nb_programs - final_nb_programs,  " were deleted"))
    
    return(robust_NMF_Programs) 
  }
  
  #' @title Computing a list that contains all NMF programs from different samples where all cluster values are set to zero 
  #' @name set_values_to_zero
  #' @description Sets all the cluster from all NMF programs from different samples to zero. 
  #' @param NMFs a list containing all the NMF Programs computed for all different k and samples.  
  #' @return A list containing all the NMF Programs for all different k and samples with the cluster values set to zero.  
  
  set_values_to_zero <- function(NMFs){ 
    for(i in 1:length(NMFs)){   #i browses different samples
      for(j in 1:length(NMFs[[i]])){   #j browses different k 
        for(s in 1:length(NMFs[[i]][[j]])){   #s browses different clusters 
          NMFs[[i]][[j]][[s]] <- 0    #set cluster value to 0
        } #close s loop
      } #close j loop
    } #close i loop
    return(NMFs)
  }
  
  #' @title Computing a data frame from the overlap scores list. 
  #' @name compute_scoresdf 
  #' @description Computes a data frame from an overlap scores list. The data frame contains the patient, the k and the cluster associated to each NMF Program 
  #' @param overlap_scores a list of all NMF Programs for all different samples and k with all cluster values set to 0. This list will keep the accumulated scores for each program. 
  #' @return A 4 column data frame that contains the patient, the k, the cluster and the score associated to each NMF Program. The data frame is group by patient and sorted by score in decreasing order. 
  
  compute_scoresdf <- function(overlap_scores){
    overlap_scores_df <- data.frame(patient=character(),k=character(),cluster=character(),score=numeric()) #initializing an empty data frame with 4 columns: patient, k, cluster and score.
    for(i in 1:length(overlap_scores)){ # i browses all different patients 
      for(j in 1:length(overlap_scores[[i]])){ # j browses all different k's of i 
        for(s in 1:length(overlap_scores[[i]][[j]])){ # s browses all different clusters of j
          overlap_scores_df <- overlap_scores_df %>% add_row(patient=names(overlap_scores)[i], k=names(overlap_scores[[i]])[j], cluster=names(overlap_scores[[i]][[j]])[s], score = overlap_scores[[i]][[j]][[s]]) #add row into the data frame with all the NMF Program information(patient, k cluster and score)
        } #end s loop
      } #end j loop
    }#end i loop 
    overlap_scores_df <- overlap_scores_df  %>% group_by(patient)  %>%  arrange(desc(score),.by_group = TRUE) %>% filter(score > 0)
    return(overlap_scores_df)
  }
  
  
  #' @title Compute a robust between samples set of NMF Programs 
  #' @name compute_robust_programs2
  #' @description Compares each program of each sample to all the programs from the rest of the samples and keeps the ones that have a jaccard index greater than a given threshold. Keeps the accumulated scores by accumulating the jaccard indexes for each program. A data frame called overlap_scores_df is saved in the environment with the score associated to each NMF Program. 
  #' @param NMFs a list containing all the NMF Programs computed for all different samples and k.
  #' @return A reduced list of robust between samples NMF programs
  
  compute_robust_programs2 <- function(NMFs,threshold){
    jaccard_index_list <- set_values_to_zero(NMFs) #copy data structure of the NMF_Programs list to associate a jaccard index to each NMF Program
    initial_nb_programs <- compute_num_programs(NMFs) #initial number of NMF Programs
    overlap_scores <- set_values_to_zero(NMFs) #list of all NMF Programs for all different samples and k with all named vectors containing genes defined by the cluster set to 0
    robust_NMF_Programs <- list() #list that will store the robust between sample NMF Programs 
    for (i in 1:(length(NMFs)-1)){  # i browses all the different samples 
      for (j in (i+1):length(NMFs)){ # once a sample chosen with i, j browses the rest of the samples. We do not compare programs within the same sample.  
        l1 = length(NMFs[[i]])      # l1 keeps the number of k in i 
        l2 = length(NMFs[[j]])      # l2 keeps the number of k in j
        for(m in 1:l1){             # m browses all different k's of i
          for(n in 1:l2){           # n browses all different k's of j 
            l3 = length(NMFs[[i]][[m]])  #l3 keeps the number of clusters in m
            l4 = length(NMFs[[j]][[n]])  #l4 keeps the number of clusters in n 
            for(q in 1:l3){         # q browses all different clusters of m
              for(s in 1:l4){       # s browses all different clusters of n
                #overlap_score <- length(intersect(names(NMFs[[i]][[m]][[q]]),names(NMFs[[j]][[n]][[s]])))
                overlap_score <- jaccard((NMFs[[i]][[m]][[q]]),(NMFs[[j]][[n]][[s]])) #computing the jaccard index between two programs 
                #Condition to keep an NMF Program between samples : Jaccard index > 0.2
                jaccard_index_list[[i]][[m]][[q]] <- max(jaccard_index_list[[i]][[m]][[q]],overlap_score)
                jaccard_index_list[[j]][[n]][[s]] <- max(jaccard_index_list[[j]][[n]][[s]],overlap_score)
                if(overlap_score > threshold){  
                  robust_NMF_Programs[[names(NMFs)[i]]][[names(NMFs[[i]])[m]]][[names(NMFs[[i]][[m]])[q]]] <- NMFs[[i]][[m]][[q]] #updating the list of robust programs, using names to keep the right k and cluster names
                  robust_NMF_Programs[[names(NMFs)[j]]][[names(NMFs[[j]])[n]]][[names(NMFs[[j]][[n]])[s]]] <- NMFs[[j]][[n]][[s]] #updating the list of robust programs, using names to keep the right k and cluster names
                  overlap_scores[[i]][[m]][[q]] <- sum(overlap_score,overlap_scores[[i]][[m]][[q]]) #updating the list of overlap scores by accumulating the scores
                  overlap_scores[[j]][[n]][[s]] <- sum(overlap_score,overlap_scores[[j]][[n]][[s]]) #updating the list of overlap scores by accumulating the scores 
                } #close if statement
              } #close s loop
            } #close q loop
          }  #close n loop
        } #close m loop
      } #close j loop
    } #close i loop
    overlap_scores_df <<- compute_scoresdf(overlap_scores) # data frame containing all scores grouped by patient and sorted by scores in decreasing order.  
    diagnostic(jaccard_index_list) #diagnostic function that computes the histogram of jaccard indexes   
    final_nb_programs <- compute_num_programs(robust_NMF_Programs)  # final number of NMF Programs
    print(paste0("The initial number of NMF programs was ", initial_nb_programs, " , there are ", final_nb_programs, " NMF programs left, " , initial_nb_programs - final_nb_programs,  " were deleted"))
    
    return(robust_NMF_Programs)
  }
  
  
  #' @title Compute a non_redundant within the sample set of NMF Programs 
  #' @name compute_robust_programs3
  #' @description For each sample, compares each of the sample programs to all the programs from the rest of the sample and removes the ones that have a jaccard index of at least 0.2.
  #' @param overlap_scores_df a data frame with all NMF Programs scores grouped by patient and arranged by score in decreasing order. 
  #' @param robust_NMFs a list of NMF Programs robust between and within samples
  #' @return A list of non-redundant within the sample NMF Programs 
  
  compute_robust_programs3 <- function(overlap_scores_df,robust_NMFs,threshold){
    jaccard_index_list <- set_values_to_zero(robust_NMFs) #copy data structure of the NMF_Programs list to associate a jaccard index to each NMF Program
    initial_nb_programs <- compute_num_programs(robust_NMFs) #initial number of NMF Programs
    l <- list()
    updated_robust_NMFs <- list()
    overlap_scores_df <- overlap_scores_df %>% mutate(keep=TRUE)
    
    for( i in unlist(unique(overlap_scores_df["patient"]))){
      temp_df <- overlap_scores_df %>% filter(patient == i)
      for(j in 1:(nrow(temp_df)-1)){
        for(s in (j+1):nrow(temp_df)){
          current_nmf <- robust_NMFs[[i]][[unlist(temp_df[j,'k'])]][[unlist(temp_df[j,'cluster'])]]
          iterating_nmf <- robust_NMFs[[i]][[unlist(temp_df[s,'k'])]][[unlist(temp_df[s,'cluster'])]]
          #jaccard_index <- length(intersect(names(current_nmf),names(iterating_nmf)))
          jaccard_index <- jaccard(current_nmf,iterating_nmf)
          if(jaccard_index > threshold){
            temp_df[s,'keep'] <- FALSE
          }
          # if(temp_df[j,'keep'] == FALSE){
          #   jaccard_index_list[[i]][[unlist(temp_df[j,'k'])]][[unlist(temp_df[j,'cluster'])]] <- max(jaccard_index_list[[i]][[unlist(temp_df[j,'k'])]][[unlist(temp_df[j,'cluster'])]],jaccard_index)
          # }
          jaccard_index_list[[i]][[unlist(temp_df[j,'k'])]][[unlist(temp_df[j,'cluster'])]] <- max(jaccard_index_list[[i]][[unlist(temp_df[j,'k'])]][[unlist(temp_df[j,'cluster'])]],jaccard_index)
          jaccard_index_list[[i]][[unlist(temp_df[s,'k'])]][[unlist(temp_df[s,'cluster'])]] <- max(jaccard_index_list[[i]][[unlist(temp_df[s,'k'])]][[unlist(temp_df[s,'cluster'])]],jaccard_index)
          
        }
  
      }
      temp_df <-temp_df %>% filter(keep == TRUE)
      for(l in 1:nrow(temp_df)){
        updated_robust_NMFs[[unlist(temp_df[l,'patient'])]][[unlist(temp_df[l,'k'])]][[unlist(temp_df[l,'cluster'])]] <- robust_NMFs[[unlist(temp_df[l,'patient'])]][[unlist(temp_df[l,'k'])]][[unlist(temp_df[l,'cluster'])]]
        jaccard_index_list[[unlist(temp_df[l,'patient'])]][[unlist(temp_df[l,'k'])]][[unlist(temp_df[l,'cluster'])]]  <- (threshold - 1)  
      }
    }
    
    final_nb_programs <- compute_num_programs(updated_robust_NMFs)  #final number of NMF Programs 
    diagnostic(jaccard_index_list) #diagnostic function that computes the histogram of jaccard indexes 
    print(paste0("The initial number of NMF programs was ", initial_nb_programs, " , there are ", final_nb_programs, " NMF programs left, " , initial_nb_programs - final_nb_programs,  " were deleted"))
    return(updated_robust_NMFs)
  }
  
  
  #' @title Compute a reduce list of robust NMF Programs for each sample. 
  #' @name reduce_NMF_programs
  #' @description Reduces a list of NMF Programs for all samples and all possible k values to a list of NMF Programs where only the sample value is kept. 
  #' @param robust_NMFs a list of NMF Programs robust between and within samples
  #' @return A list of robust within the sample, between samples and non-redundant within the sample NMF Programs.  
  
  reduce_NMF_programs <- function(robust_NMFs){
    num_program = 1  # initial counter of number of robust NMF Programs 
    robust_NMF_clustering <- list()
    for(i in 1:length(names(robust_NMFs))){
      for(j in 1:length(names(robust_NMFs[[i]]))){
        for(s in 1:length(names(robust_NMFs[[i]][[j]]))){
          robust_NMF_clustering[[paste0(str_sub(names(robust_NMFs)[[i]],-2),"_",names(robust_NMFs[[i]])[[j]],"_",names(robust_NMFs[[i]][[j]])[[s]])]] <- robust_NMFs[[i]][[j]][[s]] #return one list with all robust NMF Programs associated to their sample the number of robust NMF Program. 
          num_program = num_program +1
        }
      }
    }
    return(robust_NMF_clustering)
  }
  
  #' @title Compute the number of genes shared by two list 
  #' @name gene_overlap
  #' @description Computes the number of genes overlap two named vectors
  #' @param a a named  vector of genes
  #' @param b a named vector of genes
  #' @return A value corresponding to the number of shared genes
  
  gene_overlap <- function(a, b) {
    return (length(intersect(names(a), names(b))))  #jaccard index defined by intersection/union
  }
  
  
  #' @title Compute a pairwise data frame with jaccard indexes between all robust programs. 
  #' @name compute_pairwise_df
  #' @description Computes a pairwise data frame from a list of robust NMF programs, the values of the dataframe are the jaccard indexes between all programs.
  #' @param robust_NMFs a list of NMF Programs robust between and within samples
  #' @return A list of robust within the sample, between samples and non-redundant within the sample NMF Programs.  
  
  compute_pairwise_df <- function(robust_NMFs,threshold){
    c <- vector()   
    for(i in names(robust_NMFs)){
      c <- append(c,i) #get all robust NMF Programs names and keep them in a vector
    }
    
    df <- data.frame(matrix(ncol = length(c), nrow = length(c)))  #create an original pairwise data frame
    rownames(df) <- c 
    colnames(df) <- c
    
    for(i in 1:ncol(df)){ 
      df[i,] <- sapply(robust_NMFs,jaccard,robust_NMFs[[i]]) #set all elements of the matrix to the jaccard index value between row and column value 
      df[i,i] <- 0 #setting diagonal values to 0 
    }
    temp_df <- df   #temporary data frame 
    temp_df[temp_df < threshold] <- 0 #set all values smaller than the threshold to 0
    temp_df <- temp_df %>% add_row(summarise(., across(where(is.numeric), sum)))  #add row that adds up all row values per column 
    row_to_add <- tail(temp_df,n=1) #store the recently added row into an object 
    df <- rbind(df,row_to_add)  # bind the row into the original data frame 
    rownames(df)[nrow(df)] <-  "Total overlapping genes"   #name row
    
    temp_row <- apply(apply(df[-nrow(df),],2,function(x) x>threshold),2,sum)  #check and count how many row values are greater than a threshold by column 
    df <- rbind(df,temp_row)   #add row to original pairwise data frame
    rownames(df)[nrow(df)] <-  "Total overlapping events"  #name row 
    
    return(df)
  }
  
  
  #' @title Compute all Meta Programs from a set of robust NMF Programs. 
  #' @name compute_metaprograms
  #' @description Computes all MP based on a defined clustering algorithm done by the function compute_metaprogram
  #' @param robust_NMFs a list of NMF Programs robust between and within samples
  #' @return A list containing all Meta Programs.  
  
  compute_metaprograms <- function(robust_NMFs,threshold1,threshold2,max_jac){
    MP_list <- list()   #Initialize list that will contain all MP 
    robust_NMFs_m <- robust_NMFs  #Copy list of robust NMF programs to modify it
    pass <- TRUE  #variable pass will confirm that the clustering condition is true 
    #pairwisedf <<- compute_pairwise_df(robust_NMFs_m,threshold1)
    
    while(pass == TRUE && length(robust_NMFs_m)>=threshold2){  #Verify that there are at least 4 nmf programs to be compared. 
      pairwise_df <- compute_pairwise_df(robust_NMFs_m,threshold1)  # Initial pairwise matrix to find cluster founder 
      values <- pairwise_df["Total overlapping events",]  #Get all columns that correspond to all overlapping events
      max_overlapping_events_p <- which.max(values) #Get column number of program that has max overlapping events
      max_overlapping_events <- pairwise_df["Total overlapping events", max_overlapping_events_p]  # Get value of max overlapping events 
      
      if(max_overlapping_events >= threshold2){   
        cl_founder_name <- colnames(pairwise_df)[max_overlapping_events_p] #Get cluster founder name 
        cl_founder <- robust_NMFs[[cl_founder_name]] #Get cluster founder from original robust NMF list
      }else{
        pass <- FALSE
        break
      }
      
      values2 <- pairwise_df[cl_founder_name,] #Get all jaccard indexes between cluster founder and the rest of NMF programs 
      max_jaccard_index_p <- which.max(values2) #Get column number of program that shares the most genes with the cluster founder
      max_jaccard_index <- pairwise_df[cl_founder_name,max_jaccard_index_p] #Get the value of the max jaccard index between cluster founder and the nmf program 
      
      cl_2_name <- colnames(pairwise_df)[max_jaccard_index_p] #Get cluster member name
      cl_2 <- robust_NMFs[[cl_2_name]]  #Get cluster member from original robust NMFs 
      
      all_programs_mp.df <- data.frame(c(names(cl_founder),names(cl_2)),c(cl_founder,cl_2),row.names = NULL) #Data frame that will keep all genes from the MP(cluster founder and nmf program)
      names(all_programs_mp.df) <- c("gene","score")
      
      cluster.df <- all_programs_mp.df
      cluster.df <- cluster.df  %>% group_by(gene)  %>% mutate(count=n())  %>%  ungroup() %>%  #group by gene and count the number of occurrences
                    arrange(desc(count),desc(score))  %>%  #sort by score and number of occurences
                    distinct(gene,.keep_all = TRUE) #remove duplicates and keep genes with the highest score
      
      MP <- setNames(cluster.df$score[1:50],cluster.df$gene[1:50]) #Defining a MP as a named vector containing the top 50 genes of the data frame 
      
      robust_NMFs_m[["MP"]] <- MP #Add MP as a named vector into our copy of NMF list
      robust_NMFs_m[[cl_founder_name]] <- NULL #Remove cluster founder from the copy of robust NMF list
      robust_NMFs_m[[cl_2_name]] <- NULL #Remove cluster member from the copy of robust NMF list
      
      MPname <- paste0(cl_founder_name,".",cl_2_name)  #Keep the names of the robust NMF programs that contribute to the MP
      
      pairwise_df_m <- compute_pairwise_df(robust_NMFs_m,threshold1) #Define new pairwise matrix with the MP included 
      values <- pairwise_df_m["MP",] #Get the jaccard indexes between all NMF programs and the MP
      max_jaccard_index_p <- which.max(values) #Get column number of program that shares the most genes with the MP
      max_jaccard_index <- pairwise_df_m["MP",max_jaccard_index_p] #Get the value of the max jaccard index between MP and the nmf program
      
      while(max_jaccard_index>max_jac && length(robust_NMFs_m)>1){
        cl_m_name <- colnames(pairwise_df_m)[max_jaccard_index_p]  #Get NMF program name
        cl_m <- robust_NMFs[[cl_m_name]] #Get NMF program gene list
        
        # Update data frame with both MP and program with max jaccard index 
        
        df1 <- data.frame(names(cl_m),cl_m, row.names = NULL) #Create data frame with the NMF program's 50 genes
        names(df1) <- c("gene","score")
        all_programs_mp.df <- rbind(all_programs_mp.df,df1) #Add data frame with the NMF program's to the data frame that keeps record of all genes 
        
        cluster.df <- all_programs_mp.df #Copy data frame with all genes to modify it
        cluster.df <- cluster.df %>% group_by(gene) %>% mutate(count = n())  %>% ungroup() %>% #group by gene and count the number of occurrences
          arrange(desc(count),desc(score)) %>% # sort by score and number of occurences
          distinct(gene,.keep_all= TRUE) # remove duplicates and keep the genes with the highest score
        MPname <- paste0(MPname,".",cl_m_name) #Update the name of the MP by adding the name of the nmf program that is contributing to the MP 
        MP <- setNames(cluster.df$score[1:50], cluster.df$gene[1:50]) # updated MP formated as a named vector 
        
        robust_NMFs_m[["MP"]] <- MP   #Updating the MP 
        robust_NMFs_m[[cl_m_name]] <- NULL # Remove cluster member 
        
        
        if(length(robust_NMFs_m)>1){   #Check that once the list of nmf programs has been updated, the length of it is still at least 2 to avoid an inifinite loop
          pairwise_df_m <- compute_pairwise_df(robust_NMFs_m,threshold1) # Update pairwise matrix with the new MP and without the nmf program that contributed to this MP 
          values <-  pairwise_df_m["MP",]   #Get list of jaccard index between the new MP and all the left nmf programs 
          max_jaccard_index_p <- which.max(values)  #Get NMF program that shares the most genes with the new MP
          max_jaccard_index <- pairwise_df_m["MP",max_jaccard_index_p]
        } #end if statement
      } #end while loop
      
      robust_NMFs_m[["MP"]] <- NULL
      MP_list[[MPname]] <- MP
      table_df <- as.data.frame(table(cluster.df$count[1:50]),row.names=NULL)  #Keep the number of appearances for different genes
      colnames(table_df) <- c("Y times", "X genes")
      table_df <- table_df  %>% arrange(desc(`Y times`)) #Sort them by number of appearances
      print(paste0("In ", MPname, " X genes appear Y times"))
      print(table_df)
      
      
    }
    #list_pdf <<- list_pdf
    return(MP_list)
  }
  
  #' @title Compute a heat map of all meta programs.  
  #' @name compute_heatmap
  #' @description Computes all MP based on a defined clustering algorithm done by the function compute_metaprogram
  #' @param robust_NMFs a list of NMF Programs robust between andn   within samples
  #' @return A data frame containing all Meta Programs.  
  
  compute_heatmap <- function(metaprograms){
    programs <- character()  #initialize character vector
    sorted_robust_NMF <- list()  #list that will keep the sorted NMF Programs
    for(i in names(metaprograms)){   # i browses all meta programs names (each MP name contains the names of all NMF programs)
      for(j in strsplit(i,".",fixed=TRUE)){  # j browses all NMF programs contributing to a given MP
        programs <- c(programs, j)     #add nmf program name to the programs vector
        for (s in programs){
          sorted_robust_NMF[[s]] <- robust_NMFs[[s]]  #update the sorted NMF Programs list 
        }
      }
    }
    
    sorted_df <- compute_pairwise_df(sorted_robust_NMF,0.2) #compute pairwise data frame from the sorted list
    sorted_df <- sorted_df[1:(nrow(sorted_df)-2),] #select all columns except the last two("nb")
   
    for(i in 1:nrow(sorted_df)){
      sorted_df[i,i] <- 50
    }
    
    data.matrix(sorted_df)
    s_df <<- sorted_df
    
    ggcorrplot(sorted_df, show.legend = TRUE,legend.title = "Jaccard index") + scale_fill_gradientn(colours = c("#6D9EC1", "white", "#E46726"), name="Jaccard", limits = c(0,0.5), na.value = "red") #+
      # theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),
      #       axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
                                                                                                                                                    
  
  }

  
  #' @title Computes the functional enrichment of a list of meta programs.  
  #' @name compute_functional_enrichment
  #' @description Computes the functional enrichment of multiples meta programs using an universal enrichment analyzer . A list called annotations is saved into the environment, it contains all the information for each MP annotation. 
  #' @param metaprograms a list of meta programs and the genes associated with each one. 
  #' @param ga_database input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene 
  #' @return A data frame containing for each meta program, the statistically significant signatures.   
  
  compute_functional_enrichment <- function(metaprograms, ga_database){
    annotations <- list() #initialize list that will contain the statistically significant signatures 
    for(i in 1:length(metaprograms)){ # i browses all meta programs 
      #Select genes from MP:
      genes <- names(metaprograms[[i]])  # get the genes for each meta program
      set.seed(42)    # set seed to be able to reproduce results 
      ego <- enricher(genes, pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE=ga_database,minGSSize=5) #function enricher that does the gene annotation, by setting pvalueCutoff = 0.05 we know that results are statistically significant. 
      annotations[[paste0("MP",i)]] <- ego@result # populate annotations list with the enricher function results 
      
    }
    
    annotations <<- annotations # list with all results will be saved into the environment
    
    pathways_df <- do.call(cbind,lapply(annotations, function(x){y <- x$Description[1:200] ; return(y)}))  # create dataframe based on annotations list. Choose top 200 statistically significant pthways  
    
    return(pathways_df)
  }
  
  #' @title Computes the functional enrichment of a list of metaprograms.  
  #' @name compute_mp_abundance
  #' @description Computes for each meta program, the abundance of nmf programs that belong to each of the two smoking status category
  #' @param counters_patients List containing for each MP, the number of nmf programs that come from each patient. ( How many nmf programs each patient contributed to each MP)   
  #' @return A data frame containing for each meta program, the abundance of nmf programs from each ss category and an associated p value.  
  
  
  ### function uses functions compute_pvalue and compute_A ###
  
  compute_mp_abondance <- function(counters_patients){
    counter_mp <- 1
    l_mp <- lapply(counters_patients,function(x){y <- data.frame(unlist(x$summary)); y['MP'] <- rep(counter_mp,nrow(y)) ; counter_mp <<- counter_mp +1 ; return(y)})
    df_mp <- do.call(rbind.data.frame, l_mp) %>% mutate(SS = case_when(patient %in% smokers ~ "smoker" , !patient %in% smokers ~ "non smoker")) %>% arrange(MP,SS)
    rownames(df_mp) <- NULL
    df_mp <- df_mp[c("MP", "patient", "SS", "Freq")]
    df_patients <<- df_mp %>% group_by(patient,SS)  %>% summarise(ff= sum(Freq))  %>% arrange(desc(ff))
    
    
    j <- df_mp %>% group_by(SS) %>% summarise(ff= sum(Freq))
    
    df_mp <- df_mp %>% mutate(Total_nb = sum(df_mp$Freq)) %>% 
             mutate(Total_nb_SS= case_when(SS == "smoker" ~ as.numeric(j[2,2]), SS =="non smoker" ~ as.numeric(j[1,2]))) %>% 
             group_by(MP,SS,Total_nb_SS,Total_nb) %>% summarise(SS_MP = sum(Freq), .groups = "drop") %>% 
             group_by(MP) %>% mutate(Total_MP= sum(SS_MP)) %>% ungroup() %>%
             pivot_wider(names_from="SS", values_from = "SS_MP") %>% filter(!is.na(smoker))
    df_mp['non smoker'] = c(df_mp['Total_MP']-df_mp['smoker'])
    df_mp  <- df_mp %>% mutate(Total_smokers = Total_nb - Total_nb_SS)
    
    df_mp$p_value <-apply(df_mp,MARGIN = 1,compute_pvalue)
    #df_mp$A <-apply(df_mp,MARGIN = 1,compute_A)
    df_mp <- df_mp  %>% mutate(adj_pvalue = p.adjust(df_mp$p_value,method="BH")) %>% 
             rename(smokers_total = Total_nb_SS, non_smokers_total = Total_smokers, all_mp = Total_nb, total_mp = Total_MP )
    
    return(df_mp)
  }
  
 
  table_patient <- function(robust_NMFs){
    patients <- numeric()
      for(i in names(robust_NMFs)){
        counter <- numeric()
        patient <- str_sub(i,-2,-1)
        for(j in names(robust_NMFs[[i]])){
          counter <- c(counter, rep(patient,length(names(robust_NMFs[[i]][[j]]))))
        }
      patients <- c(patients,counter)
      }
    print(table(patients))
  }
}

##############################################################
######### LOADING SCE AND RUNNING NMFs #######################
##############################################################

#Patient's smoking status 

smokers <- c("22","24","29","36","38")
non_smokers <- c("31","25","33","39","42")

#Define k range
k_range <- 5:11

#Compute NMF Programs for each patient
NMF_22 <- compute_programs(sce_22,k_range)
NMF_25 <- compute_programs(sce_25,k_range)
NMF_36 <- compute_programs(sce_36,k_range) 
NMF_29 <- compute_programs(sce_29,k_range)
NMF_31 <- compute_programs(sce_31,k_range)
NMF_33 <- compute_programs(sce_33,k_range)
NMF_38 <- compute_programs(sce_38,k_range)
NMF_39 <- compute_programs(sce_39,k_range)
NMF_42 <- compute_programs(sce_42,k_range)
NMF_24 <- compute_programs(sce_24,k_range)

#Create list that contains all NMF programs for all patients
NMFs <- list()
NMFs[[paste0("NMF_",22)]] <- NMF_22
NMFs[[paste0("NMF_",25)]] <- NMF_25
NMFs[[paste0("NMF_",36)]] <- NMF_36
NMFs[[paste0("NMF_",24)]] <- NMF_24
NMFs[[paste0("NMF_",29)]] <- NMF_29
NMFs[[paste0("NMF_",31)]] <- NMF_31
NMFs[[paste0("NMF_",33)]] <- NMF_33
NMFs[[paste0("NMF_",38)]] <- NMF_38
NMFs[[paste0("NMF_",39)]] <- NMF_39
NMFs[[paste0("NMF_",42)]] <- NMF_42


#Defining a set of robust NMF Programs : 
#1. Robust within the sample

threshold_within_sample <- 0.6
robust_NMFs <- compute_robust_programs1(NMFs,k_range,threshold_within_sample)
table_patient(robust_NMFs)
#2. Robust between samples

threshold_between_samples <- 0.4
robust_NMFs <- compute_robust_programs2(robust_NMFs,threshold_between_samples)
table_patient(robust_NMFs)

#3. Non redundant within the sample 
#threshold2_within_sample <- 25
threshold2_within_sample <- 0.4
robust_NMFs <- compute_robust_programs3(overlap_scores_df,robust_NMFs,threshold2_within_sample)

#Change data structure for defining the MP !!!! Necessary to compute the metaprograms. 
robust_NMFs <- reduce_NMF_programs(robust_NMFs)

#Create pairwise data frame to get gene overlapping between programs.
#No need to compute this separately, the function compute metaprogram does it but it's not shown in the environment
threshold_overlapping_event <- 0.4 
pairwise_df <- compute_pairwise_df(robust_NMFs,threshold_overlapping_event)  



##############################################################
######### COMPUTING METAPROGRAMS #############################
##############################################################

threshold_overlapping_event <- 0.4
threshold_mp_founder <- 3
max_jac <- 0.4
metaprograms <- compute_metaprograms(robust_NMFs,threshold_overlapping_event,threshold_mp_founder,max_jac)


### Figure 2 Rapport.###
compute_heatmap(metaprograms)




##############################################################
######### BIOLOGICAL ANNOTATION OF MPs #######################
##############################################################

#Loading MP annotation database:  
ga_database <-  msigdbr(species = "Homo sapiens")
categories <- c("H","C5","C6")
subcategories_C5 <- c("HPO")
ga_database <- ga_database %>% filter(gs_cat %in% categories & !gs_subcat %in% subcategories_C5) %>% select(gs_name, gene_symbol)
pathways_df <- compute_functional_enrichment(metaprograms, ga_database)


####
counters_patients <- list()
for(i in 1:length(metaprograms)){
  patients <- strsplit(names(metaprograms)[i],"[.]")
  for(j in patients){
    patient <- substr(j,1,2)
    counters_patients[[paste0("MP",i)]][["summary"]] <- table(patient)
  }
}

counters_genes <- list()
for(i in 1:length(metaprograms)){
  patients <- unlist(strsplit(names(metaprograms)[i],"[.]"))
  for(j in names(metaprograms[[i]])){
    counter = 0
    for(k in patients){
      if(j %in% names(robust_NMFs[[k]])){
        counter = counter +1 
      }
      counters_genes[[paste0("MP",i)]][[j]] <- counter
    }
  }
}

compute_pvalue <- function(lin){
  pval_dat <- data.frame("MP" = c(as.integer(lin['non smoker']),as.integer(lin['smoker'])),"Total" = c(as.integer(lin['Total_smokers']),as.integer(lin['Total_nb_SS'])),row.names = c("Non Smoker", "Smoker"),stringsAsFactors = FALSE)
  p_value <- fisher.test(pval_dat)$p.value
  return(p_value)
}

mp_abundance <- compute_mp_abondance(counters_patients)


#List of NMF programs that are included in at least 1 MP
#Get list of unused nmf programs: 

patients <- character()
for(i in 1:length(metaprograms)){
  patients <- c(patients,unlist(strsplit(names(metaprograms)[i],"[.]")))
}

##############################################################
######### RESIDUAL NMFS ######################################
##############################################################

unused_nmf_programs <- robust_NMFs
for(i in patients){
  unused_nmf_programs[[i]] <- NULL
}

counter <- numeric()
for(i in names(unused_nmf_programs)){
  patient <- substr(i,1,2)
  counter <- c(patient, counter)
}

threshold_overlapping_event <- 0.3
threshold_mp_founder <- 1
max_jac <- 0.3
metaprograms_2 <- compute_metaprograms(unused_nmf_programs,threshold_overlapping_event,threshold_mp_founder,max_jac)
compute_heatmap(metaprograms_2)
pathways_df_2 <- compute_functional_enrichment(metaprograms_2, ga_database)


#Heatmap construction: 

# #Get for each patient, a list of its nmf programs included in at least one MP 
# summary_patients <- list()
# for(patient in patients){
#   patient_number <- substr(patient,1,2)
#   summary_patients[[patient_number]][[patient]] <- robust_NMFs[[patient]]
# }
# 
# for(i in names(summary_patients)){  #browse through patients 
#   gene_list <- character() # initialize gene list for each patient
#   for(j in names(summary_patients[[i]])){  
#     gene_list <- c(gene_list, names(summary_patients[[i]][[j]]))
#   }
#   table_genes <- as.data.frame(table(gene_list))
#   table_genes <- subset(table_genes,Freq == 1)
#   gene_list <- table_genes$gene_list
# }


#### Comparison between two MP: function jaccard needs to be computed differently #### 

jaccard2 <- function(a, b) {
  intersection = length(intersect(a,b)) # compute the number of shared genes between a and b
  union = length(a) + length(b) - intersection #compute genes from a and b minus the intersection
  return (intersection/union)  #jaccard index defined by intersection/union
}

MP_genes2 <- read_excel("//clb.loc/Recherche/NAS/EQUIPE SAINTIGNY/Sonia CANJURA/excels/MP_genes2.xlsx")

#Compare all 6 MP with Puram and Gavish datasets :
list_index <- list()
for(i in 1: length(metaprograms)){
  jaccard_index <- numeric()
  for(j in 1:length(MP_genes2)){
    jaccard_index <- c(jaccard_index, jaccard2(names(metaprograms[[i]]), MP_genes2[[j]]))
  }
  list_index[[i]] <- jaccard_index
}

df_comparison <- as.data.frame(list_index, col.names= c(names(counters_patients)),row.names = colnames(MP_genes2)) #View comparison table. 


#############################################
############### FIGURES #####################
#############################################

######## Patients distribution across MP. Figure 3A Rapport.  ############
my.cols <- piratepal(palette = "basel")
names(my.cols) <- c("CLB-HN22","CLB-HN24","CLB-HN25","CLB-HN29","CLB-HN31","CLB-HN33","CLB-HN36","CLB-HN38","CLB-HN39","CLB-HN42")


counters_patients <- list()
for(i in 1:length(metaprograms)){
  patients <- strsplit(names(metaprograms)[i],"[.]")
  for(j in patients){
    patient <- substr(j,1,2)
    counters_patients[[paste0("MP",i)]][["summary"]] <- table(patient)
  }
}

patients <- c("22","24","25","29","31","33","36","38","39","42")
df <- data.frame(matrix(nrow = length(NMFs)))
colnames(df) <- c("init")
rownames(df) <- c("22","24","29","36","38","25","31","33","39","42")
for(i in 1:length(counters_patients)){
  df[paste0("MP", i)] <- ifelse(rownames(df) %in% names(counters_patients[[i]][["summary"]]),1,0)
}
df <- subset(df, select=-c(init))
rownames(df) <- c("CLB-HN22","CLB-HN24","CLB-HN29","CLB-HN36","CLB-HN38","CLB-HN25","CLB-HN31","CLB-HN33","CLB-HN39","CLB-HN42")
colnames(df) <- c("EpiDiff1","EpiDiff2","Stress","Cell cycle", "pEMT","EpiSen")
df <- reshape2::melt(df %>% rownames_to_column("patient"))
df <- df %>% mutate(bin_value = case_when(value == '1' ~ 'present', value == '0' ~ 'absent'))
df$value <- as.numeric(df$value)
df$order_patient <- factor(df$patient,levels = c("CLB-HN22","CLB-HN24","CLB-HN29","CLB-HN36","CLB-HN38","CLB-HN25","CLB-HN31","CLB-HN33","CLB-HN39","CLB-HN42"))
ggplot(df, aes(x=variable, y=order_patient))+ geom_tile(aes(fill=bin_value), color="gray", size=0.25) + scale_fill_manual(values = c("white", "steelblue"))+ guides(fill = guide_legend(title = NULL)) + xlab("Metaprogram (MP)") + ylab(NULL)



###### Smokers/Non smokers distribution across MP. Figure 3B Rapport #####

my.cols <- piratepal(palette = "basel")
names(my.cols) <- c("CLB-HN22","CLB-HN24","CLB-HN25","CLB-HN29","CLB-HN31","CLB-HN33","CLB-HN36","CLB-HN38","CLB-HN39","CLB-HN42")


counters_patients <- list()
for(i in 1:length(metaprograms)){
  patients <- strsplit(names(metaprograms)[i],"[.]")
  for(j in patients){
    patient <- substr(j,1,2)
    counters_patients[[paste0("MP",i)]][["summary"]] <- table(patient)
  }
}

counter_mp <- 1
l_mp <- lapply(counters_patients,function(x){y <- data.frame(unlist(x$summary)); y['MP'] <- rep(counter_mp,nrow(y)) ; counter_mp <<- counter_mp +1 ; return(y)})
df_mp <- do.call(rbind.data.frame, l_mp) %>% mutate(SS = case_when(patient %in% smokers ~ "smoker" , !patient %in% smokers ~ "non smoker")) %>% arrange(MP,SS)
rownames(df_mp) <- NULL
df_mp <- df_mp[c("MP", "patient", "SS", "Freq")]
final_df <- df_mp %>% select(-c("patient")) %>% group_by(MP,SS) %>% summarise(NMFs= sum(Freq))
ggplot(final_df, aes(fill=SS, y=NMFs, x=MP)) + geom_bar(position='stack', stat='identity') + theme_minimal() + scale_fill_manual(values = c("darksalmon", "darkslateblue"))+
  scale_x_discrete(limits = c("MP1",  "MP2", "MP3", "MP4", "MP5", "MP6"),labels = c("EpiDiff1","EpiDiff2","Stress","Cell cycle", "pEMT","EpiSen")) + scale_y_continuous(breaks=seq(11))+
  guides(fill = guide_legend(title = NULL))+ xlab("Metaprogram (MP)")+ ylab("Number of NMF Programs")


##### Cell cycle correlation.  Figure S1 Rapport######

pvalue_ttest <- function(lin){
  lin <- unname(unlist(lin))
  return(-log10(t.test(lin,mu=0)$p.value))
}


list_genes <- lapply(metaprograms,function(x){y <- names(x); return(y)})
names(list_genes) <- c("MP1","MP2","MP3","MP4","MP5","MP6")


corr_list <- list()
list_sce <- list(sce_22,sce_24,sce_25,sce_29,sce_31,sce_33,sce_36,sce_38,sce_39,sce_42)
names(list_sce) <- c("sce_22","sce_24","sce_25","sce_29","sce_31","sce_33","sce_36","sce_38","sce_39","sce_42")
list_counts <- list()

for(i in names(list_sce)){
  log_counts <- SummarizedExperiment::assay(list_sce[[i]],"logcounts")
  colnames(log_counts) <- make.unique(colnames(log_counts))
  log_counts <- na.omit(log_counts)
  
  listdf <- list()
  
  for(j in 1:length(list_genes)){
    temp_counts <- log_counts[rownames(log_counts) %in% unlist(list_genes[[j]]),]
    temp_counts <- as.data.frame(apply(temp_counts,2,mean))
    colnames(temp_counts) <- c(paste0("MP", j))
    listdf[paste0("MP", j)] <- temp_counts
  }
  
  cell <- listdf[[4]]
  listdf[[4]] <- NULL
  
  for(k in 1:length(listdf)){
    corr_list[[str_extract(i,"[:digit:].")]][[k]] <- cor(cell,listdf[[k]])
    #corr_list[[str_extract(deparse(substitute(i)),"[:digit:].")]][[k]] <- cor(cell,listdf[[k]])
  }
}

corr_df <- as.data.frame(do.call("cbind",corr_list)) 
colnames(corr_df) <- c("p22" ,  "p24" ,  "p25"  ,"p29", "p31" ,  "p33"  , "p36" ,  "p38" ,  "p39"  , "p42")
corr_df <- corr_df %>% mutate(pval = apply(corr_df,1,pvalue_ttest))
corr_df <- corr_df %>% rowwise() %>% mutate(average = mean(c(p22,p24,p25,p31,p33,p36,p38,p39,p42)))
rownames(corr_df) <- c("EpiDiff1", "EpiDiff2", "Stress", "pEMT", "EpiSen")


my.cols <- piratepal(palette = "basel")
names(my.cols) <- c("CLB-HN22","CLB-HN24","CLB-HN25","CLB-HN29","CLB-HN31","CLB-HN33","CLB-HN36","CLB-HN38","CLB-HN39","CLB-HN42")


plot_df <- corr_df %>% select(pval,average)
ggplot(plot_df, aes(x=average, y=pval)) + theme()+ geom_point() + geom_text_repel(label=rownames(corr_df)) + xlab("Mean correlation with cell cycle MP") + ylab("Significance (-log10(pvalue))") + geom_hline(yintercept= -log10(0.05), color= "blue", linetype="dashed")


##### Density and violin plot. Figure 4 and S2 rapport##### 
#There is one figure for each MP. So in log_counts <- log_counts[rownames(log_counts) %in% unlist(list_genes[6]),], change the list_genes to each MP ####

list_genes <- lapply(metaprograms,function(x){y <- names(x); return(y)})
names(list_genes) <- c("MP1","MP2","MP3","MP4","MP5","MP6")

list_sce <- list(sce_22,sce_24,sce_25,sce_29,sce_31,sce_33,sce_36,sce_38,sce_39,sce_42)
names(list_sce) <- c("sce_22","sce_24","sce_25","sce_29","sce_31","sce_33","sce_36","sce_38","sce_39","sce_42")

df_all <- list()
for(i in names(list_sce)){
  initial_log_counts <- SummarizedExperiment::assay(list_sce[[i]],"logcounts") 
  log_counts <- initial_log_counts
  log_counts <- log_counts[rownames(log_counts) %in% unlist(list_genes[6]),]
  log_counts <- apply(log_counts,2,mean)

  df_patient <- as.data.frame(log_counts)
  colnames(df_patient) <- c("values")
  df_patient <- df_patient %>% mutate(patient = paste0("CLB-HN",str_extract(i,"[:digit:].")))
  
  df_all <- rbind(df_all,df_patient)
}
rownames(df_all) <- NULL

my.cols <- piratepal(palette = "basel")
names(my.cols) <- c("CLB-HN22","CLB-HN24","CLB-HN25","CLB-HN29","CLB-HN31","CLB-HN33","CLB-HN36","CLB-HN38","CLB-HN39","CLB-HN42")


order_ss <- c("CLB-HN22","CLB-HN24","CLB-HN29","CLB-HN36","CLB-HN38","CLB-HN25","CLB-HN31","CLB-HN33","CLB-HN39","CLB-HN42")
order_ss <- c("CLB-HN42","CLB-HN39","CLB-HN33","CLB-HN31","CLB-HN25","CLB-HN38","CLB-HN36","CLB-HN29","CLB-HN24","CLB-HN22")
df_all$order_patient <- factor(df_all$patient, levels = order_ss)

#Violin plot
plt <- ggbetweenstats(
  data = df_all,
  plot.type = "violin",
  x = order_patient,
  y = values,
  pairwise.comparisons = FALSE,
  bf.message = FALSE, 
  centrality.plotting = F,
  results.subtitle = FALSE,
  var.equal = FALSE,
  package = "yarrr",
  palette = "basel", 
  xlab = "Patient",
  ylab = "MP expression level"
)

plt + scale_colour_manual(values = my.cols, breaks = order_ss) + theme_minimal() +theme(axis.text = element_text(size = 12))

#Density plot

ggplot(df_all, aes(x= values, y= order_patient, fill = order_patient)) + theme_minimal()+ geom_density_ridges() + labs(x= " MP expression level") + scale_fill_manual(values = my.cols, breaks = order_ss) + ylab(NULL)
 
