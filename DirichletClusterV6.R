library(seqinr)

## arguments
args<- commandArgs(TRUE)
#args<-c("5_test.txt",
#        "Paired_5_test.txt",
#        "../Result/5_test",
#        "1",
#        "Overlap_5_test.txt")


## functions
sorted_overlap<- function(overlap_info) {
  seq_cluster_index<- lapply(lapply(strsplit(readLines(overlap_info), "\t"), function(x){x<-x[x!=""]}), as.integer)
  l<- unlist(lapply(seq_cluster_index,length))
  sorted_list<- seq_cluster_index[order(l, decreasing = T)]
  sorted_list_length<- lapply(sorted_list,length)
  return(list("sorted_overlap_index" = unlist(sorted_list), 
              "sorted_overlap_length" = unlist(sorted_list_length)))
}

index_map<- function(sorted_index) {
  res<- c()
  for(l in 1:length(sorted_index)) {
    res<- c(res, 2 * sorted_index[l] - 1)
    res<- c(res, 2 * sorted_index[l])
  }
  return(res)
}

args_to_params<- function(args) {
  # order of the markov chain
  m<- as.numeric(args[4])
  
  output_keyword<- args[3]
  
  seqs_raw<- read.table(args[1], colClasses = "character")
  seq_length<- nchar(seqs_raw[1, 1])
  transition_counts_raw<- read.table(args[2])
  
  # reorder raw sequences
  sorted_index_info<- sorted_overlap(args[5])
  sorted_overlap_index<- sorted_index_info$sorted_overlap_index + 1
  sorted_overlap_length<- sorted_index_info$sorted_overlap_length
  seqs_index<- index_map(sorted_overlap_index)
  seqs<- as.matrix(seqs_raw[seqs_index, ], colClasses = "character")
  seqs_details<- apply(seqs, 2, strsplit, split = "")[[1]]
  transition_counts<- transition_counts_raw[sorted_overlap_index, ]
  
  return(list("output_keyword" = output_keyword,
              "m" = m,
              "transition_counts" = transition_counts,
              "sorted_overlap_index" = sorted_overlap_index,
              "sorted_overlap_length" = sorted_overlap_length,
              "seqs" = seqs,
              "seqs_details" = seqs_details))
}

calc_gc_content<- function(transition_counts, m) {
  return(sum(transition_counts[(4^m+1):(3*4^m)])/sum(transition_counts))
}

counts_to_prob<- function(counts) {
  if(sum(counts) == 0) {
    return(0)
  } else {
    return(counts/sum(counts))
  }
}

# calculate starting point probability using counts
start_p_counts<- function(counts, m) {
  start_counts<- unlist(tapply(counts, rep(1:(2*2^m), each = 4), sum))
  startp<- start_counts/sum(start_counts)
  return(startp)
}

# calculate transition probability using counts
trans_p_counts<- function(counts, m) {
  transp<- unlist(tapply(counts, rep(1:(2*2^m), each = 4), counts_to_prob))
  return(transp)
}
  
# estimate new group starting point and transition probablity using kNN method on gc content
new_group<- function(transition_counts, m, n_seqs, neighbor) {
  gc<- apply(transition_counts, 1, calc_gc_content, m = m)
  new_group_transp<- matrix(0, n_seqs, 4^(m+1))
  new_group_startp<- matrix(0, n_seqs, 4)
  
  gc_order<- order(gc)
  
  for(i in 1:length(gc_order)) {
    # evenly split neighbors
    subgroup_index<- seq(ifelse(i-neighbor/2 > 0, i-neighbor/2, 1),
                         ifelse(i+neighbor/2< length(gc_order), i+neighbor/2, length(gc_order)))
    subgroup_index<- gc_order[subgroup_index]
    subgroup<- apply(transition_counts[subgroup_index, ], 2, sum)
    new_group_startp[gc_order[i], ]<- start_p_counts(subgroup, m)
    new_group_transp[gc_order[i], ]<- trans_p_counts(subgroup, m)
  }
  return(list("startp" = new_group_startp, "transp" = new_group_transp))
}

prior_p_temp<- function(seq_index, params, z, alpha) {
  p<- vector(mode = "list", length = params$n_particles)
  group<- max(z)
  for(k in 1:params$n_particles) {
    group_col<- unique(z[z[, k] != 0, k])
    p[[k]]<- matrix(alpha/(alpha + seq_index - 1)/(group + 1 - length(group_col)), group + 1)
    for(j in group_col) {
      p[[k]][j]<- sum(z[, k] == j)/(seq_index - 1 + alpha)
    }
  }
  return(p)
}

prior_p<- function(seq_index, params, w, z, alpha) {
  group<- max(z)
  p<- matrix(0, group + 1, params$n_particles)
  temp_prior<- prior_p_temp(seq_index, params, z, alpha)
  for(k in 1:params$n_particles) {
    p[, k]<- w[k] * temp_prior[[k]]
  }
  p<- rowSums(p)
  return(p)
}

trans_p<- function(seq_index, params, new_group_params, z) {
  p<- vector(mode = "list", length = params$n_particles)
  group<- max(z)
  new_trans_p<- new_group_params$transp[seq_index, ]
  for(k in 1:params$n_particles) {
    p[[k]]<- array(0, dim = c(group + 1, 4^(params$m + 1)))
    for(j in 1:group) {
      cluster_index<- which(z[, k] == j)
      if(length(cluster_index) > 1) {
        p[[k]][j, ]<- trans_p_counts(colSums(params$transition_counts[cluster_index, ]), params$m)
      } else if(length(cluster_index) == 1) {
        p[[k]][j, ]<- trans_p_counts(as.numeric(params$transition_counts[cluster_index, ]), params$m)
      } 
      #else if(j == (group + 1)) {
      #  p[[k]][j, ]<- new_trans_p
      #} 
      else {
        p[[k]][j, ]<- 0.0 * new_trans_p
      }
    }
    p[[k]][(group + 1), ]<- new_trans_p
  }
  # res<- lapply(p, function(x) {ifelse(x == 0.0, min(x[x>0])/1000, x)})
  # return(res)
  return(p)
}

startwhere_original<-function(seq_index, seqs_details){
  
  if (seqs_details[[seq_index]][1]=="A"){
    where = 1
  } else if (seqs_details[[seq_index]][1]=="C"){
    where = 2
  }else if (seqs_details[[seq_index]][1]=="G"){
    where = 3
  }else {
    where = 4
  }
  return(where)  
}

startwhere_compreverse<-function(seq_index, seqs_details){
  
  if (seqs_details[[seq_index]][100]=="A"){
    where = 4
  } else if (seqs_details[[seq_index]][100]=="C"){
    where = 3
  }else if (seqs_details[[seq_index]][100]=="G"){
    where = 2
  }else {
    where = 1
  }
  return(where)  
}

start_p<- function(seq_index, params, new_group_params, z) {
  p<- vector(mode = "list", length = params$n_particles)
  group<- max(z)
  new_start_p<- new_group_params$startp[seq_index, ]
  for(k in 1:params$n_particles) {
    p[[k]]<- array(0, dim = c(group + 1, 4))
    for(j in 1:group) {
      cluster_index<- which(z[, k] == j)
      if(length(cluster_index) > 1) {
        p[[k]][j, ]<- start_p_counts(colSums(params$transition_counts[cluster_index, ]), params$m)
      } else if(length(cluster_index) == 1) {
        p[[k]][j, ]<- start_p_counts(as.numeric(params$transition_counts[cluster_index, ]), params$m)
      }
      #else if(j == (group + 1)) {
      #  p[[k]][j, ]<- new_start_p
      #} 
      else {
        p[[k]][j, ]<- 0.0 * new_start_p
      }
    }
    p[[k]][(group + 1), ]<- new_start_p 
  }
  # res<- lapply(p, function(x) {ifelse(x == 0.0, min(x[x>0])/1000, x)})
  # return(res)
  return(p)
}

post_p_temp<- function(seq_index, params, new_group_params, w, z) {
  group<- max(z)
  p<- vector(mode = "list", length = params$n_particles)
  
  temp_trans_p<- trans_p(seq_index, params, new_group_params, z)
  
  start_o_1<- startwhere_original(2*seq_index - 1, params$seqs_details)
  start_r_1<- startwhere_compreverse(2*seq_index - 1, params$seqs_details)
  start_o_2<- startwhere_original(2*seq_index, params$seqs_details)
  start_r_2<- startwhere_compreverse(2*seq_index, params$seqs_details)
  temp_start_p<- start_p(seq_index, params, new_group_params, z)
  
  for(k in 1:params$n_particles) {
    p[[k]]<- matrix(0, group + 1)
    for(j in 1:(group + 1)) {
      if(sum(temp_start_p[[k]][j, ]) > 0 & sum(temp_trans_p[[k]][j, ]) > 0){
        p[[k]][j]<- w[k]*exp(
          log(temp_start_p[[k]][j, start_o_1]) + 
          log(temp_start_p[[k]][j, start_r_1]) + 
          log(temp_start_p[[k]][j, start_o_2]) + 
          log(temp_start_p[[k]][j, start_r_2])
        ) * exp(sum(params$transition_counts[seq_index, ] * # log(temp_trans_p[[k]][j, ])
                    log(if(0 %in% temp_trans_p[[k]][j, ]) unlist(lapply(temp_trans_p[[k]][j, ], function(x) {ifelse(x == 0, 0.00000001, x)})) else temp_trans_p[[k]][j, ])
                    )
                )
      } else {
        p[[k]][j]<- 0
      }
    }
  }
  p<- Reduce("+", p)
  return(p)
}

post_p<- function(seq_index, params, new_group_params, w, z, alpha) {
  group<- max(z)
  p<- rep(0, group + 1)
  
  temp_post_p<- post_p_temp(seq_index, params, new_group_params, w, z)
  temp_prior_p<- prior_p(seq_index, params, w, z, alpha)
  
  for(j in 1:(group + 1)) {
    p[j]<- temp_prior_p[j] * temp_post_p[j, ]
  }
  return(p/sum(p))
}

cluster_one_seq<- function(seq_index, params, new_group_params, w, z, alpha) {
  posterior<- post_p(seq_index, params, new_group_params, w, z, alpha)
  print(posterior)
  return(sample.int(length(posterior), params$n_particles, replace = T, prob = posterior))
}

cluster_overlap_seqs<- function(seqs_list, params, new_group_params, w, z, alpha) {
  temp_z<- matrix(0, ncol = params$n_particles, nrow = length(seqs_list))
  counter<- 1
  for(seq_index in seqs_list) {
    temp_z[counter, ]<- cluster_one_seq(seq_index, params, new_group_params, w, z, alpha)
    counter<- counter + 1
  }
  return(temp_z)
}

get_mode<- function(x) {
  ux<- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

check_majority_vote<- function(z, p_majority) {
  vote<- get_mode(apply(z, 1, get_mode))
  if(length(which(apply(z, 1, get_mode) == vote)) < p_majority * length(apply(z, 1, get_mode))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

update_weights<- function(seq_index, params, new_group_params, w, z) {
  p<- rep(0, params$n_particles)
  temp_trans_p<- trans_p(seq_index, params, new_group_params, z)
  temp_start_p<- start_p(seq_index, params, new_group_params, z)
  start_o_1<- startwhere_original(2*seq_index - 1, params$seqs_details)
  start_r_1<- startwhere_compreverse(2*seq_index - 1, params$seqs_details)
  start_o_2<- startwhere_original(2*seq_index, params$seqs_details)
  start_r_2<- startwhere_compreverse(2*seq_index, params$seqs_details)
  for(k in 1:params$n_particles) {
    p[k]<- w[k]*exp(
      log(temp_start_p[[k]][z[seq_index, k], start_o_1]) + 
      log(temp_start_p[[k]][z[seq_index, k], start_r_1]) + 
      log(temp_start_p[[k]][z[seq_index, k], start_o_2]) + 
      log(temp_start_p[[k]][z[seq_index, k], start_r_2])
    ) * exp(sum(params$transition_counts[seq_index, ] * # log(temp_trans_p[[k]][j, ])
                  log(if(0 %in% temp_trans_p[[k]][z[seq_index, k], ]) unlist(lapply(temp_trans_p[[k]][z[seq_index, k], ], function(x) {ifelse(x == 0, 0.00000001, x)})) else temp_trans_p[[k]][z[seq_index, k], ])
                )
            )
  }
  res<- p/sum(p)
  return(res)
}

resample_particles<- function(seq_index, params, w, z) {
  eff_n<- 1/sum(w^2)
  # print(paste(eff_n, params$threshold, params$n_particles))
  if(eff_n < params$threshold * params$n_particles) {
    resample_index<- sample.int(params$n_particles, replace = T, prob = w)
    w<- rep(1/params$n_particles, params$n_particles)
    return(list("z" = z[seq_index, resample_index], "w" = w))
  } else {
    return(list("z" = z[seq_index, ], "w" = w))
  }
}

# implement dirichlet clustering to a single group of sequences
dirichlet_cluster_single<- function(params, alpha, z_mode_lower) {
  # initialize hidden variable
  z<- matrix(0, ncol = params$n_particles, nrow = params$n_seqs)
  w<- rep(1/params$n_particles, params$n_particles)
  
  actual_index<- cumsum(params$sorted_overlap_length)
  z[(1:params$sorted_overlap_length[1]), ]<- 1
  new_group_params<- new_group(params$transition_counts,
                               params$m,
                               params$n_seqs,
                               params$neighbor)
  # container for removed index
  rm_index<- list()
  
  for(i in 2:length(actual_index)) {
    seqs_list<- (actual_index[i-1]+1):actual_index[i]
    # if overlapped seqs
    if(length(seqs_list) > 1) {
      # hack: check two alpha values to see if results are stable
      # alternative: set alpha/1000, alpha*1000 to two contant values
      temp_z_lower<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, alpha/100)
      temp_z_upper<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, min(alpha*(10^6),1))
      # temp_z_lower<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, 10^(-6))
      # temp_z_upper<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, 1)
      is_marjority_unified_lower<- check_majority_vote(temp_z_lower, params$p_marjority)
      is_marjority_unified_upper<- check_majority_vote(temp_z_lower, params$p_marjority)
      majority_vote_lower<- get_mode(apply(temp_z_lower, 1, get_mode))
      majority_vote_upper<- get_mode(apply(temp_z_upper, 1, get_mode))
      if(is_marjority_unified_lower && is_marjority_unified_upper && (majority_vote_upper == majority_vote_lower)) {
#         temp_z<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, alpha)
#         majority_vote<- get_mode(apply(temp_z, 1, get_mode))
        for(seq_index in seqs_list) {
          z[seq_index, ]<- majority_vote_lower
        }
       # write.table(z, paste(params$output_keyword, "_z.txt", sep=""))
      } else {
        rm_index[[length(rm_index) + 1]]<- seqs_list
        file.remove(paste(params$output_keyword, "_rm_index.txt", sep=""))
        lapply(rm_index, function(x) {write(x, append = T, ncolumns = length(x),
                                            file = paste(params$output_keyword, "_rm_index.txt", sep=""))})
      }
    } else {
      # if a single seq
      z[seqs_list, ]<- cluster_one_seq(seqs_list, params, new_group_params, w, z, alpha)
      w<- update_weights(seqs_list, params, new_group_params, w, z)
      resample_status<- resample_particles(seqs_list, params, w, z)
      z[seqs_list, ]<- resample_status$z
      w<- resample_status$w
      #print(w)
    }
    write.table(z, paste(params$output_keyword, "_z.txt", sep=""))
    print(i)
    print(z[seqs_list,])
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  }
  
  # process previously removed seqs
  if(length(rm_index) > 0) {
    for(i in 1:length(rm_index)) {
      seqs_list<- rm_index[[i]]
      temp_z<- cluster_overlap_seqs(seqs_list, params, new_group_params, w, z, alpha)
      majority_vote<- get_mode(apply(temp_z, 1, get_mode))
      for(seq_index in seqs_list) {
        z[seq_index, ]<- majority_vote
      }
      write.table(z, paste(params$output_keyword, "_z.txt", sep=""))
    }
    
    
  }
  print(apply(z, 1, get_mode) + z_mode_lower)
  return(apply(z, 1, get_mode) + z_mode_lower)
}

reverse_index_map<- function(sorted_overlap_length, index) {
  return(unique(unlist(lapply(index, function(x) {return(min(which(x <= cumsum(sorted_overlap_length))))}))))
}

# iteratively implement dirichlet clustering to sub-groups of sequences
dirichlet_cluster<- function(params, alpha, z_mode) {
  updated_z_mode<- z_mode
  temp_z_mode_lower<- 0
  for(c in unique(z_mode)) {
    temp_index<- which(z_mode == c)
    temp_params<- params
    temp_params$transition_counts<- temp_params$transition_counts[temp_index, ]
    temp_params$sorted_overlap_length<- temp_params$sorted_overlap_length[reverse_index_map(params$sorted_overlap_length,
                                                                                            temp_index)]
    temp_params$seqs<- temp_params$seqs[index_map(temp_index), ]
    temp_params$seqs_details<- temp_params$seqs_details[index_map(temp_index)]
    temp_params$n_seqs<- dim(temp_params$transition_counts)[1]
    
    temp_z_mode<- dirichlet_cluster_single(temp_params, alpha, temp_z_mode_lower)
    updated_z_mode[temp_index]<- temp_z_mode 
    temp_z_mode_lower<- max(temp_z_mode)
  }
  return(updated_z_mode)
}

## main function
main<- function(args, alphas) {
	params<- args_to_params(args)
  
  # set additional parameters
  params$n_particles<- 100
  params$threshold<- 0.9 # used in resampling
  params$neighbor<- 1000
  params$p_marjority<- 0.5
  params$n_seqs<- dim(params$transition_counts)[1]

  # initial clustering
  z_mode<- matrix(1, ncol = 1, nrow = params$n_seqs)
  old_z_mode_max<- max(z_mode)
  
  for(alpha in alphas) {
    z_mode<- dirichlet_cluster(params, alpha, z_mode)
    write.table(z_mode, paste(params$output_keyword, alpha, "z_mode.txt", sep="_"))
    if(max(z_mode) == old_z_mode_max) {
      break
    }
    old_z_mode_max<- max(z_mode)
  }
  write.table(z_mode, paste(params$output_keyword, "_z_mode.txt", sep=""))
}

# run algorithm
main(args, alphas = c(10^(-8), 10^(-5)))
