# POSTERIOR PREDICTIVE CHECKS


# generate simulated data (yrep) using parameters estimated by the model
# The output, yrep_list, is the set of data that would have been observed if 
# the model (with the sampled parameters) were the true data generating machine
simulate_yrep <- function(pi_samples, phi_samples, N, P, V, S = 1000) {
  yrep_list <- vector("list", S)
  zrep_list <- vector("list", S)
  for (s in 1:S) {
    pi_s  <- pi_samples[[s]]
    phi_s <- phi_samples[[s]]
    yrep <- matrix(NA, nrow = N, ncol = P)
    z_s <- integer(N)
    for (i in 1:N) {
      z_i <- sample(1:length(pi_s), size = 1, prob = pi_s)
      z_s[i] <- z_i
      for (d in 1:P) {
        yrep[i, d] <- sample(1:V[d], size = 1, prob = phi_s[[z_i]][[d]])
      }
    }
    yrep_list[[s]] <- yrep
    zrep_list[[s]] <- z_s
  }
  return(list(yrep_list = yrep_list, zrep_list = zrep_list))
}
yrep_list <- simulate_yrep(pi_samples, phi_samples, N = nrow(X), 
                           P = ncol(X), V = V, S = 1460)





barplot_ppc_cluster <- function(yrep_list, zrep_list, X, z_samples, V, var_index, cluster_index) {
  S <- length(yrep_list)
  freq_mat <- matrix(NA, nrow = S, ncol = V[var_index])
  for (s in 1:S) {
    z_rep <- zrep_list[[s]]
    y_rep <- yrep_list[[s]]
    
    idx_cluster <- which(z_rep == cluster_index)
    
    if (length(idx_cluster) > 0) {
      tab <- table(factor(y_rep[idx_cluster, var_index], levels = 1:V[var_index]))
      freq_mat[s, ] <- tab / sum(tab)
    }
  }
  z_last <- z_samples[[length(z_samples)]]
  idx_obs <- which(z_last == cluster_index)
  tab_obs <- table(factor(X[idx_obs, var_index], levels = 1:V[var_index]))
  freq_obs <- tab_obs / sum(tab_obs)
  freq_mean <- colMeans(freq_mat, na.rm = TRUE)
  freq_sd   <- apply(freq_mat, 2, sd, na.rm = TRUE)
  bar_centers <- barplot(freq_mean,
                         ylim = c(0, 1),
                         beside = TRUE,
                         col = "#BDD7E7",
                         names.arg = paste("Lev", 1:V[var_index]),
                         main = paste("Var", var_index, "Cluster", cluster_index),
                         ylab = "frequency")
  arrows(bar_centers, freq_mean - freq_sd,
         bar_centers, freq_mean + freq_sd,
         angle = 90, code = 3, length = 0.05)
  points(bar_centers, freq_obs, col = "#a50f15", pch = 19)
}
sim_results <- simulate_yrep(pi_samples, phi_samples, N = nrow(X), P = ncol(X), V = V, S = 1000)
yrep_list <- sim_results$yrep_list
zrep_list <- sim_results$zrep_list



par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 1, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 1, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 1, cluster_index = 3)
par(mfrow = c(1, 1))


par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 2, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 2, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 2, cluster_index = 3)
par(mfrow = c(1, 1))


par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 3, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 3, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 3, cluster_index = 3)
par(mfrow = c(1, 1))


par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 4, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 4, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 4, cluster_index = 3)
par(mfrow = c(1, 1))


par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 5, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 5, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 5, cluster_index = 3)
par(mfrow = c(1, 1))ù


par(mfrow = c(1, 3))
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 6, cluster_index = 1)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 6, cluster_index = 2)
barplot_ppc_cluster(yrep_list = yrep_list, zrep_list = zrep_list, X = X, z_samples = z_samples, 
                    V = V, var_index = 6, cluster_index = 3)
par(mfrow = c(1, 1))



# more precise implementatio of plots:
barplot_ppc_cluster <- function(yrep_list, zrep_list, X, z_samples, V, var_index, cluster_index) {
  S <- length(yrep_list)
  freq_mat <- matrix(NA, nrow = S, ncol = V[var_index])
  
  # Frequenze simulate
  for (s in 1:S) {
    z_rep <- zrep_list[[s]]
    y_rep <- yrep_list[[s]]
    
    idx_cluster <- which(z_rep == cluster_index)
    
    if (length(idx_cluster) > 0) {
      tab <- table(factor(y_rep[idx_cluster, var_index], levels = 1:V[var_index]))
      freq_mat[s, ] <- tab / sum(tab)
    }
  }
  
  # Frequenze osservate
  z_last <- z_samples[[length(z_samples)]]
  idx_obs <- which(z_last == cluster_index)
  tab_obs <- table(factor(X[idx_obs, var_index], levels = 1:V[var_index]))
  freq_obs <- tab_obs / sum(tab_obs)
  
  # Calcola media e intervallo
  freq_mean <- colMeans(freq_mat, na.rm = TRUE)
  freq_sd   <- apply(freq_mat, 2, sd, na.rm = TRUE)
  
  # Plot
  bar_centers <- barplot(freq_mean,
                         ylim = c(0, 1),
                         beside = TRUE,
                         col = "lightblue",
                         names.arg = paste("Lev", 1:V[var_index]),
                         main = paste("Var", var_index, "Cluster", cluster_index),
                         ylab = "frequency")
  
  # Aggiungi errore standard
  arrows(bar_centers, freq_mean - freq_sd,
         bar_centers, freq_mean + freq_sd,
         angle = 90, code = 3, length = 0.05)
  
  # Aggiungi osservato
  points(bar_centers, freq_obs, col = "red", pch = 19)
}
# Supponiamo tu abbia eseguito:
sim_results <- simulate_yrep_with_z(pi_samples, phi_samples, N = nrow(X), P = ncol(X), V = V, S = 1000)
yrep_list <- sim_results$yrep_list
zrep_list <- sim_results$zrep_list

# Ora puoi fare il plot:
barplot_ppc_cluster(
  yrep_list = yrep_list,
  zrep_list = zrep_list,
  X = X,
  z_samples = z_samples,
  V = V,
  var_index = 3,        # cambia con la variabile desiderata
  cluster_index = 2     # cambia con il cluster desiderato