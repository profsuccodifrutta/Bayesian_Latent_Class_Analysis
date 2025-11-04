library(coda)



# GIBBS SAMPLER PARAMETERS INITIALIZATION
X = df
K = 3
V = c(2, 2, 5, 5, 4, 3)
alpha = rep(1, K)
beta_list <- list(rep(1, V[1]), rep(1, V[2]), rep(1, V[3]), rep(1, V[4]),
                  rep(1, V[5]), rep(1, V[6]))



# GIBBS SAMPLER
gibbs_sampler_time <- function(X, K, V, alpha, beta_list, n_iter = 1000) {
  N <- nrow(X)  # numero di osservazioni
  P <- ncol(X)  # numero di variabili categoriche
  
  z <- sample(1:K, N, replace = TRUE)
  pi <- rep(1 / K, K)
  phi <- full_cond_phi(X, z, beta_list, K, V)
  
  samples <- list(pi = vector("list", n_iter),
                  phi = vector("list", n_iter),
                  z = vector("list", n_iter))
  
  start_time <- Sys.time()  # Inizio temporizzazione
  
  for (t in 1:n_iter) {
    pi <- full_cond_pi(z, alpha)
    phi <- full_cond_phi(X, z, beta_list, K, V)
    z <- full_cond_z(X, pi, phi, K)
    
    samples$pi[[t]] <- pi
    samples$phi[[t]] <- phi
    samples$z[[t]] <- z
    
    # Stampa ogni 100 iterazioni
    if (t %% 100 == 0) {
      elapsed <- Sys.time() - start_time
      cat(sprintf("Iterazione %d completata - Tempo trascorso: %s\n", t, elapsed))
    }
  }
  
  return(samples)
}


# thinning and burn in of the chain
thin <- 10
burn_in <- 400
pi_samples  <- samples$pi[seq(burn_in + 1, length(samples$pi), by = thin)]
phi_samples <- samples$phi[seq(burn_in + 1, length(samples$phi), by = thin)]
z_samples   <- samples$z[seq(burn_in + 1, length(samples$z), by = thin)]




# POSTERIOR ANALYSIS OF PI
pi_mean <- Reduce("+", pi_samples) / length(pi_samples) #posterior mean pi

# traceplots and acf
par(mfrow = c(3,2))
for (k in 1:K) {
  traceplot <- sapply(pi_samples, function(p) p[k])
  plot(traceplot, type = "l", main = paste("Traceplot pi[", k, "]"), ylab = "pi_k")
  acf(traceplot, main = paste("ACF pi[", k, "]"))
}

# convert pi sample into a matrices (n_iter x K)
pi_matrix <- do.call(rbind, pi_samples)

# Convert in mcmc 
pi_mcmc <- mcmc(pi_matrix)

# compute ESS for pi
effectiveSize(pi_mcmc)




# Credible intervals (95%) for each pi
ci_pi <- apply(pi_matrix, 2, quantile, probs = c(0.025, 0.975))

# pi histogram and credible intervals 
par(mfrow = c(1, K))

for (k in 1:K) {
  hist(pi_matrix[, k], breaks = 30, probability = TRUE,
       main = paste0("Posterior π[", k, "]"),
       xlab = expression(pi), col = "lightblue", border = "white")
  
  abline(v = pi_mean[k], col = "red", lwd = 2, lty = 2)  # Posterior mean
  abline(v = ci_pi[1, k], col = "darkgreen", lwd = 2, lty = 3)  # 2.5%
  abline(v = ci_pi[2, k], col = "darkgreen", lwd = 2, lty = 3)  # 97.5%
  
}
legend("topright", legend = c("Mean", "95% CI"),
       col = c("red", "darkgreen"), lty = c(2, 3), lwd = 2)




#  POSTERIOR ANALYSIS OF PHI

# mean of phi per each class k and variable d
phi_mean <- vector("list", K)
for (k in 1:K) {
  phi_mean[[k]] <- vector("list", length(V))  # for each variable
  for (d in 1:length(V)) {
    # sum over all iterations
    sum_phi_kd <- Reduce("+", lapply(phi_samples, function(phi) phi[[k]][[d]]))
    phi_mean[[k]][[d]] <- sum_phi_kd / length(phi_samples)
  }
}

# example: what is the frequency of the levels of variable 4 in the first cluster?
# these results are shown in the barplots
print(phi_mean[[1]][[4]])




par(mfrow = c(2,3))

# barplots shows the distribution of each level of each variable among all cluster
# to produce barplot for each variable modify the following: colors = number of levels of the 
# variable taken under consideration and change the title (could be done with a for loop)

# then run once for each variable

barplot_matrix <- do.call(rbind, lapply(1:K, function(k) {
  phi_mean[[k]][[5]]  # variable number
}))

barplot(t(barplot_matrix),
        beside = TRUE,
        col = c("lightblue","blue", "darkgreen", "purple"), # number of colors = number of levels 
        legend.text = TRUE,
        names.arg = paste("Cluster", 1:K),
        main = "Distribuzione Businnes")



# traceplots  and acf of phi

# intractable to monitor all traceplot and acf of phi, select just few variables, 
# (the once that discriminates better the clusters) to check the overall convergence
num_clusters <- length(1:3)     # number of cluster 
num_vars <- length(1)          # select 1 variable
num_levels <- V[2]             # number of levels of the variable selected
total_panels <- num_clusters * num_vars * num_levels

par(mfrow = c(total_panels, 2)) 

plot_trace_phi_with_acf <- function(phi_samples, clusters = c(1, 2,3), vars = c(2), V) {
  for (k in clusters) {
    for (d in vars) {
      for (v in 1:V[d]) {
        trace <- sapply(phi_samples, function(phi) phi[[k]][[d]][v])
        
        # Traceplot
        plot(trace, type = "l",
             main = paste0("Trace: φ[", k, "][", d, "][", v, "]"),
             ylab = "Value", xlab = "Iteration")
        
        # Autocorrelation plot
        acf(trace, main = paste0("ACF: φ[", k, "][", d, "][", v, "]"))
      }
    }
  }
}

#esempio : 2 grafici: traceplot e acf, ripetuti per tutti i livelli della variabile.
# variabile 2 ha 2 livelli, 2x2xcluster2 = 4 x 2 (traceplot e acf) = 8grafici

# example: traceplot and acf will be produced for each level of the variable, if the 
# variable has 2 levels: 2 x 2 x number of cluster = 12 graphs produced
plot_trace_phi_with_acf(phi_samples, clusters = 1:3, vars = 1, V = V) 




# POSTERIOR ANALYSIS OF Z
z_samples <- samples$z[burn_in:length(samples$z)] 

# traceplots z

# diagnostics for z are tricky: check the overall size of cluster assigment, and checks
# for possible problems of labels switching

par(mfrow = c(1,1))
cluster_freq <- sapply(z_samples, function(z) table(factor(z, levels = 1:K)))
matplot(t(cluster_freq), type = "l", lty = 1, col = 1:K,
        ylab = "Numero osservazioni", xlab = "Iterazione",
        main = "Dimensioni dei cluster nel tempo")
legend("topright", legend = paste("Cluster", 1:K), col = 1:K, lty = 1)
