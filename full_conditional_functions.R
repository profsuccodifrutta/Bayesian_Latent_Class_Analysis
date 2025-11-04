# defining the full conditional functions 
library(MCMCpack)
full_cond_pi <- function(z, alpha) {
  # z: vettore lungo N con valori in 1,...,K
  # alpha: vettore lungo K (parametri della Dirichlet prior)
  
  K <- length(alpha)
  N_k <- tabulate(z, nbins = K)  # conto le osservazioni per ciascuna classe
  return(rdirichlet(1, alpha + N_k))  # restituisce un vettore di dimensione K
}


#distribuzione dei livelli della variabile d per ogni classe k 

full_cond_phi <- function(X, z, beta_list, K, V) {
  # X: matrice N x P di osservazioni categoriali
  # z: vettore lungo N con classi latenti
  # beta_list: lista di P vettori di dimensione V_d (priori per ogni variabile)
  # K: numero di classi
  # V: vettore lungo P, con V[d] = num categorie della variabile d
  
  P <- ncol(X)
  phi <- vector("list", K)  # lista lunga K
  
  for (k in 1:K) {
    phi_k <- list()  # lista lunga P, Lista per le distribuzioni delle variabili nella classe k
    for (d in 1:P) {
      v_d <- V[d] #numero categorie per la variabile d
      x_kd <- X[z == k, d]  # valori della variabile d nella classe k
      counts <- tabulate(x_kd, nbins = v_d)  # vettore lungo V_d, frequenze osservate per ciascun livello
      beta_d <- beta_list[[d]]
      phi_k[[d]] <- rdirichlet(1, counts + beta_d)  # vettore lungo V_d
    }
    phi[[k]] <- phi_k
  }
  
  return(phi)  # lista K x P con vettori di dimensione V_d
}



full_cond_z <- function(X, pi, phi, K) {
  N <- nrow(X)
  P <- ncol(X)
  z_new <- integer(N)
  
  for (n in 1:N) {
    probs <- numeric(K)  # Vettore per salvare le probabilità P(z_n = k | ...)
    for (k in 1:K) {
      p_k <- pi[k] # Inizializza con la prior π_k
      for (d in 1:P) {
        v <- X[n, d] # Valore osservato per la variabile d nell'osservazione n
        p_k <- p_k * phi[[k]][[d]][v]  # Probabilità congiunta: pi_k * prodotto delle prob categorie
      }
      probs[k] <- p_k  # Salva la probabilità congiunta (non normalizzata)
    }
    probs <- probs / sum(probs) # Normalizzazione
    z_new[n] <- sample(1:K, size = 1, prob = probs)
  }
  
  return(z_new)
}


