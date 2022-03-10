# Model for repeated binary tests as described in Wang & Hanson (2019) in JAGS
# based on code and priors from Wang et al. (2020)
# including calculations for parallel testing

model_str <- "
model{
  for (i in 1:N) {
    for (k in 1:K) {
      d1[i,k] <- se[k] ^ sum(x[i,k,]) * (1 - se[k]) ^ (J - sum(x[i,k,]))
      d2[i,k] <- sp[k] ^ (J - sum(x[i,k,])) * (1 - sp[k]) ^ sum(x[i,k,])
    }
    dp[i] <- prod(d1[i,])
    dn[i] <- prod(d2[i,])
    
    for (j in 1:J) {
      for (u in 1:(K-1)) {
        for (v in (u+1):K) {
          t1[i,j,u,v] <- cop[u,v] * ((se[u] - 1) / se[u])^x[i,u,j] * ((se[v] - 1) / se[v])^x[i,v,j] / ((1 - se[u]) * (1 - se[v]))
          t2[i,j,u,v] <- con[u,v] * (sp[u] / (sp[u] - 1))^x[i,u,j] * (sp[v] / (sp[v] - 1))^x[i,v,j] / (sp[u] * sp[v])
        }
      }
      for (v in 1:K) {
        for (u in v:K) {
          t1[i,j,u,v] <- 0
          t2[i,j,u,v] <- 0
        }
      }
    }
    tcp_tilde[i] <- sum(t1[i,,,])
    tcn_tilde[i] <- sum(t2[i,,,])
    
    for (k in 1:K) {
      for (u in 1:(J-1)) {
        for (v in (u+1):J) {
          r1[i,k,u,v] <- ((se[k] - 1) / se[k])^(x[i,k,u] + x[i,k,v]) / (1 - se[k])^2
          r2[i,k,u,v] <- (sp[k] / (sp[k] - 1))^(x[i,k,u] + x[i,k,v]) / sp[k]^2
        }
      }
      for (v in 1:J) {
        for (u in v:J) {
          r1[i,k,u,v] <- 0
          r2[i,k,u,v] <- 0
        }
      }
      r3[i,k] <- rp[k] * sum(r1[i,k,,])
      r4[i,k] <- rn[k] * sum(r2[i,k,,])
    }
    rcp_tilde[i] <- sum(r3[i,])
    rcn_tilde[i] <- sum(r4[i,])
    
    eta[i] <- dp[i] * (1 + rcp_tilde[i] + tcp_tilde[i])
    theta[i] <- dn[i] * (1 + rcn_tilde[i] + tcn_tilde[i])
    prob[i] <- (step(ref[i]) * ref[i] + (1 - step(ref[i])) * pi) * eta[i] + (step(ref[i]) * (1 - ref[i]) + (1 - step(ref[i])) * (1 - pi)) * theta[i]
    z[i] ~ dpois(-log(prob[i]))
  }
  
  for (k in 1:K) {
    sp[k] ~ dbeta(omega2 * (kappa2 - 2) + 1, (1 - omega2) * (kappa2 - 2) + 1)
    se[k] ~ dbeta(omega1 * (kappa1 - 2) + 1, (1 - omega1) * (kappa1 - 2) + 1)T(1 - sp[k], 1)
    rp[k] ~ dunif(max(-se[k]^2, -(1 - se[k])^2), (1-se[k]) * se[k])
    rn[k] ~ dunif(-(1 - sp[k])^2, (1 - sp[k]) * sp[k])
  }

  for (l in 1:(K - 1)) {
    for (h in (l + 1):K) {
      cop[l, h] ~ dunif((1 - se[l]) * (se[h] - 1), (min(se[l], se[h]) - se[l] * se[h]))
      con[l, h] ~ dunif((1 - sp[l]) * (sp[h] - 1), (min(sp[l], sp[h]) - sp[l] * sp[h]))
    }
  }
  for (h in 1:K) {
    for (l in h:K) {
      cop[l, h] <- 0
      con[l, h] <- 0
    }
  }
  
  for (w in 1:ncomb) {
    for (j in 1:J) {
      for (k in 1:K) {
        par_d1[w, j, k] <- comb[w, k]*(1 - se[k])^j + 1 - comb[w, k]
        par_d2[w, j, k] <- comb[w, k]*sp[k]^j + 1 - comb[w, k]
      }
      par_dp[w, j] <- prod(par_d1[w, j,])
      par_dn[w, j] <- prod(par_d2[w, j,])
    }
    
    for (u in 1:(K - 1)) {
      for (v in (u + 1):K) {
        par_t1[w, u, v] <- comb[w, u] * comb[w, v] * cop[u, v] * (1 - se[u])^(-1)*(1 - se[v])^(-1)
        par_t2[w, u, v] <- comb[w, u] * comb[w, v] * con[u, v] * sp[u]^(-1)*sp[v]^(-1)
      }
    }
    for (v in 1:K) {
      for (u in v:K) {
        par_t1[w, u, v] <- 0
        par_t2[w, u, v] <- 0
      }
    }
    par_tcp[w] <- sum(par_t1[w,,])
    par_tcn[w] <- sum(par_t2[w,,])
    
    for (k in 1:K) {
      par_r1[w, k] <- comb[w, k] * rp[k] * (1 - se[k])^(-2)
      par_r2[w, k] <- comb[w, k] * rn[k] * sp[k]^(-2)
    }
    par_rcp[w] <- sum(par_r1[w,])
    par_rcn[w] <- sum(par_r2[w,])
    
    for (j in 1:J) {
      par_se[w, j] <- 1 - par_dp[w, j] * (1 + j * par_tcp[w] + 0.5 * j * (j-1) * par_rcp[w])  
      z2[w,j] ~ dpois(-log(par_se[w, j]))
      par_sp[w, j] <- par_dn[w, j] * (1 + j * par_tcn[w] + 0.5 * j * (j-1) * par_rcn[w])
      z3[w,j] ~ dpois(-log(par_sp[w, j]))  
    }
  }
  
  omega1 ~ dbeta(1, 1)
  omega2 ~ dbeta(1, 1)T(0.5, )
  kappa1 <- kappaMinusTwo1 + 2
  kappaMinusTwo1 ~ dgamma(0.01, 0.01)
  kappa2 <- kappaMinusTwo2 + 2
  kappaMinusTwo2 ~ dgamma(0.01, 0.01)
  pi ~ dbeta(1, 1)
}"
writeLines(model_str, con="repeated_measurements_parallel_tests.bug")
