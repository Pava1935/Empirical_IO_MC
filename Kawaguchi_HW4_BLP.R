library(tidyverse)
library(evd)

# set the seed
set.seed(1)
# number of products
J <- 10
# dimension of product characteristics including the intercept
K <- 3
# number of markets
T <- 100
# number of consumers per market
N <- 500
# number of Monte Carlo
L <- 500

# set parameters of interests
beta <- rnorm(K); 
beta[1] <- 4
beta

sigma <- abs(rnorm(K))

mu <- 0.5
omega <- 1

# set auxiliary parameters
price_xi <- 1
sd_x <- 2
sd_xi <- 0.1
sd_c <- 0.5
sd_p <- 0.01

# make product characteristics data
X <- 
  matrix(
    sd_x * rnorm(J * (K - 1)), 
    nrow = J
  )
X <- 
  cbind(
    rep(1, J), 
    X
  )
colnames(X) <- paste("x", 1:K, sep = "_")
X <- 
  data.frame(j = 1:J, X) %>%
  tibble::as_tibble()
# add outside option
X <- 
  rbind(
    rep(0, dim(X)[2]),
    X
  ) 

# make market-product data
M <- 
  expand.grid(
    j = 1:J, 
    t = 1:T
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    xi = sd_xi * rnorm(J * T),
    c = exp(sd_c * rnorm(J * T)),
    p = exp(price_xi * xi + sd_p * rnorm(J * T)) + c
  ) 
M <- 
  M %>%
  dplyr::group_by(t) %>%
  dplyr::sample_frac(size = purrr::rdunif(1, J) / J) %>%
  dplyr::ungroup()

# add outside option
outside <- 
  data.frame(
    j = 0, 
    t = 1:T, 
    xi = 0, 
    c = 0, 
    p = 0
  )
M <- 
  rbind(
    M,
    outside
  ) %>%
  dplyr::arrange(
    t, 
    j
  )

# make consumer-market data
V <- 
  matrix(
    rnorm(N * T * (K + 1)), 
    nrow = N * T
  ) 
colnames(V) <- 
  c(
    paste("v_x", 1:K, sep = "_"), 
    "v_p"
  )
V <- 
  data.frame(
    expand.grid(
      i = 1:N, 
      t = 1:T
    ),
    V
  ) %>%
  tibble::as_tibble()

# make choice data
df <- 
  expand.grid(
    t = 1:T, 
    i = 1:N, 
    j = 0:J
  ) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(
    V, 
    by = c("i", "t")
  ) %>%
  dplyr::left_join(
    X, 
    by = c("j")
  ) %>%
  dplyr::left_join(
    M, 
    by = c("j", "t")
  ) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(
    t, 
    i, 
    j
  )

# draw idiosyncratic shocks
e <- evd::rgev(dim(df)[1])


# compute indirect utility
compute_indirect_utility <- 
  function(
    df, 
    beta, 
    sigma, 
    mu, 
    omega
  ) {
    # extract matrices
    X <- as.matrix(dplyr::select(df, dplyr::starts_with("x_")))
    p <- as.matrix(dplyr::select(df, p)) 
    v_x <- as.matrix(dplyr::select(df, dplyr::starts_with("v_x")))
    v_p <- as.matrix(dplyr::select(df, v_p))
    xi <- as.matrix(dplyr::select(df, xi))
    # random coefficients
    beta_i <- as.matrix(rep(1, dim(v_x)[1])) %*% t(as.matrix(beta)) + v_x %*% diag(sigma) 
    alpha_i <- - exp(mu + omega * v_p)
    # conditional mean indirect utility
    value <- as.matrix(rowSums(beta_i * X) + p * alpha_i + xi) 
    colnames(value) <- "u"
    return(value)
  }


# compute indirect utility
u <- 
  compute_indirect_utility(
    df = df, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )


# compute choice
compute_choice_smooth <-
  function(
    X,
    M, 
    V, 
    beta, 
    sigma, 
    mu, 
    omega
  ) {
    # constants
    T <- max(M$t)
    N <- max(V$i)
    J <- max(X$j)
    # make choice data
    df <- 
      expand.grid(
        t = 1:T, 
        i = 1:N, 
        j = 0:J
      ) %>%
      tibble::as_tibble() %>%
      dplyr::left_join(
        V, 
        by = c("i", "t")
      ) %>%
      dplyr::left_join(
        X, 
        by = c("j")
      ) %>%
      dplyr::left_join(
        M, 
        by = c("j", "t")
      ) %>%
      dplyr::filter(!is.na(p)) %>%
      dplyr::arrange(
        t, 
        i, 
        j
      )
    # compute indirect utility
    u <- 
      compute_indirect_utility(
        df, 
        beta, 
        sigma, 
        mu, 
        omega
      )
    # add u 
    df_choice <- 
      data.frame(
        df, 
        u
      ) %>%
      tibble::as_tibble()
    # make choice
    df_choice <- 
      df_choice %>%
      dplyr::group_by(
        t, 
        i
      ) %>%
      dplyr::mutate(q = exp(u)/sum(exp(u))) %>%
      dplyr::ungroup()
    # return
    return(df_choice)
  }


df_choice_smooth <-
  compute_choice_smooth(
    X = X, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )



# compute share
compute_share_smooth <-
  function(
    X, 
    M, 
    V, 
    beta, 
    sigma, 
    mu, 
    omega
  ) {
    # constants
    T <- max(M$t)
    N <- max(V$i)
    J <- max(X$j)
    # compute choice
    df_choice <- 
      compute_choice_smooth(
        X, 
        M, 
        V, 
        beta, 
        sigma,
        mu, 
        omega
      )
    # make share data
    df_share_smooth <- 
      df_choice %>%
      dplyr::select(
        -dplyr::starts_with("v_"), 
        -u, 
        -i
      ) %>%
      dplyr::group_by(
        t, 
        j
      ) %>%
      dplyr::mutate(q = sum(q)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
        t, 
        j, 
        .keep_all = TRUE
      ) %>%
      dplyr::group_by(t) %>%
      dplyr::mutate(s = q/sum(q)) %>%
      dplyr::ungroup()
    # log share difference
    df_share_smooth <- 
      df_share_smooth %>%
      dplyr::group_by(t) %>%
      dplyr::mutate(y = log(s/sum(s * (j == 0)))) %>%
      dplyr::ungroup()
    return(df_share_smooth)
  }

df_share_smooth <- 
  compute_share_smooth(
    X = X, 
    M = M, 
    V = V, 
    beta = beta, 
    sigma = sigma, 
    mu = mu, 
    omega = omega
  )
summary(df_share_smooth)

# mixed logit estimation
## draw mcmc V
V_mcmc <- 
  matrix(
    rnorm(L * T * (K + 1)), 
    nrow = L * T
  ) 
colnames(V_mcmc) <- 
  c(
    paste("v_x", 1:K, sep = "_"), 
    "v_p"
  )
V_mcmc <- 
  data.frame(
    expand.grid(
      i = 1:L, 
      t = 1:T
    ),
    V_mcmc
  ) %>%
  tibble::as_tibble() 

## draw mcmc e
df_mcmc <- 
  expand.grid(
    t = 1:T, 
    i = 1:L, 
    j = 0:J
  ) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(
    V_mcmc, 
    by = c("i", "t")
  ) %>%
  dplyr::left_join(
    X, 
    by = c("j")
  ) %>%
  dplyr::left_join(
    M, 
    by = c("j", "t")
  ) %>%
  dplyr::filter(!is.na(p)) %>%
  dplyr::arrange(
    t, 
    i, 
    j
  )
# draw idiosyncratic shocks
e_mcmc <- evd::rgev(dim(df_mcmc)[1])

# set parameters
theta <- 
  c(
    beta, 
    sigma, 
    mu, 
    omega
  )
theta

M_no <- 
  M %>%
  dplyr::mutate(xi = 0)

# compute indirect utility from delta
compute_indirect_utility_delta <- 
  function(
    df,
    delta, 
    sigma, 
    mu,
    omega
  ) {
    # extract matrices
    X <- as.matrix(dplyr::select(df, dplyr::starts_with("x_")))
    p <- as.matrix(dplyr::select(df, p)) 
    v_x <- as.matrix(dplyr::select(df, dplyr::starts_with("v_x")))
    v_p <- as.matrix(dplyr::select(df, v_p))
    # expand delta
    delta_ijt <- 
      df %>%
      dplyr::left_join(
        delta, 
        by = c("t", "j")
      ) %>%
      dplyr::select(delta) %>%
      as.matrix()
    # random coefficients
    beta_i <- v_x %*% diag(sigma) 
    alpha_i <- - exp(mu + omega * v_p) - (- exp(mu + omega^2/2))
    # conditional mean indirect utility
    value <- as.matrix(delta_ijt + rowSums(beta_i * X) + p * alpha_i) 
    colnames(value) <- "u"
    return(value)
  }


delta <- df_share_smooth |>
  mutate(delta = theta[1] * x_1 + theta[2] * x_2 + theta[3] * x_3 - exp(theta[7] +
                                                                          (theta[8] ^ 2 / 2)) * p + xi) |>
  select(t, j, delta)



# compute indirect utility from delta
u_delta <-
  compute_indirect_utility_delta(
    df, 
    delta, 
    sigma,
    mu, 
    omega
  )


# compute indirect utility from delta
compute_indirect_utility_delta <- 
  function(
    df,
    delta, 
    sigma, 
    mu,
    omega
  ) {
    # extract matrices
    X <- as.matrix(dplyr::select(df, dplyr::starts_with("x_")))
    p <- as.matrix(dplyr::select(df, p)) 
    v_x <- as.matrix(dplyr::select(df, dplyr::starts_with("v_x")))
    v_p <- as.matrix(dplyr::select(df, v_p))
    # expand delta
    delta_ijt <- 
      df %>%
      dplyr::left_join(
        delta, 
        by = c("t", "j")
      ) %>%
      dplyr::select(delta) %>%
      as.matrix()
    # random coefficients
    beta_i <- v_x %*% diag(sigma) 
    alpha_i <- - exp(mu + omega * v_p) - (- exp(mu + omega^2/2))
    # conditional mean indirect utility
    value <- as.matrix(delta_ijt + rowSums(beta_i * X) + p * alpha_i) 
    colnames(value) <- "u"
    return(value)
  }

# compute indirect utility from delta
u_delta <-
  compute_indirect_utility_delta(
    df, 
    delta, 
    sigma,
    mu, 
    omega
  )
head(u_delta)

# compute choice from delta
compute_choice_smooth_delta <-
  function(
    X, 
    M, 
    V, 
    delta, 
    sigma, 
    mu, 
    omega) {
    # constants
    T <- max(M$t)
    N <- max(V$i)
    J <- max(X$j)
    # make choice data
    df <- 
      expand.grid(
        t = 1:T, 
        i = 1:N, 
        j = 0:J
      ) %>%
      tibble::as_tibble() %>%
      dplyr::left_join(
        V, 
        by = c("i", "t")
      ) %>%
      dplyr::left_join(
        X, 
        by = c("j")
      ) %>%
      dplyr::left_join(
        M, 
        by = c("j", "t")
      ) %>%
      dplyr::filter(!is.na(p)) %>%
      dplyr::arrange(
        t, 
        i, 
        j
      )
    # compute indirect utility
    u <- 
      compute_indirect_utility_delta(
        df, 
        delta, 
        sigma, 
        mu, 
        omega
      )
    # add u 
    df_choice <- 
      data.frame(
        df, 
        u
      ) %>%
      tibble::as_tibble()
    # make choice
    df_choice <- 
      df_choice %>%
      dplyr::group_by(
        t, 
        i
      ) %>%
      dplyr::mutate(q = exp(u)/sum(exp(u))) %>%
      dplyr::ungroup()
    # return
    return(df_choice)
  }

# compute choice
df_choice_smooth_delta <- 
  compute_choice_smooth_delta(
    X, 
    M, 
    V, 
    delta, 
    sigma, 
    mu, 
    omega
  )
df_choice_smooth_delta

summary(df_choice_smooth_delta)

# compute share from delta
compute_share_smooth_delta <-
  function(
    X, 
    M, 
    V, 
    delta, 
    sigma, 
    mu, 
    omega
  ) {
    # constants
    T <- max(M$t)
    N <- max(V$i)
    J <- max(X$j)
    # compute choice
    df_choice <- 
      compute_choice_smooth_delta(
        X, 
        M, 
        V, 
        delta, 
        sigma,
        mu, 
        omega
      )
    # make share data
    df_share_smooth <-
      df_choice %>%
      dplyr::select(
        -dplyr::starts_with("v_"),
        -u,
        -i
      ) %>%
      dplyr::group_by(
        t, 
        j
      ) %>%
      dplyr::mutate(q = sum(q)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
        t, 
        j, 
        .keep_all = TRUE
      ) %>%
      dplyr::group_by(t) %>%
      dplyr::mutate(s = q/sum(q)) %>%
      dplyr::ungroup()
    # log share difference
    df_share_smooth <- 
      df_share_smooth %>%
      dplyr::group_by(t) %>%
      dplyr::mutate(y = log(s/sum(s * (j == 0)))) %>%
      dplyr::ungroup()
    return(df_share_smooth)
  }

# compute share
df_share_smooth_delta <-
  compute_share_smooth_delta(
    X, 
    M, 
    V, 
    delta, 
    sigma, 
    mu, 
    omega
  ) 

kappa <- 1
lambda <- 1e-4



# solve delta by the fixed-point algorithm
solve_delta <-
  function(
    df_share_smooth, 
    X, 
    M, 
    V, 
    delta, 
    sigma, 
    mu, 
    omega, 
    kappa, 
    lambda
  ) {
    # initial distance
    distance <- 10000
    # fixed-point algorithm
    delta_old <- delta
    while (distance > lambda) {
      # save the old delta
      delta_old$delta <- delta$delta
      # compute the share with the old delta
      df_share_smooth_predicted <- 
        compute_share_smooth_delta(
          X, 
          M, 
          V, 
          delta_old, 
          sigma, 
          mu, 
          omega
        )  
      # update the delta
      delta$delta <- 
        delta_old$delta + 
        (log(df_share_smooth$s) - log(df_share_smooth_predicted$s)) * kappa
      delta <- 
        delta %>%
        dplyr::mutate(delta = ifelse(j == 0, 0, delta))
      # update the distance
      distance <- max(abs(delta$delta - delta_old$delta))
      # print(distance)
    }
    return(delta)
  }


delta_new <-
  solve_delta(
    df_share_smooth, 
    X, 
    M, 
    V_mcmc, 
    delta, 
    sigma, 
    mu, 
    omega, 
    kappa, 
    lambda
  )

Psi <- diag(length(beta) + 1)

# compute the optimal linear parameters
compute_theta_linear <-
  function(
    df_share_smooth, 
    delta, 
    mu, 
    omega, 
    Psi
  ) {
    # extract matrices
    X <- 
      df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(dplyr::starts_with("x_")) %>%
      as.matrix()
    p <- 
      df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(p) %>%
      as.matrix()
    W <- 
      df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(dplyr::starts_with("x_"), c) %>%
      as.matrix()
    delta_m <- 
      delta %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(delta) %>%
      as.matrix()
    alpha <- - exp(mu + omega^2/2)
    # compute the optimal linear parameters
    theta_linear_1 <-
      crossprod(
        X, 
        W
      ) %*% 
      solve(
        Psi, 
        crossprod(
          W, 
          X
        )
      )
    theta_linear_2 <-
      crossprod(
        X, 
        W
      ) %*% 
      solve(
        Psi, 
        crossprod(
          W, 
          delta_m - alpha * p
        )
      )
    theta_linear <- 
      solve(
        theta_linear_1, 
        theta_linear_2
      )
    return(theta_linear)
  }

theta_linear <-
  compute_theta_linear(
    df_share_smooth, 
    delta, 
    mu, 
    omega, 
    Psi
  ) 
cbind(
  theta_linear, 
  beta
)

# solve xi associated with delta and linear parameters
solve_xi <-
  function(
    df_share_smooth, 
    delta, 
    beta, 
    mu, 
    omega
  ) {
    # extract matrices
    X1 <- 
      df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(
        dplyr::starts_with("x_"), 
        p
      ) %>%
      as.matrix()
    delta_m <- 
      delta %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(delta) %>%
      as.matrix()
    alpha <- - exp(mu + omega^2/2)
    theta_linear <- 
      c(
        beta, 
        alpha
      )
    # compute xi
    xi <- delta_m - X1 %*% theta_linear
    colnames(xi) <- "xi"
    # return
    return(xi)
  }

xi_new <- 
  solve_xi(
    df_share_smooth, 
    delta, 
    beta, 
    mu, 
    omega
  )

xi_true <-
  df_share_smooth %>%
  dplyr::filter(j != 0) %>%
  dplyr::select(xi)
summary(xi_true - xi_new)

# non-linear parmaeters
theta_nonlinear <- 
  c(
    mu, 
    omega, 
    sigma
  )

# compute GMM objective function
compute_gmm_objective_a4 <-
  function(
    theta_nonlinear, 
    delta, 
    df_share_smooth, 
    Psi,
    X, 
    M, 
    V_mcmc, 
    kappa, 
    lambda
  ) {
    # exctract parameters
    mu <- theta_nonlinear[1]
    omega <- theta_nonlinear[2]
    sigma <- theta_nonlinear[3:length(theta_nonlinear)]
    # extract matrix
    W <- 
      df_share_smooth %>%
      dplyr::filter(j != 0) %>%
      dplyr::select(
        dplyr::starts_with("x_"), 
        c
      ) %>%
      as.matrix()
    # compute the delta that equates the actual and predicted shares
    delta <- 
      solve_delta(
        df_share_smooth, 
        X, 
        M, 
        V_mcmc, 
        delta, 
        sigma, 
        mu, 
        omega, 
        kappa, 
        lambda
      )
    # compute the optimal linear parameters
    beta <-
      compute_theta_linear(
        df_share_smooth, 
        delta, 
        mu, 
        omega, 
        Psi
      ) 
    # compute associated xi
    xi <- 
      solve_xi(
        df_share_smooth, 
        delta, 
        beta, 
        mu, 
        omega
      )
    # compute objective
    objective <- crossprod(xi, W) %*% solve(Psi, crossprod(W, xi))
    print(objective)
    # return
    return(objective)
  }

objective <-
  compute_gmm_objective_a4(
    theta_nonlinear, 
    delta, 
    df_share_smooth, 
    Psi, 
    X, 
    M, 
    V_mcmc, 
    kappa, 
    lambda
  ) 


optim_param <- optim(
  theta_nonlinear,
  compute_gmm_objective_a4,
  delta = delta,
  df_share_smooth = df_share_smooth,
  Psi = Psi,
  X = X,
  M = M,
  V_mcmc = V,
  kappa = kappa,
  lambda = lambda,
  method = "BFGS"
)
