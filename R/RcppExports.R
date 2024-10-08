# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

g_cpp <- function(x, n0 = 1, a = 1, b = 1, eps = 0) {
    .Call('_twbfn_g_cpp', PACKAGE = 'twbfn', x, n0, a, b, eps)
}

lambda_obj_cpp <- function(l, n0, a, b, eps, mx = 1000L) {
    .Call('_twbfn_lambda_obj_cpp', PACKAGE = 'twbfn', l, n0, a, b, eps, mx)
}

lambda_obj_adj_cpp <- function(l, n0, a, b, eps, mx = 1000L) {
    .Call('_twbfn_lambda_obj_adj_cpp', PACKAGE = 'twbfn', l, n0, a, b, eps, mx)
}

find_lambda_cpp <- function(n0, a, b, eps) {
    .Call('_twbfn_find_lambda_cpp', PACKAGE = 'twbfn', n0, a, b, eps)
}

dgpa_cpp <- function(x, n0, a, b, eps, lambda, logd = TRUE) {
    .Call('_twbfn_dgpa_cpp', PACKAGE = 'twbfn', x, n0, a, b, eps, lambda, logd)
}

llh_cpp <- function(data, lambda, a, b, n0, eps = 0) {
    .Call('_twbfn_llh_cpp', PACKAGE = 'twbfn', data, lambda, a, b, n0, eps)
}

joint_prior_cpp <- function(a, b, n0) {
    .Call('_twbfn_joint_prior_cpp', PACKAGE = 'twbfn', a, b, n0)
}

joint_post_cpp <- function(data, lambda, a, b, n0) {
    .Call('_twbfn_joint_post_cpp', PACKAGE = 'twbfn', data, lambda, a, b, n0)
}

ppa_mcmc <- function(data, init_n0, init_a, init_b, init_lambda, n_iter = 10000L, a_sd = 0.1, b_sd = 0.1, n0_sd = 1, adapt_interval = 1000L, burn_in = 10000L) {
    .Call('_twbfn_ppa_mcmc', PACKAGE = 'twbfn', data, init_n0, init_a, init_b, init_lambda, n_iter, a_sd, b_sd, n0_sd, adapt_interval, burn_in)
}

