##
##  a u g l a g 2 . R  Augmented Lagrangian
##


# auglag2 <-
# function(x0, fn, gr = NULL, lower = NULL, upper = NULL,
#             hin = NULL, hinjac = NULL, heq = NULL, heqjac = NULL,
#             local_solver = c("LBFGS", "COBYLA"), eq_penalty = FALSE,
#             nl.info = FALSE, control = list(), ...)
# {
#     lsolver = match.arg(local_solver)
#     if (lsolver == "LBFGS") {
#         gsolver <- "NLOPT_LD_AUGLAG"
#         lsolver <- "NLOPT_LD_LBFGS"
#         need_gr <- TRUE
#     } else if (lsolver == "COBYLA") {
#         gsolver <- "NLOPT_LN_AUGLAG"
#         lsolver <- "NLOPT_LN_COBYLA"
#         need_gr <- FALSE
#     }
#     if (eq_penalty)
#         gsolver <- paste(gsolver, "EQ", sep="_")
# 
#     # Function and gradient, if needed
#     fun <- match.fun(fn)
#     fn  <- function(x) fun(x, ...)
#     
#     if (need_gr && is.null(gr)) {
#         gr <- function(x) nl.grad(x, fn)
#     }
# 
#     # Global and local options
#     opts <- nl.opts(control)
#     opts$algorithm <- gsolver
#     local_opts <- list(algorithm = lsolver,
#                         xtol_rel = 1e-6,
#                         eval_grad_f = gr)
#     opts$local_opts <- local_opts
#         
#     # Inequality constraints
#     if (!is.null(hin)) {
#         .hin <- match.fun(hin)
#         hin <- function(x) (-1) * .hin(x)   # change  hin >= 0  to  hin <= 0 !
#     }
# 
#     if (need_gr) {
#         if (is.null(hinjac)) {
#             hinjac <- function(x) nl.jacobian(x, hin)
#         } else {
#             .hinjac <- match.fun(hinjac)
#             hinjac <- function(x) (-1) * .hinjac(x)
#         }
#     }
# 
#     # Equality constraints
#     if (!is.null(heq)) {
#         .heq <- match.fun(heq)
#         heq <- function(x) .heq(x)
#     }
# 
#     if (need_gr) {
#         if (is.null(heqjac)) {
#             heqjac <- function(x) nl.jacobian(x, heq)
#         } else {
#             .heqjac <- match.fun(heqjac)
#             heqjac <- function(x) .heqjac(x)
#         }
#     }
# 
#     S0 <- nloptr(x0,
#                 eval_f = fn,
#                 eval_grad_f = gr,
#                 lb = lower,
#                 ub = upper,
#                 eval_g_ineq = hin,
#                 eval_jac_g_ineq = hinjac,
#                 eval_g_eq = heq,
#                 eval_jac_g_eq = heqjac,
#                 opts = opts)
# 
#     if (nl.info) print(S0)
#     S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
#                 convergence = S0$status, message = S0$message)
#     return(S1)
# }
# 
# 
# ##  Simple Example: HS#14
# fn <- function(x) (x[1]-2)^2 + (x[2]-1)^2
# hin <- function(x) -0.25*x[1]^2 - x[2]^2 + 1    # hin >= 0
# him <- function(x) -hin(x)                      # him <= 0
# heq <- function(x) x[1] - 2*x[2] + 1            # heq == 0
# 
# gr <- function(x) grad(fn, x)
# himjac <- function(x) jacobian(him, x)
# hinjac <- function(x) jacobian(hin, x)
# heqjac <- function(x) jacobian(heq, x)
# 
# local_opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-4)
# opts = list(algorithm = "NLOPT_LN_AUGLAG",
#             xtol_rel = 1e-4, maxeval = 100)
# opts$local_opts <- local_opts
# 
# x0 <- c(1, 1)                           # heq(x0) == 0 *must* be satisfied
# nloptr( x0 = x0,                        #              for the initial value
#         eval_f = fn,
#         eval_grad_f = NULL,
#         lb = NULL,
#         ub = NULL,
#         eval_g_ineq = him,
#         eval_jac_g_ineq = NULL,
#         eval_g_eq = heq,
#         eval_jac_g_eq = NULL,
#         opts = opts)
# 
# local_opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-4)
# opts = list(algorithm = "NLOPT_LN_AUGLAG_EQ",
#             xtol_rel = 1e-4, maxeval = 100)
# opts$local_opts <- local_opts
# 
# x0 <- c(1, 1)                           # WRONG ! even for true minimum xs !
# nloptr( x0 = x0,                        # seems hin(x) >= 0 is not considered
#         eval_f = fn,
#         eval_grad_f = NULL,
#         lb = NULL,
#         ub = NULL,
#         eval_g_ineq = him,
#         eval_jac_g_ineq = NULL,
#         eval_g_eq = heq,
#         eval_jac_g_eq = NULL,
#         opts = opts)
# 
# local_opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-4)
# opts = list(algorithm = "NLOPT_LD_AUGLAG",
#             xtol_rel = 1e-4, maxeval = 100)
# opts$local_opts <- local_opts
# 
# x0 <- c(1, 1)
# nloptr( x0 = x0,
#         eval_f = fn,
#         eval_grad_f = gr,
#         lb = NULL,
#         ub = NULL,
#         eval_g_ineq = him,
#         eval_jac_g_ineq = himjac,
#         eval_g_eq = heq,
#         eval_jac_g_eq = heqjac,
#         opts = opts)




