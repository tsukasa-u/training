// 1.4. Notation
//     (x, pi, s)  primal, sual, and slack variables for problem (GNP)
//     (x_star, pi_star, s_star)   oprimal variables for problem (GNP)
//     (x_k, pi_k, s_k)    the kth estimate of (s_star, pi_star, s_star)
//     f_m g_k, c_k, J_k   function and gradients evaluated at x_k
//     (x_k_hat, pi_k_hat, s_k_hat)    optimal variables for QP subproblem (GQP_k)

// 2. The SQP Iteration
//     (GNP)
//         minimize(x) f(x)
//         subject to c(x) >= 0
 
//     (2.1)
//         c(x_star) >= 0
//         pi_star >= 0
//         c(x_star)^T pi_star = 0
//         J(x_star)^T pi_star = g(x_star)

//     (2.2)
//         L(x, x_k, pi_k) = f(x) - pi_k^T d_L(x, x_k)

//         c_L(x, x_k) = c_k + J_k(x - x_k)
//         d_L(x, x_k) = c(x) - c_L(x, x_k)

//         nabla(L(x, x_k, pi_k)) = g(x) - (J(x) - J_k)^T pi_k
//         nabla^2 (L(x, x_k, pi_k)) = nabla^2 (f(x)) - SUM_i((pi_k)_i nabla^2(c_i(x)))

//         if nablra^2(L) is independent of x_k (and is the same as the Hessian of the conventional Lagrangian)
//             L(x, x_k, pi_k) = f_k
//             nabla(L(x, x_k, pi_k)) = g_k
//                 At x = x_k

//     
//     //     L_q(x, x_k, pi_k) = f_k + g_k^T(x - x_k) + 1/2 (x - x_k)^T nabla^2(L(x_k, x_k, pi_k)) (x - x_k)

//     (GQP_star)
//         minimaize_x L_q(x, x_k, pi_k)
//         subject to linearized constraints c_l(x, x_k) >= 0

//     (GQP_k)
//         minimize_x f_k + g_k^T (x - x_k) + 1/2 (x-x_k)^T H_k (x -x_k)
//         subjects to c_k + J_k(x - x_k) >= 0

//         (optimize) -> 
//             c_k + J_k(x_k^hat - x_k) = s_k^hat
//             pi_k^hat >= 0
//             s_k^hat >= 0
//             g_k + H_k(x_k^hat - x_k) = J_k^T pi_k^hat
//             (pi_k^hat)^T s_k^hat = 0

//                 Where c(x_star) - s_star = 0 and s_star >= 0

//     // 2.5. The Working-Set Matix W_k
//     // 2.6. THe Null-Space Matrix Z_k

//     2.7 The merit Function and Line Search
//     (2.3)
//         M_rho(x, pi, s) = f(x) - pi^T (x(x) - s) + 1/2 SUM_{i=1}^m (rho_i(c_i(x) - s_i)^2)

//         v(alpha) = (x_k, pi_k, s_k) + alpha(x_k^hat - x_k, pi_k^hat - pi_k, s_k^hat - s_k)
//         phi_rho(alpha) = M_who(v(alpha))
//             and 0 < alpha <= 1

//     (LSP_rho)
//         minimize_rho ||rho||_2^2
//         subject to phi'_rho(0) = -1/2 p_k^T H_k p_k and rho >= 0
//               Where p_k = x_k^hat - x_k

//     (2.4)
//         rho_i^bar = max{rho_i^star, rho_i^hat}
//             Where
//                 if rho_i < 4(who_i^star + Delta_rho)
//                     rho_i^hat = rho_i
//                 else
//                     eho_i^hat = (who_i(rho_i^hat + Delta_rho))^(1/2)

//             and if k = 0
//                 Delta_rho = 1