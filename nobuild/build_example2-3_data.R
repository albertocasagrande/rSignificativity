devtools::load_all()

set.seed(1)

sigma_list <- list(list("sigma" = cohen_kappa, "name" = "Cohen's kappa"),
                   list("sigma" = IA, "name" = "IA"))

M <- matrix(c(8, 0, 3, 9), nrow = 2, ncol = 2)

cat("Example 2")

for (sigma in sigma_list) {
  s_value <- sigma$sigma(M)
  n <- nrow(M)
  m <- sum(M)
  m_value <- significativity(sigma$sigma, s_value, n, m)

  cat(paste0("\n", sigma["name"], "-significativity of ", s_value, " in M_{",
             n, ",", m, "}: ", m_value, "\n"))
}

cat("\nExample 3")

for (sigma in sigma_list) {
  s_value <- sigma$sigma(M)
  n <- nrow(M)
  p_value <- significativity(sigma$sigma, s_value, n)

  cat(paste0("\n", sigma["name"], "-significativity of ", s_value, " in P_{",
             n, "}: ", p_value, "\n"))
}