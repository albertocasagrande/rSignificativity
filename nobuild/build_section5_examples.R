devtools::load_all()

M5_1000 <- matrix(c(84, 12, 26, 36, 15,
                    11, 106, 5, 22, 0,
                    6, 7, 240, 6, 52,
                    5, 43, 43, 71, 14,
                    3, 15, 42, 26, 110), nrow = 5, ncol = 5)

M2_1000 <- matrix(c(265, 129,
                    101, 505), nrow = 2, ncol = 2)

M2_1000_2 <- matrix(c(260, 132,
                      103, 505), nrow = 2, ncol = 2)

M2_12 <- matrix(c(6, 3,
                  0, 3), nrow = 2, ncol = 2)

M2_10 <- matrix(c(3, 0,
                  2, 5), nrow = 2, ncol = 2)

Ms <- list(list(name = "M2_10", matrix = M2_10),
           list(name = "M2_1000", matrix = M2_1000),
           list(name = "M5_1000", matrix = M5_1000))

sigma <- cohen_kappa

for (i in seq_len(length(Ms))) {
  M <- Ms[[i]]$matrix

  value <- sigma(M)
  s <- significativity(sigma, value, nrow(M), sum(M), 1000000)

  cat(paste0("Cohen's kappa of ", Ms[[i]]$name, " is ", value,
             ". Its significativity is ", s, ".\n"))
}
