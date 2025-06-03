devtools::load_all()

produce_delta_data <- function(sigma, sigma_name, sigma_min, sigma_max,
                               c_steps = 20, Ns = c(2, 3),
                               Ms = c(10, 100, 1000),
                               repetitions = 10, number_of_samples = 10000) {

  df <- data.frame(matrix(ncol = 1, nrow = c_steps + 1))
  colnames(df) <- sigma_name
  df[sigma_name] <- ((0:c_steps) * ((sigma_max - sigma_min)
                                    / c_steps) + sigma_min)

  for (n in Ns) {
    cat(paste0("Producing data for n=", n, "..."))

    set.seed(1)
    rho <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
      significativity(sigma, c, n,
                      number_of_samples = number_of_samples)
    })

    df[paste0("rho", n)] <- rho

    for (m in Ms) {
      set.seed(1)
      sum_images <- rep(0, c_steps + 1)
      for (i in 1:repetitions) {
        images <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
          significativity(sigma, c, n, m,
                          number_of_samples = number_of_samples)
        })
        sum_images <- sum_images + images
      }

      df[paste0("varrho", n, "_", m)] <- sum_images / repetitions
      df[paste0("delta", n, "_", m)] <- rho - (sum_images / repetitions)
    }

    cat("done\n")
  }

  filename <- paste0(sigma_name, "_delta.csv")
  write.table(df, filename, sep = " ", row.names = FALSE, quote = FALSE)

  print(paste0('\"', filename, '\" has been saved'))
}

produce_delta_data(cohen_kappa, "kappa", -1, 1)
produce_delta_data(IA, "IA", 0, 1)