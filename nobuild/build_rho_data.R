devtools::load_all()

produce_rho_data <- function(sigma, sigma_name, sigma_min, sigma_max,
                             c_steps = 20, Ns = c(2, 3, 4, 5),
                             repetitions = 100, number_of_samples = 10000) {

  df <- data.frame(matrix(ncol = 1, nrow = c_steps + 1))
  colnames(df) <- sigma_name
  df[sigma_name] <- ((0:c_steps) * ((sigma_max - sigma_min)
                                    / c_steps) + sigma_min)

  for (n in Ns) {
    cat(paste0("Producing data for n=", n, "..."))

    set.seed(1)
    sum_images <- rep(0, c_steps + 1)
    for (i in 1:repetitions) {
      images <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
        significativity(sigma, c, n,
                        number_of_samples = number_of_samples)
      })
      sum_images <- sum_images + images
    }
    mean_images <- sum_images / repetitions

    df[paste0("rho", n)] <- mean_images

    set.seed(1)
    sum_images <- rep(0, c_steps + 1)
    for (i in 1:repetitions) {
      images <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
        significativity(sigma, c, n,
                        number_of_samples = number_of_samples)
      })
      sum_images <- sum_images + (images - mean_images)^2
    }

    df[paste0("rho", n, "_err")] <- sqrt(sum_images / repetitions)

    cat("done\n")
  }

  filename <- paste0(sigma_name, "_rho.csv")
  write.table(df, filename, sep = " ", row.names = FALSE, quote = FALSE)

  print(paste0('\"', filename, '\" has been saved'))
}


produce_rho_data(cohen_kappa, "kappa", -1, 1)
produce_rho_data(IA, "IA", 0, 1)