devtools::load_all()

Ns <- c(2, 3)
Ms <- c(2, 4, 8, 16)
k <- 20
repetitions <- 100
sigma <- cohen_kappa
sigma_name <- "kappa"
min_sigma <- -1
max_sigma <- 1

for (n in Ns) {
  df <- data.frame(matrix(ncol = 1, nrow = k + 1))
  colnames(df) <- sigma_name
  df[sigma_name] <- (0:k) * ((max_sigma - min_sigma) / k) + min_sigma

  for (m in Ms) {
    cat(paste0("Dealing with n=", n, " and m=", m, "..."))
    test_name <- paste0("R_", sigma_name, "_", n, "_", m)

    df[test_name] <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
      significativity(sigma, c, n, m, number_of_samples = NULL)
    })

    number_of_samples <- ceiling(sqrt(choose(n^2 + m - 1, m)))

    df[paste0("num_of_samples_", test_name)] <- number_of_samples

    set.seed(1)
    sum_images <- rep(0, k + 1)
    delta_sum_images <- rep(0, k + 1)
    for (i in 1:repetitions) {
      images <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
        significativity(sigma, c, n, m,
                        number_of_samples = number_of_samples)
      })

      sum_images <- sum_images + images
      delta_sum_images <- delta_sum_images + (images - df[test_name])
    }
    mean_images <- sum_images / repetitions
    delta_mean_images <- delta_sum_images / repetitions

    df[paste0("mean_sampled_", test_name)] <- mean_images

    df[paste0("mean_delta_", test_name)] <- delta_mean_images

    set.seed(1)
    sum_images <- rep(0, k + 1)
    delta_sum_images <- rep(0, k + 1)
    for (i in 1:repetitions) {
      images <- apply(df[, sigma_name, drop = FALSE], 1, function(c) {
        significativity(sigma, c, n, m,
                        number_of_samples = number_of_samples)
      })

      sum_images <- sum_images + (images - mean_images)^2
      delta_sum_images <- (delta_sum_images + ((images - df[test_name])
                                               - delta_mean_images)^2)
    }

    df[paste0("stdev_sampled_", test_name)] <- sqrt(sum_images / repetitions)

    df[paste0("stdev_delta_", test_name)] <- sqrt(delta_sum_images
                                                  / repetitions)
    cat("done\n")
  }

  filename <- paste0(sigma_name, "_", n, ".csv")
  write.table(df, filename, sep = " ", row.names = FALSE, quote = FALSE)

  cat(paste0('\"', filename, '\" has been saved\n'))
}