devtools::load_all()

n <- 2
num_of_samples <- 2000
c <- 0.3
sigma <- IA

k <- n * n
filenames <- c("in_set.txt", "out_set.txt")

for (filename in filenames) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
}

set.seed(1)
for (i in 1:num_of_samples) {
  sample <- sample_prob_simplex(k)

  matrix <- matrix(sample, nrow = n, ncol = n)

  if (sigma(matrix) < c) {
    filename <- filenames[1]
  } else {
    filename <- filenames[2]
  }

  cat(sample, file = filename, append = TRUE, sep = " ")
  cat("\n", file = filename, append = TRUE)
}
