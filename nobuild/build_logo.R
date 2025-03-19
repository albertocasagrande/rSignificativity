devtools::load_all()
library(hexSticker)
library(ggplot2)
#install.packages("hexSticker")

kappa <- function(conf_matrix) {
  p_o <- sum(diag(conf_matrix)) / sum(conf_matrix)

  # Calculate the expected agreement (P_e)
  row_totals <- rowSums(conf_matrix)
  col_totals <- colSums(conf_matrix)
  total <- sum(conf_matrix)
  p_e <- sum((row_totals * col_totals) / total^2)

  # Calculate Cohen's kappa
  kappa <- (p_o - p_e) / (1 - p_e)

  return(kappa)
}

k <- 10;

f <- function(n) {
  return(function(c) {
    return(significativity(kappa, c, n, 100))
  })
}


v_in <- (-k:k) / k
ns <- (2:4)

v_out2 <- sapply(v_in, f(2))
v_out3 <- sapply(v_in, f(3))
v_out4 <- sapply(v_in, f(4))
v_out5 <- sapply(v_in, f(5))

data <- data.frame(
  x = rep(v_in, 4),
  y = c(v_out2, v_out3, v_out4, v_out5),
  Function = rep(c("2", "3", "4", "5"), each = length(v_in)) # Function labels
)

graph <- ggplot(data, aes(x, y, color = Function)) +
  geom_line() + geom_point() +
  theme_minimal() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

logo_dir <- "../man/figures"

dir.create(logo_dir, recursive = TRUE)

s <- sticker(graph, p_y = 1.47, p_color = "#000000",
             package = "rSignificativity", p_size = 30, s_x = 0.98, s_y = .8,
             s_width = 1.2, s_height = 1.05,
             h_fill = "#ffffff", h_color = "#cccccc",
             filename = file.path(logo_dir, "logo.png"), dpi = 600)

file.remove("Rplots.pdf")
