is_a_confusion_matrix <- function(confusion_matrix) {
  # Check if the input is a matrix
  if (!is.matrix(confusion_matrix)) {
    stop("Input must be a confusion matrix (a matrix).")
  }

  # Check if the confusion matrix is square
  if (nrow(confusion_matrix) != ncol(confusion_matrix)) {
    stop("Confusion matrix must be square.")
  }
}

#' @name Cohen's kappa
#' @title Cohen's \eqn{\kappa}
#' @description This method computes Cohen's \eqn{\kappa} of a confusion matrix.
#' @param conf_matrix The confusion matrix.
#' @return The Cohen's \eqn{\kappa} of `conf_matrix`.
#' @examples
#' # create a confusion matrix
#' conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)
#'
#' # evaluate Cohen's kappa of `conf_matrix`
#' cohen_kappa(conf_matrix)

cohen_kappa <- function(conf_matrix) {
  is_a_confusion_matrix(conf_matrix)

  row_totals <- rowSums(conf_matrix)
  col_totals <- colSums(conf_matrix)
  total <- sum(conf_matrix)

  p_o <- sum(diag(conf_matrix)) / total
  p_e <- sum((row_totals * col_totals) / total^2)

  # Calculate Cohen's kappa
  return((p_o - p_e) / (1 - p_e))
}

entropy <- function(param) {
  sum <- 0
  if (is.vector(param)) {
    for (i in seq_along(param)) {
      if (param[i] != 0) {
        sum <- sum + param[i] * log2(param[i])
      }
    }
  } else {
    if (is.matrix(param)) {
      for (i in 1:nrow(param)) {
        for (j in 1:ncol(param)) {
          if (param[i, j] != 0) {
            sum <- sum + param[i, j] * log2(param[i, j])
          }
        }
      }
    }
  }

  return(-sum)
}

#' @name Information Agreement
#' @title Information Agreement
#' @description This method computes the Information Agreement \eqn{\text{IA}}
#'   of a confusion matrix.
#' @param conf_matrix The confusion matrix.
#' @return The \eqn{\text{IA}} of `conf_matrix`.
#' @examples
#' # create a confusion matrix
#' conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)
#'
#' # evaluate the Information Agreement of `conf_matrix`
#' IA(conf_matrix)
IA <- function(conf_matrix) {
  is_a_confusion_matrix(conf_matrix)

  matrix_sum <- sum(conf_matrix)
  if (matrix_sum == 0) {
    return(0)
  }

  prob_matrix <- conf_matrix / matrix_sum

  prob_x <- colSums(prob_matrix)
  prob_y <- rowSums(prob_matrix)

  row_zeros <- sum(prob_y == 0)
  col_zeros <- sum(prob_x == 0)

  if (row_zeros == 1) {
    return(1 - (col_zeros / nrow(prob_matrix)))
  }

  if (col_zeros == 1) {
    return(1 - (row_zeros / nrow(prob_matrix)))
  }

  hx <- entropy(prob_x)
  hy <- entropy(prob_y)

  return((hx + hy - entropy(prob_matrix)) / min(hx, hy))
}


#' @name Scott's pi
#' @title Scott's \eqn{\pi}
#' @description This method computes Scott's \eqn{\pi} of a confusion matrix.
#' @param conf_matrix The confusion matrix.
#' @return The Scott's \eqn{\pi} of `conf_matrix`.
#' @examples
#' # create a confusion matrix
#' conf_matrix <- matrix(1:9, nrow = 3, ncol = 3)
#'
#' # evaluate Scott's pi of `conf_matrix`
#' scott_pi(conf_matrix)

scott_pi <- function(conf_matrix) {
  is_a_confusion_matrix(conf_matrix)

  row_totals <- rowSums(conf_matrix)
  col_totals <- colSums(conf_matrix)
  total <- sum(conf_matrix)

  p_e <- sum(((row_totals + col_totals) / (2 * total))^2)
  p_o <- sum(diag(conf_matrix)) / total

  return((p_o - p_e) / (1 - p_e))
}
