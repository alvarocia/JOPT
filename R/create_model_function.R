#' Create a Model Function from Expressions
#'
#' This function takes a vector of mathematical expressions (as character strings)
#' and generates a function that, given an input vector `x`, computes the specified
#' expressions and returns the results as a column matrix.
#'
#' @param expressions A character vector of mathematical expressions to define the model.
#'   Each expression should be valid R code and reference elements of `x` (e.g., `"x[1]"`, `"x[2]^2"`).
#'
#' @return A function that takes an input vector `x` and evaluates the functions in
#'   `expressions`, returning a column matrix of the results.
#'
#' @examples
#' # Define the model expressions
#' expressions <- c("1", "x[1]", "x[1]*x[2]^2")
#'
#' # Create the model function
#' model_function <- create_model_function(expressions)
#'
#' # Test the model function with an input vector
#' input_vector <- c(2, 3) # x[1] = 2, x[2] = 3
#' result <- model_function(input_vector)
#' print(result)
#'
#' @export
create_model_function <- function(expressions) {
  # Validate the input
  if (!is.character(expressions)) {
    stop("The 'expressions' parameter must be a character vector.")
  }

  # Parse each string in 'expressions' into an R expression
  parsed_expressions <- lapply(expressions, function(expr) parse(text = expr)[[1]])

  # Create and return a function that evaluates the expressions for an input vector 'x'
  model_function <- function(x) {
    # Validate the input vector
    if (!is.numeric(x)) {
      stop("The input vector 'x' must be numeric.")
    }
    # Evaluate each expression in the context of 'x'
    results <- sapply(parsed_expressions, function(expr) eval(expr, envir = list(x = x)))
    # Convert the results to a column matrix
    return(matrix(results, ncol = 1))
  }

  return(model_function)
}
