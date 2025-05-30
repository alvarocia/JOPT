% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_model_function.R
\name{create_model_function}
\alias{create_model_function}
\title{Create a Model Function from Expressions}
\usage{
create_model_function(expressions)
}
\arguments{
\item{expressions}{A character vector of mathematical expressions to define the model.
Each expression should be valid R code and reference elements of \code{x} (e.g., \code{"x[1]"}, \code{"x[2]^2"}).}
}
\value{
A function that takes an input vector \code{x} and evaluates the functions in
\code{expressions}, returning a column matrix of the results.
}
\description{
This function takes a vector of mathematical expressions (as character strings)
and generates a function that, given an input vector \code{x}, computes the specified
expressions and returns the results as a column matrix.
}
\examples{
# Define the model expressions
expressions <- c("1", "x[1]", "x[1]*x[2]^2")

# Create the model function
model_function <- create_model_function(expressions)

# Test the model function with an input vector
input_vector <- c(2, 3) # x[1] = 2, x[2] = 3
result <- model_function(input_vector)
print(result)

}
