# Package Name: JOPT

JOPT is an R package designed to implement J-optimal subsample selection for regression models. The package provides efficient tools for selecting subsets of data to optimize statistical efficiency in modeling. It includes functions for model creation, subsample selection, and running comparison examples. The methodology is based on Cia-Mina et al. (2025).

## Installation

### From Source

To install the package from source, first clone or download the repository, then use the following commands in R:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install the package
devtools::install("path/to/JOPT")
```

## Usage

### 1. Creating a Model

#### Definition
The function `create_model_function()` is used to define a regression model for subsample selection. This function takes user-defined mathematical expressions and constructs a function that evaluates them based on an input vector `x`.

#### Inputs
- `expressions`: A character vector of mathematical expressions to define the model. Each expression should be valid R code and reference elements of `x` (e.g., "x[1]", "x[2]^2").

#### Outputs
- A function that takes an input vector `x` and evaluates the expressions in `expressions`, returning a column matrix of the results.

#### Example Usage
```r
# Define the model expressions
expressions <- c("1", "x[1]", "x[1]*x[2]^2")

# Create the model function
model_function <- create_model_function(expressions)

# Test the model function with an input vector
input_vector <- c(2, 3) # x[1] = 2, x[2] = 3
result <- model_function(input_vector)
print(result)
```

### 2. Selecting J-Optimal Subsamples

#### Definition
The function `jseq()` implements the J-optimal subsample selection algorithm. It selects subsets of data that maximize model efficiency based on statistical criteria.

#### Inputs
- `x`: A dataset (data frame) containing the covariates for the regression model.
- `alpha`: A numeric value between 0 and 1 specifying the subsample proportion.
- `model_vec`: A character vector defining the regression model. Each element should represent a term in the model.
- `k0`: Initial size of the subsample (default: `5 * length(model_vec)`).
- `q`: A numeric value between `0.5` and `1` (default: `5/8`).
- `gamma`: A numeric value between `0` and `q-0.5` (default: `1/10`).
- `eps1`: A small positive value (default: `0`).

#### Outputs
- A list with the following components:
  - `x_j`: A subsample of `x` containing the selected observations.
  - `idx`: A vector of indices corresponding to the selected rows of `x`.

#### Example Usage
```r
# Example 1: Bivariate regression
set.seed(123)
x1 <- runif(1e3, min = -1, max = 1)
x2 <- runif(1e3, min = -1, max = 1)
x <- data.frame(x1 = x1, x2 = x2)
model_vec <- c("1", "x[1]", "x[2]", "x[1]*x[2]", "x[1]^2", "x[2]^2")
result <- jseq(x, 0.3, model_vec)

# Plot the full dataset and the selected subsample
plot(x$x1, x$x2, col = "black", pch = 16, cex = 0.7, xlab = "x1", ylab = "x2")
points(result$x_j$x1, result$x_j$x2, col = "red", pch = 16, cex = 0.7)
title(main = "J-OPT", line = 1)
```

### 3. Running Efficiency Comparisons

#### Definition
The function `run_efficiency_comparison_example()` provides an example of how to compare different subsample selection strategies using the package.

#### Inputs
- None.

#### Outputs
- An output demonstrating efficiency comparisons of different subsampling methods.

#### Example Usage
```r
# To run the example:
run_efficiency_comparison_example()
```

## References

Cia-Mina et al. (2025). *J-Optimal Subsample Selection in Regression Models: A Computational Approach*.

## License

JOPT is licensed under the MIT License.

## Contact

For issues or feature requests, please submit them via the package's GitHub repository.

