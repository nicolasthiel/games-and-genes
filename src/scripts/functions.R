boolean_expression_matrix <- function(A_SD, A_SR, lower_percentile, upper_percentile) {
    # Builds a boolean matrix describing which gene is differentially expressed in which sample
    # based on percentiles.
    # INPUT:
    #   A_SD - Expression matrix of diseased samples (data frame)
    #   A_SR - Expression matrix of reference samples (data frame)
    #   lower_percentile - Lower percentile threshold
    #   upper_percentile - Upper percentile threshold

    # Ensure A_SD and A_SR are matrices for numerical operations
    A_SD <- as.matrix(A_SD)
    A_SR <- as.matrix(A_SR)

    # Calculate the lower and upper percentiles for each row (gene)
    p_lower <- apply(A_SR, 1, quantile, probs = lower_percentile / 100)
    p_upper <- apply(A_SR, 1, quantile, probs = upper_percentile / 100)

    # Create a boolean matrix (same size as A_SD) initialized to FALSE
    B <- matrix(FALSE, nrow = nrow(A_SD), ncol = ncol(A_SD))

    # Identify differentially expressed genes
    B <- (A_SD <= p_lower) | (A_SD >= p_upper)

    # Convert to a data frame to match input type
    return(as.data.frame(B))
}

support_of <- function(B) {
    # Finds the support of a boolean matrix B using row names.
    # INPUT:
    #   B - a boolean data frame
    # OUTPUT:
    #   sp_B - the support of B, a list where each element contains the row names of TRUE values in each column

    # Ensure B is a matrix for proper boolean indexing
    B <- as.matrix(B)

    # Get row names (if they exist)
    row_names <- rownames(B)

    # Apply function to each column: find row names where the value is TRUE
    sp_B <- lapply(seq_len(ncol(B)), function(col) {
        true_indices <- which(B[, col]) # Find indices where TRUE
        if (!is.null(row_names)) {
            return(row_names[true_indices]) # Return row names
        } else {
            return(true_indices) # Fall back to indices if no row names
        }
    })

    # Name the list elements with column names (optional)
    names(sp_B) <- colnames(B)

    return(sp_B)
}

find_coalitions <- function(sp_B) {
    # Finds the distinct sets in sp_B for all elements.
    # INPUT:
    #   sp_B - a list containing subsets (equivalent to MATLAB's cell array)
    # OUTPUT:
    #   coalitions - a list containing the distinct subsets / coalitions

    # Initialize coalitions list
    coalitions <- list()

    for (i in seq_along(sp_B)) {
        if (length(sp_B[[i]]) == 0) { # Skip if empty set
            next
        }

        to_append <- TRUE
        for (j in seq_along(coalitions)) {
            if (identical(coalitions[[j]], sp_B[[i]])) {
                to_append <- FALSE
                break
            }
        }

        if (to_append) {
            coalitions <- append(coalitions, list(sp_B[[i]]))
        }
    }

    return(coalitions)
}
