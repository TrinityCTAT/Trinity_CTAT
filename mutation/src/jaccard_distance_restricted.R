dist_matrix = function(x, by.row = FALSE) {
    
    # Create a matrix and transpose if by row
    x = as.matrix(x)
    if (by.row == FALSE) {
        x = t(x)
    }
    
    # Make distance matrix
    m = matrix(nrow=nrow(x), ncol=nrow(x))
    diag(m) = 0
    colnames(m) = rownames(x)
    rownames(m) = rownames(x)
    
    # Look and see if rna is in the file name
    vf_row_is_rna = grepl( "rna", rownames(m))
    
    # Traverse rows
    for (i in 1:(nrow(x)-1)) {

        message("row (",i,")")
        
        # Traverse rows excluding rows already done
        for ( j in ( i + 1 ):nrow( x ) ) {
            # Check to see if one of the files is DNA and one is RNA
            if( !xor( vf_row_is_rna[i], vf_row_is_rna[j] ) ){
              m[i,j] = 1
              m[j,i] = 1
            } else {
              # How many are the same excluding NA
              num_same = sum(na.omit(x[i,] == x[j,]))
              # How many are different exlcuding NA
              num_diff = sum(na.omit(x[i,] != x[j,]))
              # Total calls where there were no NA
              total = num_same + num_diff
              # Can not divide by zero so handle special case
              jaccard = ifelse (total > 0, 1 - (num_same / (num_same + num_diff)), 1)
              # Store as square matrix
              m[i,j] = jaccard
              m[j,i] = jaccard
            }
        }
    }
    return(m)
}
