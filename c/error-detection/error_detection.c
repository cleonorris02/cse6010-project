 /*         
 * Outline and logic generated with ChatGPT (OpenAI), Oct 2025.         
 * Reviewed and modified by Cleo Norris.         
 * */

 /*
 This file builds a parity-augmented DNA matrix, validates it, 
 detects a single mutation, and corrects it if possible. 
 Parity is sum of digits mod 4 with A=0, T=1, G=2, C=3.
 */

#include "error_detection.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// mappings to move between bases and modular arithmetic
static int base_to_digit(char base) {
    switch (base) {
        case DNA_BASE_A:
            return 0;
        case DNA_BASE_T:
            return 1;
        case DNA_BASE_G:
            return 2;
        case DNA_BASE_C:
            return 3;
        default:
            return -1;
    }
}
static char digit_to_base(int digit) {
    static const char mapping[] = {DNA_BASE_A, DNA_BASE_T, DNA_BASE_G, DNA_BASE_C};
    if (digit < 0 || digit > 3) {
        return '?';
    }
    return mapping[digit];
}

// helper to free a partially allocated block in case of error
static void free_partial_block(char **block, size_t allocated_rows) {
    if (!block) {
        return;
    }
    for (size_t i = 0; i < allocated_rows; ++i) {
        free(block[i]);
    }
    free(block);
}

// constructs a parity-augmented block from input DNA rows
char **build_parity_block(const char **rows,
                          size_t row_count,
                          size_t *out_rows,
                          size_t *out_cols) {
    // validate inputs
    if (!rows || row_count == 0 || !out_rows || !out_cols) {
        return NULL;
    }

    size_t row_length = strlen(rows[0]);
    if (row_length == 0) {
        return NULL;
    }

    for (size_t i = 1; i < row_count; ++i) {
        if (!rows[i] || strlen(rows[i]) != row_length) {
            return NULL;
        }
    }
    
    // compute dimensions of augmented block
    size_t total_rows = row_count + 1;
    size_t total_cols = row_length + 1;

    // block is an array of total_rows pointers
    // each row allocates total_cols+1 chars
    // the extra +1 is to allow null-termination for printing
    char **block = calloc(total_rows, sizeof(char *));
    if (!block) {
        return NULL;
    }

    for (size_t i = 0; i < total_rows; ++i) {
        block[i] = calloc(total_cols + 1, sizeof(char));
        if (!block[i]) {
            free_partial_block(block, i);
            return NULL;
        }
    }

    size_t total_sum = 0;

    // copy data and compute row parities
    for (size_t i = 0; i < row_count; ++i) {
        size_t row_sum = 0;
        for (size_t j = 0; j < row_length; ++j) {
            if (!rows[i]) {
                free_partial_block(block, total_rows);
                return NULL;
            }
            // convert to uppercase and validate base
            char base = (char)toupper((unsigned char)rows[i][j]);
            if (base_to_digit(base) < 0) {
                free_partial_block(block, total_rows);
                return NULL;
            }
            block[i][j] = base;
            // accumulate row sum
            row_sum += (size_t)base_to_digit(base);
        }
        // accumulate total sum
        total_sum += row_sum;
        // set row parity nucleotide
        block[i][row_length] = digit_to_base((int)(row_sum % 4));
    }

    // compute column parities
    for (size_t j = 0; j < row_length; ++j) {
        size_t column_sum = 0;
        for (size_t i = 0; i < row_count; ++i) {
            // for each data column j, sum digits down rows i=0..row_count-1
            column_sum += (size_t)base_to_digit(block[i][j]);
        }
        block[row_count][j] = digit_to_base((int)(column_sum % 4));
    }

    block[row_count][row_length] = digit_to_base((int)(total_sum % 4));

    // set output dimensions
    *out_rows = total_rows;
    *out_cols = total_cols;
    return block;
}

// deallocates the memory returned by by build_parity_block
void free_parity_block(char **block, size_t total_rows) {
    if (!block) {
        return;
    }
    for (size_t i = 0; i < total_rows; ++i) {
        free(block[i]);
    }
    free(block);
}

// utility helper that prints a parity block to stdout
static int validate_block(char **block, size_t total_rows, size_t total_cols) {
    // check pointers and minimal sizes
    if (!block || total_rows < 2 || total_cols < 2) {
        return 0;
    }
    for (size_t i = 0; i < total_rows; ++i) {
        // ensure each row pointer is non-null
        if (!block[i]) {
            return 0;
        }
        // ensure rows are null-terminated for printing convenience
        block[i][total_cols] = '\0';
    }
    return 1;
}

// detects and corrects a single mutation in a parity-augmented block
int detect_and_correct_parity_block(char **block,
                                    size_t total_rows,
                                    size_t total_cols,
                                    size_t *corrected_row,
                                    size_t *corrected_col) {
    if (!validate_block(block, total_rows, total_cols)) {
        return PARITY_INVALID_INPUT;
    }

    // dimensions of data portion (excluding parity row/column)
    size_t data_rows = total_rows - 1;
    size_t data_cols = total_cols - 1;

    // arrays to hold sums, expected parities, and mismatches
    size_t *row_sums = calloc(total_rows, sizeof(size_t));
    size_t *col_sums = calloc(total_cols, sizeof(size_t));
    int *row_expected = calloc(total_rows, sizeof(int));
    int *col_expected = calloc(total_cols, sizeof(int));
    int *row_parity = calloc(total_rows, sizeof(int));
    int *col_parity = calloc(total_cols, sizeof(int));
    size_t *row_mismatches = calloc(total_rows, sizeof(size_t));
    size_t *col_mismatches = calloc(total_cols, sizeof(size_t));

    // check allocations
    if (!row_sums || !col_sums || !row_expected || !col_expected ||
        !row_parity || !col_parity || !row_mismatches || !col_mismatches) {
        free(row_sums);
        free(col_sums);
        free(row_expected);
        free(col_expected);
        free(row_parity);
        free(col_parity);
        free(row_mismatches);
        free(col_mismatches);
        return PARITY_INVALID_INPUT;
    }

    // counters for mismatches
    size_t row_mismatch_count = 0;
    size_t col_mismatch_count = 0;

    // for each data cell (i,j), add digit to row_sums[i] and col_sums[j]
    for (size_t i = 0; i < data_rows; ++i) {
        for (size_t j = 0; j < data_cols; ++j) {
            int digit = base_to_digit(block[i][j]);
            // if any base invalid, return PARITY_INVALID_INPUT
            if (digit < 0) {
                goto cleanup_invalid;
            }
            row_sums[i] += (size_t)digit;
            col_sums[j] += (size_t)digit;
        }
        // sums rows mod 4
        row_expected[i] = (int)(row_sums[i] % 4);
        // read stored row parity from block
        int stored_row_parity = base_to_digit(block[i][data_cols]);
        if (stored_row_parity < 0) {
            goto cleanup_invalid;
        }
        // stored parity digits read from block
        row_parity[i] = stored_row_parity;
        // accumulate for bottom-right parity
        col_sums[data_cols] += (size_t)stored_row_parity;
        // record mismatches
        if (stored_row_parity != row_expected[i]) {
            row_mismatches[row_mismatch_count++] = i;
        }
    }

    for (size_t j = 0; j < data_cols; ++j) {
        // sums columns mod 4
        col_expected[j] = (int)(col_sums[j] % 4);
        // read stored column parity from block
        int stored_col_parity = base_to_digit(block[data_rows][j]);
        if (stored_col_parity < 0) {
            goto cleanup_invalid;
        }
        // stored parity digits read from block
        col_parity[j] = stored_col_parity;
        // add stored column parity into row_sums[data_rows]
        row_sums[data_rows] += (size_t)stored_col_parity;
        // record mismatches
        if (stored_col_parity != col_expected[j]) {
            col_mismatches[col_mismatch_count++] = j;
        }
    }

    // compare bottom-right parity cell to column sums
    col_expected[data_cols] = (int)(col_sums[data_cols] % 4);
    int stored_bottom_right = base_to_digit(block[data_rows][data_cols]);
    if (stored_bottom_right < 0) {
        goto cleanup_invalid;
    }
    col_parity[data_cols] = stored_bottom_right;
    // record column mismatch
    if (stored_bottom_right != col_expected[data_cols]) {
        col_mismatches[col_mismatch_count++] = data_cols;
    }
    // compare bottom-right parity cell to row sums
    row_parity[data_rows] = stored_bottom_right;
    row_expected[data_rows] = (int)(row_sums[data_rows] % 4);
    // record row mismatch
    if (stored_bottom_right != row_expected[data_rows]) {
        row_mismatches[row_mismatch_count++] = data_rows;
    }

    // analyze mismatches
    if (row_mismatch_count == 0 && col_mismatch_count == 0) {
        free(row_sums);
        free(col_sums);
        free(row_expected);
        free(col_expected);
        free(row_parity);
        free(col_parity);
        free(row_mismatches);
        free(col_mismatches);
        return PARITY_OK;
    }

    // if more than one mismatch in rows or columns, cannot correct
    if (row_mismatch_count != 1 || col_mismatch_count != 1) {
        goto cleanup_unrecoverable;
    }

    size_t row_idx = row_mismatches[0];
    size_t col_idx = col_mismatches[0];

    // single mismatch found, identify what kind
    // data cell error
    if (row_idx < data_rows && col_idx < data_cols) {
        int current_digit = base_to_digit(block[row_idx][col_idx]);
        if (current_digit < 0) {
            goto cleanup_invalid;
        }
        // recompute what the digit must be to satisfy both row and column
        size_t row_without = row_sums[row_idx] - (size_t)current_digit;
        size_t col_without = col_sums[col_idx] - (size_t)current_digit;
        int target_row_parity = row_parity[row_idx];
        int target_col_parity = col_parity[col_idx];
        int needed_row = (int)((target_row_parity - (int)(row_without % 4) + 4) % 4);
        int needed_col = (int)((target_col_parity - (int)(col_without % 4) + 4) % 4);
        if (needed_row != needed_col) {
            // conflicting requirements, cannot correct
            goto cleanup_unrecoverable;
        }
        // write corrected base and return PARITY_CORRECTED
        block[row_idx][col_idx] = digit_to_base(needed_row);
        if (corrected_row) {
            *corrected_row = row_idx;
        }
        if (corrected_col) {
            *corrected_col = col_idx;
        }
        goto cleanup_corrected;
    // row parity cell error
    } else if (row_idx < data_rows && col_idx == data_cols) {
        // replace the row parity base with the expected value
        block[row_idx][data_cols] = digit_to_base(row_expected[row_idx]);
        // adjust bottom-right parity cell accordingly
        int new_bottom_right = (int)((col_sums[data_cols] - row_parity[row_idx] + row_expected[row_idx]) % 4);
        block[data_rows][data_cols] = digit_to_base(new_bottom_right);
        if (corrected_row) {
            *corrected_row = row_idx;
        }
        if (corrected_col) {
            *corrected_col = data_cols;
        }
        goto cleanup_corrected;
    // column parity cell error
    } else if (row_idx == data_rows && col_idx < data_cols) {
        // replace the column parity base with the expected value
        block[data_rows][col_idx] = digit_to_base(col_expected[col_idx]);
        // adjust bottom-right parity cell accordingly
        int new_bottom_right = (int)((row_sums[data_rows] - col_parity[col_idx] + col_expected[col_idx]) % 4);
        block[data_rows][data_cols] = digit_to_base(new_bottom_right);
        if (corrected_row) {
            *corrected_row = data_rows;
        }
        if (corrected_col) {
            *corrected_col = col_idx;
        }
        goto cleanup_corrected;
    // bottom-right parity cell error
    } else if (row_idx == data_rows && col_idx == data_cols) {
        // replace bottom-right parity base with expected value
        int corrected_digit = col_expected[data_cols];
        block[data_rows][data_cols] = digit_to_base(corrected_digit);
        if (corrected_row) {
            *corrected_row = data_rows;
        }
        if (corrected_col) {
            *corrected_col = data_cols;
        }
        goto cleanup_corrected;
    }

cleanup_unrecoverable:
    free(row_sums);
    free(col_sums);
    free(row_expected);
    free(col_expected);
    free(row_parity);
    free(col_parity);
    free(row_mismatches);
    free(col_mismatches);
    return PARITY_UNRECOVERABLE;

cleanup_invalid:
    free(row_sums);
    free(col_sums);
    free(row_expected);
    free(col_expected);
    free(row_parity);
    free(col_parity);
    free(row_mismatches);
    free(col_mismatches);
    return PARITY_INVALID_INPUT;

cleanup_corrected:
    free(row_sums);
    free(col_sums);
    free(row_expected);
    free(col_expected);
    free(row_parity);
    free(col_parity);
    free(row_mismatches);
    free(col_mismatches);
    return PARITY_CORRECTED;
}

void print_parity_block(char **block, size_t total_rows, size_t total_cols) {
    if (!block) {
        return;
    }
    // print each row as a null-terminated string
    for (size_t i = 0; i < total_rows; ++i) {
        if (!block[i]) {
            continue;
        }
        block[i][total_cols] = '\0';
        printf("%s\n", block[i]);
    }
}