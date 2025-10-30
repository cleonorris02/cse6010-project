#ifndef ERROR_DETECTION_H
#define ERROR_DETECTION_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Mapping used throughout the module */
#define DNA_BASE_A 'A'
#define DNA_BASE_T 'T'
#define DNA_BASE_G 'G'
#define DNA_BASE_C 'C'

/* Error codes returned by detect_and_correct_parity_block */
#define PARITY_OK 0
#define PARITY_CORRECTED 1
#define PARITY_UNRECOVERABLE -1
#define PARITY_INVALID_INPUT -2

/*
 * Build a parity-augmented block from a set of fixed-length DNA sequences.
 *
 * @param rows           Array of pointers to null-terminated DNA strings.
 * @param row_count      Number of sequences in the array.
 * @param out_rows       Output parameter receiving the total number of rows in
 *                       the augmented block (includes the parity row).
 * @param out_cols       Output parameter receiving the total number of columns
 *                       in the augmented block (includes the parity column).
 * @return               Newly allocated matrix of characters with dimensions
 *                       (*out_rows) x (*out_cols). The caller owns the matrix
 *                       and must free it with free_parity_block(). Returns
 *                       NULL if allocation fails or input is invalid.
 */
char **build_parity_block(const char **rows,
                          size_t row_count,
                          size_t *out_rows,
                          size_t *out_cols);

/*
 * Free a matrix allocated by build_parity_block.
 */
void free_parity_block(char **block, size_t total_rows);

/*
 * Detect and correct a single mutation in a parity-augmented block.
 *
 * The function checks row and column parity as described in the block sum
 * method and attempts to correct a single mutated nucleotide (data or parity).
 *
 * @param block          Matrix previously produced by build_parity_block.
 * @param total_rows     Number of rows in the matrix (data rows + parity row).
 * @param total_cols     Number of columns in the matrix (data columns + parity column).
 * @param corrected_row  Optional pointer that receives the row index of a corrected
 *                       nucleotide. Ignored when NULL.
 * @param corrected_col  Optional pointer that receives the column index of a corrected
 *                       nucleotide. Ignored when NULL.
 * @return               PARITY_OK if no errors were detected, PARITY_CORRECTED
 *                       if a single error was detected and corrected, or
 *                       PARITY_UNRECOVERABLE when the corruption cannot be
 *                       corrected (e.g., multiple errors). PARITY_INVALID_INPUT
 *                       indicates malformed data.
 */
int detect_and_correct_parity_block(char **block,
                                    size_t total_rows,
                                    size_t total_cols,
                                    size_t *corrected_row,
                                    size_t *corrected_col);

/*
 * Utility helper that prints a parity block to stdout.
 */
void print_parity_block(char **block, size_t total_rows, size_t total_cols);

#ifdef __cplusplus
}
#endif

#endif /* ERROR_DETECTION_H */
