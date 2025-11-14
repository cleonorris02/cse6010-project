 /*         
 * Outline and logic generated with ChatGPT (OpenAI), Oct 2025.         
 * Reviewed and modified by Cleo Norris.         
 * */

// header guards to prevent the header from being included twice
#ifndef ERROR_DETECTION_H
#define ERROR_DETECTION_H

// for standard definitions
#include <stddef.h>

// for when compiling with a C++ compiler
#ifdef __cplusplus
extern "C" {
#endif

// create macro constant for the character literal for a nucleotide base
#define DNA_BASE_A 'A'
#define DNA_BASE_T 'T'
#define DNA_BASE_G 'G'
#define DNA_BASE_C 'C'

// error codes returned by detect_and_correct_parity_block
#define PARITY_OK 0                 // no errors detected
#define PARITY_CORRECTED 1          // single error detected and corrected
#define PARITY_UNRECOVERABLE -1     // multiple errors detected, cannot correct
#define PARITY_INVALID_INPUT -2     // bad characters, null pointers, or memory issues

/*
Build augmented DNA matrix that includes parity rows and columns from set of 
fixed-length DNA sequences.
 
@param rows: Array of input DNA strings (each string = one row of data)
@param row_count: number of DNA strings in the input array
@param out_rows: Output parameter receiving the total number of rows in 
                    the augmented block (includes the parity row).
@param out_cols: Output parameter receiving the total number of columns in 
                    the augmented block (includes the parity column).
@return: Newly allocated matrix of characters with dimensions (*out_rows) x (*out_cols). 
            Returns NULL if allocation fails or input is invalid.
 */
char **build_parity_block(const char **rows,
                          size_t row_count,
                          size_t *out_rows,
                          size_t *out_cols);

// deallocates the memory returned by by build_parity_block
void free_parity_block(char **block, size_t total_rows);

/*
Detect and correct a single mutation in a parity-augmented block.
 
The function checks row and column parity as described in the block sum
method and attempts to correct a single mutated nucleotide (data or parity).

@param block: Matrix previously produced by build_parity_block.
@param total_rows: Number of rows in the matrix (data rows + parity row).
@param total_cols: Number of columns in the matrix (data columns + parity column).
@param corrected_row: Optional pointer that receives the row index of a corrected
                      nucleotide. Ignored when NULL.
@param corrected_col: Optional pointer that receives the column index of a corrected
                      nucleotide. Ignored when NULL.
@return: PARITY_OK if no errors were detected, PARITY_CORRECTED
                      if a single error was detected and corrected, or
                      PARITY_UNRECOVERABLE when the corruption cannot be
                      corrected (e.g., multiple errors). PARITY_INVALID_INPUT
                      indicates malformed data.
*/
int detect_and_correct_parity_block(char **block,
                                    size_t total_rows,
                                    size_t total_cols,
                                    size_t *corrected_row,
                                    size_t *corrected_col);

/*
Utility helper that prints a parity block to std output for inspection.
Each row is null-terminated ('\0' at the end), so it can be printed as a C string.
*/ 
void print_parity_block(char **block, size_t total_rows, size_t total_cols);

#ifdef __cplusplus
}
#endif

#endif /* ERROR_DETECTION_H */