#include "error_detection.h"

#include <stdio.h>
#include <string.h>

static void introduce_mutation(char **block, size_t row, size_t col, char new_base) {
    printf("Introducing mutation at (%zu, %zu): %c -> %c\n",
           row, col, block[row][col], new_base);
    block[row][col] = new_base;
}

static void run_check(char **block, size_t rows, size_t cols) {
    size_t corrected_row = 0;
    size_t corrected_col = 0;
    int status = detect_and_correct_parity_block(block, rows, cols,
                                                 &corrected_row, &corrected_col);
    if (status == PARITY_OK) {
        printf("Parity check: no errors detected.\n");
    } else if (status == PARITY_CORRECTED) {
        printf("Parity check: corrected nucleotide at (%zu, %zu).\n",
               corrected_row, corrected_col);
    } else if (status == PARITY_UNRECOVERABLE) {
        printf("Parity check: unrecoverable corruption detected.\n");
    } else {
        printf("Parity check: invalid input encountered.\n");
    }
}

int main(void) {
    const char *dna_rows[] = {
        "AACGGATGA",
        "TTAGGCATA",
        "CGTATTCGG"
    };
    size_t block_rows = 0;
    size_t block_cols = 0;

    char **block = build_parity_block(dna_rows, 3, &block_rows, &block_cols);
    if (!block) {
        fprintf(stderr, "Failed to construct parity block.\n");
        return 1;
    }

    printf("Initial block with parity nucleotides:\n");
    print_parity_block(block, block_rows, block_cols);
    putchar('\n');

    introduce_mutation(block, 0, 0, DNA_BASE_T);
    printf("Block after mutation:\n");
    print_parity_block(block, block_rows, block_cols);
    run_check(block, block_rows, block_cols);
    printf("Block after correction:\n");
    print_parity_block(block, block_rows, block_cols);
    putchar('\n');

    introduce_mutation(block, 1, block_cols - 1, DNA_BASE_A);
    printf("Block after parity mutation:\n");
    print_parity_block(block, block_rows, block_cols);
    run_check(block, block_rows, block_cols);
    printf("Block after correcting parity nucleotide:\n");
    print_parity_block(block, block_rows, block_cols);
    putchar('\n');

    introduce_mutation(block, block_rows - 1, 2, DNA_BASE_G);
    printf("Block after column parity mutation:\n");
    print_parity_block(block, block_rows, block_cols);
    run_check(block, block_rows, block_cols);
    printf("Restored block:\n");
    print_parity_block(block, block_rows, block_cols);

    free_parity_block(block, block_rows);
    return 0;
}
