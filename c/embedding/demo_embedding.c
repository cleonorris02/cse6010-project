/*         
 * Outline and logic generated with ChatGPT (OpenAI), Nov 2025.         
 * Reviewed and modified by Cleo Norris.         
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "embedding.h"

int main(void) {
    const char *sequence = "ACGTACGTACGT";
    CandidateSNP candidates[] = {
        {0, 'A', NULL, 0},
        {1, 'C', NULL, 0},
        {2, 'G', NULL, 0},
        {3, 'T', NULL, 0},
        {4, 'A', NULL, 0},
        {5, 'C', NULL, 0},
        {6, 'G', NULL, 0},
        {7, 'T', NULL, 0}
    };
    const uint8_t payload[] = {0xB6};
    EmbeddingResult result;
    char *error = NULL;
    size_t i;

    if (embed_bitstream(sequence,
                        candidates,
                        sizeof(candidates) / sizeof(candidates[0]),
                        payload,
                        sizeof(payload),
                        &result,
                        &error) != 0) {
        fprintf(stderr, "Embedding failed: %s\n", error ? error : "unknown error");
        free(error);
        return EXIT_FAILURE;
    }

    printf("Original sequence: %s\n", sequence);
    printf("Embedded sequence: %s\n", result.sequence);
    printf("Encoded alleles:%s\n", result.num_alleles ? "" : " (none)");
    for (i = 0; i < result.num_alleles; ++i) {
        const EmbeddedAllele *allele = &result.alleles[i];
        printf("  pos=%zu ref=%c allele=%c bit=%d\n",
               allele->position,
               allele->reference,
               allele->allele,
               allele->bit);
    }

    free_embedding_result(&result);
    return EXIT_SUCCESS;
}

