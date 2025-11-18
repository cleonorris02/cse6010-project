/*         
 * Outline and logic generated with ChatGPT (OpenAI), Oct 2025.         
 * Reviewed and modified by Cleo Norris.         
*/

/*
 * Encode the encrypted bitstream into SNP positions within a DNA sequence, and make a mutation, pass it on to error-detection module.
*/ 

#include <ctype.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "embedding.h"

/*
 * SNP embedding utilities for hiding encrypted payloads within genomic
 * sequences. This module mirrors the behavior of the previous Python
 * implementation while providing a C interface suitable for integration with
 * lower-level systems.
 */

/* -- Internal helpers ----------------------------------------------------- */

static int set_error(char **out_error, const char *message);
static int base_index(char base);
static int select_allele(char reference,
                         int bit,
                         const CandidateSNP *candidate,
                         char *out_allele,
                         char **out_error);

/* Deterministic mapping used when alternates are not provided. */
static const char DEFAULT_ALLELE_MAP[4][2] = {
    /* A */ {'C', 'G'},
    /* C */ {'A', 'T'},
    /* G */ {'A', 'T'},
    /* T */ {'C', 'G'}
};

size_t calculate_capacity(const CandidateSNP *candidates, size_t num_candidates) {
    (void)candidates; /* Unused parameter retained for parity with Python API. */
    return num_candidates;
}

int embed_bitstream(const char *sequence,
                    const CandidateSNP *candidates,
                    size_t num_candidates,
                    const uint8_t *payload,
                    size_t payload_len,
                    EmbeddingResult *out_result,
                    char **out_error) {
    size_t seq_len;
    size_t bit_count;
    char *mutated;
    EmbeddedAllele *alleles;
    size_t i;

    if (out_error) {
        *out_error = NULL;
    }

    if (!sequence || !candidates || !out_result) {
        return set_error(out_error, "Invalid argument: NULL pointer supplied.");
    }

    if (payload_len > 0 && !payload) {
        return set_error(out_error, "Invalid argument: payload data is NULL.");
    }

    seq_len = strlen(sequence);
    bit_count = payload_len * 8;

    if (bit_count > num_candidates) {
        return set_error(out_error,
                         "Insufficient candidate SNPs for payload capacity.");
    }

    mutated = (char *)malloc(seq_len + 1);
    if (!mutated) {
        return set_error(out_error, "Failed to allocate memory for sequence.");
    }
    memcpy(mutated, sequence, seq_len + 1);

    if (bit_count > 0) {
        alleles = (EmbeddedAllele *)calloc(bit_count, sizeof(EmbeddedAllele));
        if (!alleles) {
            free(mutated);
            return set_error(out_error, "Failed to allocate memory for alleles.");
        }
    } else {
        alleles = NULL;
    }

    for (i = 0; i < bit_count; ++i) {
        const CandidateSNP *candidate = &candidates[i];
        size_t pos = candidate->position;
        char expected = (char)toupper((unsigned char)candidate->reference);
        char current;
        char allele;
        size_t byte_index = i / 8;
        size_t bit_offset = 7 - (i % 8);
        int bit = (payload[byte_index] >> bit_offset) & 1;

        if (pos >= seq_len) {
            set_error(out_error, "Candidate SNP position outside sequence bounds.");
            free(mutated);
            free(alleles);
            return -1;
        }

        current = (char)toupper((unsigned char)mutated[pos]);
        if (current != expected) {
            set_error(out_error,
                      "Reference base does not match sequence at candidate position.");
            free(mutated);
            free(alleles);
            return -1;
        }

        if (select_allele(expected, bit, candidate, &allele, out_error) != 0) {
            free(mutated);
            free(alleles);
            return -1;
        }

        mutated[pos] = allele;
        alleles[i].position = pos;
        alleles[i].reference = expected;
        alleles[i].allele = allele;
        alleles[i].bit = bit;
    }

    out_result->sequence = mutated;
    out_result->alleles = alleles;
    out_result->num_alleles = bit_count;
    return 0;
}

void free_embedding_result(EmbeddingResult *result) {
    if (!result) {
        return;
    }
    free(result->sequence);
    free(result->alleles);
    result->sequence = NULL;
    result->alleles = NULL;
    result->num_alleles = 0;
}

static int set_error(char **out_error, const char *message) {
    if (out_error) {
        size_t len = strlen(message);
        char *copy = (char *)malloc(len + 1);
        if (!copy) {
            *out_error = NULL;
            return -1;
        }
        memcpy(copy, message, len + 1);
        *out_error = copy;
    }
    return -1;
}

static int base_index(char base) {
    switch (base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

static int select_allele(char reference,
                         int bit,
                         const CandidateSNP *candidate,
                         char *out_allele,
                         char **out_error) {
    char normalized[4];
    size_t normalized_count = 0;
    size_t i;
    int ref_index;

    ref_index = base_index(reference);
    if (ref_index < 0) {
        return set_error(out_error, "Unsupported reference nucleotide.");
    }

    if (candidate && candidate->alternates && candidate->num_alternates > 0) {
        for (i = 0; i < candidate->num_alternates && normalized_count < 4; ++i) {
            char alt = (char)toupper((unsigned char)candidate->alternates[i]);
            size_t j;
            int is_duplicate = 0;

            if (alt == reference || base_index(alt) < 0) {
                continue;
            }

            for (j = 0; j < normalized_count; ++j) {
                if (normalized[j] == alt) {
                    is_duplicate = 1;
                    break;
                }
            }
            if (!is_duplicate) {
                normalized[normalized_count++] = alt;
            }
        }
    }

    if (normalized_count >= 2) {
        *out_allele = normalized[bit & 1];
        return 0;
    }

    if (normalized_count == 1) {
        if (bit == 0) {
            *out_allele = normalized[0];
            return 0;
        }
        /* Need an alternate distinct from both reference and provided allele. */
        if (DEFAULT_ALLELE_MAP[ref_index][0] != normalized[0]) {
            *out_allele = DEFAULT_ALLELE_MAP[ref_index][0];
            return 0;
        }
        if (DEFAULT_ALLELE_MAP[ref_index][1] != normalized[0]) {
            *out_allele = DEFAULT_ALLELE_MAP[ref_index][1];
            return 0;
        }
        return set_error(out_error, "Unable to determine fallback allele.");
    }

    *out_allele = DEFAULT_ALLELE_MAP[ref_index][bit & 1];
    return 0;
}

