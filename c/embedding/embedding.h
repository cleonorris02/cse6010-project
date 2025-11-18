/*         
 * Outline and logic generated with ChatGPT (OpenAI), Nov 2025.         
 * Reviewed and modified by Cleo Norris.         
*/

#ifndef EMBEDDING_H
#define EMBEDDING_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Description of a candidate SNP position within a genomic sequence. */
typedef struct {
    size_t position;        /* Zero-based index in the sequence. */
    char reference;         /* Expected reference nucleotide. */
    const char *alternates; /* Optional array of alternate alleles. */
    size_t num_alternates;  /* Number of entries in `alternates`. */
} CandidateSNP;

/* Details about an allele chosen to embed a bit. */
typedef struct {
    size_t position; /* SNP position within the sequence. */
    char reference;  /* Reference nucleotide at that position. */
    char allele;     /* Allele substituted to encode the bit. */
    int bit;         /* Encoded bit value (0 or 1). */
} EmbeddedAllele;

/* Aggregated result from embedding a payload in a sequence. */
typedef struct {
    char *sequence;             /* Mutated sequence with embedded payload. */
    EmbeddedAllele *alleles;    /* Array describing each embedded SNP. */
    size_t num_alleles;         /* Number of encoded SNPs. */
} EmbeddingResult;

/* Returns the number of SNPs available for embedding. */
size_t calculate_capacity(const CandidateSNP *candidates, size_t num_candidates);

/*
 * Embeds a payload bitstream into the provided sequence using the supplied
 * candidate SNPs. Returns 0 on success or -1 on failure. The caller owns the
 * memory in `out_result` and should release it using `free_embedding_result`.
 */
int embed_bitstream(const char *sequence,
                    const CandidateSNP *candidates,
                    size_t num_candidates,
                    const uint8_t *payload,
                    size_t payload_len,
                    EmbeddingResult *out_result,
                    char **out_error);

/* Releases memory allocated inside an EmbeddingResult. */
void free_embedding_result(EmbeddingResult *result);

#ifdef __cplusplus
}
#endif

#endif /* EMBEDDING_H */

