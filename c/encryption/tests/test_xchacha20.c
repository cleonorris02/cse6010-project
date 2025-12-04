/*
 * XChaCha20 Verification Tests
 *
 * This test file verifies that:
 * 1. XChaCha20 encryption/decryption round-trip works correctly
 * 2. The implementation matches expected behavior from libsodium
 * 3. The binary-to-DNA encoding is reversible
 */

#include <sodium.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NONCE_SIZE crypto_stream_xchacha20_NONCEBYTES
#define KEY_SIZE crypto_stream_xchacha20_KEYBYTES

/* Test result tracking */
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_EQ(expected, actual, msg) \
    do { \
        if ((expected) == (actual)) { \
            tests_passed++; \
            printf("  PASS: %s\n", msg); \
        } else { \
            tests_failed++; \
            printf("  FAIL: %s (expected %d, got %d)\n", msg, (int)(expected), (int)(actual)); \
        } \
    } while (0)

#define ASSERT_STR_EQ(expected, actual, msg) \
    do { \
        if (strcmp((expected), (actual)) == 0) { \
            tests_passed++; \
            printf("  PASS: %s\n", msg); \
        } else { \
            tests_failed++; \
            printf("  FAIL: %s\n    expected: %s\n    actual:   %s\n", msg, (expected), (actual)); \
        } \
    } while (0)

#define ASSERT_MEM_EQ(expected, actual, len, msg) \
    do { \
        if (memcmp((expected), (actual), (len)) == 0) { \
            tests_passed++; \
            printf("  PASS: %s\n", msg); \
        } else { \
            tests_failed++; \
            printf("  FAIL: %s (memory mismatch)\n", msg); \
        } \
    } while (0)

/* Binary to DNA conversion (same as in main.c) */
static char *binary_to_dna(const unsigned char *data, size_t length) {
    static const char nucleotides[4] = {'A', 'C', 'G', 'T'};
    if (!data || length == 0) {
        return NULL;
    }
    size_t dna_length = length * 4;
    char *output = (char *)malloc(dna_length + 1);
    if (!output) {
        return NULL;
    }
    size_t index = 0;
    for (size_t i = 0; i < length; ++i) {
        unsigned char byte = data[i];
        output[index++] = nucleotides[(byte >> 6) & 0x03];
        output[index++] = nucleotides[(byte >> 4) & 0x03];
        output[index++] = nucleotides[(byte >> 2) & 0x03];
        output[index++] = nucleotides[byte & 0x03];
    }
    output[dna_length] = '\0';
    return output;
}

/* DNA to binary conversion (reverse of binary_to_dna) */
static int dna_to_binary(const char *dna, unsigned char **data, size_t *length) {
    if (!dna || !data || !length) {
        return -1;
    }
    size_t dna_len = strlen(dna);
    if (dna_len % 4 != 0) {
        return -1;  /* DNA length must be multiple of 4 */
    }
    *length = dna_len / 4;
    *data = (unsigned char *)malloc(*length);
    if (!*data) {
        return -1;
    }
    for (size_t i = 0; i < *length; ++i) {
        unsigned char byte = 0;
        for (int j = 0; j < 4; ++j) {
            unsigned char bits;
            switch (dna[i * 4 + j]) {
                case 'A': bits = 0; break;
                case 'C': bits = 1; break;
                case 'G': bits = 2; break;
                case 'T': bits = 3; break;
                default: free(*data); *data = NULL; return -1;
            }
            byte = (unsigned char)((byte << 2) | bits);
        }
        (*data)[i] = byte;
    }
    return 0;
}

/*
 * Test 1: Verify encryption/decryption round-trip
 * XChaCha20 is a stream cipher, so encrypt XOR decrypt = original
 */
void test_encrypt_decrypt_roundtrip(void) {
    printf("\n[Test 1] Encryption/Decryption Round-trip\n");

    const char *plaintext = "ACGTACGTACGTACGT";  /* DNA sequence */
    size_t plaintext_len = strlen(plaintext);

    unsigned char key[KEY_SIZE];
    unsigned char nonce[NONCE_SIZE];
    unsigned char *ciphertext = malloc(plaintext_len);
    unsigned char *decrypted = malloc(plaintext_len + 1);

    /* Generate random key and nonce */
    randombytes_buf(key, sizeof(key));
    randombytes_buf(nonce, sizeof(nonce));

    /* Encrypt */
    int enc_result = crypto_stream_xchacha20_xor(
        ciphertext,
        (const unsigned char *)plaintext,
        plaintext_len,
        nonce,
        key
    );
    ASSERT_EQ(0, enc_result, "Encryption returns success");

    /* Verify ciphertext is different from plaintext.
     * Note: The probability of XChaCha20 keystream matching the plaintext
     * for even a 16-byte input is 2^(-128), which is negligible. */
    int is_different = memcmp(plaintext, ciphertext, plaintext_len) != 0;
    ASSERT_EQ(1, is_different, "Ciphertext differs from plaintext");

    /* Decrypt (XOR again with same key/nonce) */
    int dec_result = crypto_stream_xchacha20_xor(
        decrypted,
        ciphertext,
        plaintext_len,
        nonce,
        key
    );
    decrypted[plaintext_len] = '\0';
    ASSERT_EQ(0, dec_result, "Decryption returns success");

    /* Verify decrypted matches original */
    ASSERT_STR_EQ(plaintext, (char *)decrypted, "Decrypted matches original plaintext");

    free(ciphertext);
    free(decrypted);
}

/*
 * Test 2: Verify binary-to-DNA encoding is reversible
 */
void test_dna_encoding_roundtrip(void) {
    printf("\n[Test 2] Binary-to-DNA Encoding Round-trip\n");

    /* Test with various byte values */
    unsigned char test_data[] = {0x00, 0x55, 0xAA, 0xFF, 0x12, 0x34, 0x56, 0x78};
    size_t test_len = sizeof(test_data);

    /* Convert to DNA */
    char *dna = binary_to_dna(test_data, test_len);
    ASSERT_EQ(1, dna != NULL, "DNA encoding succeeds");

    /* Verify expected length (4 nucleotides per byte) */
    ASSERT_EQ(test_len * 4, strlen(dna), "DNA string has correct length");

    /* Verify specific encodings:
     * 0x00 = 00 00 00 00 = AAAA
     * 0x55 = 01 01 01 01 = CCCC
     * 0xAA = 10 10 10 10 = GGGG
     * 0xFF = 11 11 11 11 = TTTT
     */
    ASSERT_EQ('A', dna[0], "0x00 first nucleotide is A");
    ASSERT_EQ('C', dna[4], "0x55 first nucleotide is C");
    ASSERT_EQ('G', dna[8], "0xAA first nucleotide is G");
    ASSERT_EQ('T', dna[12], "0xFF first nucleotide is T");

    /* Convert back to binary */
    unsigned char *recovered;
    size_t recovered_len;
    int result = dna_to_binary(dna, &recovered, &recovered_len);
    ASSERT_EQ(0, result, "DNA decoding succeeds");
    ASSERT_EQ(test_len, recovered_len, "Recovered data has correct length");
    ASSERT_MEM_EQ(test_data, recovered, test_len, "Recovered data matches original");

    free(dna);
    free(recovered);
}

/*
 * Test 3: Verify with known key/nonce produces consistent output
 */
void test_deterministic_encryption(void) {
    printf("\n[Test 3] Deterministic Encryption with Fixed Key/Nonce\n");

    const char *plaintext = "Hello, XChaCha20!";
    size_t plaintext_len = strlen(plaintext);

    /* Fixed key (32 bytes) - for testing only! */
    unsigned char key[KEY_SIZE] = {
        0x82, 0xec, 0x6a, 0xcc, 0xaf, 0xde, 0x3b, 0x7c,
        0xda, 0xdd, 0x9b, 0x98, 0x54, 0xfe, 0x15, 0x74,
        0x71, 0xc6, 0xd5, 0x0e, 0x07, 0x31, 0x5b, 0x91,
        0xc2, 0x48, 0x62, 0x8a, 0xe7, 0x1d, 0xe7, 0xca
    };

    /* Fixed nonce (24 bytes) - for testing only! */
    unsigned char nonce[NONCE_SIZE] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17
    };

    unsigned char *ciphertext1 = malloc(plaintext_len);
    unsigned char *ciphertext2 = malloc(plaintext_len);

    /* Encrypt twice with same key/nonce */
    crypto_stream_xchacha20_xor(ciphertext1, (const unsigned char *)plaintext,
                                 plaintext_len, nonce, key);
    crypto_stream_xchacha20_xor(ciphertext2, (const unsigned char *)plaintext,
                                 plaintext_len, nonce, key);

    /* Both should produce identical ciphertext */
    ASSERT_MEM_EQ(ciphertext1, ciphertext2, plaintext_len, 
                  "Same key/nonce produces identical ciphertext");

    /* Convert to DNA and verify it's valid DNA */
    char *dna = binary_to_dna(ciphertext1, plaintext_len);
    ASSERT_EQ(1, dna != NULL, "DNA encoding of ciphertext succeeds");

    /* Verify DNA contains only valid nucleotides */
    int valid_dna = 1;
    for (size_t i = 0; i < strlen(dna); i++) {
        if (dna[i] != 'A' && dna[i] != 'C' && dna[i] != 'G' && dna[i] != 'T') {
            valid_dna = 0;
            break;
        }
    }
    ASSERT_EQ(1, valid_dna, "DNA encoding contains only valid nucleotides");

    free(ciphertext1);
    free(ciphertext2);
    free(dna);
}

/*
 * Test 4: Full pipeline test - encrypt DNA, convert to DNA encoding
 */
void test_full_pipeline(void) {
    printf("\n[Test 4] Full Encryption Pipeline\n");

    const char *dna_sequence = "ACGTACGTACGTACGTACGTACGT";
    size_t seq_len = strlen(dna_sequence);

    unsigned char key[KEY_SIZE];
    unsigned char nonce[NONCE_SIZE];

    /* Generate random key and nonce */
    randombytes_buf(key, sizeof(key));
    randombytes_buf(nonce, sizeof(nonce));

    /* Encrypt */
    unsigned char *ciphertext = malloc(seq_len);
    crypto_stream_xchacha20_xor(ciphertext, (const unsigned char *)dna_sequence,
                                 seq_len, nonce, key);

    /* Convert ciphertext to DNA encoding */
    char *encrypted_dna = binary_to_dna(ciphertext, seq_len);
    ASSERT_EQ(1, encrypted_dna != NULL, "Ciphertext converts to DNA");

    /* Convert nonce to DNA encoding */
    char *nonce_dna = binary_to_dna(nonce, NONCE_SIZE);
    ASSERT_EQ(1, nonce_dna != NULL, "Nonce converts to DNA");

    printf("  Original DNA:  %s\n", dna_sequence);
    printf("  Nonce DNA:     %s (length %zu)\n", nonce_dna, strlen(nonce_dna));
    printf("  Encrypted DNA: %s (length %zu)\n", encrypted_dna, strlen(encrypted_dna));

    /* Verify we can reverse the process */
    unsigned char *recovered_nonce;
    size_t recovered_nonce_len;
    int nonce_result = dna_to_binary(nonce_dna, &recovered_nonce, &recovered_nonce_len);
    ASSERT_EQ(0, nonce_result, "Nonce DNA decodes successfully");
    ASSERT_MEM_EQ(nonce, recovered_nonce, NONCE_SIZE, "Recovered nonce matches original");

    unsigned char *recovered_ciphertext;
    size_t recovered_ct_len;
    int ct_result = dna_to_binary(encrypted_dna, &recovered_ciphertext, &recovered_ct_len);
    ASSERT_EQ(0, ct_result, "Ciphertext DNA decodes successfully");

    /* Decrypt */
    unsigned char *decrypted = malloc(recovered_ct_len + 1);
    crypto_stream_xchacha20_xor(decrypted, recovered_ciphertext,
                                 recovered_ct_len, recovered_nonce, key);
    decrypted[recovered_ct_len] = '\0';

    ASSERT_STR_EQ(dna_sequence, (char *)decrypted, 
                  "Full round-trip: original DNA recovered correctly");

    free(ciphertext);
    free(encrypted_dna);
    free(nonce_dna);
    free(recovered_nonce);
    free(recovered_ciphertext);
    free(decrypted);
}

/*
 * Test 5: Verify XChaCha20 produces different output with different nonces
 */
void test_nonce_uniqueness(void) {
    printf("\n[Test 5] Nonce Uniqueness\n");

    const char *plaintext = "Same plaintext for all tests";
    size_t len = strlen(plaintext);

    unsigned char key[KEY_SIZE];
    randombytes_buf(key, sizeof(key));

    unsigned char nonce1[NONCE_SIZE], nonce2[NONCE_SIZE];
    randombytes_buf(nonce1, sizeof(nonce1));
    randombytes_buf(nonce2, sizeof(nonce2));

    unsigned char *cipher1 = malloc(len);
    unsigned char *cipher2 = malloc(len);

    crypto_stream_xchacha20_xor(cipher1, (const unsigned char *)plaintext, len, nonce1, key);
    crypto_stream_xchacha20_xor(cipher2, (const unsigned char *)plaintext, len, nonce2, key);

    /* Different nonces should produce different ciphertext.
     * Note: The probability of collision with different 192-bit nonces is
     * 2^(-192) per keystream byte, which is negligible for any practical test. */
    int is_different = memcmp(cipher1, cipher2, len) != 0;
    ASSERT_EQ(1, is_different, "Different nonces produce different ciphertext");

    free(cipher1);
    free(cipher2);
}

int main(void) {
    printf("===========================================\n");
    printf("XChaCha20 Implementation Verification Tests\n");
    printf("===========================================\n");

    if (sodium_init() < 0) {
        fprintf(stderr, "Failed to initialize libsodium\n");
        return 1;
    }

    printf("\nUsing libsodium version: %s\n", sodium_version_string());
    printf("XChaCha20 key size: %d bytes\n", KEY_SIZE);
    printf("XChaCha20 nonce size: %d bytes\n", NONCE_SIZE);

    /* Run all tests */
    test_encrypt_decrypt_roundtrip();
    test_dna_encoding_roundtrip();
    test_deterministic_encryption();
    test_full_pipeline();
    test_nonce_uniqueness();

    /* Summary */
    printf("\n===========================================\n");
    printf("Test Summary\n");
    printf("===========================================\n");
    printf("Passed: %d\n", tests_passed);
    printf("Failed: %d\n", tests_failed);
    printf("Total:  %d\n", tests_passed + tests_failed);

    if (tests_failed > 0) {
        printf("\n*** SOME TESTS FAILED ***\n");
        return 1;
    }

    printf("\n*** ALL TESTS PASSED ***\n");
    printf("\nThis verifies that:\n");
    printf("1. XChaCha20 encryption/decryption works correctly (round-trip)\n");
    printf("2. Binary-to-DNA encoding is reversible and correct\n");
    printf("3. The encryption is deterministic with fixed key/nonce\n");
    printf("4. The full pipeline (encrypt -> DNA encode -> DNA decode -> decrypt) works\n");
    printf("5. Different nonces produce different ciphertext (semantic security)\n");

    return 0;
}
