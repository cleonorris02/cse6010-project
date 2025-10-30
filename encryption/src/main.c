#include "hotspot.h"

#include <ctype.h>
#include <errno.h>
#include <omp.h>
#include <sodium.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#define NONCE_SIZE crypto_stream_xchacha20_NONCEBYTES
#define KEY_SIZE crypto_stream_xchacha20_KEYBYTES

static int ensure_directory(const char *path) {
    if (!path) {
        return -1;
    }
#if defined(_WIN32)
    int rc = _mkdir(path);
    if (rc != 0 && errno != EEXIST) {
        return -1;
    }
#else
    if (mkdir(path, 0700) != 0) {
        if (errno != EEXIST) {
            return -1;
        }
    }
#endif
    return 0;
}

static int hex_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static int load_key_from_hex_file(const char *path, unsigned char key[KEY_SIZE]) {
    FILE *file = fopen(path, "r");
    if (!file) {
        perror("Failed to open key file");
        return -1;
    }
    char buffer[KEY_SIZE * 2 + 4];
    size_t read_bytes = fread(buffer, 1, sizeof(buffer) - 1, file);
    fclose(file);
    if (read_bytes < KEY_SIZE * 2) {
        fprintf(stderr, "Key file must contain at least %d hexadecimal characters\n", KEY_SIZE * 2);
        return -1;
    }
    size_t hex_index = 0;
    size_t key_index = 0;
    while (key_index < KEY_SIZE && hex_index + 1 < read_bytes) {
        while (hex_index < read_bytes && isspace((unsigned char)buffer[hex_index])) {
            hex_index++;
        }
        if (hex_index + 1 >= read_bytes) {
            break;
        }
        int high = hex_value(buffer[hex_index]);
        int low = hex_value(buffer[hex_index + 1]);
        if (high < 0 || low < 0) {
            fprintf(stderr, "Invalid hex character in key file\n");
            return -1;
        }
        key[key_index++] = (unsigned char)((high << 4) | low);
        hex_index += 2;
    }
    if (key_index != KEY_SIZE) {
        fprintf(stderr, "Key file does not contain enough hex characters for a %d-byte key\n", KEY_SIZE);
        return -1;
    }
    return 0;
}

static char *build_positions_string(const HotspotRecord *record) {
    if (!record || record->position_count == 0) {
        return NULL;
    }
    size_t buffer_len = record->position_count * 21 + 1;
    char *buffer = malloc(buffer_len);
    if (!buffer) {
        return NULL;
    }
    size_t offset = 0;
    for (size_t i = 0; i < record->position_count; ++i) {
        int written = snprintf(buffer + offset, buffer_len - offset, "%zu", record->positions[i]);
        if (written < 0 || (size_t)written >= buffer_len - offset) {
            free(buffer);
            return NULL;
        }
        offset += (size_t)written;
        if (i + 1 < record->position_count) {
            if (offset + 1 >= buffer_len) {
                free(buffer);
                return NULL;
            }
            buffer[offset++] = ',';
            buffer[offset] = '\0';
        }
    }
    return buffer;
}

static char *build_plaintext_block(const HotspotRecord *record) {
    char *positions = build_positions_string(record);
    if (!positions) {
        return NULL;
    }
    size_t plaintext_len = snprintf(NULL, 0, "Hotspot Positions: %s\nReference: %s\n", positions, record->reference ? record->reference : "");
    char *plaintext = malloc(plaintext_len + 1);
    if (!plaintext) {
        free(positions);
        return NULL;
    }
    snprintf(plaintext, plaintext_len + 1, "Hotspot Positions: %s\nReference: %s\n", positions, record->reference ? record->reference : "");
    free(positions);
    return plaintext;
}

static int write_output_file(const char *directory, size_t index, const unsigned char nonce[NONCE_SIZE], const unsigned char *ciphertext, size_t ciphertext_len, const HotspotRecord *record) {
    char path[512];
    snprintf(path, sizeof(path), "%s/hotspot_%zu.bin", directory, index);
    FILE *bin = fopen(path, "wb");
    if (!bin) {
        perror("Failed to open ciphertext output file");
        return -1;
    }
    if (fwrite(nonce, 1, NONCE_SIZE, bin) != NONCE_SIZE) {
        perror("Failed to write nonce");
        fclose(bin);
        return -1;
    }
    if (ciphertext_len > 0) {
        if (fwrite(ciphertext, 1, ciphertext_len, bin) != ciphertext_len) {
            perror("Failed to write ciphertext");
            fclose(bin);
            return -1;
        }
    }
    fclose(bin);

    snprintf(path, sizeof(path), "%s/hotspot_%zu.meta", directory, index);
    FILE *meta = fopen(path, "w");
    if (!meta) {
        perror("Failed to open metadata file");
        return -1;
    }
    fprintf(meta, "Hotspot Count: %zu\n", record->position_count);
    fprintf(meta, "Reference: %s\n", record->reference ? record->reference : "");
    if (record->alternate) {
        fprintf(meta, "Alternate: %s\n", record->alternate);
    }
    fprintf(meta, "Nonce (hex): ");
    for (size_t i = 0; i < NONCE_SIZE; ++i) {
        fprintf(meta, "%02x", nonce[i]);
    }
    fprintf(meta, "\nCiphertext Length: %zu\n", ciphertext_len);
    fclose(meta);
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <dog.txt> <hex-key-file> <output-directory>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_path = argv[1];
    const char *key_path = argv[2];
    const char *output_dir = argv[3];

    if (sodium_init() < 0) {
        fprintf(stderr, "Failed to initialise libsodium\n");
        return EXIT_FAILURE;
    }

    unsigned char key[KEY_SIZE];
    if (load_key_from_hex_file(key_path, key) != 0) {
        return EXIT_FAILURE;
    }

    if (ensure_directory(output_dir) != 0) {
        perror("Failed to create output directory");
        return EXIT_FAILURE;
    }

    HotspotCollection collection;
    if (hotspot_collection_init(&collection) != 0) {
        fprintf(stderr, "Failed to initialise collection\n");
        return EXIT_FAILURE;
    }
    if (parse_hotspot_file(input_path, &collection) != 0) {
        fprintf(stderr, "Failed to parse hotspot data\n");
        hotspot_collection_free(&collection);
        return EXIT_FAILURE;
    }

    if (collection.count == 0) {
        fprintf(stderr, "No hotspot records found\n");
        hotspot_collection_free(&collection);
        return EXIT_FAILURE;
    }

    int failure = 0;

#pragma omp parallel for num_threads(7) schedule(dynamic)
    for (size_t i = 0; i < collection.count; ++i) {
        unsigned char nonce[NONCE_SIZE];
        randombytes_buf(nonce, sizeof(nonce));

        char *plaintext = build_plaintext_block(&collection.records[i]);
        if (!plaintext) {
#pragma omp critical
            {
                fprintf(stderr, "Failed to build plaintext for record %zu\n", i);
                failure = 1;
            }
            continue;
        }
        size_t plaintext_len = strlen(plaintext);
        unsigned char *ciphertext = malloc(plaintext_len);
        if (!ciphertext && plaintext_len > 0) {
#pragma omp critical
            {
                fprintf(stderr, "Failed to allocate ciphertext buffer for record %zu\n", i);
                failure = 1;
            }
            free(plaintext);
            continue;
        }

        if (crypto_stream_xchacha20_xor(ciphertext, (const unsigned char *)plaintext, plaintext_len, nonce, key) != 0) {
#pragma omp critical
            {
                fprintf(stderr, "Encryption failed for record %zu\n", i);
                failure = 1;
            }
            free(ciphertext);
            free(plaintext);
            continue;
        }

        if (write_output_file(output_dir, i, nonce, ciphertext, plaintext_len, &collection.records[i]) != 0) {
#pragma omp critical
            {
                fprintf(stderr, "Failed to write output for record %zu\n", i);
                failure = 1;
            }
        }

        free(ciphertext);
        free(plaintext);
    }

    hotspot_collection_free(&collection);
    if (failure) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

