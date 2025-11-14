#include "sequence.h"

#include <errno.h>
#include <omp.h>
#include <sodium.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NONCE_SIZE crypto_stream_xchacha20_NONCEBYTES
#define KEY_SIZE crypto_stream_xchacha20_KEYBYTES

typedef struct {
    char *nonce_dna;
    char *ciphertext_dna;
    int status;
} EncryptionResult;

static void free_result(EncryptionResult *result) {
    if (!result) {
        return;
    }
    free(result->nonce_dna);
    free(result->ciphertext_dna);
    result->nonce_dna = NULL;
    result->ciphertext_dna = NULL;
    result->status = 0;
}

static void free_results(EncryptionResult *results, size_t count) {
    if (!results) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        free_result(&results[i]);
    }
    free(results);
}

static char *duplicate_string(const char *src) {
    if (!src) {
        return NULL;
    }
    size_t len = strlen(src);
    char *copy = (char *)malloc(len + 1);
    if (!copy) {
        return NULL;
    }
    memcpy(copy, src, len + 1);
    return copy;
}

static int hex_value(char c) {
    if (c >= '0' && c <= '9') {
        return c - '0';
    }
    if (c >= 'a' && c <= 'f') {
        return 10 + (c - 'a');
    }
    if (c >= 'A' && c <= 'F') {
        return 10 + (c - 'A');
    }
    return -1;
}

static int load_key_from_hex(const char *path, unsigned char key[KEY_SIZE]) {
    FILE *file = fopen(path, "r");
    if (!file) {
        fprintf(stderr, "Failed to open key file %s: %s\n", path, strerror(errno));
        return -1;
    }
    char buffer[KEY_SIZE * 2 + 16];
    size_t read_bytes = fread(buffer, 1, sizeof(buffer) - 1, file);
    fclose(file);
    if (read_bytes < KEY_SIZE * 2) {
        fprintf(stderr, "Key file %s must contain at least %d hexadecimal characters.\n", path, KEY_SIZE * 2);
        return -1;
    }
    size_t key_index = 0;
    size_t hex_index = 0;
    while (key_index < KEY_SIZE && hex_index + 1 < read_bytes) {
        while (hex_index < read_bytes && buffer[hex_index] != '\0' && buffer[hex_index] <= ' ') {
            hex_index++;
        }
        if (hex_index + 1 >= read_bytes) {
            break;
        }
        int high = hex_value(buffer[hex_index]);
        int low = hex_value(buffer[hex_index + 1]);
        if (high < 0 || low < 0) {
            fprintf(stderr, "Encountered non-hexadecimal characters in key file %s.\n", path);
            return -1;
        }
        key[key_index++] = (unsigned char)((high << 4) | low);
        hex_index += 2;
    }
    if (key_index != KEY_SIZE) {
        fprintf(stderr, "Key file %s does not contain enough data for a %d-byte key.\n", path, KEY_SIZE);
        return -1;
    }
    return 0;
}

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

static char *build_plaintext(const SequenceRecord *record) {
    if (!record || !record->sequence) {
        return NULL;
    }
    if (record->positions || record->reference) {
        const char *positions = record->positions ? record->positions : "";
        const char *reference = record->reference ? record->reference : "";
        size_t required = snprintf(NULL, 0, "Hotspot Positions: %s\nReference: %s\nSequence: %s", positions, reference,
                                   record->sequence);
        char *buffer = (char *)malloc(required + 1);
        if (!buffer) {
            return NULL;
        }
        snprintf(buffer, required + 1, "Hotspot Positions: %s\nReference: %s\nSequence: %s", positions, reference,
                 record->sequence);
        return buffer;
    }
    return duplicate_string(record->sequence);
}

typedef struct {
    const char *input_path;
    const char *key_path;
    const char *output_path;
    int threads;
} Options;

static void print_usage(const char *program) {
    fprintf(stderr,
            "Usage: %s --input <snp_hotspots_strings.tsv> --key <key.hex> --output <encrypted.tsv> [--threads N]\n",
            program);
}

static int parse_arguments(int argc, char **argv, Options *options) {
    if (!options) {
        return -1;
    }
    options->input_path = NULL;
    options->key_path = NULL;
    options->output_path = NULL;
    options->threads = 7;

    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (strcmp(arg, "--input") == 0 && i + 1 < argc) {
            options->input_path = argv[++i];
        } else if (strcmp(arg, "--key") == 0 && i + 1 < argc) {
            options->key_path = argv[++i];
        } else if (strcmp(arg, "--output") == 0 && i + 1 < argc) {
            options->output_path = argv[++i];
        } else if (strcmp(arg, "--threads") == 0 && i + 1 < argc) {
            options->threads = atoi(argv[++i]);
        } else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
            print_usage(argv[0]);
            return 1;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", arg);
            print_usage(argv[0]);
            return -1;
        }
    }

    if (!options->input_path || !options->key_path || !options->output_path) {
        fprintf(stderr, "Missing required arguments.\n");
        print_usage(argv[0]);
        return -1;
    }
    if (options->threads <= 0) {
        options->threads = 7;
    }
    return 0;
}

int main(int argc, char **argv) {
    Options options;
    int arg_status = parse_arguments(argc, argv, &options);
    if (arg_status != 0) {
        return arg_status < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
    }

    if (sodium_init() < 0) {
        fprintf(stderr, "Failed to initialise libsodium.\n");
        return EXIT_FAILURE;
    }

    unsigned char key[KEY_SIZE];
    if (load_key_from_hex(options.key_path, key) != 0) {
        return EXIT_FAILURE;
    }

    SequenceCollection collection;
    if (sequence_collection_init(&collection) != 0) {
        fprintf(stderr, "Failed to initialise sequence collection.\n");
        return EXIT_FAILURE;
    }

    if (load_sequence_records(options.input_path, &collection) != 0) {
        sequence_collection_free(&collection);
        return EXIT_FAILURE;
    }

    if (collection.count == 0) {
        fprintf(stderr, "No sequences were loaded from %s.\n", options.input_path);
        sequence_collection_free(&collection);
        return EXIT_FAILURE;
    }

    EncryptionResult *results = (EncryptionResult *)calloc(collection.count, sizeof(EncryptionResult));
    if (!results) {
        fprintf(stderr, "Failed to allocate memory for encryption results.\n");
        sequence_collection_free(&collection);
        return EXIT_FAILURE;
    }

    omp_set_num_threads(options.threads);

    int encountered_error = 0;
    size_t total_records = collection.count;
#pragma omp parallel for schedule(dynamic)
    for (long index = 0; index < (long)total_records; ++index) {
        size_t i = (size_t)index;
        const SequenceRecord *record = &collection.records[i];
        EncryptionResult result = {0};
        size_t plaintext_length = 0;
        unsigned char *ciphertext = NULL;
        unsigned char nonce[NONCE_SIZE];
        char *plaintext = build_plaintext(record);
        if (!plaintext) {
#pragma omp critical
            {
                encountered_error = 1;
                fprintf(stderr, "Failed to build plaintext for record %s.\n", record->identifier);
            }
            continue;
        }
        plaintext_length = strlen(plaintext);
        ciphertext = (unsigned char *)malloc(plaintext_length);
        if (!ciphertext) {
#pragma omp critical
            {
                encountered_error = 1;
                fprintf(stderr, "Failed to allocate ciphertext buffer for record %s.\n", record->identifier);
            }
            free(plaintext);
            continue;
        }
        randombytes_buf(nonce, sizeof nonce);
        if (crypto_stream_xchacha20_xor(ciphertext, (const unsigned char *)plaintext, plaintext_length, nonce, key) != 0) {
#pragma omp critical
            {
                encountered_error = 1;
                fprintf(stderr, "Encryption failed for record %s.\n", record->identifier);
            }
            free(ciphertext);
            free(plaintext);
            continue;
        }
        result.nonce_dna = binary_to_dna(nonce, sizeof nonce);
        result.ciphertext_dna = binary_to_dna(ciphertext, plaintext_length);
        free(ciphertext);
        free(plaintext);
        if (!result.nonce_dna || !result.ciphertext_dna) {
#pragma omp critical
            {
                encountered_error = 1;
                fprintf(stderr, "Failed to encode encrypted data for record %s.\n", record->identifier);
            }
            free(result.nonce_dna);
            free(result.ciphertext_dna);
            continue;
        }
        result.status = 0;
        results[i] = result;
    }

    if (encountered_error) {
        fprintf(stderr, "Aborting due to errors encountered during encryption.\n");
        free_results(results, collection.count);
        sequence_collection_free(&collection);
        return EXIT_FAILURE;
    }

    FILE *output = fopen(options.output_path, "w");
    if (!output) {
        fprintf(stderr, "Failed to open output file %s: %s\n", options.output_path, strerror(errno));
        free_results(results, collection.count);
        sequence_collection_free(&collection);
        return EXIT_FAILURE;
    }

    fprintf(output, "record_id\tnonce_dna\tciphertext_dna\n");
    for (size_t i = 0; i < collection.count; ++i) {
        const SequenceRecord *record = &collection.records[i];
        EncryptionResult *result = &results[i];
        const char *identifier = record->identifier ? record->identifier : "record";
        fprintf(output, "%s\t%s\t%s\n", identifier, result->nonce_dna, result->ciphertext_dna);
    }

    fclose(output);
    free_results(results, collection.count);
    sequence_collection_free(&collection);
    return EXIT_SUCCESS;
}
