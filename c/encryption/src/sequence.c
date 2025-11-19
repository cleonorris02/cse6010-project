/*         
 * Outline and logic generated with ChatGPT (OpenAI), Oct 2025.         
 * Reviewed and modified by Viru Repalle.         
 * */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "sequence.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

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

static void free_record(SequenceRecord *record) {
    if (!record) {
        return;
    }
    free(record->identifier);
    free(record->positions);
    free(record->reference);
    free(record->sequence);
    record->identifier = NULL;
    record->positions = NULL;
    record->reference = NULL;
    record->sequence = NULL;
}

static char *trim(char *value) {
    if (!value) {
        return value;
    }
    while (isspace((unsigned char)*value)) {
        value++;
    }
    if (*value == '\0') {
        return value;
    }
    char *end = value + strlen(value) - 1;
    while (end > value && isspace((unsigned char)*end)) {
        *end-- = '\0';
    }
    return value;
}

int sequence_collection_init(SequenceCollection *collection) {
    if (!collection) {
        return -1;
    }
    collection->records = NULL;
    collection->count = 0;
    collection->capacity = 0;
    return 0;
}

void sequence_collection_free(SequenceCollection *collection) {
    if (!collection) {
        return;
    }
    for (size_t i = 0; i < collection->count; ++i) {
        free_record(&collection->records[i]);
    }
    free(collection->records);
    collection->records = NULL;
    collection->count = 0;
    collection->capacity = 0;
}

int sequence_collection_append(SequenceCollection *collection, SequenceRecord record) {
    if (!collection) {
        free_record(&record);
        return -1;
    }
    if (collection->count == collection->capacity) {
        size_t new_capacity = collection->capacity == 0 ? 16 : collection->capacity * 2;
        SequenceRecord *resized = (SequenceRecord *)realloc(collection->records, new_capacity * sizeof(SequenceRecord));
        if (!resized) {
            free_record(&record);
            return -1;
        }
        collection->records = resized;
        collection->capacity = new_capacity;
    }
    collection->records[collection->count++] = record;
    return 0;
}

static int split_columns(char *line, char ***columns_out, size_t *count_out) {
    if (!line || !columns_out || !count_out) {
        return -1;
    }
    size_t capacity = 8;
    size_t count = 0;
    char **columns = (char **)malloc(capacity * sizeof(char *));
    if (!columns) {
        return -1;
    }
    char *cursor = line;
    while (cursor && *cursor) {
        if (count == capacity) {
            size_t new_capacity = capacity * 2;
            char **resized = (char **)realloc(columns, new_capacity * sizeof(char *));
            if (!resized) {
                free(columns);
                return -1;
            }
            columns = resized;
            capacity = new_capacity;
        }
        char *next = strchr(cursor, '\t');
        if (next) {
            *next = '\0';
        }
        columns[count++] = cursor;
        if (!next) {
            break;
        }
        cursor = next + 1;
    }
    *columns_out = columns;
    *count_out = count;
    return 0;
}

static void free_columns(char **columns) {
    free(columns);
}

static ssize_t read_line(FILE *file, char **buffer, size_t *size) {
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200809L
    return getline(buffer, size, file);
#else
    if (!file || !buffer || !size) {
        return -1;
    }
    if (*buffer == NULL || *size == 0) {
        *size = 256;
        *buffer = (char *)malloc(*size);
        if (!*buffer) {
            return -1;
        }
    }
    size_t index = 0;
    int ch = 0;
    while ((ch = fgetc(file)) != EOF) {
        if (index + 1 >= *size) {
            size_t new_size = *size * 2;
            char *resized = (char *)realloc(*buffer, new_size);
            if (!resized) {
                return -1;
            }
            *buffer = resized;
            *size = new_size;
        }
        (*buffer)[index++] = (char)ch;
        if (ch == '\n') {
            break;
        }
    }
    if (index == 0 && ch == EOF) {
        return -1;
    }
    (*buffer)[index] = '\0';
    return (ssize_t)index;
#endif
}

static int locate_column(const char *name) {
    if (!name) {
        return -1;
    }
    if (strcasecmp(name, "record_id") == 0 || strcasecmp(name, "id") == 0 || strcasecmp(name, "hotspot_id") == 0) {
        return 0;
    }
    if (strcasecmp(name, "hotspot_positions") == 0 || strcasecmp(name, "positions") == 0) {
        return 1;
    }
    if (strcasecmp(name, "reference") == 0 || strcasecmp(name, "reference_sequence") == 0) {
        return 2;
    }
    if (strcasecmp(name, "hotspot_string") == 0 || strcasecmp(name, "hotspot_sequence") == 0 ||
        strcasecmp(name, "sequence") == 0 || strcasecmp(name, "dna_string") == 0) {
        return 3;
    }
    return -1;
}

int load_sequence_records(const char *path, SequenceCollection *collection) {
    if (!path || !collection) {
        return -1;
    }
    FILE *file = fopen(path, "r");
    if (!file) {
        fprintf(stderr, "Failed to open %s: %s\n", path, strerror(errno));
        return -1;
    }

    char *line = NULL;
    size_t line_size = 0;
    ssize_t length = read_line(file, &line, &line_size);
    if (length < 0) {
        fprintf(stderr, "The TSV file %s is empty or unreadable.\n", path);
        free(line);
        fclose(file);
        return -1;
    }

    char *header_line = trim(line);
    char **header_columns = NULL;
    size_t header_count = 0;
    if (split_columns(header_line, &header_columns, &header_count) != 0 || header_count == 0) {
        fprintf(stderr, "Unable to parse header columns in %s.\n", path);
        free(line);
        fclose(file);
        return -1;
    }

    int id_index = -1;
    int positions_index = -1;
    int reference_index = -1;
    int sequence_index = -1;
    for (size_t i = 0; i < header_count; ++i) {
        int column_type = locate_column(header_columns[i]);
        switch (column_type) {
            case 0:
                id_index = (int)i;
                break;
            case 1:
                positions_index = (int)i;
                break;
            case 2:
                reference_index = (int)i;
                break;
            case 3:
                sequence_index = (int)i;
                break;
            default:
                break;
        }
    }
    if (sequence_index < 0) {
        fprintf(stderr, "The TSV file %s must contain a column with DNA strings (e.g., hotspot_string).\n", path);
        free_columns(header_columns);
        free(line);
        fclose(file);
        return -1;
    }

    size_t row_index = 0;
    while ((length = read_line(file, &line, &line_size)) >= 0) {
        char *trimmed = trim(line);
        if (*trimmed == '\0' || *trimmed == '#') {
            continue;
        }
        char **columns = NULL;
        size_t column_count = 0;
        if (split_columns(trimmed, &columns, &column_count) != 0) {
            fprintf(stderr, "Failed to split columns on row %zu in %s.\n", row_index + 1, path);
            free_columns(header_columns);
            free(line);
            fclose(file);
            return -1;
        }
        SequenceRecord record = {0};
        if (id_index >= 0 && (size_t)id_index < column_count) {
            record.identifier = duplicate_string(trim(columns[id_index]));
        } else {
            char generated[32];
            snprintf(generated, sizeof(generated), "record_%zu", row_index);
            record.identifier = duplicate_string(generated);
        }
        if (!record.identifier) {
            free_columns(columns);
            free_columns(header_columns);
            free(line);
            fclose(file);
            return -1;
        }
        if (positions_index >= 0 && (size_t)positions_index < column_count) {
            record.positions = duplicate_string(trim(columns[positions_index]));
        }
        if (reference_index >= 0 && (size_t)reference_index < column_count) {
            record.reference = duplicate_string(trim(columns[reference_index]));
        }
        if ((size_t)sequence_index < column_count) {
            record.sequence = duplicate_string(trim(columns[sequence_index]));
        }
        free_columns(columns);

        if (!record.sequence || record.sequence[0] == '\0') {
            free_record(&record);
            free_columns(header_columns);
            free(line);
            fclose(file);
            fprintf(stderr, "Encountered a row without a DNA sequence at index %zu in %s.\n", row_index + 1, path);
            return -1;
        }

        if (sequence_collection_append(collection, record) != 0) {
            free_record(&record);
            free_columns(header_columns);
            free(line);
            fclose(file);
            return -1;
        }
        row_index++;
    }

    free_columns(header_columns);
    free(line);
    fclose(file);
    return 0;
}
