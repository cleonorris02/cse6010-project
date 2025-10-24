#define _POSIX_C_SOURCE 200809L

#include "hotspot.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

static void trim_newline(char *line) {
    if (!line) {
        return;
    }
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        line[--len] = '\0';
    }
}

int hotspot_collection_init(HotspotCollection *collection) {
    if (!collection) {
        return -1;
    }
    collection->records = NULL;
    collection->count = 0;
    collection->capacity = 0;
    return 0;
}

void hotspot_collection_free(HotspotCollection *collection) {
    if (!collection) {
        return;
    }
    for (size_t i = 0; i < collection->count; ++i) {
        free(collection->records[i].positions);
        free(collection->records[i].reference);
        free(collection->records[i].alternate);
    }
    free(collection->records);
    collection->records = NULL;
    collection->count = 0;
    collection->capacity = 0;
}

int hotspot_collection_append(HotspotCollection *collection, HotspotRecord record) {
    if (!collection) {
        return -1;
    }
    if (collection->count == collection->capacity) {
        size_t new_capacity = collection->capacity == 0 ? 4 : collection->capacity * 2;
        HotspotRecord *new_records = realloc(collection->records, new_capacity * sizeof(HotspotRecord));
        if (!new_records) {
            return -1;
        }
        collection->records = new_records;
        collection->capacity = new_capacity;
    }
    collection->records[collection->count++] = record;
    return 0;
}

static size_t *parse_positions(const char *line, size_t *out_count) {
    const char *prefix = "Hotspot Positions:";
    size_t prefix_len = strlen(prefix);
    if (strncmp(line, prefix, prefix_len) != 0) {
        errno = EINVAL;
        return NULL;
    }
    const char *cursor = line + prefix_len;
    while (*cursor && isspace((unsigned char)*cursor)) {
        cursor++;
    }
    size_t capacity = 8;
    size_t count = 0;
    size_t *positions = malloc(capacity * sizeof(size_t));
    if (!positions) {
        return NULL;
    }
    char *copy = strdup(cursor);
    if (!copy) {
        free(positions);
        return NULL;
    }
    char *token = strtok(copy, ",");
    while (token) {
        while (*token && isspace((unsigned char)*token)) {
            token++;
        }
        if (*token == '\0') {
            token = strtok(NULL, ",");
            continue;
        }
        char *endptr = NULL;
        errno = 0;
        unsigned long long value = strtoull(token, &endptr, 10);
        if (errno != 0 || (endptr && *endptr != '\0' && !isspace((unsigned char)*endptr))) {
            free(copy);
            free(positions);
            errno = EINVAL;
            return NULL;
        }
        if (count == capacity) {
            capacity *= 2;
            size_t *new_positions = realloc(positions, capacity * sizeof(size_t));
            if (!new_positions) {
                free(copy);
                free(positions);
                return NULL;
            }
            positions = new_positions;
        }
        positions[count++] = (size_t)value;
        token = strtok(NULL, ",");
    }
    free(copy);
    if (count == 0) {
        free(positions);
        errno = EINVAL;
        return NULL;
    }
    size_t *shrunk = realloc(positions, count * sizeof(size_t));
    if (shrunk) {
        positions = shrunk;
    }
    *out_count = count;
    return positions;
}

static char *parse_sequence_line(const char *line, const char *prefix) {
    size_t prefix_len = strlen(prefix);
    if (strncmp(line, prefix, prefix_len) != 0) {
        errno = EINVAL;
        return NULL;
    }
    const char *cursor = line + prefix_len;
    while (*cursor && isspace((unsigned char)*cursor)) {
        cursor++;
    }
    return strdup(cursor);
}

static int read_next_content_line(FILE *file, char **lineptr, size_t *n) {
    ssize_t read;
    while ((read = getline(lineptr, n, file)) != -1) {
        (void)read;
        trim_newline(*lineptr);
        if (**lineptr == '\0') {
            continue;
        }
        return 0;
    }
    return -1;
}

int parse_hotspot_file(const char *path, HotspotCollection *collection) {
    if (!path || !collection) {
        errno = EINVAL;
        return -1;
    }
    FILE *file = fopen(path, "r");
    if (!file) {
        return -1;
    }
    char *line = NULL;
    size_t line_capacity = 0;
    int status = 0;

    while (read_next_content_line(file, &line, &line_capacity) == 0) {
        if (strncmp(line, "Hotspot Positions:", 18) != 0) {
            fprintf(stderr, "Unexpected line: %s\n", line);
            status = -1;
            break;
        }
        size_t position_count = 0;
        size_t *positions = parse_positions(line, &position_count);
        if (!positions) {
            status = -1;
            break;
        }
        if (read_next_content_line(file, &line, &line_capacity) != 0) {
            fprintf(stderr, "Missing Reference line after Hotspot Positions\n");
            free(positions);
            status = -1;
            break;
        }
        char *reference = parse_sequence_line(line, "Reference:");
        if (!reference) {
            fprintf(stderr, "Malformed Reference line: %s\n", line);
            free(positions);
            status = -1;
            break;
        }

        char *alternate = NULL;
        long file_pos = ftell(file);
        if (file_pos >= 0) {
            if (read_next_content_line(file, &line, &line_capacity) == 0) {
                if (strncmp(line, "Alternate:", 10) == 0) {
                    alternate = parse_sequence_line(line, "Alternate:");
                    if (!alternate) {
                        fprintf(stderr, "Malformed Alternate line: %s\n", line);
                        free(positions);
                        free(reference);
                        status = -1;
                        break;
                    }
                } else {
                    // Not an Alternate line, rewind to allow parsing in the next iteration.
                    if (fseek(file, file_pos, SEEK_SET) != 0) {
                        fprintf(stderr, "Failed to rewind file pointer\n");
                        free(positions);
                        free(reference);
                        status = -1;
                        break;
                    }
                }
            }
        }

        HotspotRecord record = {
            .positions = positions,
            .position_count = position_count,
            .reference = reference,
            .alternate = alternate,
        };
        if (hotspot_collection_append(collection, record) != 0) {
            fprintf(stderr, "Failed to store hotspot record\n");
            free(positions);
            free(reference);
            free(alternate);
            status = -1;
            break;
        }
    }

    free(line);
    fclose(file);
    if (status != 0) {
        hotspot_collection_free(collection);
    }
    return status;
}

