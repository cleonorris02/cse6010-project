#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stddef.h>

typedef struct {
    char *identifier;
    char *positions;
    char *reference;
    char *sequence;
} SequenceRecord;

typedef struct {
    SequenceRecord *records;
    size_t count;
    size_t capacity;
} SequenceCollection;

int sequence_collection_init(SequenceCollection *collection);
void sequence_collection_free(SequenceCollection *collection);
int sequence_collection_append(SequenceCollection *collection, SequenceRecord record);

int load_sequence_records(const char *path, SequenceCollection *collection);

#endif /* SEQUENCE_H */
