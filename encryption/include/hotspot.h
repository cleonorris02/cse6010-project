#ifndef HOTSPOT_H
#define HOTSPOT_H

#include <stddef.h>

typedef struct {
    size_t *positions;
    size_t position_count;
    char *reference;
    char *alternate;
} HotspotRecord;

typedef struct {
    HotspotRecord *records;
    size_t count;
    size_t capacity;
} HotspotCollection;

int hotspot_collection_init(HotspotCollection *collection);
void hotspot_collection_free(HotspotCollection *collection);
int hotspot_collection_append(HotspotCollection *collection, HotspotRecord record);

int parse_hotspot_file(const char *path, HotspotCollection *collection);

#endif
