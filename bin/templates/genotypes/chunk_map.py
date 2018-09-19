#!/usr/bin/env python

# MAP, CHUNK_SIZE
# chunks

chrs = {}

with open('${MAP}', mode = 'r') as gwas_map:
    for line in gwas_map.readlines():
        line = line.strip().split('\t')
        chromosome = int(line[0])
        pos = int(line[3])

        chrs.setdefault(chromosome, [])
        chrs[chromosome].append(pos)

with open('chunks', 'w') as chunks:
    for chromosome in chrs:
        pos = chrs[chromosome]
        for start in range(min(pos), max(pos), ${CHUNK_SIZE}):
            end = min(start + ${CHUNK_SIZE}, max(pos))
            if [ x for x in pos if x >= start and x <= end ]:
                chunks.write('{}\\t{}\\t{}\\n'.format(chromosome, start, end))