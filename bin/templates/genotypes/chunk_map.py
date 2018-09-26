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
        start_chr = min(pos)
        end_chr = max(pos)

        for start in range(start_chr, end_chr, ${CHUNK_SIZE}):
            end = min(start + ${CHUNK_SIZE}, end_chr)

            # if the last chunk is small, add it to the penultimate chunk
            if (end_chr - end) < (${CHUNK_SIZE}/5):
                end = end_chr

            if [ x for x in pos if x >= start and x <= end ]:
                chunks.write('{}\\t{}\\t{}\\n'.format(chromosome, start, end))

                if end == end_chr:
                    break