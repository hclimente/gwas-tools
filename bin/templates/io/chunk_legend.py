#!/usr/bin/env python

# LEGEND, CHUNK_SIZE
# chunks
import gzip

chromosome = '${LEGEND}'.split('.')[0].replace('1000GP_Phase3_chr', '')
first_pos = float('inf')
last_pos = float('-inf')

with gzip.open('${LEGEND}', mode = 'rt') as legend:
  legend.readline()
  for line in legend.readlines():
    line = line.strip().split(' ')
    pos = int(line[1])
    if pos < first_pos: first_pos = pos
    if pos > last_pos: last_pos = pos

with open('chunks', 'w') as chunks:
  for i in range(first_pos, last_pos, ${CHUNK_SIZE}):
    chunks.write('{}\\t{}\\t{}\\n'.format(chromosome, i, min(i + ${CHUNK_SIZE}, last_pos)))