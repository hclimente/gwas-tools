#!/usr/bin/env python

with open('out.bed', 'w') as bed_out, open('${MAP}', 'r') as map_in:
    for line in map_in.readlines():
        chromosome, rs, gpos, pos = line.strip().split()
        bed_out.write('chr{}\\t{}\\t{}\\t{}\\n'.format(chromosome, int(pos) - 1, pos, rs))