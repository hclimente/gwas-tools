#!/usr/bin/env python

with open('${BED}', 'r') as bed_in, open('out.map', 'w') as map_out:
    for line in bed_in.readlines():
        chromosome, _, pos, rs = line.strip().split()
        map_out.write('{}\\t{}\\t0\\t{}\\n'.format(chromosome, rs, pos))