#!/usr/bin/env python

# MAP, CHUNK_SIZE, GENOME
# chunks

chrs = {}
sex_chromosomes = {'23': 'X', '25': 'X', 'X': 'X', 
                   '24': 'Y', 'Y': 'Y'}

if "${GENOME}" ==  'GRCh38':
    par1 = {'X': (10001,2781479), 'Y': (10001,2781479)}
    par2 = {'X': (155701383,156030895), 'Y': (56887903,57217415)}
elif "${GENOME}" ==  'GRCh37':
    par1 = {'X': (60001,2699520), 'Y': (10001,2649520)}
    par2 = {'X': (154931044,155260560), 'Y': (59034050,59363566)}

with open('${MAP}', mode = 'r') as gwas_map:
    for line in gwas_map.readlines():
        line = line.strip().split('\\t')
        chromosome = sex_chromosomes.get(line[0], line[0])
        pos = int(line[3])

        chrs.setdefault(chromosome, [])
        chrs[chromosome].append(pos)

with open('chunks', 'w') as chunks:
    for chromosome in chrs:
        pos = chrs[chromosome]
        start_chr = min(pos)
        end_chr = max(pos)
        chromosome_name = chromosome

        if chromosome in sex_chromosomes.keys():
            start_chr = par1[chromosome][1] + 1
            end_chr = par2[chromosome][0] - 1
            chromosome_name = chromosome + '_NONPAR'

            chunks.write('{}\\t{}\\t{}\\n'.format(chromosome + '_PAR1', 
                                                  par1[chromosome][0], par1[chromosome][1]))
            chunks.write('{}\\t{}\\t{}\\n'.format(chromosome + '_PAR2', 
                                                  par2[chromosome][0], par2[chromosome][1]))

        for start in range(start_chr, end_chr, ${CHUNK_SIZE}):
            end = min(start + ${CHUNK_SIZE}, end_chr)

            # if the last chunk is small, add it to the penultimate chunk
            if (end_chr - end) < (${CHUNK_SIZE}/5):
                end = end_chr

            if [ x for x in pos if x >= start and x <= end ]:
                chunks.write('{}\\t{}\\t{}\\n'.format(chromosome_name, start, end))

                if end == end_chr:
                    break
