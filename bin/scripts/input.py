def INPUT():
    def TARGETOME_BUILDING():
        # SETUP
        import os
        directory = os.getcwd()
        if 'temp' not in os.listdir(directory):
            os.mkdir(''.join([directory, '/temp']))
        if 'targetome' not in os.listdir(''.join([directory, '/temp'])):
            os.mkdir(''.join([directory, '/temp/targetome']))
        for a in os.listdir(''.join([directory, '/temp/targetome/'])):
            os.remove(os.path.join(''.join([directory, '/temp/targetome/']), a))
        # READING THE SETTINGS
        settings = open(''.join([directory, '/input/settings.txt']), 'r').read().splitlines()
        settings_dict = {}
        for z in list(range(len(settings))):
            if '\t' not in settings[z] and settings[z] != []:
                settings_dict[settings[z].split(': ')[0]] = {}
                b = z + 1
                while b < len(settings) and '\t' in settings[b]:
                    settings_dict[settings[z].split(': ')[0]][settings[b].split('\t')[1].split(': ')[0]] = \
                    settings[b].split('\t')[1].split(': ')[1]
                    b = b + 1
        settings = settings_dict
        del settings_dict
        # PREPARING THE INPUT FILE
        file_U1s = ''.join([directory, '/input/query.fa'])
        file_output = ''.join([directory, '/temp/targetome/targetome'])
        with open(file_U1s, 'r') as library:
            # creating a multi-fasta file with U1 sequences
            # input sequences are the sense strand of the U1
            library = library.read().splitlines()
            for a in list(range(len(library))):
                if '>' in library[a]:
                    ID = library[a].split('>')[1]
                    # U1 sequences are reverse-complemented
                    sequence = list(library[a + 1])
                    sequence.reverse()
                    for c in list(range(len(sequence))):
                        if sequence[c] == 'A' or sequence[c] == 'a':
                            sequence[c] = 'T'
                        elif sequence[c] == 'T' or sequence[c] == 't':
                            sequence[c] = 'A'
                        elif sequence[c] == 'C' or sequence[c] == 'c':
                            sequence[c] = 'G'
                        elif sequence[c] == 'G' or sequence[c] == 'g':
                            sequence[c] = 'C'
                    sequence = ''.join(sequence)
                    # duplicate U1s are discarded
                    # each U1 creates 8 categories, based on the register:
                    # normal (norm)
                    # mismatches (mm)
                    # single bulge on the sense (BS1)
                    # double bulge on the sense (BS2)
                    # single bulge on the antisense (BA1)
                    # double bulge on the antisense (BA2)
                    # asymmetric loop with double bulge on the sense (ALS)
                    # asymmetric loop with double bulge on the antisense (ALA)
                    library[a] = [ID,
                                  {'norm': sequence,
                                   'BS1': {},
                                   'BS2': {},
                                   'BA1': {},
                                   'BA2': {},
                                   'ALS': {},
                                   'ALA': {},
                                   'mm1_norm': {},
                                   'mm1_BS1': {},
                                   'mm1_BS2': {},
                                   'mm1_BA1': {},
                                   'mm1_BA2': {},
                                   'mm1_ALS': {},
                                   'mm1_ALA': {},
                                   'mm2_norm': {},
                                   'mm2_BS1': {},
                                   'mm2_BS2': {},
                                   'mm2_BA1': {},
                                   'mm2_BA2': {},
                                   'mm2_ALS': {},
                                   'mm2_ALA': {},
                                   'mm3_norm': {},
                                   'mm3_BS1': {},
                                   'mm3_BS2': {},
                                   'mm3_BA1': {},
                                   'mm3_BA2': {},
                                   'mm3_ALS': {},
                                   'mm3_ALA': {},
                                   'mm4_norm': {},
                                   'mm4_BS1': {},
                                   'mm4_BS2': {},
                                   'mm4_BA1': {},
                                   'mm4_BA2': {},
                                   'mm4_ALS': {},
                                   'mm4_ALA': {},
                                   'mm5_norm': {},
                                   'mm5_BS1': {},
                                   'mm5_BS2': {},
                                   'mm5_BA1': {},
                                   'mm5_BA2': {},
                                   'mm5_ALS': {},
                                   'mm5_ALA': {},
                                   'mm6_norm': {},
                                   'mm6_BS1': {},
                                   'mm6_BS2': {},
                                   'mm6_BA1': {},
                                   'mm6_ALS': {},
                                   'mm7_norm': {},
                                   'mm7_BS1': {},
                                   'mm7_BS2': {}}]
                    library[a + 1] = []
            while [] in library:
                library.remove([])
            del ID
            del sequence
            # labelling U1s with duplicate IDs
            single_IDs = {}
            duplicate_IDs = {}
            for a in library:
                if a[0] not in single_IDs:
                    single_IDs[a[0]] = ''
                else:
                    duplicate_IDs[a[0]] = 1
            for a in list(range(len(library))):
                ID = library[a][0]
                if ID in duplicate_IDs:
                    library[a][0] = '_'.join([ID, str(duplicate_IDs[ID])])
                    duplicate_IDs[ID] = duplicate_IDs[ID] + 1
            del ID
            del single_IDs
            del duplicate_IDs
            # creating alternative registers for each U1
            # a targetome of all possible targets (alternative registers) is created from each U1
            # for each U1 7 additional categories of alternative registers are created (no repeated sequences within the same category)
            print('creating targetome...')
            for U1 in list(range(len(library))):
                # mm_norm
                if settings['registers']['full'] != 'n':
                    for mismatches in list(range(1, len(library[U1][1]['norm']) + 1 - int(
                            settings['registers']['full']))):  # indicates the number of mismatches
                        sequence = list(library[U1][1]['norm'])
                        # determining the combinations of sequence positions to be changed
                        # all combinations are unique regardless of order, and positions are not repeated
                        # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                        import itertools
                        combinations = []
                        for combination in itertools.combinations(list(range(len(sequence))), mismatches):
                            if list(combination) != []:
                                combinations.append(list(combination))
                        # changing bases for each combination of positions
                        for combination in combinations:
                            # each combination of positions has a series of combinations of bases
                            # each combination of bases of each combination of positions is a potential new sequence
                            # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                            for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                sequence = list(library[U1][1]['norm'])
                                for change in list(range(mismatches)):
                                    sequence[combination[change]] = list(nucleotides)[change]
                                sequence = ''.join(sequence)
                                library[U1][1][''.join(['mm', str(mismatches), '_norm'])][sequence] = ''
                # BS1
                if settings['registers']['BS1'] != 'n':
                    # adding 1 base to the original U1 sequence, but only after position 1 and before the last position
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence))):
                        for modification_1 in ['A', 'T', 'C', 'G']:
                            sequence = list(library[U1][1]['norm'])
                            sequence.insert(position, modification_1)
                            sequence = ''.join(sequence)
                            if len(library[U1][1]['norm']) >= int(settings['registers']['BS1']):
                                library[U1][1]['BS1'][sequence] = ''
                                # introducing mismatches
                                for mismatches in list(range(1, len(library[U1][1]['norm']) + 1 - int(
                                        settings['registers']['BS1']))):  # indicates the number of mismatches
                                    sequence_mm = list(sequence)
                                    # determining the combinations of sequence positions to be changed
                                    # all combinations are unique regardless of order, and positions are not repeated
                                    # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                    import itertools
                                    combinations = []
                                    for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                        if list(combination) != [] \
                                                and (
                                                position not in combination and
                                                position - 1 not in combination and
                                                position + 1 not in combination): # excluding positions corresponding to the bulge or next to it
                                            combinations.append(list(combination))
                                    # changing bases for each combination of positions
                                    for combination in combinations:
                                        # each combination of positions has a series of combinations of bases
                                        # each combination of bases of each combination of positions is a potential new sequence
                                        # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                        for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                            sequence_mm = list(sequence)
                                            for change in list(range(mismatches)):
                                                sequence_mm[combination[change]] = list(nucleotides)[change]
                                            sequence_mm = ''.join(sequence_mm)
                                            library[U1][1][''.join(['mm', str(mismatches), '_BS1'])][sequence_mm] = ''
                # BS2
                if settings['registers']['BS2'] != 'n':
                    # adding 2 bases to the original U1 sequence, but only after position 1 and before the last position
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence))):
                        for modification_1 in ['A', 'T', 'C', 'G']:
                            for modification_2 in ['A', 'T', 'C', 'G']:
                                sequence = list(library[U1][1]['norm'])
                                sequence.insert(position, ''.join([modification_1, modification_2]))
                                sequence = ''.join(sequence)
                                if len(library[U1][1]['norm']) >= int(settings['registers']['BS2']):
                                    library[U1][1]['BS2'][sequence] = ''
                                    # introducing mismatches
                                    for mismatches in list(range(1, len(library[U1][1]['norm']) + 1 - int(
                                            settings['registers']['BS2']))):  # indicates the number of mismatches
                                        sequence_mm = list(sequence)
                                        # determining the combinations of sequence positions to be changed
                                        # all combinations are unique regardless of order, and positions are not repeated
                                        # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                        import itertools
                                        combinations = []
                                        for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                            if list(combination) != [] \
                                                    and (
                                                    position not in combination and
                                                    position - 1 not in combination and
                                                    position + 1 not in combination and
                                                    position + 2 not in combination): # excluding positions corresponding to the bulge or next to it
                                                combinations.append(list(combination))
                                        # changing bases for each combination of positions
                                        for combination in combinations:
                                            # each combination of positions has a series of combinations of bases
                                            # each combination of bases of each combination of positions is a potential new sequence
                                            # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                            for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                                sequence_mm = list(sequence)
                                                for change in list(range(mismatches)):
                                                    sequence_mm[combination[change]] = list(nucleotides)[change]
                                                sequence_mm = ''.join(sequence_mm)
                                                library[U1][1][''.join(['mm', str(mismatches), '_BS2'])][sequence_mm] = ''
                # BA1
                if settings['registers']['BA1'] != 'n':
                    # removing 1 base
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence) - 1)):
                        sequence = list(library[U1][1]['norm'])
                        del sequence[position]
                        sequence = ''.join(sequence)
                        if len(library[U1][1]['norm']) - 1 >= int(settings['registers']['BA1']):
                            library[U1][1]['BA1'][sequence] = ''
                            # introducing mismatches
                            for mismatches in list(range(1, len(library[U1][1]['norm']) - int(
                                    settings['registers']['BA1']))):  # indicates the number of mismatches
                                sequence_mm = list(sequence)
                                # determining the combinations of sequence positions to be changed
                                # all combinations are unique regardless of order, and positions are not repeated
                                # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                import itertools
                                combinations = []
                                for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                    if list(combination) != [] \
                                            and (
                                            position not in combination and
                                            position - 1 not in combination): # excluding positions corresponding to the bulge or next to it
                                        combinations.append(list(combination))
                                # changing bases for each combination of positions
                                for combination in combinations:
                                    # each combination of positions has a series of combinations of bases
                                    # each combination of bases of each combination of positions is a potential new sequence
                                    # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                    for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                        sequence_mm = list(sequence)
                                        for change in list(range(mismatches)):
                                            sequence_mm[combination[change]] = list(nucleotides)[change]
                                        sequence_mm = ''.join(sequence_mm)
                                        library[U1][1][''.join(['mm', str(mismatches), '_BA1'])][sequence_mm] = ''
                # BA2
                if settings['registers']['BA2'] != 'n':
                    # removing 2 bases
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence) - 2)):
                        sequence = list(library[U1][1]['norm'])
                        del sequence[position]
                        del sequence[position]
                        sequence = ''.join(sequence)
                        if len(library[U1][1]['norm']) - 2 >= int(settings['registers']['BA2']):
                            library[U1][1]['BA2'][sequence] = ''
                            # introducing mimatches
                            for mismatches in list(range(1, len(library[U1][1]['norm']) - 1 - int(
                                    settings['registers']['BA2']))):  # indicates the number of mismatches
                                sequence_mm = list(sequence)
                                # determining the combinations of sequence positions to be changed
                                # all combinations are unique regardless of order, and positions are not repeated
                                # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                import itertools
                                combinations = []
                                for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                    if list(combination) != [] \
                                            and (
                                            position not in combination and
                                            position - 1 not in combination): # excluding positions corresponding to the bulge or next to it
                                        combinations.append(list(combination))
                                # changing bases for each combination of positions
                                for combination in combinations:
                                    # each combination of positions has a series of combinations of bases
                                    # each combination of bases of each combination of positions is a potential new sequence
                                    # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                    for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                        sequence_mm = list(sequence)
                                        for change in list(range(mismatches)):
                                            sequence_mm[combination[change]] = list(nucleotides)[change]
                                        sequence_mm = ''.join(sequence_mm)
                                        library[U1][1][''.join(['mm', str(mismatches), '_BA2'])][sequence_mm] = ''
                # ALS
                if settings['registers']['ALS'] != 'n':
                    # removing 1 base and adding 2
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence) - 1)):
                        for modification_1 in ['A', 'T', 'C', 'G']:
                            for modification_2 in ['A', 'T', 'C', 'G']:
                                sequence = list(library[U1][1]['norm'])
                                del sequence[position]
                                sequence.insert(position, ''.join([modification_1, modification_2]))
                                sequence = ''.join(sequence)
                                if len(library[U1][1]['norm']) - 1 >= int(settings['registers']['ALS']):
                                    library[U1][1]['ALS'][sequence] = ''
                                    # introducing mismatches
                                    for mismatches in list(range(1, len(library[U1][1]['norm']) - int(
                                            settings['registers']['ALS']))):  # indicates the number of mismatches
                                        sequence_mm = list(sequence)
                                        # determining the combinations of sequence positions to be changed
                                        # all combinations are unique regardless of order, and positions are not repeated
                                        # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                        import itertools
                                        combinations = []
                                        for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                            if list(combination) != [] \
                                                    and (
                                                    position not in combination and
                                                    position - 1 not in combination and
                                                    position + 1 not in combination and
                                                    position + 2 not in combination): # excluding positions corresponding to the bulge or next to it
                                                combinations.append(list(combination))
                                        # changing bases for each combination of positions
                                        for combination in combinations:
                                            # each combination of positions has a series of combinations of bases
                                            # each combination of bases of each combination of positions is a potential new sequence
                                            # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                            for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                                sequence_mm = list(sequence)
                                                for change in list(range(mismatches)):
                                                    sequence_mm[combination[change]] = list(nucleotides)[change]
                                                sequence_mm = ''.join(sequence_mm)
                                                library[U1][1][''.join(['mm', str(mismatches), '_ALS'])][sequence_mm] = ''
                # ALA
                if settings['registers']['ALA'] != 'n':
                    # removing 2 bases and adding 1
                    sequence = list(library[U1][1]['norm'])
                    for position in list(range(1, len(sequence) - 2)):
                        for modification_1 in ['A', 'T', 'C', 'G']:
                            sequence = list(library[U1][1]['norm'])
                            del sequence[position]
                            del sequence[position]
                            sequence.insert(position, modification_1)
                            sequence = ''.join(sequence)
                            if len(library[U1][1]['norm']) - 2 >= int(settings['registers']['ALA']):
                                library[U1][1]['ALA'][sequence] = ''
                                # introducing mismatches
                                for mismatches in list(range(1, len(library[U1][1]['norm']) - 1 - int(
                                        settings['registers']['ALA']))):  # indicates the number of mismatches
                                    sequence_mm = list(sequence)
                                    # determining the combinations of sequence positions to be changed
                                    # all combinations are unique regardless of order, and positions are not repeated
                                    # the number of combinations is calculated with the "combination" formula = 11! / ((11 - b)! b!)
                                    import itertools
                                    combinations = []
                                    for combination in itertools.combinations(list(range(len(sequence_mm))), mismatches):
                                        if list(combination) != [] \
                                                and (
                                                position not in combination and
                                                position - 1 not in combination and
                                                position + 1 not in combination): # excluding positions corresponding to the bulge or next to it
                                            combinations.append(list(combination))
                                    # changing bases for each combination of positions
                                    for combination in combinations:
                                        # each combination of positions has a series of combinations of bases
                                        # each combination of bases of each combination of positions is a potential new sequence
                                        # the number of combinations is calculated as = (n position combinations) * (3 ^ b)
                                        for nucleotides in itertools.product('ATCG', repeat=mismatches):
                                            sequence_mm = list(sequence)
                                            for change in list(range(mismatches)):
                                                sequence_mm[combination[change]] = list(nucleotides)[change]
                                            sequence_mm = ''.join(sequence_mm)
                                            library[U1][1][''.join(['mm', str(mismatches), '_ALA'])][sequence_mm] = ''
                # removing redundant registers between different mismatch classes of the same register
                if settings['registers']['full'] == 'n': # adjustment
                    library[U1][1]['full'] = {}
                library_filtered = {}
                for register in ['norm', 'BS1', 'BS2', 'BA1', 'BA2', 'ALS', 'ALA']:
                    for mismatches in range(1, 8):
                        if ''.join(['mm', str(mismatches), '_', register]) in library[U1][1]:
                            registers_filtered = {}
                            for sequence in library[U1][1][''.join(['mm', str(mismatches), '_', register])]:
                                if sequence not in registers_filtered: # checking the same class
                                    redundant = 'ok'
                                    for register_class in range(mismatches): # checking superior classes
                                        if register_class == 0 and register == 'norm': # checking the superior class for mm1_norm
                                            if sequence == library[U1][1]['norm']:
                                                redundant = 'redundant'
                                        elif register_class == 0 and register != 'norm': # checking the superior class for any other mm1
                                            if sequence in library[U1][1][register]:
                                                redundant = 'redundant'
                                        else: # checking the superior class for any others
                                            if sequence in library[U1][1][''.join(['mm', str(register_class), '_', register])]:
                                                redundant = 'redundant'
                                    if redundant == 'ok':
                                        registers_filtered[sequence] = ''
                            library_filtered[''.join(['mm', str(mismatches), '_', register])] = registers_filtered
                for register in library[U1][1]:
                    if register in library_filtered:
                        library[U1][1][register] = library_filtered[register]
                del library_filtered
            # outputting a different file for each category
            for a in ['norm', 'BS1', 'BS2', 'BA1', 'BA2', 'ALS', 'ALA',
                      'mm1_norm', 'mm2_norm', 'mm3_norm', 'mm4_norm', 'mm5_norm', 'mm6_norm', 'mm7_norm',
                      'mm1_BS1', 'mm2_BS1', 'mm3_BS1', 'mm4_BS1', 'mm5_BS1', 'mm6_BS1', 'mm7_BS1',
                      'mm1_BS2', 'mm2_BS2', 'mm3_BS2', 'mm4_BS2', 'mm5_BS2', 'mm6_BS2', 'mm7_BS2',
                      'mm1_BA1', 'mm2_BA1', 'mm3_BA1', 'mm4_BA1', 'mm5_BA1', 'mm6_BA1',
                      'mm1_BA2', 'mm2_BA2', 'mm3_BA2', 'mm4_BA2', 'mm5_BA2',
                      'mm1_ALS', 'mm2_ALS', 'mm3_ALS', 'mm4_ALS', 'mm5_ALS', 'mm6_ALS',
                      'mm1_ALA', 'mm2_ALA', 'mm3_ALA', 'mm4_ALA', 'mm5_ALA']:
                file = '_'.join([
                    file_output,
                    a,
                    'query.fa'])
                with open(file, 'x') as output:
                    if a == 'norm' and settings['registers']['full'] != 'n':
                        for d in list(range(len(library))):
                            if d == 0:
                                output.write('>')
                                output.write(library[d][0])
                                output.write('\n')
                                output.write(library[d][1][a])
                            else:
                                output.write('\n')
                                output.write('>')
                                output.write(library[d][0])
                                output.write('\n')
                                output.write(library[d][1][a])
                    else:
                        for d in list(range(len(library))):
                            count = 1
                            for e in library[d][1][a]:
                                if d == 0 and count == 0:
                                    output.write('>')
                                    output.write('_'.join([library[d][0], str(count)]))
                                    output.write('\n')
                                    output.write(e)
                                else:
                                    output.write('\n')
                                    output.write('>')
                                    output.write('_'.join([library[d][0], str(count)]))
                                    output.write('\n')
                                    output.write(e)
                                count = count + 1
            del count
            del file
            del output
            del library
    TARGETOME_BUILDING()
    def DATABASE_BUILDING():
        # SETUP
        import os
        directory = os.getcwd()
        if 'database' not in os.listdir(''.join([directory, '/temp'])):
            os.mkdir(''.join([directory, '/temp/database']))
        # READING THE SETTINGS
        settings = open(''.join([directory, '/input/settings.txt']), 'r').read().splitlines()
        settings_dict = {}
        for z in list(range(len(settings))):
            if '\t' not in settings[z] and settings[z] != []:
                settings_dict[settings[z].split(': ')[0]] = {}
                b = z + 1
                while b < len(settings) and '\t' in settings[b]:
                    settings_dict[settings[z].split(': ')[0]][settings[b].split('\t')[1].split(': ')[0]] = \
                        settings[b].split('\t')[1].split(': ')[1]
                    b = b + 1
        settings = settings_dict
        del settings_dict
        # PREPARING EXONS SEQUENCES FOR DATABASE BUILDING
        if settings['database']['exons'] == 'y' and \
                (os.path.exists(''.join([directory, '/temp/database/database_exons.fa'])) != True and
                 os.path.exists(''.join([directory, '/temp/database/database_exons_1.fa'])) != True):
            file_annotation = ''.join([directory, '/input/annotation.gtf'])
            file_assembly = ''.join([directory, '/input/assembly.fna'])
            file_output_exons = ''.join([directory, '/temp/database/database_exons.fa'])
            with open(file_output_exons, 'x') as output_exons:
                # MT genes are excluded
                # gene names are assigned with gene IDs in order to avoid duplicates (separated by "__", as some names contain "_")
                # then gene name will be removed because otherwise the id will be too long for BLAST
                # transcripts and exons are first stored based on the gene ID
                # then introns coordinates are extrapolated
                print('extracting exon annotation (it may take a few minutes if this is the first time running utargetome)...')
                annotation = open(file_annotation, 'r').read().splitlines()
                del annotation[0:5]
                annotation_genes = {}
                for a in annotation:  # gene annotation
                    if a.split('\t')[2] == 'gene' and a.split('\t')[0] != 'MT':
                        if '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]]) in annotation_genes:
                            print('duplicate gene found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]])] = {
                                'chromosome': a.split('\t')[0],
                                'direction': a.split('\t')[6],
                                'transcript': {}}
                for a in annotation:  # transcript annotation
                    if a.split('\t')[2] == 'transcript' and a.split('\t')[0] != 'MT':
                        if a.split('\t')[8].split('"')[5] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])]['transcript']:
                            print('duplicate transcript found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]] = {'exons': [], 'introns': []}
                for a in annotation:  # exon annotation
                    if a.split('\t')[2] == 'exon' and a.split('\t')[0] != 'MT':
                        if [int(a.split('\t')[3]), int(a.split('\t')[4])] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])]['transcript'][
                            a.split('\t')[8].split('"')[5]]['exons']:
                            print('duplicate exon found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]]['exons'].append(
                                [int(a.split('\t')[3]), int(a.split('\t')[4])])
                del annotation
                for a in annotation_genes:  # extracting intron coordinates
                    for b in annotation_genes[a]['transcript']:
                        annotation_genes[a]['transcript'][b]['exons'] = sorted(
                            annotation_genes[a]['transcript'][b]['exons'],
                            key=lambda x: x[0])  # not all exons are in order
                        for c in list(range(0, len(annotation_genes[a]['transcript'][b][
                                                       'exons']))):  # transcripts with 1 exon produce no intron
                            if c != len(annotation_genes[a]['transcript'][b]['exons']) - 1:
                                annotation_genes[a]['transcript'][b]['introns'].append(
                                    [str(annotation_genes[a]['transcript'][b]['exons'][c][1] + 1),
                                     str(annotation_genes[a]['transcript'][b]['exons'][c + 1][0] - 1)])
                                annotation_genes[a]['transcript'][b]['exons'][c] = [
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][0]),
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][1])]
                            else:
                                annotation_genes[a]['transcript'][b]['exons'][c] = [
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][0]),
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][1])]
                # arranging assembly file
                print('arranging assembly...')
                assembly = open(file_assembly, 'r').read().split('>')
                for a in list(range(len(assembly))):
                    if assembly[a][0:2] != 'NC':
                        assembly[a] = []
                assembly[len(assembly) - 1] = []
                while [] in assembly:
                    assembly.remove([])
                assembly.insert(0, [])
                for a in list(range(1, len(assembly))):
                    assembly[a] = assembly[a].splitlines()
                    assembly[a].remove(assembly[a][0])
                    assembly[a] = ''.join(assembly[a])
                # extracting exon sequences
                annotation_exons = {}
                count = 0
                for a in annotation_genes:
                    for b in annotation_genes[a]['transcript']:
                        for c in annotation_genes[a]['transcript'][b]['exons']:
                            if annotation_genes[a]['chromosome'] != 'X' and annotation_genes[a][
                                'chromosome'] != 'Y':
                                sequence = assembly[int(annotation_genes[a]['chromosome'])][int(c[0]) - 1:int(c[1])]
                            elif annotation_genes[a]['chromosome'] == 'X':
                                sequence = assembly[23][int(c[0]) - 1:int(c[1])]
                            elif annotation_genes[a]['chromosome'] == 'Y':
                                sequence = assembly[24][int(c[0]) - 1:int(c[1])]
                            # reverse and complement
                            sequence = list(sequence)
                            if annotation_genes[a]['direction'] == '+':
                                for d in list(range(len(sequence))):
                                    if sequence[d] == 'a':
                                        sequence[d] = 'A'
                                    elif sequence[d] == 't':
                                        sequence[d] = 'T'
                                    elif sequence[d] == 'c':
                                        sequence[d] = 'C'
                                    elif sequence[d] == 'g':
                                        sequence[d] = 'G'
                            elif annotation_genes[a]['direction'] == '-':
                                sequence.reverse()
                                for d in list(range(len(sequence))):
                                    if sequence[d] == 'A' or sequence[d] == 'a':
                                        sequence[d] = 'T'
                                    elif sequence[d] == 'T' or sequence[d] == 't':
                                        sequence[d] = 'A'
                                    elif sequence[d] == 'C' or sequence[d] == 'c':
                                        sequence[d] = 'G'
                                    elif sequence[d] == 'G' or sequence[d] == 'g':
                                        sequence[d] = 'C'
                            sequence = ''.join(sequence)
                            # redundant exons are counted
                            if '_'.join([a.split('__')[1], c[0], c[1]]) not in annotation_exons:
                                annotation_exons['_'.join([a.split('__')[1], c[0], c[1]])] = [sequence, 1]
                            else:
                                annotation_exons['_'.join([a.split('__')[1], c[0], c[1]])][1] += 1
                    count = count + 1
                del assembly
                del annotation_genes
                # printing exon sequences
                print('printing exon database...')
                count = 0
                for a in annotation_exons:
                    if count == 0:
                        output_exons.write('>')
                        output_exons.write('_'.join([a, str(annotation_exons[a][1])]))
                        output_exons.write('\n')
                        output_exons.write(annotation_exons[a][0])
                    else:
                        output_exons.write('\n')
                        output_exons.write('>')
                        output_exons.write('_'.join([a, str(annotation_exons[a][1])]))
                        output_exons.write('\n')
                        output_exons.write(annotation_exons[a][0])
                    count = count + 1
                del count
                del sequence
                del annotation_exons
                del output_exons
        # PREPARING INTRON SEQUENCES FOR DATABASE BUILDING
        if settings['database']['introns'] == 'y' and \
                (os.path.exists(''.join([directory, '/temp/database/database_introns.fa'])) != True and
                os.path.exists(''.join([directory, '/temp/database/database_introns_1.fa'])) != True):
            file_annotation = ''.join([directory, '/input/annotation.gtf'])
            file_assembly = ''.join([directory, '/input/assembly.fna'])
            file_output_introns = ''.join([directory, '/temp/database/database_introns.fa'])
            with open(file_output_introns, 'x') as output_introns:
                # MT genes are excluded
                # gene names are assigned with gene IDs in order to avoid duplicates (separated by "__", as some names contain "_")
                # then gene name will be removed because otherwise the id will be too long for BLAST
                # transcripts and exons are first stored based on the gene ID
                # then introns coordinates are extrapolated
                print('extracting intron annotation (it may take up to 30 minutes if this is the first time running utargetome)...')
                annotation = open(file_annotation, 'r').read().splitlines()
                del annotation[0:5]
                annotation_genes = {}
                for a in annotation:  # gene annotation
                    if a.split('\t')[2] == 'gene' and a.split('\t')[0] != 'MT':
                        if '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]]) in annotation_genes:
                            print('duplicate gene found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]])] = {
                                'chromosome': a.split('\t')[0],
                                'direction': a.split('\t')[6],
                                'transcript': {}}
                for a in annotation:  # transcript annotation
                    if a.split('\t')[2] == 'transcript' and a.split('\t')[0] != 'MT':
                        if a.split('\t')[8].split('"')[5] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])]['transcript']:
                            print('duplicate transcript found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]] = {'exons': [], 'introns': []}
                for a in annotation:  # exon annotation
                    if a.split('\t')[2] == 'exon' and a.split('\t')[0] != 'MT':
                        if [int(a.split('\t')[3]), int(a.split('\t')[4])] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])]['transcript'][
                            a.split('\t')[8].split('"')[5]]['exons']:
                            print('duplicate exon found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]]['exons'].append(
                                [int(a.split('\t')[3]), int(a.split('\t')[4])])
                del annotation
                for a in annotation_genes:  # extracting intron coordinates
                    for b in annotation_genes[a]['transcript']:
                        annotation_genes[a]['transcript'][b]['exons'] = sorted(
                            annotation_genes[a]['transcript'][b]['exons'],
                            key=lambda x: x[0])  # not all exons are in order
                        for c in list(range(0, len(annotation_genes[a]['transcript'][b][
                                                       'exons']))):  # transcripts with 1 exon produce no intron
                            if c != len(annotation_genes[a]['transcript'][b]['exons']) - 1:
                                annotation_genes[a]['transcript'][b]['introns'].append(
                                    [str(annotation_genes[a]['transcript'][b]['exons'][c][1] + 1),
                                     str(annotation_genes[a]['transcript'][b]['exons'][c + 1][0] - 1)])
                                annotation_genes[a]['transcript'][b]['exons'][c] = [
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][0]),
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][1])]
                            else:
                                annotation_genes[a]['transcript'][b]['exons'][c] = [
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][0]),
                                    str(annotation_genes[a]['transcript'][b]['exons'][c][1])]
                # arranging assembly file
                print('arranging assembly...')
                assembly = open(file_assembly, 'r').read().split('>')
                for a in list(range(len(assembly))):
                    if assembly[a][0:2] != 'NC':
                        assembly[a] = []
                assembly[len(assembly) - 1] = []
                while [] in assembly:
                    assembly.remove([])
                assembly.insert(0, [])
                for a in list(range(1, len(assembly))):
                    assembly[a] = assembly[a].splitlines()
                    assembly[a].remove(assembly[a][0])
                    assembly[a] = ''.join(assembly[a])
                # extracting intron sequences
                annotation_introns = {}
                count = 0
                for a in annotation_genes:
                    for b in annotation_genes[a]['transcript']:
                        for c in annotation_genes[a]['transcript'][b]['introns']:
                            if annotation_genes[a]['chromosome'] != 'X' and annotation_genes[a][
                                'chromosome'] != 'Y':
                                sequence = assembly[int(annotation_genes[a]['chromosome'])][int(c[0]) - 1:int(c[1])]
                            elif annotation_genes[a]['chromosome'] == 'X':
                                sequence = assembly[23][int(c[0]) - 1:int(c[1])]
                            elif annotation_genes[a]['chromosome'] == 'Y':
                                sequence = assembly[24][int(c[0]) - 1:int(c[1])]
                            # reverse and complement
                            sequence = list(sequence)
                            if annotation_genes[a]['direction'] == '+':
                                for d in list(range(len(sequence))):
                                    if sequence[d] == 'a':
                                        sequence[d] = 'A'
                                    elif sequence[d] == 't':
                                        sequence[d] = 'T'
                                    elif sequence[d] == 'c':
                                        sequence[d] = 'C'
                                    elif sequence[d] == 'g':
                                        sequence[d] = 'G'
                            elif annotation_genes[a]['direction'] == '-':
                                sequence.reverse()
                                for d in list(range(len(sequence))):
                                    if sequence[d] == 'A' or sequence[d] == 'a':
                                        sequence[d] = 'T'
                                    elif sequence[d] == 'T' or sequence[d] == 't':
                                        sequence[d] = 'A'
                                    elif sequence[d] == 'C' or sequence[d] == 'c':
                                        sequence[d] = 'G'
                                    elif sequence[d] == 'G' or sequence[d] == 'g':
                                        sequence[d] = 'C'
                            sequence = ''.join(sequence)
                            # redundant introns are counted
                            if '_'.join([a.split('__')[1], c[0], c[1]]) not in annotation_introns:
                                annotation_introns['_'.join([a.split('__')[1], c[0], c[1]])] = [sequence, 1]
                            else:
                                annotation_introns['_'.join([a.split('__')[1], c[0], c[1]])][1] += 1
                    count = count + 1
                del assembly
                del annotation_genes
                # printing intron sequences
                print('printing intron database...')
                count = 0
                for a in annotation_introns:
                    if count == 0:
                        output_introns.write('>')
                        output_introns.write('_'.join([a, str(annotation_introns[a][1])]))
                        output_introns.write('\n')
                        output_introns.write(annotation_introns[a][0])
                    else:
                        output_introns.write('\n')
                        output_introns.write('>')
                        output_introns.write('_'.join([a, str(annotation_introns[a][1])]))
                        output_introns.write('\n')
                        output_introns.write(annotation_introns[a][0])
                    count = count + 1
                del count
                del sequence
                del annotation_introns
                del output_introns
        # PREPARING 5'-SS SEQUENCES FOR DATABASE BUILDING
        if settings['database']['5SSs'] == 'y' and \
                (os.path.exists(''.join([directory, '/temp/database/database_5SSs.fa'])) != True and
                os.path.exists(''.join([directory, '/temp/database/database_5SSs_1.fa'])) != True):
            file_annotation = ''.join([directory, '/input/annotation.gtf'])
            file_assembly = ''.join([directory, '/input/assembly.fna'])
            file_output = ''.join([directory, '/temp/database/database_5SSs.fa'])
            with open(file_output, 'x') as output:
                # the sequence around 5'-SSs (-100 > + 100) is extracted
                # MT genes and last exons are excluded
                # each 5'-SS is labelled through chromosome, direction, gene ID and position
                print('extracting 5-SS annotation (it may take a few minutes if this is the first time running utargetome)...')
                annotation = open(file_annotation, 'r').read().splitlines()
                del annotation[0:5]
                annotation_genes = {}  # extracting gene/transcript information to find last exons
                for a in annotation:  # gene annotation
                    if a.split('\t')[2] == 'gene' and a.split('\t')[0] != 'MT':
                        if '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]]) in annotation_genes:
                            print('duplicate gene found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]])] = {
                                'chromosome': a.split('\t')[0],
                                'direction': a.split('\t')[6],
                                'transcript': {}}
                for a in annotation:  # transcript annotation
                    if a.split('\t')[2] == 'transcript' and a.split('\t')[0] != 'MT':
                        if a.split('\t')[8].split('"')[5] in \
                                annotation_genes[
                                    '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])][
                                    'transcript']:
                            print('duplicate transcript found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]] = {'exons': []}
                for a in annotation:  # exon annotation
                    if a.split('\t')[2] == 'exon' and a.split('\t')[0] != 'MT':
                        if [int(a.split('\t')[3]), int(a.split('\t')[4])] in \
                                annotation_genes[
                                    '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])][
                                    'transcript'][a.split('\t')[8].split('"')[5]]['exons']:
                            print('duplicate exon found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]]['exons'].append(
                                [int(a.split('\t')[3]), int(a.split('\t')[4])])
                del annotation
                annotation_exons = {}  # filtering last exons
                for a in annotation_genes:
                    chromosome = annotation_genes[a]['chromosome']
                    direction = annotation_genes[a]['direction']
                    ID = a.split('__')[1]
                    for b in annotation_genes[a]['transcript']:
                        annotation_genes[a]['transcript'][b]['exons'] = sorted(
                            annotation_genes[a]['transcript'][b]['exons'],
                            key=lambda x: x[0])  # not all exons are in order
                        if direction == '+':  # determining last exon based on gene direction
                            last_exon = annotation_genes[a]['transcript'][b]['exons'][
                                len(annotation_genes[a]['transcript'][b]['exons']) - 1]
                        else:
                            last_exon = annotation_genes[a]['transcript'][b]['exons'][0]
                        for c in annotation_genes[a]['transcript'][b]['exons']:  # selecting non-last exons
                            if c != last_exon:
                                exon = '_'.join([chromosome, direction, ID, str(c[0]), str(c[1])])
                                # redundant exons are counted
                                if exon not in annotation_exons:
                                    annotation_exons[exon] = 1
                                else:
                                    annotation_exons[exon] += 1
                del exon
                del chromosome
                del direction
                del ID
                del last_exon
                del annotation_genes
                # arranging assembly file
                print('arranging assembly...')
                assembly = open(file_assembly, 'r').read().split('>')
                for a in list(range(len(assembly))):
                    if assembly[a][0:2] != 'NC':
                        assembly[a] = []
                assembly[len(assembly) - 1] = []
                while [] in assembly:
                    assembly.remove([])
                assembly.insert(0, [])
                for a in list(range(1, len(assembly))):
                    assembly[a] = assembly[a].splitlines()
                    assembly[a].remove(assembly[a][0])
                    assembly[a] = ''.join(assembly[a])
                # extracting splice site sequences
                splice_sites = {}
                count = 0
                for a in annotation_exons:
                    if a.split('_')[1] == '+':
                        if a.split('_')[0] != 'X' and a.split('_')[0] != 'Y':
                            sequence = assembly[int(a.split('_')[0])][
                                       int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        elif a.split('_')[0] == 'X':
                            sequence = assembly[23][int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        elif a.split('_')[0] == 'Y':
                            sequence = assembly[24][int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        # uppercase
                        sequence = list(sequence)
                        for c in list(range(len(sequence))):
                            if sequence[c] == 'a':
                                sequence[c] = 'A'
                            elif sequence[c] == 't':
                                sequence[c] = 'T'
                            elif sequence[c] == 'c':
                                sequence[c] = 'C'
                            elif sequence[c] == 'g':
                                sequence[c] = 'G'
                        sequence = ''.join(sequence)
                        splice_sites['_'.join(
                            [a.split('_')[2], a.split('_')[3], a.split('_')[4], str(annotation_exons[a])])] = sequence
                    elif a.split('_')[1] == '-':
                        if a.split('_')[0] != 'X' and a.split('_')[0] != 'Y':
                            sequence = assembly[int(a.split('_')[0])][
                                       int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        elif a.split('_')[0] == 'X':
                            sequence = assembly[23][int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        elif a.split('_')[0] == 'Y':
                            sequence = assembly[24][int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        # reverse and complement
                        sequence = list(sequence)
                        sequence.reverse()
                        for c in list(range(len(sequence))):
                            if sequence[c] == 'A' or sequence[c] == 'a':
                                sequence[c] = 'T'
                            elif sequence[c] == 'T' or sequence[c] == 't':
                                sequence[c] = 'A'
                            elif sequence[c] == 'C' or sequence[c] == 'c':
                                sequence[c] = 'G'
                            elif sequence[c] == 'G' or sequence[c] == 'g':
                                sequence[c] = 'C'
                        sequence = ''.join(sequence)
                        splice_sites['_'.join(
                            [a.split('_')[2], a.split('_')[3], a.split('_')[4], str(annotation_exons[a])])] = sequence
                    count = count + 1
                del assembly
                del sequence
                del annotation_exons
                del count
                # printing splice site sequences
                print('printing 5-SS database...')
                count = 0
                for a in splice_sites:
                    if count == 0:
                        output.write('>')
                        output.write(a)
                        output.write('\n')
                        output.write(splice_sites[a])
                    else:
                        output.write('\n')
                        output.write('>')
                        output.write(a)
                        output.write('\n')
                        output.write(splice_sites[a])
                    count = count + 1
                del count
                del splice_sites
                del output
        # PREPARING 3'-SS SEQUENCES FOR DATABASE BUILDING
        if settings['database']['3SSs'] == 'y' and \
                (os.path.exists(''.join([directory, '/temp/database/database_3SSs.fa'])) != True and
                os.path.exists(''.join([directory, '/temp/database/database_3SSs_1.fa'])) != True):
            file_annotation = ''.join([directory, '/input/annotation.gtf'])
            file_assembly = ''.join([directory, '/input/assembly.fna'])
            file_output = ''.join([directory, '/temp/database/database_3SSs.fa'])
            with open(file_output, 'x') as output:
                # the sequence around 3'-SSs (-100 > + 100) is extracted
                # MT genes and last exons are excluded
                # each 3'-SS is labelled through chromosome, direction, gene ID and position
                print('extracting 3-SS annotation (it may take a few minutes if this is the first time running utargetome)...')
                annotation = open(file_annotation, 'r').read().splitlines()
                del annotation[0:5]
                annotation_genes = {}  # extracting gene/transcript information to find last exons
                for a in annotation:  # gene annotation
                    if a.split('\t')[2] == 'gene' and a.split('\t')[0] != 'MT':
                        if '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]]) in annotation_genes:
                            print('duplicate gene found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[5], a.split('\t')[8].split('"')[1]])] = {
                                'chromosome': a.split('\t')[0],
                                'direction': a.split('\t')[6],
                                'transcript': {}}
                for a in annotation:  # transcript annotation
                    if a.split('\t')[2] == 'transcript' and a.split('\t')[0] != 'MT':
                        if a.split('\t')[8].split('"')[5] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])]['transcript']:
                            print('duplicate transcript found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[9], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]] = {'exons': []}
                for a in annotation:  # exon annotation
                    if a.split('\t')[2] == 'exon' and a.split('\t')[0] != 'MT':
                        if [int(a.split('\t')[3]), int(a.split('\t')[4])] in annotation_genes[
                            '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])]['transcript'][
                            a.split('\t')[8].split('"')[5]]['exons']:
                            print('duplicate exon found')
                        else:
                            annotation_genes[
                                '__'.join([a.split('\t')[8].split('"')[11], a.split('\t')[8].split('"')[1]])][
                                'transcript'][a.split('\t')[8].split('"')[5]]['exons'].append(
                                [int(a.split('\t')[3]), int(a.split('\t')[4])])
                del annotation
                annotation_exons = {}  # filtering first exons
                for a in annotation_genes:
                    chromosome = annotation_genes[a]['chromosome']
                    direction = annotation_genes[a]['direction']
                    ID = a.split('__')[1]
                    for b in annotation_genes[a]['transcript']:
                        annotation_genes[a]['transcript'][b]['exons'] = sorted(
                            annotation_genes[a]['transcript'][b]['exons'],
                            key=lambda x: x[0])  # not all exons are in order
                        if direction == '+':  # determining first exon based on gene direction
                            first_exon = annotation_genes[a]['transcript'][b]['exons'][0]
                        else:
                            first_exon = annotation_genes[a]['transcript'][b]['exons'][
                                len(annotation_genes[a]['transcript'][b]['exons']) - 1]
                        for c in annotation_genes[a]['transcript'][b]['exons']:  # selecting non-first exons
                            if c != first_exon:
                                exon = '_'.join([chromosome, direction, ID, str(c[0]), str(c[1])])
                                # redundant exons are counted
                                if exon not in annotation_exons:
                                    annotation_exons[exon] = 1
                                else:
                                    annotation_exons[exon] += 1
                del exon
                del chromosome
                del direction
                del ID
                del first_exon
                del annotation_genes
                # arranging assembly file
                print('arranging assembly...')
                assembly = open(file_assembly, 'r').read().split('>')
                for a in list(range(len(assembly))):
                    if assembly[a][0:2] != 'NC':
                        assembly[a] = []
                assembly[len(assembly) - 1] = []
                while [] in assembly:
                    assembly.remove([])
                assembly.insert(0, [])
                for a in list(range(1, len(assembly))):
                    assembly[a] = assembly[a].splitlines()
                    assembly[a].remove(assembly[a][0])
                    assembly[a] = ''.join(assembly[a])
                # extracting splice site sequences
                splice_sites = {}
                count = 0
                for a in annotation_exons:
                    if a.split('_')[1] == '+':
                        if a.split('_')[0] != 'X' and a.split('_')[0] != 'Y':
                            sequence = assembly[int(a.split('_')[0])][
                                       int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        elif a.split('_')[0] == 'X':
                            sequence = assembly[23][int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        elif a.split('_')[0] == 'Y':
                            sequence = assembly[24][int(a.split('_')[3]) - 101:int(a.split('_')[3]) + 99]
                        # uppercase
                        sequence = list(sequence)
                        for c in list(range(len(sequence))):
                            if sequence[c] == 'a':
                                sequence[c] = 'A'
                            elif sequence[c] == 't':
                                sequence[c] = 'T'
                            elif sequence[c] == 'c':
                                sequence[c] = 'C'
                            elif sequence[c] == 'g':
                                sequence[c] = 'G'
                        sequence = ''.join(sequence)
                        splice_sites['_'.join(
                            [a.split('_')[2], a.split('_')[3], a.split('_')[4], str(annotation_exons[a])])] = sequence
                    elif a.split('_')[1] == '-':
                        if a.split('_')[0] != 'X' and a.split('_')[0] != 'Y':
                            sequence = assembly[int(a.split('_')[0])][
                                       int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        elif a.split('_')[0] == 'X':
                            sequence = assembly[23][int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        elif a.split('_')[0] == 'Y':
                            sequence = assembly[24][int(a.split('_')[4]) - 100:int(a.split('_')[4]) + 100]
                        # reverse and complement
                        sequence = list(sequence)
                        sequence.reverse()
                        for c in list(range(len(sequence))):
                            if sequence[c] == 'A' or sequence[c] == 'a':
                                sequence[c] = 'T'
                            elif sequence[c] == 'T' or sequence[c] == 't':
                                sequence[c] = 'A'
                            elif sequence[c] == 'C' or sequence[c] == 'c':
                                sequence[c] = 'G'
                            elif sequence[c] == 'G' or sequence[c] == 'g':
                                sequence[c] = 'C'
                        sequence = ''.join(sequence)
                        splice_sites['_'.join(
                            [a.split('_')[2], a.split('_')[3], a.split('_')[4], str(annotation_exons[a])])] = sequence
                    count = count + 1
                del assembly
                del sequence
                del annotation_exons
                del count
                # printing splice site sequences
                print('printing 3-SS database...')
                count = 0
                for a in splice_sites:
                    if count == 0:
                        output.write('>')
                        output.write(a)
                        output.write('\n')
                        output.write(splice_sites[a])
                    else:
                        output.write('\n')
                        output.write('>')
                        output.write(a)
                        output.write('\n')
                        output.write(splice_sites[a])
                    count = count + 1
                del count
                del splice_sites
                del output
    DATABASE_BUILDING()
    def DATABASE_CHUNKS():
        # SETUP
        import os
        directory = os.getcwd()
        # READING THE SETTINGS
        settings = open(''.join([directory, '/input/settings.txt']), 'r').read().splitlines()
        settings_dict = {}
        for z in list(range(len(settings))):
            if '\t' not in settings[z] and settings[z] != []:
                settings_dict[settings[z].split(': ')[0]] = {}
                b = z + 1
                while b < len(settings) and '\t' in settings[b]:
                    settings_dict[settings[z].split(': ')[0]][settings[b].split('\t')[1].split(': ')[0]] = \
                        settings[b].split('\t')[1].split(': ')[1]
                    b = b + 1
        settings = settings_dict
        del settings_dict
        # CREATING CHUNKS
        for z in ['exons', 'introns', '5SSs', '3SSs']:
            if settings['database'][z] == 'y' and \
                    (os.path.exists(''.join([directory, '/temp/database/database_', z, '_1.fa'])) != True):
                print('creating chunks for:', z, 'database')
                import math
                if z == 'introns':
                    chunk_numbers = int(settings['partitions']['number']) * 11
                else:
                    chunk_numbers = int(settings['partitions']['number'])
                if chunk_numbers != 1: # creating chunks only for specified database
                    file = open(''.join([directory, '/temp/database/database_', z, '.fa']), 'r').read().splitlines()
                    database = []
                    for a in range(len(file)):
                        if '>' in file[a]:
                            database.append([file[a], file[a + 1]])
                    del file
                    chunk_size = math.ceil(len(database) / chunk_numbers)
                    chunk_current = 0
                    chunks = []
                    for b in range(chunk_numbers):
                        with open(''.join([directory, '/temp/database/database_', z, '_', str(b + 1), '.fa']),
                                  'x') as output:
                            for d in range(len(database[chunk_current:chunk_current + chunk_size])):
                                if d == 0:
                                    output.write(database[chunk_current:chunk_current + chunk_size][d][0])
                                    output.write('\n')
                                    output.write(database[chunk_current:chunk_current + chunk_size][d][1])
                                else:
                                    output.write('\n')
                                    output.write(database[chunk_current:chunk_current + chunk_size][d][0])
                                    output.write('\n')
                                    output.write(database[chunk_current:chunk_current + chunk_size][d][1])
                        chunk_current += chunk_size
                    del database
                    os.remove(''.join([directory, '/temp/database/database_', z, '.fa']))
                else:
                    os.rename(''.join([directory, '/temp/database/database_', z, '.fa']),
                              ''.join([directory, '/temp/database/database_', z, '_1.fa']))
    DATABASE_CHUNKS()
    def MAKE_BLAST():
        # SETUP
        import os
        directory = os.getcwd()
        if 'BLAST' not in os.listdir(''.join([directory, '/temp'])):
            os.mkdir(''.join([directory, '/temp/BLAST']))
        for a in os.listdir(''.join([directory, '/temp/BLAST/'])):
            os.remove(os.path.join(''.join([directory, '/temp/BLAST/']), a))
        # READING THE SETTINGS
        settings = open(''.join([directory, '/input/settings.txt']), 'r').read().splitlines()
        settings_dict = {}
        for z in list(range(len(settings))):
            if '\t' not in settings[z] and settings[z] != []:
                settings_dict[settings[z].split(': ')[0]] = {}
                b = z + 1
                while b < len(settings) and '\t' in settings[b]:
                    settings_dict[settings[z].split(': ')[0]][settings[b].split('\t')[1].split(': ')[0]] = \
                    settings[b].split('\t')[1].split(': ')[1]
                    b = b + 1
        settings = settings_dict
        del settings_dict
        # DETERMINING THE ANTISENSE LENGTH
        file_U1s = ''.join([directory, '/input/query.fa'])
        library = open(file_U1s, 'r').read().splitlines()
        length = len(library[1])
        del library
        del file_U1s
        # CREATING THE SCRIPT FOR BLAST
        with open(''.join([directory, '/scripts/blast_exe']), 'w') as output:
            for a in [a for a in settings['database'] if settings['database'][a] == 'y']:
                if a == 'introns':
                    chunk_numbers = 11 * int(settings['partitions']['number'])
                else:
                    chunk_numbers = 1 * int(settings['partitions']['number'])
                for b in range(chunk_numbers): # making a script for each chunk
                    output.write('\n')
                    output.write(''.join(['makeblastdb -in temp/database/database_',
                                          a,
                                          '_',
                                          str(b + 1),
                                          '.fa -input_type fasta -dbtype nucl -parse_seqids -blastdb_version 5 -max_file_sz 4GB -out temp/database/database']))
                    for c in settings['registers']:
                        if c == 'full':
                            c = 'norm'
                            word_size = str(length)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 8)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'BS1':
                            word_size = str(length + 1)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 8)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'BS2':
                            word_size = str(length + 2)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 8)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'BA1':
                            word_size = str(length - 1)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 7)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'BA2':
                            word_size = str(length - 2)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 6)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'ALS':
                            word_size = str(length + 1)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 7)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
                        elif c == 'ALA':
                            word_size = str(length - 1)
                            output.write('\n')
                            output.write(''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                                  c,
                                                  '_query.fa -task blastn-short -word_size ',
                                                  word_size,
                                                  ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt 2>/dev/null']))
                            output.write('\n')
                            output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt >> temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '.txt']))
                            output.write('\n')
                            output.write(''.join(['rm temp/BLAST/BLAST_',
                                                  c,
                                                  '_',
                                                  a,
                                                  '_',
                                                  str(b + 1),
                                                  '.txt']))
                            for e in list(range(1, 6)):
                                output.write('\n')
                                output.write(
                                    ''.join(['blastn -db temp/database/database -query temp/targetome/targetome_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_query.fa -task blastn-short -word_size ',
                                             word_size,
                                             ' -strand plus -evalue 1000000000 -max_target_seqs 1000000000 -out temp/BLAST/BLAST_',
                                             'mm',
                                             str(e),
                                             '_',
                                             c,
                                             '_',
                                             a,
                                             '_',
                                             str(b + 1),
                                             '.txt 2>/dev/null']))
                                output.write('\n')
                                output.write(''.join(['grep -e "Query=" -e ">" -e "Sbjct" temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt >> temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '.txt']))
                                output.write('\n')
                                output.write(''.join(['rm temp/BLAST/BLAST_',
                                                      'mm',
                                                      str(e),
                                                      '_',
                                                      c,
                                                      '_',
                                                      a,
                                                      '_',
                                                      str(b + 1),
                                                      '.txt']))
    MAKE_BLAST()
INPUT()
