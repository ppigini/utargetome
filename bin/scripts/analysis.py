def ANALYSIS():
    # SETUP
    import os
    directory = os.getcwd()
    for a in os.listdir(''.join([directory, '/results/'])):
        os.remove(os.path.join(''.join([directory, '/results/']), a))
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
    # COUNTING BLASTED TARGETS
    file_BLAST = ''.join([directory, '/temp/BLAST/BLAST'])
    file_output = ''.join([directory, '/results/'])
    databases = []
    for x in settings['database']:
        if settings['database'][x] == 'yes':
            databases.append(x)
    registers = {'norm' : {}, 'mm1_norm' : {}, 'mm2_norm' : {}, 'mm3_norm' : {}, 'mm4_norm' : {}, 'mm5_norm' : {}, 'mm6_norm' : {}, 'mm7_norm' : {},
                  'BS1' : {}, 'mm1_BS1' : {}, 'mm2_BS1' : {}, 'mm3_BS1' : {}, 'mm4_BS1' : {}, 'mm5_BS1' : {}, 'mm6_BS1' : {}, 'mm7_BS1' : {},
                 'BS2' : {}, 'mm1_BS2' : {}, 'mm2_BS2' : {}, 'mm3_BS2' : {}, 'mm4_BS2' : {}, 'mm5_BS2' : {}, 'mm6_BS2' : {}, 'mm7_BS2' : {},
                 'BA1' : {}, 'mm1_BA1' : {}, 'mm2_BA1' : {}, 'mm3_BA1' : {}, 'mm4_BA1' : {}, 'mm5_BA1' : {}, 'mm6_BA1' : {},
                 'BA2' : {}, 'mm1_BA2' : {}, 'mm2_BA2' : {}, 'mm3_BA2' : {}, 'mm4_BA2' : {}, 'mm5_BA2' : {},
                 'ALS' : {}, 'mm1_ALS' : {}, 'mm2_ALS' : {}, 'mm3_ALS' : {}, 'mm4_ALS' : {}, 'mm5_ALS' : {}, 'mm6_ALS' : {},
                 'ALA' : {}, 'mm1_ALA' : {}, 'mm2_ALA' : {}, 'mm3_ALA' : {}, 'mm4_ALA' : {}, 'mm5_ALA' : {}}
    for x in databases:  # all output files are analyzed
        if x in ['5SSs', '3SSs']: # defining the positions for the splice sites (NOTE: POSITION "98" EQUALS "-3")
            if settings['positions'][x] == 'range':
                positions = []
                for a in list(range(51, 151)):
                    positions.append(str(a))
            elif settings['positions'][x] == 'cumulative':
                positions = []
                if x == '5SSs':
                    for b in list(range(91, 126)):
                        positions.append(str(b))
                elif x == '3SSs':
                    for b in list(range(91, 100)):
                        positions.append(str(b))
            elif settings['positions'][x] == 'overlapping':
                file_U1s = ''.join([directory, '/input/query.fa']) # determining the original U1 length
                library = open(file_U1s, 'r').read().splitlines()
                U1_length = len(library[1])
                del library
                del file_U1s
                positions = {'norm' : list(range(100 - U1_length + 2, 101)), # each register has a different length
                             'BS1' : list(range(100 - U1_length + 1, 101)),
                             'BS2' : list(range(100 - U1_length + 0, 101)),
                             'BA1' : list(range(100 - U1_length + 3, 101)),
                             'BA2' : list(range(100 - U1_length + 4, 101)),
                             'ALS' : list(range(100 - U1_length + 1, 101)),
                             'ALA' : list(range(100 - U1_length + 3, 101))}
            else:
                positions = []
                for a in settings['positions'][x].split(','):
                    if int(a) <= 0:
                        positions.append(str(int(a) + 101))
                    else:
                        positions.append(str(int(a) + 100))
        results = {} # storage for counts
        results_registers = {} # storage for combinations of registers
        if settings['print'][x] == 'yes': # storage for targets to be printed
            results_print = {}
        # COLLECTING THE HITS FOR EACH REGISTER
        for y in registers:
            print('analyzing:', x, y)
            if 'norm' in y:
                register = 'norm'
            elif 'BS1' in y:
                register = 'BS1'
            elif 'BS2' in y:
                register = 'BS2'
            elif 'BA1' in y:
                register = 'BA1'
            elif 'BA2' in y:
                register = 'BA2'
            elif 'ALS' in y:
                register = 'ALS'
            elif 'ALA' in y:
                register = 'ALA'
            results_filtered = {}  # the dictionary prevents duplicate hits from different registers
            results_BLAST = open(''.join([file_BLAST, '_', y, '_', x, '.txt']), 'r').read().splitlines()
            for a in list(range(len(results_BLAST))):
                # finding the results for the single query
                if 'Query=' in results_BLAST[a]:
                    ID = results_BLAST[a][7:len(results_BLAST[a])]
                    if y != 'norm':  # adjusting the ID for alternative registers (removing additional label)
                        ID = '_'.join(ID.split('_')[0:len(ID.split('_')) - 1])
                    if ID not in results_filtered:
                        results_filtered[ID] = {}
                    # navigating the best scoring matches
                    c = a + 1
                    while c < len(results_BLAST) and 'Query=' not in results_BLAST[c]:
                        if '>' in results_BLAST[c]:
                            # each match must be analyzed separately, as it might include multiple matches
                            e = c + 1
                            while e < len(results_BLAST) and '>' not in results_BLAST[e]:
                                # analysis for exons and introns
                                if 'Sbjct' in results_BLAST[e] and x in ['exons', 'introns']:
                                    hit = '_'.join(
                                        [results_BLAST[c].split('>')[1].split(' ')[0], results_BLAST[e].split('  ')[1]])
                                    if hit not in results_filtered[ID]:  # avoiding duplicate hits
                                        results_filtered[ID][hit] = ''
                                # analysis for 5SSs and 3SSs
                                elif 'Sbjct' in results_BLAST[e] and x in ['5SSs', '3SSs']:
                                    # the counting changes based on the "overlapping" option
                                    if settings['positions'][x] != 'overlapping':
                                        position = results_BLAST[e].split('  ')[1]
                                        if position in positions:
                                            if position not in results_filtered[ID]:
                                                results_filtered[ID][position] = {}
                                            hit = '_'.join([results_BLAST[c].split('>')[1].split(' ')[0],
                                                            results_BLAST[e].split('  ')[1]])
                                            if hit not in results_filtered[ID][position]:  # avoiding duplicate hits
                                                results_filtered[ID][position][hit] = ''
                                    else:
                                        position = results_BLAST[e].split('  ')[1]
                                        if int(position) in positions[register]:
                                            if position not in results_filtered[ID]:
                                                results_filtered[ID][position] = {}
                                            hit = '_'.join([results_BLAST[c].split('>')[1].split(' ')[0],
                                                            results_BLAST[e].split('  ')[1]])
                                            if hit not in results_filtered[ID][position]:  # avoiding duplicate hits
                                                results_filtered[ID][position][hit] = ''
                                e = e + 1
                        c = c + 1
            del results_BLAST
            registers[y] = results_filtered
            del results_filtered
        # COUNGTING THE HITS FOR EACH COMBINATION OF REGISTERS
        registers_combinations = [['norm'], ['mm1_norm'], ['mm2_norm'], ['mm3_norm'], ['mm4_norm'], ['mm5_norm'], ['mm6_norm'], ['mm7_norm'],
                                  ['BS1'], ['mm1_BS1'], ['mm2_BS1'], ['mm3_BS1'], ['mm4_BS1'], ['mm5_BS1'], ['mm6_BS1'], ['mm7_BS1'],
                                  ['BS2'], ['mm1_BS2'], ['mm2_BS2'], ['mm3_BS2'], ['mm4_BS2'], ['mm5_BS2'], ['mm6_BS2'], ['mm7_BS2'],
                                  ['BA1'], ['mm1_BA1'], ['mm2_BA1'], ['mm3_BA1'], ['mm4_BA1'], ['mm5_BA1'], ['mm6_BA1'],
                                  ['BA2'], ['mm1_BA2'], ['mm2_BA2'], ['mm3_BA2'], ['mm4_BA2'], ['mm5_BA2'],
                                  ['ALS'], ['mm1_ALS'], ['mm2_ALS'], ['mm3_ALS'], ['mm4_ALS'], ['mm5_ALS'], ['mm6_ALS'],
                                  ['ALA'], ['mm1_ALA'], ['mm2_ALA'], ['mm3_ALA'], ['mm4_ALA'], ['mm5_ALA'],
                                  ['norm',
                                   'BS1',
                                   'BS2'],
                                  ['norm', 'mm1_norm',
                                   'BS1', 'mm1_BS1',
                                   'BS2', 'mm1_BS2',
                                   'BA1',
                                   'ALS'],
                                  ['norm', 'mm1_norm', 'mm2_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2',
                                   'BA1', 'mm1_BA1',
                                   'BA2',
                                   'ALS', 'mm1_ALS',
                                   'ALA'],
                                  ['norm', 'mm1_norm', 'mm2_norm', 'mm3_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1', 'mm3_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2', 'mm3_BS2',
                                   'BA1', 'mm1_BA1', 'mm2_BA1',
                                   'BA2', 'mm1_BA2',
                                   'ALS', 'mm1_ALS', 'mm2_ALS',
                                   'ALA', 'mm1_ALA'],
                                  ['norm', 'mm1_norm', 'mm2_norm', 'mm3_norm', 'mm4_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1', 'mm3_BS1', 'mm4_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2', 'mm3_BS2', 'mm4_BS2',
                                   'BA1', 'mm1_BA1', 'mm2_BA1', 'mm3_BA1',
                                   'BA2', 'mm1_BA2', 'mm2_BA2',
                                   'ALS', 'mm1_ALS', 'mm2_ALS', 'mm3_ALS',
                                   'ALA', 'mm1_ALA', 'mm2_ALA'],
                                  ['norm', 'mm1_norm', 'mm2_norm', 'mm3_norm', 'mm4_norm', 'mm5_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1', 'mm3_BS1', 'mm4_BS1', 'mm5_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2', 'mm3_BS2', 'mm4_BS2', 'mm5_BS2',
                                   'BA1', 'mm1_BA1', 'mm2_BA1', 'mm3_BA1', 'mm4_BA1',
                                   'BA2', 'mm1_BA2', 'mm2_BA2', 'mm3_BA2',
                                   'ALS', 'mm1_ALS', 'mm2_ALS', 'mm3_ALS', 'mm4_ALS',
                                   'ALA', 'mm1_ALA', 'mm2_ALA', 'mm3_ALA'],
                                  ['norm', 'mm1_norm', 'mm2_norm', 'mm3_norm', 'mm4_norm', 'mm5_norm', 'mm6_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1', 'mm3_BS1', 'mm4_BS1', 'mm5_BS1', 'mm6_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2', 'mm3_BS2', 'mm4_BS2', 'mm5_BS2', 'mm6_BS2',
                                   'BA1', 'mm1_BA1', 'mm2_BA1', 'mm3_BA1', 'mm4_BA1', 'mm5_BA1',
                                   'BA2', 'mm1_BA2', 'mm2_BA2', 'mm3_BA2', 'mm4_BA2',
                                   'ALS', 'mm1_ALS', 'mm2_ALS', 'mm3_ALS', 'mm4_ALS', 'mm5_ALS',
                                   'ALA', 'mm1_ALA', 'mm2_ALA', 'mm3_ALA', 'mm4_ALA'],
                                  ['norm', 'mm1_norm', 'mm2_norm', 'mm3_norm', 'mm4_norm', 'mm5_norm', 'mm6_norm', 'mm7_norm',
                                   'BS1', 'mm1_BS1', 'mm2_BS1', 'mm3_BS1', 'mm4_BS1', 'mm5_BS1', 'mm6_BS1', 'mm7_BS1',
                                   'BS2', 'mm1_BS2', 'mm2_BS2', 'mm3_BS2', 'mm4_BS2', 'mm5_BS2', 'mm6_BS2', 'mm7_BS2',
                                   'BA1', 'mm1_BA1', 'mm2_BA1', 'mm3_BA1', 'mm4_BA1', 'mm5_BA1', 'mm6_BA1',
                                   'BA2', 'mm1_BA2', 'mm2_BA2', 'mm3_BA2', 'mm4_BA2', 'mm5_BA2',
                                   'ALS', 'mm1_ALS', 'mm2_ALS', 'mm3_ALS', 'mm4_ALS', 'mm5_ALS', 'mm6_ALS',
                                   'ALA', 'mm1_ALA', 'mm2_ALA', 'mm3_ALA', 'mm4_ALA', 'mm5_ALA']]
        for y in registers_combinations:  # counting for each combination
            print('counting:', x, '+'.join(y))
            results_registers['+'.join(y)] = []
            for z in y:
                # establishing the register
                if 'norm' in z:
                    register = 'norm'
                elif 'BS1' in z:
                    register = 'BS1'
                elif 'BS2' in z:
                    register = 'BS2'
                elif 'BA1' in z:
                    register = 'BA1'
                elif 'BA2' in z:
                    register = 'BA2'
                elif 'ALS' in z:
                    register = 'ALS'
                elif 'ALA' in z:
                    register = 'ALA'
                # translating/annotating the register
                if registers[z] != {}:
                    z_adjusted = z
                    z_adjusted = z_adjusted.split('_')
                    if len(z_adjusted) == 1:
                        if z_adjusted[0] == 'norm':
                            z_adjusted = 'full'
                        else:
                            z_adjusted = z_adjusted[0]
                    else:
                        if z_adjusted[1] == 'norm':
                            z_adjusted[1] = 'full'
                        z_adjusted = '.'.join([z_adjusted[1], z_adjusted[0]])
                    results_registers['+'.join(y)].append(z_adjusted)
            for z in y:
                if registers[z] == {}:
                    results_registers['+'.join(y)] = []
            # collecting the hits from different registers
            if results_registers['+'.join(y)] != []:
                results_filtered = {}  # the dictionary prevents duplicate hits from different registers
                for a in y:
                    for b in registers[a]:
                        if x in ['exons', 'introns']:
                            if b not in results_filtered:
                                results_filtered[b] = registers[a][b]
                            else:
                                for e in registers[a][b]:
                                    if e not in results_filtered[b]:
                                        results_filtered[b][e] = ''
                        elif x in ['5SSs', '3SSs']:
                            if b not in results_filtered:
                                results_filtered[b] = registers[a][b]
                            else:
                                for e in registers[a][b]:
                                    if e not in results_filtered[b]:
                                        results_filtered[b][e] = registers[a][b][e]
                                    else:
                                        for g in registers[a][b][e]:
                                            if g not in results_filtered[b][e]:
                                                results_filtered[b][e][g] = ''
            if results_registers['+'.join(y)] != [] and x in ['exons', 'introns']:
                # storing targets for printing
                if settings['print'][x] == 'yes':
                    for b in results_filtered:
                        if b not in results_print:
                            results_print[b] = {}
                        for c in results_filtered[b]:
                            c = '_'.join([c.split('_')[0], c.split('_')[1], c.split('_')[2], c.split('_')[4]])
                            if len(y) == 1:
                                if c not in results_print[b]:
                                    results_print[b][c] = ['+'.join(results_registers['+'.join(y)])]
                                else:
                                    if '+'.join(results_registers['+'.join(y)]) not in results_print[b][c]:
                                        results_print[b][c].append('+'.join(results_registers['+'.join(y)]))
                # counting the matches
                results_counts = {}
                for a in results_filtered:
                    count = 0
                    for b in results_filtered[a]:
                        count += int(b.split('_')[3])
                    results_counts[a] = str(count)
                del results_filtered
                # collecting for later printing
                for a in results_counts:
                    if a not in results:
                        results[a] = {}
                        results[a]['+'.join(y)] = results_counts[a]
                    else:
                        results[a]['+'.join(y)] = results_counts[a]
            elif results_registers['+'.join(y)] != [] and x in ['5SSs', '3SSs']:
                # storing targets for printing
                if settings['print'][x] == 'yes':
                    for b in results_filtered:
                        if b not in results_print:
                            results_print[b] = {}
                        for c in results_filtered[b]:
                            for d in results_filtered[b][c]:
                                d = '_'.join([d.split('_')[0], d.split('_')[1], d.split('_')[2], d.split('_')[4]])
                                if len(y) == 1:
                                    if d not in results_print[b]:
                                        results_print[b][d] = ['+'.join(results_registers['+'.join(y)])]
                                    else:
                                        if '+'.join(results_registers['+'.join(y)]) not in results_print[b][d]:
                                            results_print[b][d].append('+'.join(results_registers['+'.join(y)]))
                # counting the matches
                # different ways if "cumulative" or "overlapping" option is "on" or "off"
                if settings['positions'][x] == 'cumulative' or settings['positions'][x] == 'overlapping':
                    results_counts = {}
                    for b in results_filtered:
                        count = 0
                        results_counts[b] = {}
                        results_counts[b][settings['positions'][x]] = {}
                        for c in results_filtered[b]:
                            for d in results_filtered[b][c]:
                                count += int(d.split('_')[3])
                        results_counts[b][settings['positions'][x]] = str(count)
                    del results_filtered
                else:
                    results_counts = {}
                    for b in results_filtered:
                        results_counts[b] = {}
                        for c in positions:
                            results_counts[b][c] = '0'
                        for c in results_filtered[b]:
                            count = 0
                            for d in results_filtered[b][c]:
                                count += int(d.split('_')[3])
                            results_counts[b][c] = str(count)
                    del results_filtered
                # collecting for later printing
                for a in results_counts:
                    if a not in results:
                        results[a] = {}
                    for b in results_counts[a]:
                        if b not in results[a]:
                            results[a][b] = {}
                        results[a][b]['+'.join(y)] = results_counts[a][b]
        # output
        with open(''.join([file_output, 'counts_', x, '.txt']), 'x') as output:
            if x in ['exons', 'introns']:
                output.write('\n')
                output.write('ID')
                for a in results_registers:
                    if results_registers[a] != []:
                        output.write('\t')
                        output.write('+'.join(results_registers[a]))
                for a in results:
                    output.write('\n')
                    output.write(a)
                    for b in results_registers:
                        if results_registers[b] != []:
                            output.write('\t')
                            output.write(results[a][b])
            elif x in ['5SSs', '3SSs']:
                output.write('\n')
                output.write('ID')
                if settings['positions'][x] not in ['overlapping', 'cumulative']:
                    output.write('\t')
                    output.write('position')
                for a in results_registers:
                    if results_registers[a] != []:
                        output.write('\t')
                        output.write('+'.join(results_registers[a]))
                if settings['positions'][x] in ['overlapping', 'cumulative']:
                    for b in results:
                        output.write('\n')
                        output.write(b)
                        for c in results_registers:
                            if results_registers[c] != []:
                                output.write('\t')
                                output.write(results[b][settings['positions'][x]][c])
                else:
                    for b in results:
                        for c in positions:
                            output.write('\n')
                            output.write(b)
                            output.write('\t')
                            if int(c) <= 100:
                                output.write(str(int(c) - 101))
                            else:
                                output.write(str(int(c) - 100))
                            for d in results_registers:
                                if results_registers[d] != []:
                                    output.write('\t')
                                    output.write(results[b][c][d])
        # printing targets
        if settings['print'][x] == 'yes':
            with open(''.join([file_output, 'targets_', x, '.txt']), 'x') as output:
                if x in ['exons', 'introns']:
                    output.write('\n')
                    output.write('ID')
                    output.write('\t')
                    output.write('coordinates')
                    output.write('\t')
                    if x == 'exons':
                        output.write('position within the exon')
                    if x == 'introns':
                        output.write('position within the intron')
                    output.write('\t')
                    output.write('register category')
                    for b in results_print:
                        for c in results_print[b]:
                            output.write('\n')
                            output.write(b)
                            output.write('\t')
                            output.write('_'.join([c.split('_')[0], c.split('_')[1], c.split('_')[2]]))
                            output.write('\t')
                            output.write(c.split('_')[3])
                            output.write('\t')
                            output.write('; '.join(results_print[b][c]))
                elif x in ['5SSs', '3SSs']:
                    output.write('\n')
                    output.write('ID')
                    output.write('\t')
                    output.write('exon coordinates')
                    output.write('\t')
                    if x == '5SSs':
                        output.write('position around the 5-SS')
                    if x == '3SSs':
                        output.write('position around the 3-SS')
                    output.write('\t')
                    output.write('register category')
                    for b in results_print:
                        for c in results_print[b]:
                            output.write('\n')
                            output.write(b)
                            output.write('\t')
                            output.write('_'.join([c.split('_')[0], c.split('_')[1], c.split('_')[2]]))
                            output.write('\t')
                            if int(c.split('_')[3]) <= 100:
                                output.write(str(int(c.split('_')[3]) - 101))
                            else:
                                output.write(str(int(c.split('_')[3]) - 100))
                            output.write('\t')
                            output.write('; '.join(results_print[b][c]))
ANALYSIS()