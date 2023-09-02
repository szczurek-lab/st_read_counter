import pandas as pd
import os
import sys
def read_vcf(file_name):
    with open(file_name) as f:
        for line in f:
            if line[0] == '#' and line[1] != '#':
                col_names = line[1:].split()
                break
    return pd.read_csv(file_name, sep='\t', comment='#', header=None, names=col_names, dtype={'CHROM': str})


def read_ac(file_name):
    data = pd.read_csv(file_name, sep='\t', dtype={'refContig': str})
    data['refPos'] += 1
    return data


def filter_df(data, data_type, mutations):
    pos = [p[1] for p in mutations]
    if data_type == 'wes':
        data = data[data['POS'].isin(pos)]
    elif data_type == 'st':
        data = data[data['refPos'].isin(pos)]
    d = data.copy()
    d['save'] = [0] * len(d)
    for index, row in d.iterrows():
        if data_type == 'wes':
            if tuple(row[['CHROM', 'POS', 'REF', 'ALT']]) in mutations:
                d.at[index, 'save'] = 1
        elif data_type == 'st':
            if tuple(row[['refContig', 'refPos', 'refAllele', 'base']]) in mutations:
                d.at[index, 'save'] = 1
    return data.loc[d['save'] == 1]


def filter_st_spots(data, spots, bar_dict=None):
    if 'x' in spots[0]:
        spots = [bar_dict[spot] for spot in spots]
    return data[data['spot'].isin(spots)]


def get_bar_dict(barcodes_file):
    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None)
    bar_dict = {}
    bar_dict_reverse = {}
    for index, row in barcodes.iterrows():
        bar_dict[row[0]] = str(row[1]) + 'x' + str(row[2])
        bar_dict_reverse[str(row[1]) + 'x' + str(row[2])] = row[0]
    return bar_dict, bar_dict_reverse


def get_proper_spots(section, cells_file):
    #cells = pd.read_csv(cells_file,  header=0, names=['section', 'spot', 'n_cells'])
    cells = pd.read_csv(cells_file,sep=',')
    cells = cells.loc[cells['section']=='P'+section]
    return list(cells['coordinates'])


def find_pvals(wes_data):
    infos = wes_data.iloc[1]['INFO'].split(';')
    for i in range(len(infos)):
        if infos[i].split('=')[0] == 'SSF':
            return i
    return False


def get_mutation_sets_trspots(sections, wes_files, st_files, spots_file, bar_dict_reverse, 
                              search_string='Somatic', good_p_threshold=0.1):
    wm, cm, gpm, cg = [], [], [], []
    for k in range(len(sections)):
        section = sections[k]
        print('starting section', section)
        data_wes = read_vcf(wes_files[k])
        data_st = read_ac(st_files[k])
        spots = get_proper_spots(section, spots_file)                # you can comment these two lines 
        data_st = filter_st_spots(data_st, spots, bar_dict_reverse)  # if you don't need to select only "true" spots
        wes_mutations = set([tuple(row) for row in data_wes[['CHROM', 'POS', 'REF', 'ALT']].itertuples(index=False)])
        st_mutations = set([tuple(row) for row in data_st[['refContig', 'refPos', 'refAllele', 'base']].itertuples(index=False)])
        somatic_wes = data_wes.loc[data_wes['INFO'].str.contains(search_string)]
        somatic_mutations = set([tuple(row) for row in somatic_wes[['CHROM', 'POS', 'REF', 'ALT']].itertuples(index=False)])
        common_mutations = somatic_mutations.intersection(st_mutations)
        cut_wes = filter_df(somatic_wes, 'wes', common_mutations)
        cut_st = filter_df(data_st, 'st', common_mutations)

        good_p_mutations = set()
        i = find_pvals(cut_wes)
        for index, row in somatic_wes.iterrows():
            pval = float(row['INFO'].split(';')[i].split('=')[1])
            if pval < good_p_threshold:
                good_p_mutations.add(tuple(row[['CHROM', 'POS', 'REF', 'ALT']]))
        st_good_p = filter_df(cut_st, 'st', good_p_mutations)
        common_good_p = set([tuple(row) for row in st_good_p[['refContig', 'refPos', 'refAllele', 'base']].itertuples(index=False)])
        wm.append(wes_mutations)
        cm.append(common_mutations)
        gpm.append(good_p_mutations)
        cg.append(common_good_p)
    return wm, cm, gpm, cg


def write_wes(mutations, vcf_in_file, vcf_out_file):
    with open(vcf_in_file) as f, open(vcf_out_file, 'w') as g:
        for line in f:
            if line[0] == '#':
                g.write(line)
            else:
                spl = line.split('\t')
                mut = (spl[0], int(spl[1]), spl[3], spl[4])
                if mut in mutations:
                    g.write(line)


def write_wes_union(mutations, vcf_in_files, vcf_out_file):
    seen = set()
    with open(vcf_out_file, 'w') as g:
        first_file = True
        for in_file in vcf_in_files:
            with open(in_file) as f:
                for line in f:
                    if line[0] == '#':
                        if first_file:
                            g.write(line)
                    else:
                        first_file = False
                        spl = line.split('\t')
                        mut = (spl[0], int(spl[1]), spl[3], spl[4])
                        if mut in mutations and mut not in seen:
                            g.write(line)
                            seen.add(mut)

##### PARAMETERS
print("first")
good_p_threshold = 0.1
search_string = 'Somatic'  # for filtering WES mutations
barcodes_file = sys.argv[1]
sections = sys.argv[6:]  # the code should also work for less/more sections
print(sections)
spots_file = sys.argv[2]  # from Igor

wes_files = [sys.argv[3]+'w' + section.replace('.', '') + '_blood.vcf' for section in sections]
st_files = [sys.argv[4]+'vardict2_' + section + '.ac' for section in sections]

directory =  sys.argv[5]
dir_loose = directory + 'loose/'
dir_strict = directory + 'strict/'
#####


bar_dict, bar_dict_reverse = get_bar_dict(barcodes_file)
mutations_vardict = get_mutation_sets_trspots(sections, wes_files, st_files, spots_file, 
                                              bar_dict_reverse, search_string, good_p_threshold)

n = len(sections)
common = mutations_vardict[3]  # common mutations (in WES and ST) for each section; somatic, good p-value
u = common[0]
for i in range(1,n):
    u = u.union(common[i])

u = {mut:0 for mut in u}  # counting, in how many sections each mutation is present
for muts in common:
    for mut in muts:
        u[mut]+=1

sets = [set() for i in range(n)]
for mut in u:
    for k in range(n):
        if u[mut]>k:   # mutation present in at least k+1 sections
            sets[k].add(mut)

print('\nNumber of mutations present in at least')
for k in range(n):
    if k==0:
        print(f'   {k+1} section', len(sets[k]))
    else:
        print(f'   {k+1} sections', len(sets[k]))



if not os.path.isdir(directory):
    os.mkdir(directory)
if not os.path.isdir(dir_loose):
    os.mkdir(dir_loose)
if not os.path.isdir(dir_strict):
    os.mkdir(dir_strict)


# writing files with filtered mutations
for i in range(n):    # mutations present in at least 1,..,4 sections present in each section (not all of them need to be good*)
    section = sections[i]
    for k in range(n):
        out_file = dir_loose + str(k+1) + '_' + section + '.vcf'
        write_wes(sets[k], wes_files[i], out_file)


for k in range(n):    # mutations present in at least 1,..,4 sections
    out_file = directory + str(k+1) + '_all.vcf'
    write_wes_union(sets[k], wes_files, out_file)

for i in range(n):    # mutations present in at least 1,..,4 sections present in each section (only good*)
    section = sections[i]
    for k in range(n):
        out_file = dir_strict + str(k+1) + '_' + section + '.vcf'
        write_wes(common[i].intersection(sets[k]), wes_files[i], out_file)

# *good mutations: mutations that are somatic, common with ST, p_value < good_p_threshold
