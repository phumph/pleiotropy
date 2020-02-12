import csv
import pandas as pd
from collections import Counter
from glob import glob

def get_deletion_neighborhood(stringer):
    # returns a set of all single character deletions of a string (includes the string)
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])

def find_bc_match(b, dd):
    # finds the correct barcode by looking for overlap with the deletion network
    # this accounts for all single errors - substitutions, insertions, and deletions
    if b in dd:
        return dd[b]
    else:
        hits = set()
        for b_edit in get_deletion_neighborhood(b):
            if b_edit in dd:
                hits.add(dd[b_edit])
        if len(hits) == 1:
            return hits.pop()
        else:
            return None

def call_bc_match(f):
    try: 
        td = pd.read_csv(f)
        dbc_counter = Counter()
        for entry in td.loc[td['Type'].isin(['DBC', 'unknown'])].as_matrix(['BC', 'Reads']):
            check = find_bc_match(entry[0], del_dict)
            if check:
                dbc_counter[check] += entry[1]
        ebc_counter = Counter()
        for entry in td.loc[td['Type'].isin(['EBC', 'unknown'])].as_matrix(['BC', 'Reads']):
            check = find_bc_match(entry[0], edel_dict)
            if check:
                ebc_counter[check] += entry[1]
        dbcs, ebcs = [d for d in dbc_counter], [e for e in ebc_counter]
        return [';'.join(dbcs), ';'.join([str(dbc_counter[d]) for d in dbcs]), ';'.join(ebcs), ';'.join([str(ebc_counter[e]) for e in ebcs])]
    except (FileNotFoundError, IOError):
        print("Wrong file or file path", f)
        return ['', '', '', '']

output_file = 'bcs_extracted_compiled.csv'
bfa_file = 'PLT_all_BFA_bcs.csv'
# Reading in known barcodes
bcs_in_bfa = pd.read_csv(bfa_file)

# Generating single-deletiion neighborhoods for each barcode for error correction
del_dict = dict()
edel_dict = dict()
for bc in set(bcs_in_bfa['Diverse.BC']):
    del_dict.update({d: bc for d in get_deletion_neighborhood(bc)})
for bc in set(bcs_in_bfa['Environment.BC']):
    edel_dict.update({d: bc for d in get_deletion_neighborhood(bc)})

bc_extract_files = glob('bc_counts/*.csv')

with open(output_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['File', 'dbcs', 'dbc_counts', 'ebcs', 'ebc_counts'])
    c = 0
    for f in bc_extract_files:
        m = call_bc_match(f)
        if m[0] != '':
            c += 1
        writer.writerow([f.split('/')[-1].split('.csv')[0]] + m)
    print('Found a dbc for', c, 'out of', len(bc_extract_files))
    