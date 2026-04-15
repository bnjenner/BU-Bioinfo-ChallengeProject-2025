from Bio import Entrez
import pandas as pd
import math
import re
import time

def pull_acc_nums(search_term, count):
    iterations = math.ceil(count/10000)
    accessions = set()
    for i in range(iterations):
        handle = Entrez.esearch(db = "sra", term = search_term, retstart = i*10000, retmax = 10000, idtype = "acc")
        uids = Entrez.read(handle)['IdList']
        handle.close()
        handle = Entrez.efetch(db = "sra", id=",".join(uids), retmode = "text")
        record = handle.read().decode('utf-8')
        handle.close()
        for match in re.finditer(r'(SRR\d+)', record):
                accessions.add(match.group(1))
        #acc_list[i*10000:len(record['IdList'])] = record['IdList']
    print(len(accessions))

    return accessions

def pull_n_seqs(n, seq_in_lin):
     acc_nums = []
     artic_primers = []
     for s in seq_in_lin:
        s = s.strip()
        try:
            handle = Entrez.efetch(db = "sra", id = s, retmode = "text")
        except:
            continue
        record = handle.read().decode('utf-8').lower()
        handle.close()
        if 'artic v' in record or 'articv' in record:
             acc_nums.append(s)
             primers = re.search("artic ?v\\.?[0-9\\.]+", record).group()
             artic_primers.append(primers)
        if len(acc_nums) >= n:
             break
     return [acc_nums, artic_primers]



def main():
    Entrez.email = "markerte@bu.edu"
    Entrez.api_key = "cc5ef99a135dbc3960c3219fc70e8b951e08"
    home_dir = "C:\\Users\\exmar\\Documents\\BU_Docs\\Classes\\Fall_25\\challenge project"
    ncbi_virus = pd.read_csv(f"{home_dir}/ncbi_virus_data_VOCs.csv").sample(frac=1, random_state = 321).reset_index(drop=True)
    ncbi_virus = ncbi_virus[~ncbi_virus['SRA_Accession'].str.contains(',', na=False)]
    lineage_list = ncbi_virus["Pangolin"].unique().tolist()
    lineage_list = list(filter(lambda x: x not in ["BA.2", "BA.1", "B.1"], lineage_list))
    all_acc_nums = {'SRA_Accession':[], 'Primers': []}
    for lineage in lineage_list:
        seq_list = list(ncbi_virus[ncbi_virus["Pangolin"] == lineage]["SRA_Accession"])
        lin_acc_nums = pull_n_seqs(n=15, seq_in_lin=seq_list)
        all_acc_nums['SRA_Accession'].extend(lin_acc_nums[0])
        all_acc_nums['Primers'].extend(lin_acc_nums[1])
    to_keep = pd.DataFrame(all_acc_nums)
    final = pd.merge(to_keep, ncbi_virus, on = "SRA_Accession", how = "left")
    final["SRA_Accession"] = final["SRA_Accession"].str.split(",")
    final = final.explode('SRA_Accession').reset_index(drop=True)
    final.to_csv("ncbi_virus_w_primers.csv", index = False)
    final["SRA_Accession"].to_csv("ncbi_virus_acc.txt", index = False, header = False)
    # write SRA accessions to list to download them

if __name__ == "__main__":
    main()
