import pandas as pd
import collections
import numpy as np
import math as mt
import time

def calculate_proportion(sp, cluster_num_list):
    clusters_curr = sp.loc[:,"res.0.8"]
    unique_curr, counts_curr = np.unique(clusters_curr, return_counts=True)
    count_dict_curr = dict(zip(unique_curr, counts_curr))
    curr_total = sum(count_dict_curr.values())
    cluster_dict = {}
    for c_n in cluster_num_list:
        if c_n in count_dict_curr.keys():
            cluster_dict[c_n] = count_dict_curr[c_n] * 100/curr_total
        else:
            cluster_dict[c_n] = 0
    
    return cluster_dict


def randomized_test(clones, orig):
    clone_sample_randoms = {}
    for c_1 in clones:
        curr_clone_cells = orig.loc[df['v2.2'] == c_1]
        clone_size = len(curr_clone_cells.index)
        rep_num = int(mt.ceil(len(orig.index)/clone_size))
        curr_clone_random_percentage = {}
        for rep in range(rep_num):
            # Sample without replacement
            curr_sample = orig.sample(n=clone_size)
            curr_percentage_dict = calculate_proportion(curr_sample, cluster_ls)
            curr_clone_random_percentage[rep] = curr_percentage_dict

        clone_sample_randoms[c_1] = curr_clone_random_percentage
    
    return clone_sample_randoms


df = pd.read_table("meta.clone.clean.integrated.v1.v2.1.v2.2.txt", sep="\t")
v1 = df.loc[:,"v1.1"]
v2 = df.loc[:,"v2.1"]
v22 = df.loc[:, "v2.2"]
clusters = df.loc[:,"res.0.8"]

not_na_v1 = pd.notnull(v1)
not_na_v2 = pd.notnull(v2)
not_na_v22 = pd.notnull(v22)

v1_not_na = v1[not_na_v1]
v2_not_na = v2[not_na_v2]
v22_not_na = v22[not_na_v22]

unique, counts = np.unique(v22_not_na, return_counts=True)
count_dict = dict(zip(unique, counts))
grt_5_count_dict = {}

for key, value in count_dict.items():
    if value > 5:
        grt_5_count_dict[key] = value

clones = list(grt_5_count_dict.keys())
cluster_list, counts_cls = np.unique(clusters, return_counts=True)
cluster_dict = dict(zip(cluster_list, counts_cls))
cluster_ls = list(cluster_dict.keys())

replication_number = 50
time_vec = []
replication_dict = {}
for j in range(replication_number):
    print(j)
    start_time = time.time()
    curr_replicate_rslt = randomized_test(clones, df)
    end_time = time.time()
    replication_dict[j] = curr_replicate_rslt
    time_vec.append(end_time - start_time)
    print(end_time - start_time)

# Format: {clone id1: {0:[], 1:[], 3:[], 4:[], 6:[]}, clone id2: {0:[], 1:[], 3:[], 4:[], 6:[]}, ...}
rearrange_dict = {}
for k,v in replication_dict.items():
    for sk,sv in v.items():
        curr_c = sk
        curr_pct_dict = {}
        for ssk,ssv in sv.items():
            for sssk,sssv in ssv.items():
                if sssk not in curr_pct_dict.keys():
                    curr_pct_dict[sssk] = [sssv]
                else:
                    curr_pct_dict[sssk].append(sssv)
        #print(len(curr_pct_dict[0]))
        if curr_c not in rearrange_dict.keys():
            rearrange_dict[curr_c] = curr_pct_dict
        else:
            for pct_k,pct_v in curr_pct_dict.items():
                rearrange_dict[curr_c][pct_k].extend(pct_v)

p_val_grt_overall = {}
p_val_less_overall = {}
for key_1,val_1 in rearrange_dict.items():
    curr_cl = key_1
    grt_p = {}
    less_p = {}
    for ky,vl in val_1.items():
        grt_p[ky] = sum(i > clone_null_pct[curr_cl][ky] for i in vl)/len(vl)
        less_p[ky] = sum(j < clone_null_pct[curr_cl][ky] for j in vl)/len(vl)
    p_val_grt_overall[curr_cl] = grt_p
    p_val_less_overall[curr_cl] = less_p

grt_df = pd.DataFrame(p_val_grt_overall)
less_df = pd.DataFrame(p_val_less_overall)

clusters_all = {}
for k1,vl in rearrange_dict.items():
    for k2,vl2 in vl.items():
        if k2 in clusters_all.keys():
            clusters_all[k2].extend(vl2)
        else:
            clusters_all[k2] = vl2         

clster_df = pd.DataFrame(clusters_all)
null_df = pd.DataFrame(clone_null_pct)

null_df.to_csv("permutation_clean_null_v2_2_ca.txt", sep = "\t")
clster_df.to_csv("percentages_all_clusters_v2_2_ca.txt", sep = "\t")
grt_df.to_csv("p_value_hyper_v2_2_ca.txt", sep = "\t")
less_df.to_csv("p_value_hypo_v2_2_ca.txt", sep = "\t")

