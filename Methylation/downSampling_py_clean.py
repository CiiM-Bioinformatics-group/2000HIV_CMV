#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd


phe = pd.read_csv(".", sep="\t")

phe = phe[["Record.Id", "Cohort", "SEX_BIRTH", "AGE", "BMI_BASELINE", "DNAm_SamplePlate", "CMV_IgG_Serology", "CMV_IgG_IU.mL", "ETHNICITY", "season_sin", "season_cos"]]

dis = phe[phe.Cohort=="Discovery"]
val = phe[phe.Cohort=="Validation"]

val.CMV_IgG_Serology.value_counts()


l = []
for random_seed in range(100):
    sampled_dis = dis[dis.CMV_IgG_Serology==1.0].sample(n=89, random_state=random_seed)
    ids = sampled_dis["Record.Id"].values
    ids = list(ids)
    l.append(ids)
df = pd.DataFrame(l, columns=[f'CMV+{i}' for i in range(89)])
df.index.name = "random_seed"

df.to_csv("DownSampling_Dis_IDs.csv", sep=",", index=True)


l = []
for random_seed in range(100):
    sampled_val = val[val.CMV_IgG_Serology==1.0].sample(n=28, random_state=random_seed)
    ids = sampled_val["Record.Id"].values
    ids = list(ids)
    l.append(ids)
df = pd.DataFrame(l, columns=[f'CMV+{i}' for i in range(28)])
df.index.name = "random_seed"
df.to_csv("DownSampling_Val_IDs.csv", sep=",", index=True)
