import pandas as pd
import json

metaDF = pd.read_csv("samplesheet.tsv", sep = "\t", index_col="Sample")

df = pd.DataFrame(columns=["Sample", "total_reads","totalreads4map","mappedReads","map_perc"])
for s in metaDF.index:

    # read fastp json
    fpjs = "fastp/" + s + ".json"
    f = open(fpjs)
    fpd = json.load(f)

    # read map2HOMD json
    hojs = "map2HOMD/" + s + ".json"
    f = open(hojs)
    hopd = json.load(f)
    map_perc = hopd['general']['reads']['mapped']['0']*100/hopd['general']['reads']['total']


    mlist = [s,
              str(fpd['summary']['before_filtering']['total_reads']),
              str(hopd['general']['reads']['total']),
              str(hopd['general']['reads']['mapped']['0']),
              str(map_perc)
              ]
    result = "\t".join(mlist)
    df = df.append(pd.Series(mlist, index=df.columns), ignore_index=True)




df.to_csv("metrics.tsv", sep='\t', index=False)