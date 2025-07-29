# PlantHGT-Review
Associated files and methodological details for the analysis presented in:

#### Fridrich, A. &amp; Irwin, N. A. T. 2025. Cross-kingdom gene transfer as a driver of land plant evolution.

## Methods
We aimed approximate the proportion of the plant genes which trace their ancestry to prokaryotic sources since the streptophyte ancestor. To do this, we blasted each gene from a series of representative plant genomes against the NCBI non-redundant database, excluded streptophtye hits, and classified the origin of the gene based on the majority rule taxonomy of the top 25 hits.

#### 1. Collect genome-predicted proteomes for representative plant taxa and take only the longest protein prediction per gene
Download genomes:
| Genome  | Source | Link |
| ------------- | ------------- | ------------- |
| Arabidopsis thaliana TAIR10  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Zea mays RefGenv4.0  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Populus trichocarpa v4.1  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Azolla filiculoides v1.1  | FernBase  | https://fernbase.org/  |
| Selaginella moellendorffii v1.0  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Physcomitrium patens v6.1  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Marchantia polymorpha v7.1  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Antheroceros angustus v1.0 | Dryad  | https://datadryad.org/dataset/doi:10.5061/dryad.msbcc2ftv  |

Reformat fasta file headers and select the longest protein for gene when necessary
```
# Arabidopsis
fasta = open('Arabidopsis_thaliana.TAIR10.pep.all.fa','r').read().split('>')[1:]
out = open('Arabidopsis_thaliana.TAIR10.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID         
for seq in seq_d:
    out.write('>3702.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Zea
fasta = open('Zmays_493_RefGen_V4.protein.fa','r').read().split('>')[1:]
out = open('Zmays_493_RefGen_V4.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('_')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                 
for seq in seq_d:
    out.write('>4577.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Populus
fasta = open('Ptrichocarpa_533_v4.1.protein.fa','r').read().split('>')[1:]
out = open('Ptrichocarpa_533_v4.1.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[1]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                 
for seq in seq_d:
    out.write('>3694.'+seq+'\n'+seq_d[seq]+'\n')
out.close()


# Azolla
fasta = open('Azolla_filiculoides.protein.highconfidence_v1.1.fasta','r').read().split('>')[1:]
out = open('Azolla_filiculoides.protein.highconfidence_v1.1.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[1]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                
for seq in seq_d:
    out.write('>84609.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Selaginella
fasta = open('Smoellendorffii_91_protein.fa','r').read().split('>')[1:]
out = open('Smoellendorffii_91_protein.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                 
for seq in seq_d:
    out.write('>88036.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Physcomitrium
fasta = open('Ppatens_870_v6.1.protein.fa','r').read().split('>')[1:]
out = open('Ppatens_870_v6.1.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                 
for seq in seq_d:
    out.write('>3218.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Marchantia
fasta = open('1480154.MpTak1_v7.1.protein.fa','r').read().split('>')[1:]
out = open('1480154.MpTak1_v7.1.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('.')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                
for seq in seq_d:
    out.write('>1480154.'+seq+'\n'+seq_d[seq]+'\n')
out.close()

# Anthorocerus
fasta = open('Anthoceros.angustus.coding.gene.protein.fa','r').read().split('>')[1:]
out = open('Anthoceros.angustus.coding.gene.protein.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('\t')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                
for seq in seq_d:
    out.write('>48387.'+seq+'\n'+seq_d[seq]+'\n')
out.close()
```
#### 2. Search for non-streptophyte homologs for each plant protein
BLAST all proteomes against NCBI NR, excluding streptophtyes streptophyte (NCBI TaxaID: 33090) hits, to identify closest homologs outside of streptophytes. Use DIAMOND BLASTP (sensitive mode, E < 1E-5, percent ID > 25, max target sequences = 250, query coverage > 50%). Synthetic constructs were also excluded from the blast results (NCBI TaxaID: 2787854).
```
# search each proteome against NCBI nr using Diamond BLASTP
for i in *rename.fa; do diamond blastp --query $i --db /data/NCBI_nr/nr/nr.dmnd -o $i.nr.sensitive.blastp --threads 32 --taxon-exclude 33090,2787854 --max-target-seqs 250 --id 25 --query-cover 50 --evalue 1e-5 --sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids; done
# cut out the resulting hits
for i in *rename.fa; do cut -f 2 $i.nr.sensitive.blastp > $i.nr.sensitive.blastp.hits; done
```
#### 3. Classify the taxonomic origin of each protein based on its BLAST hits
```
# load modules
from ete3 import NCBITaxa
from collections import Counter
from glob import glob
import subprocess
ncbi = NCBITaxa()

# for each plant proteomes...
for fname in glob('*rename.fa'):
    # create dictionaries to record BLAST hits and taxonomic information
    gene_to_hit = {} # protein: {hit:[evalue,id]}
    acc_to_taxid = {} # accession:taxid
    # read the blast output file
    blast = open(fname+'.nr.sensitive.blastp','r').readlines()
    # record each blast hit and its taxonomy to a maximum of 25 per query protein
    for h in blast:
        if '.' in h.split('\t')[1]: # exclude pdb
            acc_to_taxid[h.split('\t')[1]] = [h.split('\t')[-1].strip()]
            try:
                if len(gene_to_hit[h.split('\t')[0]]) == 25:
                    pass
                else:
                    gene_to_hit[h.split('\t')[0]][h.split('\t')[1]] = [h.split('\t')[10],h.split('\t')[2]]
            except:
                gene_to_hit[h.split('\t')[0]] = {h.split('\t')[1]:[h.split('\t')[10],h.split('\t')[2]]}
    # make record keeping dictionaries for the classifications
    gene_taxonomy = {} # gene: [best_domain, n_best_domain, best_hit, pid, evalue]
    # for each query protein, check the taxonomic domain of the hits (e.g., eukaryota, prokaryota, other), assign the domain based on majority rule
    for i in gene_to_hit:
        # get a list of all of the best hits
        accessions = list(gene_to_hit[i].keys())
        # for each hit...
        for h in accessions:
            # assign a taxonomic domain
            domain = 'NA'
            taxid = acc_to_taxid[h][0]
            if taxid != '':
                try:
                    # 2759 = eukaryota
                    if 2759 in ncbi.get_lineage(taxid):
                        domain = 'Eukaryota'
                    # 2 = bacteria, 2157 = archaea
                    elif (2 in ncbi.get_lineage(taxid)) or (2157 in ncbi.get_lineage(taxid)):
                        domain = 'Prokaryota'
                    # usually would represent viruses
                    else:
                        domain = 'Other'
                except:
                    pass
            else:
                pass
            acc_to_taxid[h].append(domain)
        # check frequency of taxonomic domains in top hits
        domains = []
        for hit in gene_to_hit[i]:
            try:
                domains.append(acc_to_taxid[hit][1])
            except:
                pass
        domains = dict(Counter(domains))
        # get the most frequently observed domain
        maj_domain = max(domains, key=domains.get)
        # get best hit with majority domain
        best_hit = ''
        for hit in gene_to_hit[i]:
            try:
                if (acc_to_taxid[hit][1] == maj_domain) and (best_hit == ''):
                    best_hit = hit
            except:
                best_hit = 'NA'
        # record info for query protein
        try:
            gene_taxonomy[i] = [maj_domain,domains[maj_domain],best_hit,gene_to_hit[i][best_hit][0],gene_to_hit[i][best_hit][1]]
        except:
            gene_taxonomy[i] = ['NA','NA','NA','NA','NA']

    # output results
    out = open(fname+'.nr.blastp.taxonomy','w')
    out.write('gene\tdomain\tdomain_n\tbest_hit\tevalue\tid\n')
    for g in gene_taxonomy:
        out.write(g+'\t'+'\t'.join([str(i) for i in gene_taxonomy[g]])+'\n')
    out.close()
```
Add in proteins with no hits to create the final annotation files. Unclassify genes that had less than 5 hits to their taxonomic annotation.
```
from glob import glob

# for each plant proteome... 
for fname in glob('*rename.fa'):
    # make a dictionary to record gene info
    gene_d = {} # gene:line
    # load the original annotation file
    f1 = open(fname+'.nr.blastp.taxonomy','r').readlines()
    # make final output file - called *all.taxonomy
    out = open(fname+'.nr.blastp.all.taxonomy','w')
    out.write('taxa\t'+f1[0])
    # get info for all annotated genes
    for i in f1[1:]:
        gene_d[i.split('\t')[0]] = i
    # load fasta file with all of the query proteins and record query species taxa ID
    fasta = open(fname,'r').read().split('>')[1:]
    taxa = fasta[0].split('.')[0]
    # record info for annotated genes with greater than or equal to five hits to their majority domain
    for i in gene_d:
        if int(gene_d[i].split('\t')[2]) >= 5:
            out.write(taxa+'\t'+gene_d[i])
    # add in genes with less than five hits and record as unknown
    for i in fasta:
        gene = i.split('\n')[0]
        try:
            gene_d[gene]
            if int(gene_d[gene].split('\t')[2]) < 5:
                out.write(taxa+'\t'+gene+'\tUnknown\t'+'\t'.join(['NA']*4)+'\n')
        except:
            out.write(taxa+'\t'+gene+'\tUnknown\t'+'\t'.join(['NA']*4)+'\n')
    out.close()
```
#### 4. Plot the results
Plot the results using a stacked bar chart
```
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd
from collections import Counter
from glob import glob

data = {}
ids = ['Eukaryota','Prokaryota','Other','Unknown']
for fname in glob('*all.taxonomy'):
    df = pd.read_csv(fname,sep='\t')
    counts = Counter(df['domain'])
    taxa = list(set(df['taxa']))[0]
    data[taxa] = [counts[a] for a in ids]
df = pd.DataFrame.from_dict(data).transpose()
df.columns = ids
df['taxa'] = df.index
df = df.melt(id_vars=['taxa'])
pivot_df = df.pivot(index='taxa', columns='variable', values='value').fillna(0)
column_order = ['Eukaryota','Prokaryota', 'Other', 'NoHit']
pivot_df = pivot_df[column_order]
# x-axis order
pivot_df = pivot_df.loc[[3702,4577,3694,84609,49495,88036,3218,1480154,48387]] 
 
pivot_df.plot(kind='bar', stacked=True, colormap='tab10')

# plot and save
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

plt.savefig("PlantHGT.bar.pdf", format="pdf")
plt.close()
```
<img src="https://github.com/nickatirwin/PlantHGT-Review/blob/main/PlantHGT.bar.png" alt="PlantHGTbar" width="500" height="450">
