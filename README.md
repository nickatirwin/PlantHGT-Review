# PlantHGT-Review
Associated files and methodological details for the analysis presented in:

#### Fridrich, A. &amp; Irwin, N. A. T. 2025. Cross-kingdom gene transfer as a driver of land plant evolution.

## Methods
We aimed to approximate the proportion of plant genes which trace their ancestry to prokaryotic sources since the last common streptophyte ancestor. To do this, we blasted each protein from a series of representative plant genomes against the NCBI non-redundant database, excluded streptophtye hits, and classified the origin of the gene based on the majority rule taxonomy of the top 25 hits. Only a single protein was taken per gene from each species.

#### 0. Requirements
Install dependencies with conda:
```
conda install bioconda::diamond
conda install conda-forge::ete3
```

#### 1. Collect genome-predicted proteomes for representative plant taxa and take only the longest protein prediction per gene
Download genomes:
| Genome  | Source | Link |
| ------------- | ------------- | ------------- |
| Arabidopsis thaliana TAIR10  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Zea mays RefGenv4.0  | Phytozome  | https://phytozome-next.jgi.doe.gov/  |
| Gnetum montanum v1.0  | TreeGenesDB  | https://treegenesdb.org/FTP/Genomes/Gnmo/v1.0/annotation/  |
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

# Gnetum
fasta = open('Gnmo.1_0.pep.fa','r').read().split('>')[1:]
out = open('Gnmo.1_0.pep.rename.fa','w')
# record all sequences - take the longest protein per gene
seq_d = {}
for seq in fasta:
    gene = seq.split('\n')[0]
    sequence = seq.split('\n',1)[1].replace('\n','').replace('*','')
    try:
        if len(seq_d[gene]) < len(sequence):
            seq_d[gene] = sequence
    except:
        seq_d[gene] = sequence
# output sequences with clean headers and an NCBI taxonomy ID                 
for seq in seq_d:
    out.write('>3381.'+seq+'\n'+seq_d[seq]+'\n')
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
Note: adjust NCBI nr database path
```
# search each proteome against NCBI nr using Diamond BLASTP
for i in *rename.fa; do diamond blastp --query $i --db /data/NCBI_nr/nr/nr.dmnd -o $i.nr.sensitive.blastp --threads 32 --taxon-exclude 33090,2787854 --max-target-seqs 250 --id 25 --query-cover 50 --evalue 1e-5 --sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids; done
# cut out the resulting hits
for i in *rename.fa; do cut -f 2 $i.nr.sensitive.blastp > $i.nr.sensitive.blastp.hits; done
```
#### 3. Classify the taxonomic origin of each protein based on its BLAST hits
```
from ete3 import NCBITaxa
from collections import Counter
from glob import glob
import subprocess
ncbi = NCBITaxa()


for fname in glob('*rename.fa'):
    # look at taxonomy of the best hits for every sequence
    # get the 25 best hits per gene
    gene_to_hit = {} # gene: {hit:[evalue,id]}
    acc_to_taxid = {} # acc:taxid
    blast = open(fname+'.nr.sensitive.blastp','r').readlines()
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
    # make record dictionaries
    gene_taxonomy = {} # gene: [best_domain, n_best_domain, perc_prok, best_hit, pid, evalue]
    # tranlsate to taxon ids and record domains
    for i in gene_to_hit:
        # get domains for the best hits
        accessions = list(gene_to_hit[i].keys())
        for h in accessions:
            # assign domains
            domain = 'NA'
            taxid = acc_to_taxid[h][0]
            if taxid != '':
                try:
                    if 2759 in ncbi.get_lineage(taxid):
                        domain = 'Eukaryota'
                    elif (2 in ncbi.get_lineage(taxid)) or (2157 in ncbi.get_lineage(taxid)):
                        domain = 'Prokaryota'
                    else:
                        domain = 'Other'
                except:
                    pass
            else:
                pass
            acc_to_taxid[h].append(domain)
        # check frequency of domains in top hits
        domains = []
        for hit in gene_to_hit[i]:
            try:
                domains.append(acc_to_taxid[hit][1])
            except:
                pass
        domains = dict(Counter(domains))
        maj_domain = max(domains, key=domains.get)
        # get the proportion of hits that were prokaryotic
        total = sum([int(domains[c]) for c in domains])
        try:
            prok = domains['Prokaryota']
        except:
            prok = 0
        perc_prok = round(float(prok/total)*100,2)
        # get best hit with majority domain
        best_hit = ''
        for hit in gene_to_hit[i]:
            try:
                if (acc_to_taxid[hit][1] == maj_domain) and (best_hit == ''):
                    best_hit = hit
            except:
                best_hit = 'NA'
        # record info for gene
        try:
            gene_taxonomy[i] = [maj_domain,domains[maj_domain],str(perc_prok),best_hit,gene_to_hit[i][best_hit][0],gene_to_hit[i][best_hit][1]]
        except:
            gene_taxonomy[i] = ['NA','NA','NA','NA','NA','NA']
    
    out = open(fname+'.nr.blastp.taxonomy','w')
    out.write('gene\tdomain\tdomain_n\tperc_prok\tbest_hit\tevalue\tid\n')
    for g in gene_taxonomy:
        out.write(g+'\t'+'\t'.join([str(i) for i in gene_taxonomy[g]])+'\n')
    out.close()
```
Extract proteins where the top 25 best hits (out of 250, minimum of 5) were prokaryotic
```
for i in *taxonomy; do awk '$2 == "Prokaryota" && $3 >= 5 && $4 == 100' $i > $i.filt; done
# extract data for plotting
wc -l *blastp.taxonomy.filt > Prokaryotic.genes.txt
awk '{print $2, $1}' Prokaryotic.genes.txt
```
<img src="https://github.com/nickatirwin/PlantHGT-Review/blob/main/PlantHGT.bar.png" alt="PlantHGTbar" width="500" height="450">
