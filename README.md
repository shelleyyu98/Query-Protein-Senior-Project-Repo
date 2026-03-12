# SPTY2D1 Project Repository
---
# Lab 3


## Create a new working directory


Creating a new working directory for this analysis
```bash
mkdir ~/lab03-$MYGIT/SPTY2D1
````
Move into the working directory
```bash
cd ~/lab03-$MYGIT/SPTY2D1
```
Confirm your current location
```bash
pwd     # Prints working directory path
```


## Retrieve your Query Protein

Extracting the SPTY2D1 human protein sequence from the large FASTA
```bash
samtools faidx ~/lab03-$MYGIT/allprotein.fas 'Hsap|NP_919261.2|SPTY2D1'
```
Saving that extracted sequence into a new FASTA file inside this directory
```bash
samtools faidx ~/lab03-$MYGIT/allprotein.fas 'Hsap|NP_919261.2|SPTY2D1' > ~/lab03-$MYGIT/SPTY2D1/NP_919261.2.1.fa
```
Viewing the FASTA file to confirm sequence retrieval
```bash
less NP_919261.2.1.fa
```
List the file to make sure it exists
```bash
ls NP_919261.2.1.fa
```


## Perform the BLAST Search (Typical)

Performing a blast search using the query protein:
```bash
blastp -db ~/lab03-$MYGIT/allprotein.fas \
       -query ~/lab03-$MYGIT/SPTY2D1/NP_919261.2.1.fa \
       -outfmt 0 \
       -max_hsps 1 \
       -out ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.typical.out
```


## Perform BLAST Search (tabular output)

Running BLASTP and request tabular output for a more detailed and easier-to-process output of the same analysis
```bash
blastp -db ~/lab03-$MYGIT/allprotein.fas \
       -query ~/lab03-$MYGIT/SPTY2D1/NP_919261.2.1.fa \
       -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" \
       -max_hsps 1 \
       -out ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.out
```
Counting how many hits belong to Homo sapiens
```bash
grep -c Hsap ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.out
```
To print all Homo sapiens hits
```bash
grep Hsap ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.out
```

## Filtering the BLAST output for high-scoring putative homologs

To keep only lines with E-value < 1e-30
```bash
awk '{if ($6 < 1e-30) print $1 }' \
    ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.out \
    > ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out
```
Counts number of filtered homologs
```bash
wc -l ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out
```


## How many paralogs are in each species?

Number of hits in each species
```bash
cut -d '|' -f 1 ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out | sort | uniq -c
```

---
# Lab 4


## Create a new working directory

Create a new directory
```bash
mkdir ~/lab04-$MYGIT/SPTY2D1
````
Move into the working directory
```bash
cd ~/lab04-$MYGIT/SPTY2D1
```
Confirm current working directory
```bash
pwd     # Prints the directory path
```

## Retrieving the Homolog Sequences from Lab 3

Extracting the sequences matching filtered homolog seuqneces from Lab 3 output
```bash

seqkit grep --pattern-file ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out \
            ~/lab03-$MYGIT/allprotein.fas \
    | seqkit grep -v -p "carpio" \
    > ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas
```


## Perform a global multiple sequence alignment in muscle

Running global multiple sequence alignment with MUSCLE
```bash
muscle -align ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas \
       -output ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas
```

Visualizing the alignment using alv
```bash
alv -kli ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas | less -RS
```
Visualize the alignment using the "majority" option
```bash
alv -kli --majority ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas | less -RS
```
Plotting the alignment using Rscript
```bash
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R \
        ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas
```


## Get some information about the alignment

Calculating the width (length) of the alignment:
```bash
alignbuddy -al ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas
```
Calculating the length of the alignment after removing any column with gaps
```bash
alignbuddy -trm all ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas \
  | alignbuddy -al
```
Calculating the length of the alignment after removing any column with gaps
```bash
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas \
  | alignbuddy -al
```


## Calculate the Average Percent Identity

Calculating the average percent identity using t_coffee
```bash
t_coffee -other_pg seq_reformat \
         -in ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas \
         -output sim
```

Repeat calculating the average percent identity using alignbuddy
```bash
alignbuddy -pi ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas \
  | awk '(NR>2){ for(i=2;i<=NF;i++){ sum+=$i; num++ } } \
         END{ print(100*sum/num) }'
```
---
# Lab 5


## Create a new working directory

Create new directory
```bash
mkdir ~/lab05-$MYGIT/SPTY2D1
````
Move into the directory
```bash
cd ~/lab05-$MYGIT/SPTY2D1
```


## Building a Maximum Likelihood Tree Using IQ-TREE

Run an IQ-TREE with 1000 ultrafast bootstrap replicates using 2 threads
```bash
iqtree -s ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.fas \
       -bb 1000 \
       -nt 2 \
       -pre ~/lab05-$MYGIT/SPTY2D1/SPTY2D1
```


## Looking at the unrooted tree

Displaying the unrooted tree in terminal
```bash
nw_display ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.treefile
```

Plotting the unrooted tree using custom R script
```bash
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R \
        ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.treefile \
        ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.treefile.pdf \
        0.4 15
```


## Midpoint Rooting

Generates a midpoint-root the IQ-TREE as an output
```bash
gotree reroot midpoint \
       -i ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.treefile \
       -o ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile
```


## Visualizing and Exporting Rooted Phylogeny

Displaying a midpoint-rooted tree at the command line
```bash
nw_order -c n ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile \
  | nw_display -
```
Exporting the midpoint-rooted phylogram as a PDF
```bash
nw_order -c n ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile \
  | nw_display -s -w 1000 \
       -l 'font-size:10px;font-family:sans' \
       -i 'font-size:8px;font-family:sans' \
       -b 'visibility:hidden' \
       -I l -n -4 -W 4.0 -v 24 - \
  | convert -density 300 svg:- \
       ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile.pdf
```

Exporting the midpoint-rooted cladogram as a PDF
```bash
nw_order -c n ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile \
  | nw_topology - \
  | nw_display -s -w 1000 \
       -l 'font-size:10px;font-family:sans' \
       -i 'font-size:8px;font-family:sans' \
       -b 'visibility:hidden' \
       -S -I l -n -4 -W 4.0 -v 24 - \
  | convert -density 300 svg:- \
       ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.cladogram.pdf
```
---

# Lab 6


## Create a new working directory

```bash
mkdir ~/lab06-$MYGIT/SPTY2D1
cd ~/lab06-$MYGIT/SPTY2D1
```


## Viewing the Raw Unaligned Sequences

```bash
less ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas
```


## Run RPS-BLAST

```bash
rpsblast -query ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas \
  -db ~/data/Pfam/Pfam \
  -out ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out \
  -outfmt "6 qseqid qlen qstart qend evalue stitle" \
  -evalue .0000000001
```


## Plotting Predicted Pfam Domains on Phylogeny

Using final gene tree from Lab 5, generate script indicating:
* **`Rscript`**: The Rscript program allows you to run an R script from the command line, without opening up the R console.
* **`--vanilla`**: This command flag indicates to R not save or restore a workspace or previous settings.
* **`plotTreeAndDomains.r`**: The script that Dr. Rest wrote. You can take a look at it using less, for example.
* Now provide the names of three files, in order: tree file, rps-blast output, and the name of the output file that you want. For example:
*  **`~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile`**: Our midpoint rooted tree file for the SPTY2D1 proteins.
*  **`~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out`**: Our predicted pfram domains from rps-blast.  
* **`~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas`**: Unaligned gene sequences.
*  **`~/lab06-$MYGIT/SPTY2D1/SPTY2D1.tree.rps.pdf`**: The name of the output pdf file we will create. 
```bash
less ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile

Rscript --vanilla ~/lab06-$MYGIT/plotTreeAndDomains2.r \
  ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile \
  ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out \
  ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas \
  ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.tree.rps.pdf
```


## Looking at the predicted domains in more detail

Command to see the annotations in SPTY2D1.rps-blast.out
```bash
mamba activate base

mlr --inidx --ifs "\t" --opprint cat ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | tail -n +2 | less -S
```

Counting hits per gene
```bash
cut -f 1 ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | sort | uniq -c
```
Counting domain types
```bash
cut -f 6 ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | sort | uniq -c
```
Sorting by alignment length
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | sort -k2nr
```

Viewing the E-values of Protein Domain Hits

```bash
cut -f 1,5 -d $'\t' ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out
```


## Map RefGENBANK REFSEQ IDs to Uniprot

Extracting unique RefSeq accessions from your file
```bash
awk -F'|' '{print $2}' ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out \
  | sort -u > ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_ids.txt
```
Mapping to UniProt
```bash
awk 'NR==FNR{m[$1]=$2; next} {if($0 in m) print $0"\t"m[$0]; else print $0"\tMISSING"}' \
  ~/lab06-$MYGIT/refseq_uniprot.tsv \
  ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_ids.txt \
  > ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv

less ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv
```


## Retrieving Predicted Protein Structures (PDBs) from Alphafold or Swissprot

```bash
chmod +x ~/lab06-$MYGIT/afetch_by_uniprot.sh

~/lab06-$MYGIT/afetch_by_uniprot.sh \
  ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv \
  ~/lab06-$MYGIT/SPTY2D1
```


## Looking at one of the PDB Files

```bash
less ~/lab06-$MYGIT/SPTY2D1/XP_068416897.1__AF-P02177-F1-model_v4.pdb
```


## Running Charge Calculation on the Downloaded PDB Structures

Create and run the charge calculation script
```bash

#!/usr/bin/env bash
set -euo pipefail

PH=${PH:-6.5}  # allow overriding with: PH=7.0 ./script.sh

for pdb in *.pdb; do
  # skip already-protein-only files
  [[ "$pdb" == *.protein.pdb ]] && continue

  base="${pdb%.pdb}"                                  # e.g., XP_...-model_v4
  prot="${base}.protein.pdb"                           # protein-only PDB
  pqr="${base}.pqr"                                    # output PQR

  # 1) Remove heteroatoms (make protein-only PDB)
  pdb_delhetatm "$pdb" > "$prot"

  # 2) Protonate & assign charges at chosen pH
  pdb2pqr --ff=AMBER --with-ph="$PH" --drop-water --keep-chain "$prot" "$pqr" \
    > "${base}.pdb2pqr.log" 2>&1

  # 3) Sum charges (col 9 in PQR ATOM/HETATM lines)
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9} END{printf "%.3f", s}' "$pqr")

  printf "%-40s pH=%s  Z=%s\n" "$base" "$PH" "$Z"
done

echo
echo "Original PDB: $(ls ~/lab06-$MYGIT/SPTY2D1/*.pdb 2>/dev/null | grep -v '.protein.pdb' | wc -l), \
Protein-only PDB: $(ls ~/lab06-$MYGIT/SPTY2D1/*.protein.pdb 2>/dev/null | wc -l), \
PQR: $(ls ~/lab06-$MYGIT/SPTY2D1/*.pqr 2>/dev/null | wc -l)"
```


## Computing Net Charge for Each PQR

```bash
out=~/lab06-$MYGIT/SPTY2D1/net_charges.tsv
: > "$out"   # truncate/create

for f in ~/lab06-$MYGIT/SPTY2D1/*.pqr; do
  base=$(basename "$f")
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9} END{printf "%.3f", s}' "$f")
  printf "%s\tNetCharge=%s\n" "$base" "$Z" >> "$out"
done

less ~/lab06-$MYGIT/SPTY2D1/net_charges.tsv
```


## Visualizing the Structures with Charges (static PNGs)

```bash
mamba activate bio312
mamba install -y -c conda-forge pymol-open-source

python ~/lab06-$MYGIT/pqr_to_charges_pdb.py ~/lab06-$MYGIT/SPTY2D1

python3 ~/lab06-$MYGIT/render_charge_pngs.py ~/lab06-"$MYGIT"/SPTY2D1
```


## Helper Files

Provides help figuring out which images belong to which species, and which are aquatic and terrestrial.
First, go to the images directory then build the helper maps.
```bash
cd ~/lab06-"$MYGIT"/SPTY2D1

# 1) Build helper maps:
#    - abbreviation -> status (aquatic/terrestrial)
awk -F, 'NR>1 {print $2"\t"$5}' ~/lab06-"$MYGIT"/species_key.csv \
  | sort -u > abbr_status.tsv

#    - refseq -> abbreviation (from lab03 list like Hsap|NP_...|MB)
awk -F'|' '{print $2"\t"$1}' ~/lab03-"$MYGIT"/SPTY2D1/SPTY2D1.blastp.detail.filtered.out \
  | sort -u > refseq_abbr.tsv

# 2) Map refseq -> png filename (derives refseq from the part before "__")
ls *.png 2>/dev/null \
  | sed 's/\.png$//' \
  | awk -F'__' '{print $1"\t"$0".png"}' \
  | sort -u > refseq_png.tsv

# 3) Join to get refseq -> abbr -> status
#    (join by abbreviation, so flip refseq_abbr first)
awk -F'\t' '{print $2"\t"$1}' refseq_abbr.tsv | sort -u > abbr_refseq.tsv
join -t $'\t' -1 1 -2 1 abbr_refseq.tsv abbr_status.tsv > abbr_refseq_status.tsv
awk -F'\t' '{print $2"\t"$1"\t"$3}' abbr_refseq_status.tsv | sort -u > refseq_abbr_status.tsv
# Columns now: RefSeq \t Abbr \t Status

# 4) Attach PNG filename (join on RefSeq)
join -t $'\t' -1 1 -2 1 refseq_png.tsv refseq_abbr_status.tsv > png_refseq_abbr_status.tsv
# Columns: RefSeq \t PNG \t Abbr \t Status

# 5) Split into lists and make a pretty guide
awk -F'\t' '$4=="aquatic"{print $2}' png_refseq_abbr_status.tsv > aquatic_pngs.txt
awk -F'\t' '$4=="terrestrial"{print $2}' png_refseq_abbr_status.tsv > terrestrial_pngs.txt

# Human-readable table
echo
echo "PNG GUIDE (RefSeq | PNG | Abbr | Status):"
column -t -s $'\t' png_refseq_abbr_status.tsv | tee PNG_GUIDE.txt

echo
echo "Counts -> Aquatic: $(wc -l < aquatic_pngs.txt 2>/dev/null || echo 0), Terrestrial: $(wc -l < terrestrial_pngs.txt 2>/dev/null || echo 0)"

less PNG_GUIDE.txt
```


## Compare Aquatic vs Terrestrial Proteins

```bash
python3 - <<'EOF'
import os, re, pandas as pd
import matplotlib
matplotlib.use('Agg')     # headless
import matplotlib.pyplot as plt

# --- env & inputs ---
MYGIT = os.environ.get("MYGIT")
if not MYGIT:
    raise SystemExit("ERROR: $MYGIT is not set")

charges_path   = "net_charges.tsv"
species_key    = os.path.expanduser(f"~/lab06-{MYGIT}/species_key.csv")
lab03_list     = os.path.expanduser(f"~/lab03-{MYGIT}/SPTY2D1/SPTY2D1.blastp.detail.filtered.out")

# --- 1) read net_charges.tsv (File \t NetCharge=VAL) ---
rows = []
with open(charges_path) as fh:
    for line in fh:
        line=line.strip()
        if not line: continue
        try:
            fname, ztok = line.split('\t', 1)
        except ValueError:
            continue
        if not ztok.startswith("NetCharge="): continue
        try:
            z = float(ztok.split("=",1)[1])
        except ValueError:
            continue
        rows.append((fname, z))

charges = pd.DataFrame(rows, columns=["File","NetCharge"])

# --- 2) RefSeq accession from filename (before the __) ---
charges["RefSeq"] = charges["File"].str.split("__", n=1, expand=True)[0]

# --- 3) Build RefSeq -> species abbr map from lab03 list (Hsap|NP_...|MB) ---
acc2abbr = {}
with open(lab03_list) as fh:
    for line in fh:
        line=line.strip()
        if not line or line.startswith("#"): continue
        parts = line.split("|")
        if len(parts) >= 2:
            abbr, refseq = parts[0], parts[1]
            acc2abbr[refseq] = abbr

charges["Abbr"]   = charges["RefSeq"].map(acc2abbr)

# --- 4) Abbr -> status (aquatic/terrestrial) from species_key.csv ---
sk = pd.read_csv(species_key)  # expects headers: species_name,abbreviation,taxid,common_name,aquatic
abbr2status = dict(zip(sk["abbreviation"], sk["aquatic"]))
charges["Status"] = charges["Abbr"].map(abbr2status).fillna("unknown")

# --- 5) Save labeled table ---
charges.to_csv("net_charges_labeled.tsv", sep="\t", index=False)

# --- 6) Plot boxplot (aquatic vs terrestrial) ---
df = charges[charges["Status"].isin(["aquatic","terrestrial"])].copy()
if df.empty:
    print("WARN: no aquatic/terrestrial rows found after mapping; check species_key.csv and lab03 list.")
else:
    ax = df.boxplot(column="NetCharge", by="Status", grid=False)
    plt.title("Net charge by habitat")
    plt.suptitle("")
    plt.ylabel("Net charge (from .pqr at your chosen pH)")
    plt.tight_layout()
    plt.savefig("netcharge_boxplot.png", dpi=150)
    print("Saved: net_charges_labeled.tsv, netcharge_boxplot.png")

# quick counts
print("\nCounts by Status:")
print(charges["Status"].value_counts(dropna=False))
EOF
```

## Running DSSP Batch Summary

First, install the software
```bash
mamba activate bio312
mamba install -y -c conda-forge mdtraj pandas matplotlib
```
Complete the analysis of secondary structure
```bash
python3 ~/lab06-$MYGIT/dssp_batch_summary_mdtraj.py \
  --pdb-dir   ~/lab06-$MYGIT/SPTY2D1 \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out \
  --out-csv   ~/lab06-$MYGIT/SPTY2D1/dssp_summary.csv \
  --plots
```
---

# Lab 7


## Create a new working directory

```bash
mkdir ~/lab07-$MYGIT/SPTY2D1
````


## Start R Session

```bash
R --vanilla
```
```bash
MYGIT <- system("echo $MYGIT", intern = TRUE)
```
Load phylogenetics library
```bash
library(ape)

tree <- read.tree(paste0("~/lab07-",MYGIT,"/species_abr.tre"))
tree <- compute.brlen(tree, method = "Grafen")
```

Load species key to get species info and habitat (aquatic/terrestrial)
```bash
species_key <- read.csv(paste0("~/lab06-",MYGIT,"/species_key.csv"), stringsAsFactors = TRUE)
head(species_key)
```


## Load Net Charge Data from Lab 6

```
charges <- read.csv(paste0("~/lab06-",MYGIT,"/SPTY2D1/net_charge_summary.csv"), header = TRUE)
colnames(charges)[3] <- "Abbr"

charges <- read.table(paste0("~/lab06-",MYGIT,"/SPTY2D1/net_charges_labeled.tsv"),
                      sep = "\t", header = TRUE)

charges_summary <- aggregate(NetCharge ~ Abbr, data = charges, FUN = mean)
head(charges_summary)
```

## Load Secondary Structure Summary Data from Lab 6

```
dssp <- read.csv(paste0("~/lab06-",MYGIT,"/SPTY2D1/dssp_summary.csv"), header = TRUE)

dssp_summary <- aggregate(cbind(frac_helix, frac_sheet, frac_coil, mean_ASA) ~ abbr,
                          data = dssp, FUN = mean)
head(dssp_summary)
```


## RPS-BLAST Domain Content

Read RPS-BLAST output to get domain counts and lengths per species
Script extracts: 
- DomainCountPerGene – how many domains are there? (average per gene, per species)
This is the domain count (how many domains were found in the myoglobin sequence of each species). 
- CommonDomainLengthPerGene – How long is the most common domain? (average per gene, per species). This iss the most common domain length (if a domain is repeated, we take the length of the domain that occurs most frequently; for myoglobin, we expect one major globin domain length).
Note: Myoglobin is typically a single-domain protein (the globin domain ~150 amino acids long), so we expect each species to have 1 domain. But to illustrate the process (and to handle other gene families with multiple domains), we will still parse and summarize the domain data.

```bash
rps <- read.table(paste0("~/lab06-",MYGIT,"/SPTY2D1/SPTY2D1.rps-blast.out"),
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
```

```
names(rps) <- c("qseqid","qlen","qstart","qend","evalue","stitle")

parts <- strsplit(rps$qseqid, "\\|")
rps$Abbr <- sapply(parts, `[`, 1)
rps$RefSeq <- sapply(parts, `[`, 2)
```
Compute the hit length for each domain
```bash
rps$AlignLength <- abs(rps$qend - rps$qstart) + 1
```


## Per-Gene Domain Summary

```
gene_summary <- aggregate(
  AlignLength ~ Abbr + RefSeq,
  data = rps,
  FUN = function(x) {
    c(
      count = length(x),
      common_length = as.integer(names(sort(table(x), decreasing = TRUE)[1]))
    )
  }
)

gene_summary$DomainCountPerGene        <- gene_summary$AlignLength[, "count"]
gene_summary$CommonDomainLengthPerGene <- gene_summary$AlignLength[, "common_length"]
gene_summary$AlignLength <- NULL
```

```
species_means <- aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
                           data = gene_summary, FUN = mean)

n_gene_copies <- aggregate(RefSeq ~ Abbr, data = gene_summary,
                           FUN = function(x) length(unique(x)))
names(n_gene_copies)[2] <- "NumGeneCopies"
```


## Final Per-species Summary: mean domains per gene and mean domain length

```
domain_summary2 <- merge(
  aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
            data = gene_summary, FUN = mean),
  aggregate(RefSeq ~ Abbr, data = gene_summary,
            FUN = function(x) length(unique(x))),
  by = "Abbr"
)

names(domain_summary2)[names(domain_summary2) == "RefSeq"] <- "NumGeneCopies"
```


Combining all trait Data
```bash

traits0 <- merge(charges_summary, dssp_summary,
                 by.x = "Abbr", by.y = "abbr", all = TRUE)

traits <- merge(traits0, domain_summary2, by = "Abbr", all = TRUE)

traits$Status <- species_key$aquatic[ match(traits$Abbr, species_key$abbreviation) ]

traits
```

Matching Species in Tree and Data
```bash

species_in_data <- traits$Abbr
tree <- drop.tip(tree, setdiff(tree$tip.label, species_in_data))

row.names(traits) <- traits$Abbr

length(tree$tip.label)
nrow(traits)
```


## Phylogenetic Regression in R

Test if habitat (aquatic vs terrestrial) is a significant predictor of a trait, while accounting for phylogeny.
Use the phylolm package to fit a phylogenetic linear model.
```bash
library(phylolm)

traits$Status <- factor(traits$Status, levels = c("terrestrial", "aquatic"))

model_netcharge <- phylolm(NetCharge ~ Status, data = traits, phy = tree, model = "BM")
summary(model_netcharge)
```


## Plotting the Results

```
traits$aquatic01 <- ifelse(traits$Status == "aquatic", 1, 0)

png(paste0("~/lab07-",MYGIT,"/NetCharge_vs_Habitat.png"), width = 600, height = 400)
plot(jitter(traits$aquatic01, amount = 0.1), traits$NetCharge,
     xlab = "Habitat (0 = terrestrial, 1 = aquatic)",
     ylab = "SPTY2D1 Net Charge (pH 6.5)",
     main = "SPTY2D1 Net Charge in Terrestrial vs Aquatic Mammals",
     pch = 19, col = "blue")
abline(a = coef(model_netcharge)[1], b = coef(model_netcharge)[2], col = "red", lwd = 2)
dev.off()
```


## Running Regression for Helix Fraction

```
model_helix <- phylolm(frac_helix ~ Status, data = traits, phy = tree, model = "BM")
summary(model_helix)
```



## Running Regression for All Additional Traits

```
vars <- c("NetCharge","frac_helix","frac_sheet","frac_coil","mean_ASA",
          "DomainCountPerGene","CommonDomainLengthPerGene","NumGeneCopies")
```

```
fit_trait <- function(varname) {

  sub <- traits[!is.na(traits[[varname]]) & !is.na(traits$Status),
                c("Abbr","Status", varname)]
  names(sub)[3] <- "y"

  if (nrow(sub) < 5) {
    return(data.frame(
      Trait = varname, N = nrow(sub),
      Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
      Note = "Too few species or no variance"
    ))
  }

  tr <- drop.tip(tree, setdiff(tree$tip.label, sub$Abbr))
  sub <- sub[match(tr$tip.label, sub$Abbr), ]

  if (sd(sub$y) < .Machine$double.eps) {
    return(data.frame(
      Trait = varname, N = nrow(sub),
      Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
      Note = "No variation in response"
    ))
  }

  mod <- phylolm(y ~ Status, data = sub, phy = tr, model = "BM")
  sm  <- summary(mod)

  sub$aquatic01 <- ifelse(sub$Status == "aquatic", 1, 0)

  png(paste0(varname, "_vs_Habitat.png"), width = 650, height = 420)
  plot(jitter(sub$aquatic01, 0.08), sub$y,
       xlab = "Habitat (0 = terrestrial, 1 = aquatic)",
       ylab = varname, pch = 19)
  abline(a = coef(mod)[1], b = coef(mod)[2], col = "red", lwd = 2)
  title(paste0("Phylogenetic regression: ", varname, " ~ Status (BM)"))
  dev.off()

  data.frame(
    Trait = varname,
    N = nrow(sub),
    Beta_aquatic = unname(coef(mod)["Statusaquatic"]),
    SE = sm$coefficients["Statusaquatic","StdErr"],
    t = sm$coefficients["Statusaquatic","t.value"],
    p = sm$coefficients["Statusaquatic","p.value"],
    R2 = unname(sm$r.squared),
    Note = ""
  )
}
```

```
res_list <- lapply(vars, fit_trait)
res <- do.call(rbind, res_list)

write.csv(res,
          paste0("~/lab07-",MYGIT,"/SPTY2D1/SPTY2D1.trait_results.csv"),
          row.names = FALSE)

res
```
