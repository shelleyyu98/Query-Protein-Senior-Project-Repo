# Lab 6: Protein Domain and Structure Analysis in Terrestrial and Diving Mammals

## Start the lab, make sure your instance is running on EC2 and log in via ssh.
<img src="img/Fig1.png" alt= “” width="150" height="150"> <img src="img/Fig2.png" alt= “” width="150" height="150">

## Clone Lab 6
**You only need to do this ONE time.**
On the command line, clone the lab 6 repository.

```bash
git clone git@github.com:Bio312/lab06-$MYGIT.git
```

Git will now clone a copy of today's lab into a folder called lab06-$MYGIT (where \$MYGIT is your GitHub username). Go there:

```bash
cd lab06-$MYGIT
```

Take a look at what is in the folder using the `ls` command. At the end of the lab, you will push the changes and files back into the online repository.

**Do all the work for this lab in this new folder you have cloned.** Check this at any time by typing ```pwd```.
You should also edit the README.md within the instance. 
![](img/VSCode.png)

**We are continuing to work on primary analysis for your research paper.** Today, we are working on protein domain prediction and protein structure analysis, which will be a critical addition to your research paper.

# 2. Domain Identification for myoglobin Proteins 

We will be using [RPS-BLAST](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBWhat) to idenfity [Pfam](https://pfam.xfam.org/) domains within our protein sequences. Please read about RPS-BLAST, and especially about Pfam.

## Input Protein Sequences
We need to use unaligned protein sequence. 
We will use `myoglobin.homologs.fas`, which were our original, unaligned sequences from lab 4. 

A directory for the myoglobin sequences already exists, so you'll get an error here (but don't forget this step when you repeat the lab for your own gene family):
```bash
mkdir ~/lab06-$MYGIT/myoglobin
```

Change into that directory:
```bash
cd ~/lab06-$MYGIT/myoglobin
```

We will use our raw unaligned sequence. Make sure the file is there and take a look at it:

```bash
less ~/lab04-$MYGIT/myoglobin/myoglobin.homologs.fas
```

## The Pfam database is already downloaded on your instance. 

The [Pfam](https://www.ebi.ac.uk/interpro/) database has alrady been downloaded to `~/data/Pfam/`. It contains models that enable "functional analysis of proteins by classifying them into families and predicting domains and important sites. "  

## Run RPS-BLAST
Here is the command to run rps blast:

```bash
rpsblast -query ~/lab04-$MYGIT/myoglobin/myoglobin.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
><img src="img/github.png" alt= “” width="20" height="20">[Github] Please take a look at the output. What do you see? Describe it here.   

```answer/notes
I see a very low e-value presented after ue.
ue .0000000001

```

><img src="img/github.png" alt= “” width="20" height="20">[Github] Note the e-value. How does adjusting the e-value to 0.1 affect the domains that are outputted?   

```answer/notes
ue .0000000001
Adjusting the e-valie to 0.1 would increase the number of domains reported (more hits), which could increase the likelihood of there being false positives.
```

*You should use an evalue of .0000000001 to answer the Brightspace questions about the myoglobin genes. However, you may need to adjust this for your own gene family.*  

## Plot the predicted Pfam domains on the phylogeny.

We will use our final gene tree from lab 5. Make sure it is there and take a look at it:

```bash
less ~/lab05-$MYGIT/myoglobin/myoglobin.homologs.al.mid.treefile
```

We are going to use a script that Josh wrote in [R](https://www.r-project.org/) to plot the pfam domain predictions from rps-blast next to their cognate protein on the phylogeny.  

On the left hand side of the plot, the package [ggtree](https://guangchuangyu.github.io/software/ggtree/) will plot the phylogeny.  

On the right hand side of the plot, the package [drawProteins](https://www.bioconductor.org/packages/devel/bioc/vignettes/drawProteins/inst/doc/drawProteins_BiocStyle.html) will plot the pfam domains we detected with rps-blast.  

The script that Josh wrote, called plotTreeAndDomains2.r, was included when you cloned your lab.

Here are the parts of the script. These must be provided in order:
* **`Rscript`**: The Rscript program allows you to run an R script from the command line, without opening up the R console.
* **`--vanilla`**: This command flag indicates to R not save or restore a workspace or previous settings.
* **`plotTreeAndDomains.r`**: The script that Dr. Rest wrote. You can take a look at it using less, for example.
* Now provide the names of three files, in order: tree file, rps-blast output, and the name of the output file that you want. For example:
*  **`~/lab05-$MYGIT/myoglobin/myoglobin.homologs.al.mid.treefile`**: Our midpoint rooted tree file for the myoglobin proteins.
*  **`~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out`**: Our predicted pfram domains from rps-blast.  
* **`~/lab04-$MYGIT/myoglobin/myoglobin.homologs.fas`**: Unaligned gene sequences.
*  **`~/lab06-$MYGIT/myoglobin/myoglobin.tree.rps.pdf`**: The name of the output pdf file we will create.  

Now, run the script:

```bash
Rscript --vanilla ~/lab06-$MYGIT/plotTreeAndDomains2.r ~/lab05-$MYGIT/myoglobin/myoglobin.homologs.al.mid.treefile ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out ~/lab04-$MYGIT/myoglobin/myoglobin.homologs.fas ~/lab06-$MYGIT/myoglobin/myoglobin.tree.rps.pdf
```

In the output pdf, you should see the predicted domains annotated by the tips of your tree. You should also see a legend with the names of all the PFAM domains. 

For your own gene family, if there are formatting issues, let us know and we may be able to help you adjust the script so that the graphic looks better.  

<img src="img/github.png" alt= “” width="20" height="20">[Github] This file should be part of your remote repository. 

## Looking at the predicted domains in more detail

We will be focusing on predictions in the Pfam database. You can look up these Pfam domain descriptions on the [Pfam website](http://pfam.xfam.org/).

There are two type of files you could look at now to try and make sense of these domains:
1. The tab delimited annotations in myoglobin.rps-blast.out. Use less to look at it, put it in a spreadsheet program, or here is a nice way to look at it:  

```bash
mamba activate base
mlr --inidx --ifs "\t" --opprint  cat ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out | tail -n +2 | less -S
```
2. The domain-on-tree graphic that you just produced using Rests's R script.  

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 1. What is the order of the columns in the file ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out? Hint: see the command that was used to produce these results.   

### Examine the distribution of Pfam domains across proteins. 
      
><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 2. Which proteins have more than one annotation?

Here's an easy command to do this at the command line, if you don't feel like counting by hand:

```bash
cut -f 1 ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out | sort | uniq -c
```

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 3. Which Pfam domain annotation is most commonly found?

Here's an easy command to do this at the command line, if you don't feel like counting by hand:

```bash
cut -f 6 ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out | sort | uniq -c
```

><img src="img/brightspace.png" alt= “” width="20" height="20">
[Brightspace] 4. Which protein has the longest annotated protein domain?  

Here's an easy command to do this at the command line, if you don't feel like calculating by hand:
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out |  sort  -k2nr
```
#### The `sort` program.  
We have already used the sort command in the above examples. Let's learn a bit more about it.  
Use `man sort` to read the manual for this command. (or use google)

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 5. Match the command flag with what it does to the sort program.

#### E-values of protein domain hits
><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 6. Which protein has a domain annotation with the best e-value?   
Here's a command to pull out just the e-values:

```bash
cut -f 1,5 -d $'\t' ~/lab06-$MYGIT/myoglobin/myoglobin.rps-blast.out 
```

## Map RefGENBANK REFSEQ IDs to Uniprot
The goal of the next part of this lab is to retrieve and analyze tertiary structure information using Alphafold predicted structures. However, Alphafold uses UNIPROT IDs, rather than GENBANK REFSEQ IDs. So, our first job is to convert our GENBANK REFSEQ IDs into UNIPROT IDs.

Extract unique RefSeq accessions from your file
```bash
awk -F'|' '{print $2}' /home/bio312-user/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
 | sort -u > /home/bio312-user/lab06-$MYGIT/myoglobin/myoglobin_refseq_ids.txt
```

### Map RefGENBANK REFSEQ IDs to Uniprot IDs, make note of any missing IDs
```bash
awk 'NR==FNR{m[$1]=$2; next} {if($0 in m) print $0"\t"m[$0]; else print $0"\tMISSING"}' ~/lab06-$MYGIT/refseq_uniprot.tsv ~/lab06-$MYGIT/myoglobin/myoglobin_refseq_ids.txt > ~/lab06-$MYGIT/myoglobin/myoglobin_refseq_uniprot_subset.tsv
```

Taka a look at the output. 
```bash
less ~/lab06-$MYGIT/myoglobin/myoglobin_refseq_uniprot_subset.tsv
```
Note how many Genbank IDs map to Uniprot IDs, and how many are "MISSING".

In this exercise, we are retrieving protein structures to compare across species.  
Because **GenBank (RefSeq) accessions** do not have consistent structural resources available, we will first **map RefSeq IDs to UniProt IDs**, which are the standard identifiers used by both AlphaFold and SWISS-MODEL.  

- The mapping was performed by Josh using UniParc/UniProt cross-references to ensure that each RefSeq sequence could be linked to its corresponding UniProt entry.  
- Josh provided you with this **RefSeq → UniProt mapping file**, so you do not need to run the mapping yourself, but it’s important to know the reason for the switch: structural resources like AlphaFold and SWISS-MODEL are organized by UniProt, not GenBank.  
- The above command looks up each of your Genbank Refseq proteins in the mapping file to find the corresponding uniprot ID.
- Occasionally, a Genbank Refseq protein doesn't have an equivalent in Uniprot.

## afetch_by_uniprot.sh: Retreive predicted protein structures (PDBs) from Alphafold or Swissprot

Make sure you have write permissions on the retrieval script by running this command:
```bash
chmod +x ~/lab06-$MYGIT/afetch_by_uniprot.sh
```

Then, to retrieve predicted protein structures (PDBs) from Alphafold or Swissprot, run this command:

```bash
~/lab06-$MYGIT/afetch_by_uniprot.sh ~/lab06-$MYGIT/myoglobin/myoglobin_refseq_uniprot_subset.tsv ~/lab06-$MYGIT/myoglobin
 ```

Take a look: Within your  ~/lab06-$MYGIT/myoglobin folder, you should now have a .pdb file and a .json file (which contains predicted error) for each UNIPROT accession.

### What the afetch_by_uniprot.sh script does

When you run `afetch_by_uniprot.sh`, the script attempts to download **predicted protein structures** for each protein ID in the provided mapping file.

1. **AlphaFold (preferred source)**  
   - The script first checks the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).  
   - If a model is available, it downloads:
     - `*.pdb` → the 3D coordinates of the protein structure.  
     - `*.json` → predicted aligned error values (confidence estimates for each region).  

2. **SWISS-MODEL (fallback source)**  
   - If no AlphaFold structure is found, the script queries the [SWISS-MODEL Repository](https://swissmodel.expasy.org/repository/).  
   - If available, it downloads one or more `*.pdb` structure models.  

3. **Missing cases**  
   - If no UniProt mapping exists or no structural models are available in either database, nothing is downloaded for that protein.

### Interpreting the script output

Each line of the script reports what happened for a protein:

- **`AF OK`** → AlphaFold structure found and saved.  
- **`AF NONE`** → No AlphaFold model exists for this protein.  
- **`SWM OK`** → SWISS-MODEL structure found and saved (only tried if AlphaFold failed).  
- **`SWM NA`** → No SWISS-MODEL structure available.  
- **`MISS`** → No UniProt ID mapping was available, so no structure could be retrieved.  

### AlphaFold vs SWISS-MODEL: what to know

Both databases provide predicted structures, but they differ:

- **AlphaFold** uses a machine learning model trained on evolutionary and structural data.  
  - Provides per-residue confidence scores (highly useful for deciding which parts of the model are reliable).  
  - Covers the entire proteome for many organisms.  

- **SWISS-MODEL** uses **comparative (homology) modeling**: it builds structures based on experimental templates from similar proteins.  
  - Quality depends on the availability of a suitable template.  
  - Often available for proteins or species not yet modeled in AlphaFold.  

For **comparative structure analysis**, AlphaFold is generally preferred when available, because of its proteome-wide consistency and built-in confidence measures.  
However, when AlphaFold models are missing, it is reasonable to include SWISS-MODEL structures — especially since they are sequence-identical to your protein of interest. Just be aware that confidence scores are reported differently, and mixing models from different pipelines may introduce small differences in quality or coverage.

### Take a look at one of your .pdb files

For example:  
```bash
less ~/lab06-$MYGIT/myoglobin/XP_068416897.1__AF-P02177-F1-model_v4.pdb
```

The file starts with a header that includes important identifying information about the structure.  
But the most important lines to notice are the ones that begin with **`ATOM`**.  

Each `ATOM` line describes a single atom in the protein and follows a standard column format.  
For AlphaFold models, after the `MODEL` line you’ll see many `ATOM` entries like:

```mathematica
ATOM      1  N   MET A   1      -6.981  -5.221  16.808  1.00 80.51           N  
ATOM      2  CA  MET A   1      -5.635  -5.695  16.455  1.00 80.51           C  
ATOM      3  C   MET A   1      -4.681  -4.614  16.898  1.00 80.51           C  
ATOM      4  CB  MET A   1      -5.525  -5.950  14.950  1.00 80.51           C  
```

Breaking down the columns:  

- **ATOM** – record type (this is an atom in the structure).  
- **1, 2, 3, …** – atom serial number.  
- **N, CA, C, CB…** – atom name (nitrogen, alpha carbon, carbon, beta carbon, etc.).  
- **MET** – residue (amino acid) name.  
- **A** – chain identifier.  
- **1** – residue sequence number.  
- **-6.981 -5.221 16.808** – x, y, z 3D coordinates of the atom.  
- **1.00** – occupancy (here always 1.00 in AlphaFold).  
- **80.51** – B-factor field. In AlphaFold files this encodes the **pLDDT confidence score** (not thermal motion as in experimental PDBs).  
- **N / C / O / …** – element symbol.  

👉 For AlphaFold PDBs, the **pLDDT score** is especially important:  
- High pLDDT (near 100) = very confident prediction.  
- Low pLDDT = more uncertain or flexible/disordered region.  

So when you “look inside” one of these PDBs, you’re seeing the entire protein atom by atom, including its amino acid sequence (through the residues) and the AlphaFold confidence for each part of the structure.  


## Comparative Analysis of the Protein Structures

We have been using the 'base' conda environment, which contains much of the software that we use. Activate a different conda environment installed on your system 'bio312', which Josh set up and contains software for analyzing PDB structures.

 ```bash
 mamba activate bio312
 ```

We are missing on important piece of software in this environment, so install it now: (only needs to be done once)
 ```bash
 mamba install -y -c conda-forge pdb2pqr
```

We already have PDB files, 
Why install pdb2pqr? It adds hydrogens, predicts protonation states (via built-in pKa heuristics), and writes a PQR file where the 9th column is the per-atom partial charge.
Summing those charges ≈ net charge of the protein at the chosen pH.
We’ll use pH 6.5, matching the myoglobin net-charge analysis.

### Run the charge calculation on the downloaded PDB structures

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

```

You’ve now generated .pqr files with pdb2pqr.
Check to see that you've generated cleaned .protein.pdb and .pqr files, and that you have one for each .pdb file:
```bash
echo; echo "Original PDB: $(ls ~/lab06-$MYGIT/myoglobin/*.pdb 2>/dev/null | grep -v '.protein.pdb' | wc -l), \
Protein-only PDB: $(ls ~/lab06-$MYGIT/myoglobin/*.protein.pdb 2>/dev/null | wc -l), \
PQR: $(ls ~/lab06-$MYGIT/myoglobin/*.pqr 2>/dev/null | wc -l)"
```
For myoglobin, you should see:
Original PDB: 17, Protein-only PDB: 17, PQR: 17

# Two ways to visualize and analyze the charges of the protein structures
You have a choice here: you can use the cool Jupyter notebook approach to visualize and analyze the charges (Option 1, below). If you run into technical troubles, Option 2 can be run in your bash shell and will accomplish the goal, but without some of the same cool features.

## Option 1: Using a Jupyter Notebook

We will use Jupyter Notebooks for some structure visualization and analysis.
You have two options:
- Run Jupyter in your web browser (classic)
- Run Jupyter directly inside VS Code (install Python and Jupyter extensions))
We will use the web browser, which is a more typical experience.

This should already be acitvated:
```bash
mamba activate bio312
```

```bash
jupyter lab --no-browser --ip=0.0.0.0 --port=8888 --ServerApp.token=''
```
*A message should appear; click "Open in Browswer"*
If you miss this, or it doesn't work, pointing your browswer to the address indicated in the terminal, such as:
http://localhost:8888/lab
http://127.0.0.1:8888/lab

Double click on "BIO312_Myoglobin_Visualization.ipynb" document within myoglobin.
*Note: when you run this on your OWN gene family, you should make a copy of the notebook: `cp ~/lab06-$MYGIT/myoglobin/BIO312_Myoglobin_Visualization.ipynb  ~/lab06-$MYGIT/YOURGENE/BIO312_YOURGENE_Visualization.ipynb`

Now, run the code from top to bottom in the *BIO312: Myoglobin Structure Visualization* notebook.

## Option 2: Using the BASH shell

In this part of the lab, we want to **visualize myoglobin structures** and compare **electrostatic charge distributions** between aquatic and terrestrial mammals.  
Previously, we generated `.pqr` files with **pdb2pqr**. Now we’ll use simple bash + python steps to make **static PNGs** of the protein surfaces.

## Compute net charge for each PQR
Charges are already embedded in your `.pqr` files (last two columns). This script will sum them:

```bash
out=~/lab06-$MYGIT/myoglobin/net_charges.tsv
: > "$out"   # truncate/create

for f in ~/lab06-$MYGIT/myoglobin/*.pqr; do
  base=$(basename "$f")
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9} END{printf "%.3f", s}' "$f")
  printf "%s\tNetCharge=%s\n" "$base" "$Z" >> "$out"
done
```

Look at the resulting tab-separated file `net_charges.tsv` with each protein and its net charge.  
```bash
less ~/lab06-$MYGIT/myoglobin/net_charges.tsv
```

## 4. Visualize structures with charges (static PNGs)
We can’t use interactive viewers outside Jupyter, but we can generate PNG snapshots with `py3Dmol` (Python).  
First, install pymol (one time only)
```bash
mamba activate bio312
mamba install -y -c conda-forge pymol-open-source
```

First, we need to create pdbs that have our pqr charges on them:
```bash
python ~/lab06-$MYGIT/pqr_to_charges_pdb.py ~/lab06-$MYGIT/myoglobin
```
If this worked it should say "Done: created #/# charges PDBs."

Next, run the following loop to render the PDB structures colored by these charges:

```bash
python3 ~/lab06-$MYGIT/render_charge_pngs.py ~/lab06-"$MYGIT"/myoglobin
```

**Output:**  
A `.png` image for each `.pqr` file in the same directory. These images show protein surfaces colored by charge:  
- **Blue = positive charge**  
- **Red = negative charge**  
- **White/grey = neutral**  

---

Take a look at the images. Need help figuring out which images belong to which species, and which are aquatic and terrestrial? This will build some helper files for you:
```bash
# Go to your images directory
cd ~/lab06-"$MYGIT"/myoglobin

# 1) Build helper maps:
#    - abbreviation -> status (aquatic/terrestrial)
awk -F, 'NR>1 {print $2"\t"$5}' ~/lab06-"$MYGIT"/species_key.csv \
  | sort -u > abbr_status.tsv

#    - refseq -> abbreviation (from lab03 list like Hsap|NP_...|MB)
awk -F'|' '{print $2"\t"$1}' ~/lab03-"$MYGIT"/myoglobin/myoglobin.blastp.detail.filtered.out \
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
```

Now, take a look at PNG_GUIDE.txt to figure out which PNGs are terrestrial and which are aquatic.
```
less PNG_GUIDE.txt
```
* What to pay attention to in the images
Color meaning (from .pqr charges):
- Blue: positive (Lys/Arg/His side chains are often contributors)
- Red: negative (Asp/Glu)
- White/grey: near-neutral
* Overall trend you’re testing:
Diving/semi-aquatic mammals should show more positive (blue) surface patches than terrestrial mammals. This is tied to aggregation avoidance at high intracellular myoglobin concentration.

Where the color shows up:
- Focus on the outer surface rather than the cartoon backbone (cartoon color is just a rainbow by residue index unless you changed it).
- Look for broad blue patches vs mixed/neutral surface.

Caveats:
- Our coloring reflects total per-atom charge assigned by pdb2pqr at the chosen pH; it’s a useful proxy for surface charge distribution, but not the exact “net surface charge” metric in Mirceta et al.
- Structures are not guaranteed to be aligned; orientation may differ between species.



###  Compare aquatic vs terrestrial proteins
Now combine your net charge results with the `species_key.csv` file to make a simple boxplot:

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
lab03_list     = os.path.expanduser(f"~/lab03-{MYGIT}/myoglobin/myoglobin.blastp.detail.filtered.out")

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

**Output:**  
- A PNG file `netcharge_boxplot.png` showing net charge distributions for aquatic vs terrestrial species.  

Now, take a look at this png.


#### What you should find
- **Terrestrial species** often have lower net charge (closer to neutral).  
- **Aquatic/diving species** typically have more positive net charges.  
- The PNG images give a **qualitative** view of surface charge patches.  
- The boxplot provides a **quantitative comparison** between habitats.  

Is this what you see?

## After you have finished Option 1 or Option 2
Once you have finished the Jupyter Notebook, you should have:
- net_charge_boxplot.png
- net_charge_summary.csv

><img src="img/github.png" alt= “” width="20" height="20">[Github] These two files should be part of your GitHub repository.

Now, continue with using DSSP to extract secondary structure characteristics of myoglobin.

First, install some more software:
```
mamba activate bio312
mamba install -y -c conda-forge mdtraj pandas matplotlib
```

Josh made this script to complete the analysis of secondary structure:
```bash
python3 ~/lab06-$MYGIT/dssp_batch_summary_mdtraj.py \
  --pdb-dir   ~/lab06-$MYGIT/myoglobin \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-csv   ~/lab06-$MYGIT/myoglobin/dssp_summary.csv \
  --plots
```

This previous python script:
- Generated .dssp files
- Parse the generated .dssp files
- Counted fraction helix (H/G/I), sheet (E/B), coil (everything else)
- Reported average backbone ASA (exposedness proxy)
- Saved a CSV and make simple plots grouped by aquatic/terrestrial using species_key.csv

Examine each of these.


# 3. Work on your own gene family: Motif and Domain Identification Using RPS-BLAST and Pfam, and Python Analysis
Repeat the above lab, using your own gene family.

><img src="img/github.png" alt= “” width="20" height="20">[Github] **Keep track of *all* the commands that you used here! (i.e. cut and paste them here, so you could recreate your work!)**.
1. mkdir ~/lab06-$MYGIT/SPTY2D1
2. cd ~/lab06-$MYGIT/SPTY2D1
3. less ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas
4. rpsblast -query ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
5. less ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile
6. Rscript --vanilla ~/lab06-$MYGIT/plotTreeAndDomains2.r ~/lab05-$MYGIT/SPTY2D1/SPTY2D1.homologs.al.mid.treefile ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out ~/lab04-$MYGIT/SPTY2D1/SPTY2D1.homologs.fas ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.tree.rps.pdf
7. mamba activate base
mlr --inidx --ifs "\t" --opprint  cat ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | tail -n +2 | less -S
8. cut -f 1 ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | sort | uniq -c
9. cut -f 6 ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out | sort | uniq -c
10. awk '{a=$4-$3;print $1,'\t',a;}' ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out |  sort  -k2nr
11. cut -f 1,5 -d $'\t' ~/lab06-$MYGIT/SPTY2D1/SPTY2D1.rps-blast.out 
12. awk -F'|' '{print $2}' /home/bio312-user/lab03-$MYGIT/SPTY2D1/SPTY2D1.blastp.detail.filtered.out \
 | sort -u > /home/bio312-user/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_ids.txt
13. awk 'NR==FNR{m[$1]=$2; next} {if($0 in m) print $0"\t"m[$0]; else print $0"\tMISSING"}' ~/lab06-$MYGIT/refseq_uniprot.tsv ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_ids.txt > ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv
14. less ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv
15. chmod +x ~/lab06-$MYGIT/afetch_by_uniprot.sh
16. ~/lab06-$MYGIT/afetch_by_uniprot.sh ~/lab06-$MYGIT/SPTY2D1/SPTY2D1_refseq_uniprot_subset.tsv ~/lab06-$MYGIT/SPTY2D1
17. #!/usr/bin/env bash
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
18. echo; echo "Original PDB: $(ls ~/lab06-$MYGIT/SPTY2D1/*.pdb 2>/dev/null | grep -v '.protein.pdb' | wc -l), \
Protein-only PDB: $(ls ~/lab06-$MYGIT/SPTY2D1/*.protein.pdb 2>/dev/null | wc -l), \
PQR: $(ls ~/lab06-$MYGIT/SPTY2D1/*.pqr 2>/dev/null | wc -l)"
19. mamba activate bio312
20. jupyter lab --no-browser --ip=0.0.0.0 --port=8888 --ServerApp.token=''
21. cp ~/lab06-$MYGIT/myoglobin/BIO312_Myoglobin_Visualization.ipynb  ~/lab06-$MYGIT/SPTY2D1/BIO312_SPTY2D1_Visualization.ipynb
22. out=~/lab06-$MYGIT/SPTY2D1/net_charges.tsv
: > "$out"   # truncate/create

for f in ~/lab06-$MYGIT/SPTY2D1/*.pqr; do
  base=$(basename "$f")
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9} END{printf "%.3f", s}' "$f")
  printf "%s\tNetCharge=%s\n" "$base" "$Z" >> "$out"
done
23. less ~/lab06-$MYGIT/SPTY2D1/net_charges.tsv
24. mamba activate bio312
mamba install -y -c conda-forge pymol-open-source
25. python ~/lab06-$MYGIT/pqr_to_charges_pdb.py ~/lab06-$MYGIT/SPTY2D1
26. python3 ~/lab06-$MYGIT/render_charge_pngs.py ~/lab06-"$MYGIT"/SPTY2D1
27. # Go to your images directory
cd ~/lab06-"$MYGIT"/SPTY2D1

# 1) Build helper maps:
#    - abbreviation -> status (aquatic/terrestrial)
awk -F, 'NR>1 {print $2"\t"$5}' ~/lab06-"$MYGIT"/species_key.csv \
  | sort -u > abbr_status.tsv

#    - refseq -> abbreviation (from lab03 list like Hsap|NP_...|MB)
awk -F'|' '{print $2"\t"$1}' ~/lab03-"$MYGIT"/myoglobin/myoglobin.blastp.detail.filtered.out \
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
28. less PNG_GUIDE.txt
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
lab03_list     = os.path.expanduser(f"~/lab03-{MYGIT}/myoglobin/myoglobin.blastp.detail.filtered.out")

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
29. mamba activate bio312
mamba install -y -c conda-forge mdtraj pandas matplotlib
30. python3 ~/lab06-$MYGIT/dssp_batch_summary_mdtraj.py \
  --pdb-dir   ~/lab06-$MYGIT/myoglobin \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-csv   ~/lab06-$MYGIT/myoglobin/dssp_summary.csv \
  --plots


><img src="img/github.png" alt= “” width="20" height="20">[Github] Provide at least one citation that you used to read about the domains that you expect to be present in your gene family. What domains does the citation indicate will be there? Did you detect them?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Provide at least one citation that you used to read about the secondary or tertiary structure of one or more proteins in your gene family. [May be the same or different from above].

><img src="img/github.png" alt= “” width="20" height="20">[Github] Examine the pdf plot of your gene tree and domain annotations. This pdf should be part of your repository. 

><img src="img/github.png" alt= “” width="20" height="20">[Github] Which pfam domains are found in the most genes?
18 pfam08243, SPT2, SPT2 chromatin protein.  This family includes the Saccharomyces cerevisiae protein SPT2 which is a chromatin protein involved in transcriptional regulation.

><img src="img/github.png" alt= “” width="20" height="20">[Github] Are any pfam domains found multiple times in a single gene?
No

><img src="img/github.png" alt= “” width="20" height="20">[Github] Examine the  plots of surface charge and secondary structure traits summarized by aquatic/terrestrial (mean_ASA,frac_coil,frac_sheet, frac_helix, ) for your gene family. These graphics should be part of your repository. 

><img src="img/github.png" alt= “” width="20" height="20">[Github] Describe whether the plots mentioned above show interesting differences between aquatic and terrestrial mammals. Were you expecting any differences based on your hypothesis and predictions?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Examine csv file: dssp_summary.csv. This should be part of your repository.

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 8. Show your TA the five graphics (pngs) and two tables (csv) that have been generated FOR YOUR GENE FAMILY. They will give you a code word. 

### Optional questions to answer about PFAM domains for your gene family:
><img src="img/github.png" alt= “” width="20" height="20">[Github] Which gene has the most pfam domains?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Which gene has the longest pfam domain annotation? When domain is this?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Which domain annotations has the best e-value?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Compare the domains in different lineages of interest for your paper. Are the predictions different or similar?

><img src="img/github.png" alt= “” width="20" height="20">[Github] 
Are there some gene copies that have all the domains, and others that are missing some domains? What does this suggest about the function of your gene family in some lineages we are studying?

><img src="img/github.png" alt= “” width="20" height="20">[Github] Is there any evidence that some gene copies have gained or lost function in some lineages?


# 4. Save your history and push your files into the repository.

Save your history:

```bash
history > lab6.commandhistory.txt
```
  
## Push any new files to the remote Git repository.

Add all of your results to the repository, commit them, and push them to the remote repository at GitHub.com

First, ensure you are in the main lab5 directory:
```bash
cd ~/lab06-$MYGIT
```
Now, commit your changes. **Make sure each line works before you continue with the next line.**

```bash
find . -size +5M | sed 's|^\./||g' | cat >> .gitignore; awk '!NF || !seen[$0]++' .gitignore
```
```bash
git add .
```
Take a look to see which changes have been staged, which haven’t, and if any files aren’t being tracked by Git. 
```bash
git status
```
```bash
git commit -a -m "Adding all new data files I generated in AWS to the repository."
```
```bash
git pull --no-edit
```
Did all of the above work? **Make sure!** If so, then run the following command:
```bash
git push 
```
