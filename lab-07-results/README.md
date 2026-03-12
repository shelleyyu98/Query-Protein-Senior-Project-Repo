# Lab 7: Phylogenetic Regression between Structural Traits and Lifestyle Transitions

## Background: Why Use Phylogenetic Regression?
Comparative studies must account for the shared evolutionary history of species. Phylogenetic regression is a statistical method that incorporates phylogenies into regression analyses of trait data. In traditional analyses (e.g., t-tests or ANOVA comparing group means), we assume each data point (species) is independent. However, closely related species often have similar traits due to common ancestry, violating the independence assumption. This can lead to misleading results – for example, two traits with no true relationship can show a high correlation if the species sharing certain trait values are phylogenetically clustered. Phylogenetic regression (also known as phylogenetic generalized least squares, or PGLS) addresses this by adjusting for the covariance among species induced by their phylogeny. In essence, it asks: are trait differences associated with an ecological factor (like habitat) above and beyond what we'd expect from phylogenetic relatedness?

In this lab, we test whether mammals that transitioned from terrestrial to aquatic environments evolved distinct myoglobin traits (surface charge, secondary structure, domain features). **In lab 6, we simply compared the mean trait values of aquatic vs. terrestrial mammals**. In Lab 7, we will use phylogenetic regression to account for the fact that aquatic mammals (e.g., whales, seals) are more closely related to some terrestrial species than to others. This method provides a more rigorous test of our hypothesis by controlling for shared ancestry. If a trait is truly associated with aquatic adaptation, a phylogenetic regression can detect that signal even when relatedness is taken into account. Conversely, if a simple group difference disappears under phylogenetic control, it suggests the difference might have arisen through phylogeny alone (e.g., all aquatic species happen to be from one lineage) rather than independent adaptation.

## Why R? 
We will use the R programming environment for this analysis because it has excellent support for phylogenetic comparative methods (packages like ape, phylolm, nlme, etc.) and is well-suited for statistical analysis. You can run R in our AWS environment by typing R at the shell prompt. This will launch the interactive R prompt (denoted by >). You can then enter the R commands provided below. (To exit R later, use q() or press Ctrl+D.)

**Reminder: R is interactive** – you can run one command at a time and inspect the output. Use this to your advantage to understand each step. 

**Data Setup: Reading Trees and Trait Data**
We will use data generated in previous labs:  
- Species tree: A phylogenetic tree of the mammal species (file species_abr.tre). This tree contains the evolutionary relationships among species, identified by their abbreviations (the same abbreviations used in our data files). 
- Species key: A CSV file (species_key.csv) mapping each species’ full name, abbreviation, and whether it’s aquatic or terrestrial. 
- **Trait data files:**
- net_charge_summary.csv (OPTION 1) and/or net_charges_labeled.tsv (OPTION 2) – net surface charge  (at pH 6.5) for each species (calculated in Lab 6). 
- dssp_summary.csv – secondary structure summary (e.g., fraction helix, sheet, coil, and average solvent accessibility) for each species (from Lab 6 DSSP analysis). 
- rps-blast.out – results of RPS-BLAST (domain search) (from Lab 6), which we’ll parse to get the domain content (number of domains) and the typical domain length for each species’ gene.

## Start the lab, make sure your instance is running on EC2 and log in via ssh.
<img src="img/Fig1.png" alt= “” width="150" height="150"> <img src="img/Fig2.png" alt= “” width="150" height="150">

## Clone Lab 7
**You only need to do this ONE time.**
On the command line, clone the lab 6 repository.

```bash
git clone git@github.com:Bio312/lab07-$MYGIT.git
```

Git will now clone a copy of today's lab into a folder called lab06-$MYGIT (where \$MYGIT is your GitHub username). Go there:

```bash
cd lab06-$MYGIT
```

Take a look at what is in the folder using the `ls` command. At the end of the lab, you will push the changes and files back into the online repository.



Create a directory
```bash
mkdir ~/lab07-$MYGIT/myoglobin
```

Open R by typing the following at the command line:
```bash
R --vanilla
```

```grab the $MYGIT system variable
MYGIT <- system("echo $MYGIT", intern = TRUE)
```

Let’s begin by loading the species phylogeny and ensuring it’s ready for analysis:

Load phylogenetics library
```R
library(ape)
```

Read the species tree; we have a slightly updated version in the current repository that ensures the branching patterns are interpreted accurately by all software.
```R
tree <- read.tree(paste0("~/lab07-",MYGIT,"/species_abr.tre"))
```

The species tree has no branch lenghts, but we need branch lengths to conduct our analyses. We can use the following approach to assign some arbirtray branch lengths based on the size of each clade.
```R
tree <- compute.brlen(tree, method = "Grafen")
```

Next, load the species key and trait datasets:

Load species key to get species info and habitat (aquatic/terrestrial)
```R
species_key <- read.csv(paste0("~/lab06-",MYGIT,"/species_key.csv"), stringsAsFactors = TRUE)
```
quick peek at the beginning of the table 
```R
head(species_key)  
```
### Load the Net Charge calculations from Lab 6
Next, we'll load the net protein charge data that we calculated in Lab 6. You may have either calculated the charge data with Option 1 and/or Option 2.

**If you used Option 1, use the following code:**
```R
charges <- read.csv(paste0("~/lab06-",MYGIT,"/myoglobin/net_charge_summary.csv"), header = TRUE)
colnames(charges)[3] <- "Abbr"
```

**If you used Option 2, use the following code:**
```R
charges <- read.table(paste0("~/lab06-",MYGIT,"/myoglobin/net_charges_labeled.tsv"), sep = "\t", header = TRUE)
```

**Continue here for either option 1 or option 2**
Summarize the average net charge per species:
```R
charges_summary <- aggregate(NetCharge ~ Abbr, data = charges, FUN = mean)
```

quick peek at the beginning of the table 
```R
head(charges_summary)
```
or at the FULL table:
```R
charges_summary
```

### Load secondary structure summary data from Lab 6
```R
dssp <- read.csv(paste0("~/lab06-",MYGIT,"/myoglobin/dssp_summary.csv"), header = TRUE)
```

Summarize the data per species; this will be  asummary for the following calculations:
- frac_helix — Fraction of residues in α-helix secondary structures.
- frac_sheet — Fraction of residues forming β-strands or β-bridges (sheet structures).
- frac_coil — Fraction of residues not in helix or sheet (loops, turns, irregular regions).
- mean_ASA — Mean solvent-accessible surface area per residue, averaged across the protein.

```R
dssp_summary <- aggregate(cbind(frac_helix, frac_sheet, frac_coil, mean_ASA) ~ abbr, data = dssp, FUN = mean)
```

quick peek at the beginning of the table 
```R
head(dssp_summary)
```
### RPS-BLAST domain content
Now, we need to parse the RPS-BLAST output to get domain info per species. The RPS-BLAST results (rps-blast.out) contain the domains detected in each protein sequence. We will create a small script to extract: 
- DomainCountPerGene – how many domains are there? (average per gene, per species)
This is the domain count (how many domains were found in the myoglobin sequence of each species). 
- CommonDomainLengthPerGene – How long is the most common domain? (average per gene, per species). This iss the most common domain length (if a domain is repeated, we take the length of the domain that occurs most frequently; for myoglobin, we expect one major globin domain length).
Note: Myoglobin is typically a single-domain protein (the globin domain ~150 amino acids long), so we expect each species to have 1 domain. But to illustrate the process (and to handle other gene families with multiple domains), we will still parse and summarize the domain data.
 
Read RPS-BLAST output to get domain counts and lengths per species
```R
rps <- read.table(paste0("~/lab06-",MYGIT,"/myoglobin/myoglobin.rps-blast.out"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
```

take a quick peek
```R
rps
```

Notice that there are no column names. We can assign them:
```R
names(rps) <- c("qseqid","qlen","qstart","qend","evalue","stitle")
```
now, take a quick peek again.
```R
rps
```

The following lines of code will parse species abbreviation and RefSeq ID from qseqid 
- Example qseqid: "Hsap|NP_001349775.1|MB"
If you have trouble running the entire code block, then run it line by line.
The first field = abbreviation
second field = RefSeq accession
```R
parts <- strsplit(rps$qseqid, "\\|")
rps$Abbr <- sapply(parts, `[`, 1)
rps$RefSeq <- sapply(parts, `[`, 2)
```

Next, compute the hit length for each domain
```R
rps$AlignLength <- abs(rps$qend - rps$qstart) + 1
```
><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 1. Which species has the longest predicted myoglobin domain?

Per-gene summary: how many domains in each gene copy, and its typical domain length
If you have trouble running the entire code block, then run it line by line.
```R
gene_summary <- aggregate(AlignLength ~ Abbr + RefSeq, data = rps,
                          FUN = function(x) {
                            c(count = length(x),
                              common_length = as.integer(names(sort(table(x), decreasing = TRUE)[1])))
                          })

gene_summary$DomainCountPerGene        <- gene_summary$AlignLength[, "count"]
gene_summary$CommonDomainLengthPerGene <- gene_summary$AlignLength[, "common_length"]
gene_summary$AlignLength <- NULL
```

Per-species summary from per-gene values
Mean domains per gene (should be ~1 for myoglobin) and mean domain length
```R
species_means <- aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
                           data = gene_summary, FUN = mean)
```

Also report how many gene copies we observed per species.
If you have trouble running the entire code block, then run it line by line.
```R
n_gene_copies <- aggregate(RefSeq ~ Abbr, data = gene_summary,
                           FUN = function(x) length(unique(x)))
                           
names(n_gene_copies)[2] <- "NumGeneCopies"
```

Take a look
```R
n_gene_copies
```

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 2. Which species has the most myogobin gene copies?

**Final per-species summary: mean domains per gene and mean domain length**
```R
domain_summary2 <- merge(
  aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
            data = gene_summary, FUN = mean),
  aggregate(RefSeq ~ Abbr, data = gene_summary,
            FUN = function(x) length(unique(x))),
  by = "Abbr"
)
names(domain_summary2)[names(domain_summary2) == "RefSeq"] <- "NumGeneCopies"
```

Take a look:
```R
domain_summary2
```
### Combining all trait data
Now, let’s combine all trait data into one data frame for analysis.
If you have trouble running the entire code block, then run it line by line.

Merge trait data frames by species abbreviation and add habitat status (aquatic/terrestrial) from species_key

```R
traits0 <- merge(charges_summary,dssp_summary,by.x = "Abbr",by.y="abbr",all=TRUE)
traits <- merge(traits0,domain_summary2, by = "Abbr",all=TRUE)
traits$Status <- species_key$aquatic[ match(traits$Abbr, species_key$abbreviation) ]
```

Take a look at the combined data:
```R
traits
```

At this point, the traits data frame should have one row per species, and columns including: 
- Abbr (species code) 
- NetCharge (net surface charge of myoglobin) 
- frac_helix, frac_sheet, frac_coil (fractions of secondary structure types) 
- mean_ASA (mean solvent accessible surface area of backbone, an indicator of how exposed the protein is) 
- DomainCountPerGene (average number of domains in  protein) - CommonDomainLengthPerGene (average length of the most common domain) 
- NumGeneCopies  
- Status (aquatic or terrestrial)

><img src="img/github.png" alt= “” width="20" height="20">[Github] Paste the output of head(traits) here (the first ~6 lines of the traits data frame). Ensure it shows the columns described above. This confirms that your data merging was successful.


```
 Abbr NetCharge frac_helix frac_sheet frac_coil   mean_ASA DomainCountPerGene
1 Acar      -1.0  0.7483871          0 0.2516129 0.03510972                  1
2 Bacu       2.0  0.7467532          0 0.2532468 0.03584368                  1
3 Btau      -1.0  0.7402597          0 0.2597403 0.03627062                  1
4 Casa       1.0  0.7581699          0 0.2418301 0.03629714                  1
5 Caur       1.5  0.7482993          0 0.2517007 0.03614842                  1
6 Chir      -1.0  0.7402597          0 0.2597403 0.03663393                  1
  CommonDomainLengthPerGene NumGeneCopies      Status
1                       112             1 terrestrial
2                       112             1     aquatic
3                       116             1 terrestrial
4                       116             1 terrestrial
5                       113             2     aquatic
6                       116             1 terrestrial

```

Before proceeding, we should ensure that our trait data and phylogeny have the same set of species. We should drop any species from the tree that we have no trait data for (and vice versa):

Prune the tree to match the trait data species, and vice versa, Drop tips from the tree that are not in our data, Make sure row names of the data match species names for some functions. 

```R
species_in_data <- traits$Abbr
tree <- drop.tip(tree, setdiff(tree$tip.label, species_in_data))
row.names(traits) <- traits$Abbr
```

Quick check:
number of species in tree
```R
length(tree$tip.label) 
```

number of species in trait data
```R
nrow(traits)            
```

These two numbers should match now (each representing the count of species we are analyzing).

## Phylogenetic Regression in R
With data prepared, we can perform a phylogenetic regression. We will test if habitat (aquatic vs terrestrial) is a significant predictor of a trait, while accounting for phylogeny. Specifically, *for myoglobin*, we’ll start with the hypothesis from earlier labs: aquatic mammals have higher myoglobin net surface charge than terrestrial mammals (thought to prevent protein aggregation at high myoglobin concentrations).

We'll use the `phylolm` package to fit a phylogenetic linear model. This works similarly to a standard linear model (lm), but requires a tree and assumes a model of trait evolution (we will use a Brownian motion model by default). 

First, load the phylolm library:
```R
library(phylolm)
```

Now, set up the regression model. We treat Status (aquatic vs terrestrial) as a categorical predictor. In R, if Status is a factor with two levels, the regression will effectively estimate an intercept (for the baseline level) and a coefficient for the other level. Let’s ensure Status is a factor and set the baseline to "terrestrial":
```R
traits$Status <- factor(traits$Status, levels = c("terrestrial", "aquatic"))
```

Now fit the phylogenetic regression of net charge on habitat by fitting the phylogenetic linear model: NetCharge as response, Status (habitat) as predictor
```R
model_netcharge <- phylolm(NetCharge ~ Status, data = traits, phy = tree, model = "BM")
```

Show a summary of the model
```R
summary(model_netcharge)
```

The summary output will show the estimated coefficients for the regression. Because we encoded Status as a factor with "terrestrial" as the base level, the output will include: 
- (Intercept) – this is the modeled mean net charge for a terrestrial species (on the phylogeny). 
- Statusaquatic – this is the difference in net charge for aquatic species relative to terrestrial. In other words, add this value to the intercept to get the aquatic species’ mean. A positive coefficient means aquatic species tend to have higher net charge than terrestrial species, consistent with our hypothesis.
- If p < 0.05 for Statusaquatic, we conclude there is a statistically significant difference in myoglobin net charge between aquatic and terrestrial mammals even after controlling for phylogeny. 

><img src="img/github.png" alt= “” width="20" height="20">[Github]  Paste the regression summary output here. Include the coefficient for Statusaquatic, its standard error, and the test statistic/p-value. This output is evidence of whether habitat is a significant predictor of net charge.

```
Call:
phylolm(formula = NetCharge ~ Status, data = traits, phy = tree, 
    model = "BM")

   AIC logLik 
25.151 -9.575 

Raw residuals:
    Min      1Q  Median      3Q     Max 
-0.8077 -0.1828  0.1923  1.0048  1.1923 

Mean tip height: 1
Parameter estimate(s) using ML:
sigma2: 0.5316329 

Coefficients:
              Estimate   StdErr t.value   p.value    
(Intercept)   -0.19233  0.40323 -0.4770    0.6407    
Statusaquatic  2.25010  0.24789  9.0771 3.056e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-squared: 0.8548       Adjusted R-squared: 0.8444 

```

**Interpreting the Regression Results**
Now, interpret what you see: 
- The intercept (terrestrial baseline) gives the estimated net charge for a terrestrial mammal’s myoglobin (at pH 6.5). 
- The Statusaquatic coefficient tells how much higher (or lower) the aquatic species’ myoglobin net charge is, on average, compared to terrestrial. We expect this to be positive (since previous analysis and literature indicated higher net charge in diving species). 
- The p-value (if provided) for the Statusaquatic effect tests the null hypothesis that this coefficient is zero (no difference). A small p-value (e.g., p < 0.05) means we reject the null and conclude there is a significant difference in net charge between aquatic and terrestrial species beyond what would be expected from phylogenetic relatedness alone. If the p-value is large, it suggests any difference in means could be explained by chance given the phylogeny (i.e., not a significant adaptive signal).
For example, if the summary output showed Statusaquatic = 5.0 ± 1.5 (SE), t = 3.3, p = 0.003, we would interpret that aquatic species have, on average, a net charge ~5 units higher than terrestrial species, and this difference is statistically significant. (This is just an illustrative example – use the actual numbers from your output.)

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 3. Review the observed p-value of the Statusaquatic effect for net surface charge of myoglobin. What is the interpretation of the result?

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 4. In the regression output, you should see an estimate of the coefficient labeled Statusaquatic with a value of 2.25. What does this estimate represent in biological terms?

### Plotting the Results
Let's visualize the relationship between habitat and net charge, along with the regression result. We will create a scatter plot of net charge vs. habitat and add the regression line. Since habitat is categorical, we can represent terrestrial vs aquatic on the x-axis (as 0 and 1, or simply two categories).

Plot NetCharge vs habitat (aquatic status)
- We'll use 0 for terrestrial and 1 for aquatic for plotting purposes

```R
traits$aquatic01 <- ifelse(traits$Status == "aquatic", 1, 0)
```

Save plot to a png file
```R
png(paste0("~/lab07-",MYGIT,"/NetCharge_vs_Habitat.png"), width = 600, height = 400)
plot(jitter(traits$aquatic01, amount = 0.1), traits$NetCharge,
     xlab = "Habitat (0 = terrestrial, 1 = aquatic)",
     ylab = "Myoglobin Net Charge (pH 6.5)",
     main = "Myoglobin Net Charge in Terrestrial vs Aquatic Mammals",
     pch = 19, col = "blue")
abline(a = coef(model_netcharge)[1], b = coef(model_netcharge)[2], col = "red", lwd = 2)
dev.off()
```

><img src="img/github.png" alt= “” width="20" height="20">[Github]  This file will be part of your GitHub repository.

A few notes on the plotting: 
- We jittered the x-values (aquatic01) by a small amount (0.1) just to spread the points for visibility. This is helpful if multiple species overlap at exactly 0 or 1. 
- Points are colored blue and plotted with solid circles (pch = 19). 
- The red line is drawn using the model coefficients. This line will essentially connect the average terrestrial value at x=0 to the average aquatic value at x=1. 

After running the code, open the PNG to ensure it rendered correctly. You should see two clusters of points (left for terrestrial, right for aquatic) and a line between them. Does the line slope upward? That would indicate aquatic species have higher net charges on average, as expected. The spread of points within each group shows the variation among species; the regression considers this variation in combination with the phylogeny.

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 5. You plotted Myoglobin Net Charge versus Habitat (0 = terrestrial, 1 = aquatic) and added the regression line from your model. What does an upward-sloping red line in your plot indicate?

**Testing Other Traits**
We hypothesized multiple structural characteristics might differ with habitat. Repeat the above regression analysis for other traits: 
- Secondary structure fractions: e.g., fraction helix (frac_helix) vs habitat. 
- Solvent accessibility: e.g., mean_ASA vs habitat. 
- Domain count or length: For myoglobin, all have one domain, so we expect no difference here (and no variation means no regression possible in that case). For a multi-domain gene, one could test DomainCount ~ Status.

For practice, let's try one more: Does the fraction of helical structure in myoglobin differ between aquatic and terrestrial mammals? We’ll run a similar model:

```R
model_helix <- phylolm(frac_helix ~ Status, data = traits, phy = tree, model = "BM")
summary(model_helix)
```

Examine the output. It’s likely that frac_helix (the proportion of amino acids in α-helix) will not show a significant habitat difference – myoglobin’s secondary structure is highly conserved across mammals, so we expect no systematic change for divers. A non-significant result here is a good reminder that not all traits, or genes, will show differences related to habitat. In fact, many genes will not exhibit any significant adaptive signal for aquatic vs terrestrial lifestyle. Myoglobin’s net charge is a special case where a clear adaptive pattern has been documented. It’s important to acknowledge this so we don’t fall prey to confirmation bias or assume every observed difference is meaningful.

Limitation: In our analysis, we used one trait value per species. If a species had multiple myoglobin gene copies, we would have to decide how to handle that (e.g., take the average trait value, or choose one representative sequence). Averaging across gene copies can blur the data and assume each copy contributes equally – this is a simplification. Also, our phylogenetic regression assumes the trait differences evolved along the species tree; it doesn’t account for within-species variation or gene duplications. These factors are worth mentioning when interpreting results.

#### All additional traits
Traits to test (we already did NetCharge and frac_helix as examples, but we will include them here for completeness)
```R
vars <- c("NetCharge","frac_helix","frac_sheet","frac_coil","mean_ASA",
          "DomainCountPerGene","CommonDomainLengthPerGene","NumGeneCopies")
```

 Defome a new function called fit_trait
 This function will fit the model safely, prune tree to data, and return a tidy row with the results for each variable you apply it to.
```R
fit_trait <- function(varname) {
# Keep only rows with non-missing response & Status
sub <- traits[!is.na(traits[[varname]]) & !is.na(traits$Status),
              c("Abbr","Status", varname)]
names(sub)[3] <- "y"

# If there are too few species, return an empty row and stop here
if (nrow(sub) < 5) {
  return(data.frame(
    Trait = varname, N = nrow(sub),
    Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
    Note = "Too few species or no variance"
  ))
}

# --- continue only if there are enough species ---
# Prune tree to species in sub
tr <- drop.tip(tree, setdiff(tree$tip.label, sub$Abbr))
sub <- sub[match(tr$tip.label, sub$Abbr), ]  # align rows to tree tip order

  # Guard: if response has (near) zero variance, skip
  if (sd(sub$y) < .Machine$double.eps) {
    return(data.frame(
      Trait = varname, N = nrow(sub),
      Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
      Note = "No variation in response"
    ))
  }

  # Fit phylogenetic linear model (BM)
  mod <- phylolm(y ~ Status, data = sub, phy = tr, model = "BM")
  sm  <- summary(mod)

  # Save a quick plot
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

Run the fit_trait function across all traits in vars and save a results table. 
```R
res_list <- lapply(vars, fit_trait)
res <- do.call(rbind, res_list)
write.csv(res, paste0("~/lab07-",MYGIT,"/myoglobin/myoglobin.trait_results.csv"), row.names = FALSE)
res
```

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 6. You ran phylogenetic regressions for several myoglobin traits comparing aquatic and terrestrial mammals (see table res). Which statement best summarizes the results?

### Think about it - have we learned something new? Or have we just generated a new hypothesis? 
* Learning something often leads to more questions than answers. *
In this analysis of myoglobin, CommonDomainLengthPerGene also showed a statistically significant difference between aquatic and terrestrial mammals. Unlike NetCharge, this result was not previously reported or predicted. In comparative bioinformatics, such unexpected patterns don’t necessarily mean they’re wrong—but they also don’t mean they’re already “known.” They point to possible new hypotheses that could be explored in future work. This is similar to what you will be trying to do for your own gene family.

><img src="img/brightspace.png" alt= “” width="20" height="20">[Brightspace] 7. The difference in CommonDomainLengthPerGene between aquatic and terrestrial mammals was not previously reported or predicted. What is the best way to interpret this kind of result?


# Work on your own gene family: Phylogenetic Regression between Structural Traits and Lifestyle Transitions
Repeat the above lab, using your own gene family.

><img src="img/github.png" alt= “” width="20" height="20">[Github] **Keep track of *all* the commands that you used here! (i.e. cut and paste them here, so you could recreate your work!)**.
mkdir ~/lab07-$MYGIT/SPTY2D1
R --vanilla
MYGIT <- system("echo $MYGIT", intern = TRUE)
library(ape)
tree <- read.tree(paste0("~/lab07-",MYGIT,"/species_abr.tre"))
tree <- compute.brlen(tree, method = "Grafen")
species_key <- read.csv(paste0("~/lab06-",MYGIT,"/species_key.csv"), stringsAsFactors = TRUE)
head(species_key)  
Load the Net Charge Calculations from Lab 6
charges <- read.csv(paste0("~/lab06-",MYGIT,"/SPTY2D1/net_charge_summary.csv"), header = TRUE)
colnames(charges)[3] <- "Abbr"
charges <- read.csv(paste0("~/lab06-",MYGIT,"/SPTY2D1/net_charge_summary.csv"), header = TRUE)
colnames(charges)[3] <- "Abbr"
charges <- read.table(paste0("~/lab06-",MYGIT,"/SPTY2D1/net_charges_labeled.tsv"), sep = "\t", header = TRUE)
charges_summary <- aggregate(NetCharge ~ Abbr, data = charges, FUN = mean)
head(charges_summary)
Charges_summary
Load secondary structure summary data from Lab 6
dssp <- read.csv(paste0("~/lab06-",MYGIT,"/SPTY2D1/dssp_summary.csv"), header = TRUE)
dssp_summary <- aggregate(cbind(frac_helix, frac_sheet, frac_coil, mean_ASA) ~ abbr, data = dssp, FUN = mean)
head(dssp_summary)
RPS-BLAST domain content
rps <- read.table(paste0("~/lab06-",MYGIT,"/SPTY2D1/SPTY2D1.rps-blast.out"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Rps
names(rps) <- c("qseqid","qlen","qstart","qend","evalue","stitle")
Rps
parts <- strsplit(rps$qseqid, "\\|")
rps$Abbr <- sapply(parts, `[`, 1)
rps$RefSeq <- sapply(parts, `[`, 2)
rps$AlignLength <- abs(rps$qend - rps$qstart) + 1
gene_summary <- aggregate(AlignLength ~ Abbr + RefSeq, data = rps,
                          FUN = function(x) {
                            c(count = length(x),
                              common_length = as.integer(names(sort(table(x), decreasing = TRUE)[1])))
                          })

gene_summary$DomainCountPerGene        <- gene_summary$AlignLength[, "count"]
gene_summary$CommonDomainLengthPerGene <- gene_summary$AlignLength[, "common_length"]
gene_summary$AlignLength <- NULL
species_means <- aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
                           data = gene_summary, FUN = mean)
n_gene_copies <- aggregate(RefSeq ~ Abbr, data = gene_summary,
                           FUN = function(x) length(unique(x)))
                           
names(n_gene_copies)[2] <- "NumGeneCopies"
n_gene_copies
Final per-species summary: mean domains per gene and mean domain length
domain_summary2 <- merge(
  aggregate(cbind(DomainCountPerGene, CommonDomainLengthPerGene) ~ Abbr,
            data = gene_summary, FUN = mean),
  aggregate(RefSeq ~ Abbr, data = gene_summary,
            FUN = function(x) length(unique(x))),
  by = "Abbr"
)
names(domain_summary2)[names(domain_summary2) == "RefSeq"] <- "NumGeneCopies"
domain_summary2
Combining all trait data
traits0 <- merge(charges_summary,dssp_summary,by.x = "Abbr",by.y="abbr",all=TRUE)
traits <- merge(traits0,domain_summary2, by = "Abbr",all=TRUE)
traits$Status <- species_key$aquatic[ match(traits$Abbr, species_key$abbreviation) ]
Traits
species_in_data <- traits$Abbr
tree <- drop.tip(tree, setdiff(tree$tip.label, species_in_data))
row.names(traits) <- traits$Abbr
length(tree$tip.label) 
nrow(traits)            
Phylogenetic Regression in R
library(phylolm)
traits$Status <- factor(traits$Status, levels = c("terrestrial", "aquatic"))
model_netcharge <- phylolm(NetCharge ~ Status, data = traits, phy = tree, model = "BM")
summary(model_netcharge)
Plotting the Results
traits$aquatic01 <- ifelse(traits$Status == "aquatic", 1, 0)
png(paste0("~/lab07-",MYGIT,"/NetCharge_vs_Habitat.png"), width = 600, height = 400)
plot(jitter(traits$aquatic01, amount = 0.1), traits$NetCharge,
     xlab = "Habitat (0 = terrestrial, 1 = aquatic)",
     ylab = "SPTY2D1 Net Charge (pH 6.5)",
     main = "SPTY2D1 Net Charge in Terrestrial vs Aquatic Mammals",
     pch = 19, col = "blue")
abline(a = coef(model_netcharge)[1], b = coef(model_netcharge)[2], col = "red", lwd = 2)
dev.off()
model_helix <- phylolm(frac_helix ~ Status, data = traits, phy = tree, model = "BM")
summary(model_helix)
vars <- c("NetCharge","frac_helix","frac_sheet","frac_coil","mean_ASA",
          "DomainCountPerGene","CommonDomainLengthPerGene","NumGeneCopies")
fit_trait <- function(varname) {
# Keep only rows with non-missing response & Status
sub <- traits[!is.na(traits[[varname]]) & !is.na(traits$Status),
              c("Abbr","Status", varname)]
names(sub)[3] <- "y"

# If there are too few species, return an empty row and stop here
if (nrow(sub) < 5) {
  return(data.frame(
    Trait = varname, N = nrow(sub),
    Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
    Note = "Too few species or no variance"
  ))
}

# --- continue only if there are enough species —
# Prune tree to species in sub
tr <- drop.tip(tree, setdiff(tree$tip.label, sub$Abbr))
sub <- sub[match(tr$tip.label, sub$Abbr), ]  # align rows to tree tip order

  # Guard: if response has (near) zero variance, skip
  if (sd(sub$y) < .Machine$double.eps) {
    return(data.frame(
      Trait = varname, N = nrow(sub),
      Beta_aquatic = NA, SE = NA, t = NA, p = NA, R2 = NA,
      Note = "No variation in response"
    ))
  }

  # Fit phylogenetic linear model (BM)
  mod <- phylolm(y ~ Status, data = sub, phy = tr, model = "BM")
  sm  <- summary(mod)

  # Save a quick plot
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
res_list <- lapply(vars, fit_trait)
res <- do.call(rbind, res_list)
write.csv(res, paste0("~/lab07-",MYGIT,"/SPTY2D1/SPTY2D1.trait_results.csv"), row.names = FALSE)
Res
 
><img src="img/github.png" alt= “” width="20" height="20">[Github] Did you achieve any statistically significant results? How will you follow up on these?
Number of gene copies of Caur is 2. The rest are 1.

Phylogenetic Regression in R:
Coefficients:
                Estimate     StdErr t.value   p.value    
(Intercept)    0.2275569  0.0106664 21.3340 6.917e-07 ***
Statusaquatic -0.0073411  0.0095538 -0.7684    .4714    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-squared: 0.08959      Adjusted R-squared: -0.06214 
                      Trait  N  Beta_aquatic           SE          t         p
1                 NetCharge  8 -2.317824e+04 2.554495e+04 -0.9073514 0.3991950
2                frac_helix  8 -7.341138e-03 9.553806e-03 -0.7683993 0.4714014
3                frac_sheet  8  5.262047e-04 1.488383e-03  0.3535413 0.7357735
4                 frac_coil  8  6.814933e-03 1.0269090e-02  0.6636357 0.5315802
5                  mean_ASA  8  1.110397e-04 7.119947e-04  0.1559557 0.8811828
6        DomainCountPerGene 17            NA           NA         NA        NA
7 CommonDomainLengthPerGene 17 -1.381157e+00 4.452580e+00 -0.3101926 0.7606826
8             NumGeneCopies 17  5.174899e-02 6.982987e-02  0.7410723 0.4700931
           R2                     Note
1 0.120658370                         
2 0.089590034                         
3 0.020406794                         
4 0.068382639                         
5 0.004037333                         
6          NA No variation in response
7 0.006373744                         
8 0.035319411                         

The p-value is larger than 0.05, so there is no statistical difference in the net charge between aquatic and terrestrial species beyond what would be expected from phylogenetic relatedness alone. This means that the net charge differences are by phylogenetic chance. 

><img src="img/github.png" alt= “” width="20" height="20">[Github] Which figures from today's lab might you include in your paper? Why?
I would include the plot comapring the net charges between aquatic and terrestrial mammals to show the negative correlation between  the net charge and habitat. It is evident that there is a weak correlation not only because of hte low R^2 value, but also because of how scattered the data points are. 

# Save your history and push your files into the repository.

Save your history:

```bash
history > lab7.commandhistory.txt
```
  
## Push any new files to the remote Git repository.

Add all of your results to the repository, commit them, and push them to the remote repository at GitHub.com

First, ensure you are in the main lab7 directory:
```bash
cd ~/lab07-$MYGIT
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
