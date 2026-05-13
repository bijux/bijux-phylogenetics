# Evidence 001

This is a registered cross-reference example for the full lecture script
[`PCM1_plots_signal.R`](external:lund/pcm1-plots-signal/script).

Identity:

- study id: `primate-longevity-signal`
- evidence id: `evidence-001`
- reference tools: `ape, geiger, phytools, treeio, tidytree`

Evidence contract:

- one checked-in R reference script
- one checked-in Python script using `bijux-phylogenetics`
- one machine-readable comparison index plus per-block payloads
- explicit block verdicts that report agreement or deviation without hiding gaps

Files used in this report:

- [R reference checks](../reference/primate_lifespan_signal_reference_r.R)
- [Python `bijux-phylogenetics` checks](../reference/primate_lifespan_signal_reference_bijux.py)
- [Example manifest](manifest.json)
- [R results JSON](r_reference_results.json)
- [Python results JSON](bijux_reference_results.json)
- [Comparison JSON](comparison.json)
- [Per-block payloads](./results/block-payloads/)
- [Comparative validation suite](./results/comparative_reference_validation_suite.json)
- [R ecosystem comparison](./results/r_ecosystem_comparison.json)
- [Trusted examples gallery](./results/trusted_examples_gallery.json)
- [Reviewer audit checklist](./results/reviewer_audit_checklist.json)
- [Reproducibility package](./results/reproducibility_package.json)
- [Method maturity registry](./results/method_maturity_registry.json)
- [Scientific debt register](./results/scientific_debt_register.json)

Coverage summary:

- executable script lines tracked: `185` / `185`
- uncovered executable lines: `[]`

Status summary:

- `artifact_only`: 2
- `plot_only`: 10
- `seeded_input_only`: 1
- `verified`: 10
- `verified_with_tolerance`: 6
- `workflow_only`: 1

## Block-by-Block Ledger

| Block | Script lines | Status | What was checked |
| --- | --- | --- | --- |
| environment-and-package-contract | 1-22 | workflow_only | Environment, package loading, and citation workflow |
| primate-data-preprocessing | 23-79 | verified | Raw primate preprocessing to checked-in analysis table |
| tree-import-and-pruning | 80-120 | verified | Tree import, checking, node labels, and pruning |
| processed-analysis-artifacts | 122-134 | artifact_only | Save processed files for later analysis |
| ape-plotting-basics | 141-155 | plot_only | APE plotting basics |
| ape-alternate-layouts | 157-158,160 | plot_only | APE alternate layouts: cladogram and fan |
| unrooted-tree-demo | 159 | verified | Unrooted tree demo |
| phytools-tree-plotting | 165-166 | plot_only | Phytools tree plotting exploration |
| extract-clade-node-77 | 168-170 | verified | Extract clade descended from R node 77 |
| rotate-nodes-behavior | 172-192 | verified | Rotate-nodes teaching demo |
| ggtree-tree-visualization | 197-220 | plot_only | Ggtree tree-visualization exploration |
| tip-order-alignment | 222-240 | verified | Tip-order alignment for joining data to the tree |
| ape-longevity-overlay | 242-245 | plot_only | APE tip overlay with longevity |
| treeio-node-mapping-and-join | 248-263 | verified | Treeio node mapping and joined tree-data object |
| joined-ggtree-trait-plotting | 265-276 | plot_only | Joined ggtree trait plotting |
| random-simulation-inputs | 287-288,294-295,309,315 | seeded_input_only | Random simulation scenarios |
| random-simulation-plotting | 290-318 | plot_only | Random trait plotting surfaces |
| random-signal-lambda-fits | 324-331 | verified_with_tolerance | Random-data lambda fits |
| primate-longevity-visual-inspection | 337,343-345 | plot_only | Primate longevity histogram and tip-size plot |
| primate-longevity-vector-assembly | 340-341 | verified | Primate longevity vector assembly |
| primate-lambda-fit | 347-354 | verified_with_tolerance | Primate longevity lambda fit |
| lambda-zero-visual-comparison | 357-371 | plot_only | Lambda-zero visual tree comparison |
| lambda-zero-covariance-and-lrt | 375-388 | verified_with_tolerance | Lambda-zero covariance and likelihood-ratio test |
| continuous-ancestral-point-estimates | 395-399 | verified_with_tolerance | Continuous ancestral point estimates |
| continuous-ancestral-intervals | 400 | verified_with_tolerance | Continuous ancestral 95% intervals |
| ancestral-table-assembly | 404-412 | verified | Assemble ancestral table and node mapping |
| ancestral-colored-tree-plot | 414-419 | plot_only | Ancestral colored tree plot |
| bonobo-gibbon-mrca-estimate | 421-429 | verified_with_tolerance | Bonobo/Gibbon MRCA estimate |
| lifespan-increase-counts | 431-448 | verified | How many times lifespan increased across primates |
| final-workspace-artifact | 450-451 | artifact_only | Save final analysis workspace |

## Block Sections

### environment-and-package-contract Environment, package loading, and citation workflow

- status: `workflow_only`
- script lines: `1-22`
- verdict: This setup block is documented for reproducibility, but it is not an analysis-equivalence target.

**R Lecture Block**

```r
#***********************************************************************************
#BIOR90: Phylogenetic Comparative Analyses ----
#Trait evolution by Charlie Cornwallis
#PCM1 Script to accompany lecture
#***********************************************************************************

#***********************************************************************************
#Packages ----
#***********************************************************************************
#pacman::p_load makes loading packages easier
pacman::p_load(openxlsx,tidyverse,RColorBrewer,ape,MCMCglmm,picante,geiger,phytools,ggtree, treeio,ggimage)

#If packages don't install, google for solutions. 
#ggtree often has this problem that can usually be solved with
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree",force = TRUE)

#If you need to citing R and R packages in your work you can check how to do this
citation() #R
citation("ape") #Specific package
```

**Python `bijux-phylogenetics` Block**

```python
# Environment/setup block.
# The credibility report records the exact R package versions used on the
# reference side and runs the Python checks in the repository environment.
```

**Comparison Evidence**

- R package versions recorded: `ape, geiger, openxlsx, phytools, tidytree, treeio`
- Python reference script: `evidence-book/studies/primate-longevity-signal/reference/primate_lifespan_signal_reference_bijux.py`

Raw comparison payload: [`block-payloads/environment-and-package-contract.json`](./results/block-payloads/environment-and-package-contract.json)

### primate-data-preprocessing Raw primate preprocessing to checked-in analysis table

- status: `verified`
- script lines: `23-79`
- verdict: The R reference reconstruction from the workbook matches the checked-in processed CSV.

**R Lecture Block**

```r
#***********************************************************************************
#1. Data importing, checking and manipulation ----
#***********************************************************************************
#If you started R studio using the Rproject file then the correct working directory should be set automatically and work on any computer

#Data on Primates - read in from excel
primate_raw <- readWorkbook("Data/primate_raw.xlsx")

#***************************
#Checking how the data looks and variable classifications
#***************************
#What are the column names?
names(primate_raw)

#How does the dataset look?
View(primate_raw)

#How are variables classified?
str(primate_raw) #Oh dear many are incorrectly classified

#Solution: - Reclassify manually e.g.
primate_raw$sex_dimorphism<-as.factor(primate_raw$sex_dimorphism)

#Or can change several columns at once
primate_raw <- primate_raw %>% mutate(across(body_mass:social_group_size, ~as.numeric(.)))

str(primate_raw) #Has it worked?

#***************************
#Check for missing data (NAs) and multiple measurements per species?
#***************************
#How many NAs are there per variable?
primate_raw %>% summarise(across(order:mating_system, ~sum(is.na(.))))

#Show rows with missing data (not (!) complete cases)
primate_raw %>% filter(!complete.cases(.))

#Are their multiple measurements for species?
primate_raw %>% group_by(species) %>% 
                summarise(n=n())  %>% 
                dplyr::filter(n>1)

#Lets check out this species with multiple measurements and missing values
primate_raw %>% filter(species == "Galago_alleni")

#Solution - summarise data and remove NAs at the same time
primate <- primate_raw %>% group_by(family, species) %>% 
                           summarise(across(body_mass:social_group_size, 
                                      .fns = ~mean(.x,na.rm=T), #fns stands for function which in this case is the "mean" and "na.rm=T" removes missing values
                                      .names = "{col}"),#This could be used to rename columns if you wanted
                                      sex_dimorphism=sex_dimorphism[!is.na(sex_dimorphism)], #Sexual dimorphism stays the same but removes NAs (!is.na())
                                      mating_system=mating_system[!is.na(mating_system)]) #Mating system stays the same but removes NAs

primate %>% group_by(species) %>% 
  summarise(n=n())  %>% 
  filter(n>1)
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path
import csv

traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")
with traits_path.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

species = [row["species"] for row in rows]
assert len(rows) == 75
assert len(set(species)) == 75
```

**Comparison Evidence**

- processed rows: `75`
- processed species: `75`
- checked-in `primate.csv` matches workbook-derived reference: `True`

Raw comparison payload: [`block-payloads/primate-data-preprocessing.json`](./results/block-payloads/primate-data-preprocessing.json)

### tree-import-and-pruning Tree import, checking, node labels, and pruning

- status: `verified`
- script lines: `80-120`
- verdict: The checked-in trimmed tree matches the R-pruned reference tree.

**R Lecture Block**

```r
#***********************************************************************************
#2. Tree importing, checking and manipulation ----
#***********************************************************************************
#Read in Tree using Package ape, read.tree multipurpose function, for nexus, newick and most tree formats
primatetree<-read.tree("Data/primatetree.nex")

#See what the file contains
primatetree

#You can access these with the variable operator $
primatetree$tip.label

#Is the tree rooted?
is.rooted(primatetree)

#Does it have polytomies
is.binary(primatetree)

#Is the tree ultrametric? Ultrametric means that all tips have same distance to the root (genetic distances between all species are relative rather than absolute), which is important for some analyses
is.ultrametric(primatetree)

#Does it have node labels?
primatetree$node.label

#Lets label the nodes so we can keep track of ancestors
primatetree <- makeNodeLabel(primatetree)

#Check correspondence between data and tree

#Are all species in the dataset in the tree?
table(primate$species %in% primatetree$tip.label)

#Are all tree tips in the dataset?
table(primatetree$tip.label %in% primate$species)

#Which are missing?
missingtips<-primatetree$tip.label[primatetree$tip.label %in% primate$species == "FALSE"]
missingtips

#Remove these species from the tree using drop.tip
primatetree<-ape::drop.tip(primatetree,missingtips)
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path
import csv

from bijux_phylogenetics import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa

original_tree = Path("PCM1_plots_signal/Lecture/R/Data/primatetree.nex")
trimmed_tree = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_csv = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

inspection = inspect_tree_path(trimmed_tree)
validation = validate_tree_path(trimmed_tree)
requested_taxa = [row["species"] for row in csv.DictReader(traits_csv.open())]
_, pruning_report = prune_tree_to_requested_taxa(original_tree, requested_taxa=requested_taxa)
```

**Comparison Evidence**

- original tree tips: `77`
- trimmed tree tips: `75`
- removed taxa: `Lagothrix_lagotricha, Eulemur_fulvus`
- checked-in trimmed tree matches R-trimmed reference: `True`

Raw comparison payload: [`block-payloads/tree-import-and-pruning.json`](./results/block-payloads/tree-import-and-pruning.json)

### processed-analysis-artifacts Save processed files for later analysis

- status: `artifact_only`
- script lines: `122-134`
- verdict: This block writes the processed CSV and trimmed tree that the Python evidence pass later consumes.

**R Lecture Block**

```r
#***********************************************************************************
#* Saving processed files ----
#***********************************************************************************
#After spending time sorting out the files you may want to save the processed files for using / importing later on for analyses

#We can easily write & read data and trees to specific files, for example:
write.csv(primate,"Data/primate.csv")
write.tree(primatetree,"Data/trimmed_primatetree.nex")

#Or you can save packages of files in one RData object
save(primate,primatetree, file = "Data/primate.RData") # I prefer this option as it creates less files, less coding and is simpler to keep track of

save.image("Data/primate.RData") #Or all files present in your environment. Simple but trash can accummulate.
```

**Python `bijux-phylogenetics` Block**

```python
# Artifact block.
# The report consumes the checked-in processed files written by the R workflow:
# `Data/primate.csv` and `Data/trimmed_primatetree.nex`.
```

**Comparison Evidence**

- processed CSV written by R: `external:lund/pcm1-plots-signal/data/primate.csv`
- trimmed tree written by R: `external:lund/pcm1-plots-signal/data/trimmed_primatetree.nex`

Raw comparison payload: [`block-payloads/processed-analysis-artifacts.json`](./results/block-payloads/processed-analysis-artifacts.json)

### ape-plotting-basics APE plotting basics

- status: `plot_only`
- script lines: `141-155`
- verdict: This is a visual exploration block; the current report does not claim figure-equivalence for base `ape` plots.

**R Lecture Block**

```r
par(mar=c(0,0,0,0))
plot(primatetree)

#Change size of tip labels
plot(primatetree,cex=0.5)

#Add node numbers
nodelabels(cex=0.5) 

#Add tip numbers
tiplabels(cex=0.5) 

#Offset labels
plot(primatetree,cex=0.5,label.offset = 0.1)
tiplabels(cex=0.5)
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# `bijux-phylogenetics` is not yet claiming rendered-figure equivalence for the
# base `ape` tree plotting surface.
```

**Comparison Evidence**

- reason: This is a visual exploration block; the current report does not claim figure-equivalence for base `ape` plots.

Raw comparison payload: [`block-payloads/ape-plotting-basics.json`](./results/block-payloads/ape-plotting-basics.json)

### ape-alternate-layouts APE alternate layouts: cladogram and fan

- status: `plot_only`
- script lines: `157-158,160`
- verdict: These layout variants are tracked as visual surfaces only.

**R Lecture Block**

```r
#Can also plot in different ways
plot(primatetree,cex=0.5,type="cladogram",main="cladogram")

# ...

plot(primatetree,cex=0.5,type="fan",main="fan")
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# The report tracks alternate `ape` layouts separately from numerical checks.
```

**Comparison Evidence**

- reason: These layout variants are tracked as visual surfaces only.

Raw comparison payload: [`block-payloads/ape-alternate-layouts.json`](./results/block-payloads/ape-alternate-layouts.json)

### unrooted-tree-demo Unrooted tree demo

- status: `verified`
- script lines: `159`
- verdict: Both sides produce an unrooted representation with the same tip set and 75 tips.

**R Lecture Block**

```r
plot(unroot(primatetree),cex=0.5,type="unrooted",main="unrooted")
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics.core.topology import unroot_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
unrooted_tree, report = unroot_tree(tree_path)
assert unrooted_tree.tip_count == 75
assert len(unrooted_tree.root.children) == 3
```

**Comparison Evidence**

- R unrooted tip count: `75`
- bijux unrooted tip count: `75`
- bijux root child count after unrooting: `3`

Raw comparison payload: [`block-payloads/unrooted-tree-demo.json`](./results/block-payloads/unrooted-tree-demo.json)

### phytools-tree-plotting Phytools tree plotting exploration

- status: `plot_only`
- script lines: `165-166`
- verdict: The `phytools::plotTree()` surface is tracked, but no rendered-figure claim is made here.

**R Lecture Block**

```r
plotTree(primatetree,fsize=0.5)
nodelabels()
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# `phytools::plotTree()` exploration is tracked here, but figure rendering is
# not part of the current equivalence claim.
```

**Comparison Evidence**

- reason: The `phytools::plotTree()` surface is tracked, but no rendered-figure claim is made here.

Raw comparison payload: [`block-payloads/phytools-tree-plotting.json`](./results/block-payloads/phytools-tree-plotting.json)

### extract-clade-node-77 Extract clade descended from R node 77

- status: `verified`
- script lines: `168-170`
- verdict: The descendant tip set matches exactly when compared by stable taxon signature rather than recycled node labels.

**R Lecture Block**

```r
#now extract the clade descended from node #77
tt77<-extract.clade(primatetree,77)
plotTree(tt77,fsize=0.5)
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
tree = load_tree(tree_path)
internal_clades = {
    "|".join(sorted(node_descendant_taxa(node))): sorted(node_descendant_taxa(node))
    for node in tree.iter_nodes()
    if not node.is_leaf()
}

# Compared by descendant-taxon identity rather than relying on a recycled label.
assert any(len(taxa) == 54 for taxa in internal_clades.values())
```

**Comparison Evidence**

- R source node: numeric `77`, label `Node2`
- R extracted tip count: `54`
- bijux extracted tip count: `54`
- matched bijux internal node label: `Node2`
- exact descendant taxon set match: `True`

Raw comparison payload: [`block-payloads/extract-clade-node-77.json`](./results/block-payloads/extract-clade-node-77.json)

### rotate-nodes-behavior Rotate-nodes teaching demo

- status: `verified`
- script lines: `172-192`
- verdict: The child-order rotation results match the R `rotateNodes` tip order for both the single-node and all-node variants.

**R Lecture Block**

```r
#rotate the tree around node #130
rt.130<-rotateNodes(primatetree,130)
plotTree(rt.130,fsize=0.5)

#Compare trees
par(mfrow=c(1,2))
plotTree(primatetree,fsize=0.5)
nodelabels()
plotTree(rt.130,fsize=0.5)
nodelabels()

# now rotate all nodes
rt.all<-rotateNodes(primatetree,"all")
plotTree(rt.all,fsize=0.5)
plotTree(primatetree,fsize=0.5)

# check if tree & rt.all are equal
all.equal(primatetree,rt.130)

# check if tree & tt77 are equal
all.equal(primatetree,tt77)
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics.core.topology import rotate_all_internal_nodes, rotate_named_node

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
rotated_once, rotate_once_report = rotate_named_node(tree_path, clade_name="Node56")
rotated_all, rotate_all_report = rotate_all_internal_nodes(tree_path)

assert rotated_once.tip_count == 75
assert rotated_all.tip_count == 75
```

**Comparison Evidence**

- single-node rotation label in R: `Node56`
- single-node rotation label in bijux: `Node56`
- single-node tip order match: `True`
- all-node tip order match: `True`

Raw comparison payload: [`block-payloads/rotate-nodes-behavior.json`](./results/block-payloads/rotate-nodes-behavior.json)

### ggtree-tree-visualization Ggtree tree-visualization exploration

- status: `plot_only`
- script lines: `197-220`
- verdict: These `ggtree` examples are tracked as visual surfaces only.

**R Lecture Block**

```r
ggtree(primatetree)
ggtree(primatetree,layout="circular")
ggtree(primatetree, branch.length="none") #You can easily turn your tree into a cladogram with the branch.length = “none” parameter.

#add tiplabels
ggtree(primatetree) + geom_tiplab()
ggtree(primatetree,layout="circular") + geom_tiplab()

#add scale
ggtree(primatetree) + geom_treescale(offset=-2,linesize=1)

#add nodelabels
ggtree(primatetree) + geom_text(aes(label=node),size=3)

#add titles
ggtree(primatetree) + ggtitle("A great title")

#annotate selected clades in trees.
ggtree(primatetree) + xlim(NA, 2)+ geom_cladelabel(node=130, label="test label", align=T, color='red') +
geom_cladelabel(node=78, label="another clade", align=T, color='blue') 

#Can add layers similar to ggplot
p <- ggtree(primatetree) + geom_text(aes(label=node),size=3)
p + geom_cladelabel(node=130, label="test label", align=T, color='red') +geom_cladelabel(node=78, label="another clade", align=T, color='blue')
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# `ggtree` rendering examples are tracked here as visual surfaces, not as
# numerical equivalence claims.
```

**Comparison Evidence**

- reason: These `ggtree` examples are tracked as visual surfaces only.

Raw comparison payload: [`block-payloads/ggtree-tree-visualization.json`](./results/block-payloads/ggtree-tree-visualization.json)

### tip-order-alignment Tip-order alignment for joining data to the tree

- status: `verified`
- script lines: `222-240`
- verdict: The aligned species order matches the trimmed tree tip order.

**R Lecture Block**

```r
#************************************
#* Combine data and trees ----
#************************************
#** ape ----
#Lets do this with ape - need to make sure data is in the same order as tree
#Make a dataset with order of tips
tip_order <- data.frame(tip=primatetree$tip.label,
                        order=1:length(primatetree$tip.label))
head(tip_order)

#Add this to the primate dataframe
primate$tip_order <- tip_order$order[match(primate$species,tip_order$tip)]

#Now you can sort data frame according to tip order
primate <- primate %>% arrange(tip_order)

#Has it worked?
head(primate)
primatetree$tip.label[1:6]
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
assert readiness.analysis_taxa[:6] == readiness.analysis_taxa[:6]
```

**Comparison Evidence**

- R aligned species exactly follow tip order: `True`
- bijux aligned species exactly follow tip order: `True`
- first six aligned species match across tools: `True`

Raw comparison payload: [`block-payloads/tip-order-alignment.json`](./results/block-payloads/tip-order-alignment.json)

### ape-longevity-overlay APE tip overlay with longevity

- status: `plot_only`
- script lines: `242-245`
- verdict: This is a rendered trait-overlay surface and is tracked separately from ordering correctness.

**R Lecture Block**

```r
#Now we can plot the data on the tree
par(mar=c(0,0,0,0))
plot(primatetree, show.tip.label=T, label.offset=0.05,cex=0.6)
tiplabels(pch=21, bg="grey38",cex=primate$longevity/400) # cex = size in plot function. need to /300 otherwise symbols are too big
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# The `ape` tip overlay is kept separate from the data/tree ordering proof.
```

**Comparison Evidence**

- reason: This is a rendered trait-overlay surface and is tracked separately from ordering correctness.

Raw comparison payload: [`block-payloads/ape-longevity-overlay.json`](./results/block-payloads/ape-longevity-overlay.json)

### treeio-node-mapping-and-join Treeio node mapping and joined tree-data object

- status: `verified`
- script lines: `248-263`
- verdict: Representative node ids align and the joined object size is consistent with the 75-taxon dataset.

**R Lecture Block**

```r
#** ggtree ----
#Lets do this with ggtree and the treeio package
primate$node<-nodeid(primatetree,primate$species)

#Lets see what this has done - in one step it has done what we did above
View(primate)

#can now use this variable to combine tree and data into one object
p_tree_data <- treeio::full_join(primatetree, primate, by = 'node')

#Create an S4 object = different slots for different types of data
p_tree_data

#You can still access the data and tree e.g.
p_tree_data@phylo$tip.label
View(p_tree_data@extraInfo)
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait, summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
summary = summarize_numeric_trait(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)

assert readiness.analysis_taxa[:6] == summary.taxa[:6]
```

**Comparison Evidence**

- representative node ids in R: `{"Pan_paniscus": 33, "Hylobates_lar": 29, "Node32": 107}`
- representative node ids in bijux: `{"Pan_paniscus": 33, "Hylobates_lar": 29, "Node32": 107}`
- R joined object tip count: `75`
- Python analysis taxon count: `75`

Raw comparison payload: [`block-payloads/treeio-node-mapping-and-join.json`](./results/block-payloads/treeio-node-mapping-and-join.json)

### joined-ggtree-trait-plotting Joined ggtree trait plotting

- status: `plot_only`
- script lines: `265-276`
- verdict: Joined-data `ggtree` figures are tracked here as visual outputs only.

**R Lecture Block**

```r
#Now you have this combined object you can do more complicated plots
ggtree(p_tree_data) +
  geom_tiplab(aes(colour="red"),offset=0.1) +
  geom_tippoint(aes(x=x+0.,size=longevity),fill="grey",colour="black",shape=21)

#Plot area needs adjusting and perhaps some other adjustments could make it nicer
ggtree(p_tree_data) +
  geom_tiplab(colour="red",offset=0.1,size=2) +
  geom_tippoint(aes(x=x+0.05,size=longevity),fill="grey",colour="black",shape=21)+
  xlim(0,1.5)+
  ylim(0,75)+
  guides(size=guide_legend(title="Longevity (months)"),color="none")
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# Joined `ggtree` trait figures are tracked here without claiming pixel-level
# rendering equivalence yet.
```

**Comparison Evidence**

- reason: Joined-data `ggtree` figures are tracked here as visual outputs only.

Raw comparison payload: [`block-payloads/joined-ggtree-trait-plotting.json`](./results/block-payloads/joined-ggtree-trait-plotting.json)

### random-simulation-inputs Random simulation scenarios

- status: `seeded_input_only`
- script lines: `287-288,294-295,309,315`
- verdict: For credibility, the report freezes these R-generated simulation inputs with `set.seed(1)` and reuses the resulting artifacts on both sides.

**R Lecture Block**

```r
random_tree<-rcoal(30)
random_data<-rTraitCont(random_tree, model="BM", sigma = 0.5,root.value=1)

# ...

random_data2<-rTraitCont(random_tree, model="BM", sigma = 0.5,root.value=10)
random_data3<-rTraitCont(random_tree, model="BM", sigma = 5,root.value=1)

# ...

random_data4<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=0)#alpha is the strength of attraction to an optimum. 0 = no attractionso should be the same as the BM model.

# ...

random_data5<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=5)
```

**Python `bijux-phylogenetics` Block**

```python
# Seeded-input block.
# The report freezes the random tree and trait tables on the R side and then
# reuses those exact artifacts for cross-tool comparison.
```

**Comparison Evidence**

- shared random tree artifact: `evidence-book/studies/primate-longevity-signal/evidence-001/results/random_tree_seed1.nwk`
- random examples frozen from R: `random_data, random_data2, random_data3, random_data4, random_data5`
- random tree tip count: `30`

Raw comparison payload: [`block-payloads/random-simulation-inputs.json`](./results/block-payloads/random-simulation-inputs.json)

### random-simulation-plotting Random trait plotting surfaces

- status: `plot_only`
- script lines: `290-318`
- verdict: These are visual teaching plots and are tracked separately from the simulation inputs and signal fits.

**R Lecture Block**

```r
plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data*2, bg="orange")

#Have a play around with the root value and sigma (standarde deviation in trait)
random_data2<-rTraitCont(random_tree, model="BM", sigma = 0.5,root.value=10)
random_data3<-rTraitCont(random_tree, model="BM", sigma = 5,root.value=1)

par(mfrow=c(1,3))
plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data, bg="orange")

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data2, bg="orange")

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data3, bg="orange")

#Different traits may evolve in different ways aka 'the mode of evolution may be different' 
par(mfrow=c(1,2))
random_data4<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=0)#alpha is the strength of attraction to an optimum. 0 = no attractionso should be the same as the BM model.

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data4, bg="orange")

#Lets increase the alpha and see what happens
random_data5<-rTraitCont(random_tree, model="OU", sigma = 5,root.value=1,alpha=5)

plot(random_tree, show.tip.label=F)
tiplabels(pch=21, cex=random_data5, bg="orange")
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# Random trait plotting is separated from the simulation and fit checks.
```

**Comparison Evidence**

- reason: These are visual teaching plots and are tracked separately from the simulation inputs and signal fits.

Raw comparison payload: [`block-payloads/random-simulation-plotting.json`](./results/block-payloads/random-simulation-plotting.json)

### random-signal-lambda-fits Random-data lambda fits

- status: `verified_with_tolerance`
- script lines: `324-331`
- verdict: The report checks the explicit random-data fit calls and extends the implied checks to the other generated examples.

**R Lecture Block**

```r
Lambda_random_data<-fitContinuous(phy=random_tree, dat = random_data,model = "lambda")
Lambda_random_data

#Try with random_data2 and random_data3. What do you find?

#What happens if you do it for random_data4 and random_data5?
Lambda_random_data5<-fitContinuous(phy=random_tree, dat = random_data5,model = "lambda",control = list(niter = 1000))
Lambda_random_data5
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda

out_dir = Path("evidence-book/studies/primate-longevity-signal/evidence-001")
random_tree = out_dir / "random_tree_seed1.nwk"

for name in ["random_data", "random_data2", "random_data3", "random_data4", "random_data5"]:
    report = estimate_pagels_lambda(
        random_tree,
        out_dir / f"{name}.csv",
        trait="value",
        taxon_column="species",
        fine_step=0.001,
    )
```

**Comparison Evidence**

- `random_data` lambda: R `0.9999681477851246` vs bijux `1.0`
- `random_data2` lambda: R `1` vs bijux `1.0`
- `random_data3` lambda: R `1` vs bijux `1.0`
- `random_data4` lambda: R `1` vs bijux `1.0`
- `random_data5` lambda: R `1` vs bijux `1.0`

Raw comparison payload: [`block-payloads/random-signal-lambda-fits.json`](./results/block-payloads/random-signal-lambda-fits.json)

### primate-longevity-visual-inspection Primate longevity histogram and tip-size plot

- status: `plot_only`
- script lines: `337,343-345`
- verdict: These are visual inspection surfaces only.

**R Lecture Block**

```r
hist(primate$longevity)

# ...

par(mfrow=c(1,1))
plot(primatetree, show.tip.label=T, label.offset=0.05,cex=0.6)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# The histogram and primate tip-size plot are tracked as visual outputs only.
```

**Comparison Evidence**

- reason: These are visual inspection surfaces only.

Raw comparison payload: [`block-payloads/primate-longevity-visual-inspection.json`](./results/block-payloads/primate-longevity-visual-inspection.json)

### primate-longevity-vector-assembly Primate longevity vector assembly

- status: `verified`
- script lines: `340-341`
- verdict: The longevity vector is aligned to the trimmed tree tip order.

**R Lecture Block**

```r
primateLL<-as.numeric(primate$longevity[match(primatetree$tip.label,primate$species)])
names(primateLL)<-primatetree$tip.label
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
assert len(readiness.analysis_taxa) == 75
```

**Comparison Evidence**

- R vector length: `75`
- bijux vector length: `75`
- R names match tip order: `True`
- bijux names match tip order: `True`

Raw comparison payload: [`block-payloads/primate-longevity-vector-assembly.json`](./results/block-payloads/primate-longevity-vector-assembly.json)

### primate-lambda-fit Primate longevity lambda fit

- status: `verified_with_tolerance`
- script lines: `347-354`
- verdict: The `bijux-phylogenetics` lambda estimate is within a small numerical tolerance of the R fit.

**R Lecture Block**

```r
#*******************************************
#* Pagels Lamda using fitContinuous in Geiger ----
#*******************************************

#Estimate Lambda
Lambda_LL<-fitContinuous(phy= primatetree, dat = primateLL,model = "lambda")
Lambda_LL
Lambda_LL$opt$lnL
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

lambda_report = estimate_pagels_lambda(
    tree_path, traits_path, trait="longevity", taxon_column="species", fine_step=0.001
)
```

**Comparison Evidence**

- lambda absolute difference: `0.00026987974695991124`
- log-likelihood absolute difference: `5.460539682644594e-06`

Raw comparison payload: [`block-payloads/primate-lambda-fit.json`](./results/block-payloads/primate-lambda-fit.json)

### lambda-zero-visual-comparison Lambda-zero visual tree comparison

- status: `plot_only`
- script lines: `357-371`
- verdict: The side-by-side real-tree versus lambda=0 plots are tracked as visual outputs only.

**R Lecture Block**

```r
#*******************************************
#* Testing if Pagels Lamda  is significantly different from 0 (start phylogeny) ----
#*******************************************
##First transform the tree to Lambda = 0 
primatetreeL0<-rescale(primatetree, "lambda", 0)

##See how they look
par(mfrow=c(1,2),mar=c(0,0,2,0))
plot(primatetree,cex=0.4,label.offset=0.1)
title("Real Phylogeny", cex=28)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)

plot(primatetreeL0, cex=0.4,label.offset=0.1)
title("Lambda = 0", cex=28)
tiplabels(pch=21, bg="grey38",cex=primateLL/300)
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# The real-tree versus lambda=0 tree rendering is tracked separately from the
# covariance and likelihood-ratio checks.
```

**Comparison Evidence**

- reason: The side-by-side real-tree versus lambda=0 plots are tracked as visual outputs only.

Raw comparison payload: [`block-payloads/lambda-zero-visual-comparison.json`](./results/block-payloads/lambda-zero-visual-comparison.json)

### lambda-zero-covariance-and-lrt Lambda-zero covariance and likelihood-ratio test

- status: `verified_with_tolerance`
- script lines: `375-388`
- verdict: The covariance surface and lambda-vs-zero test agree within numerical tolerance.

**R Lecture Block**

```r
vcv.phylo(primatetreeL0)[1:3,1:3]
vcv.phylo(primatetree)[1:3,1:3]

#Calculate the loglikelihood when Lambda is set to 0
Lambda_LL0<-fitContinuous(phy= primatetreeL0, dat = primateLL,model = "lambda")
Lambda_LL0$opt$lnL

#Use a log-likelihood ratio test to calculate if Lambda is different from 0
##Calculate the  log likelihood ratio
LLDiff0 <- -2*(Lambda_LL0$opt$lnL - Lambda_LL$opt$lnL)
LLDiff0

##Perform Chi sq test
pchisq(LLDiff0, df=1,lower.tail = FALSE)
```

**Python `bijux-phylogenetics` Block**

```python
from math import erfc, sqrt
from pathlib import Path

from bijux_phylogenetics import estimate_pagels_lambda
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    lambda_transform_covariance,
    load_comparative_dataset,
)

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

lambda_report = estimate_pagels_lambda(
    tree_path, traits_path, trait="longevity", taxon_column="species", fine_step=0.001
)
ll_diff0 = -2.0 * (lambda_report.null_log_likelihood - lambda_report.log_likelihood)
p_value = erfc(sqrt(ll_diff0 / 2.0))
dataset = load_comparative_dataset(tree_path, traits_path, trait="longevity", taxon_column="species")
covariance = build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
lambda0_covariance = lambda_transform_covariance(covariance, 0.0)
```

**Comparison Evidence**

- likelihood-ratio absolute difference: `1.092107925160235e-05`
- p-value absolute difference: `8.374938977420018e-16`
- lambda=0 covariance top-left 3x3 max diff: `9.99999993922529e-09`
- real covariance top-left 3x3 max diff: `0.0`

Raw comparison payload: [`block-payloads/lambda-zero-covariance-and-lrt.json`](./results/block-payloads/lambda-zero-covariance-and-lrt.json)

### continuous-ancestral-point-estimates Continuous ancestral point estimates

- status: `verified_with_tolerance`
- script lines: `395-399`
- verdict: Clade-aligned ancestral point estimates match to floating-point noise.

**R Lecture Block**

```r
long_ace<-ace(primate$longevity,primatetree,type="continuous",method="pic")

#Lets look at the output
names(long_ace$ace) #nodes names
long_ace$ace #ancestral estimates
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
```

**Comparison Evidence**

- shared internal clades compared: `74`
- max point-estimate absolute difference: `5.115907697472721e-13`

Raw comparison payload: [`block-payloads/continuous-ancestral-point-estimates.json`](./results/block-payloads/continuous-ancestral-point-estimates.json)

### continuous-ancestral-intervals Continuous ancestral 95% intervals

- status: `verified_with_tolerance`
- script lines: `400`
- verdict: The Brownian/PIC confidence-interval surface now matches the R reference to floating-point noise.

**R Lecture Block**

```r
long_ace$CI95 #)95% CIs
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
# Point estimates agree with R; the remaining open question is interval equivalence.
```

**Comparison Evidence**

- shared internal clades compared: `74`
- max lower-95 difference: `5.684341886080801e-13`
- max upper-95 difference: `5.115907697472721e-13`

Raw comparison payload: [`block-payloads/continuous-ancestral-intervals.json`](./results/block-payloads/continuous-ancestral-intervals.json)

### ancestral-table-assembly Assemble ancestral table and node mapping

- status: `verified`
- script lines: `404-412`
- verdict: The ancestral table assembly is consistent: 75 tip rows plus 74 internal rows on both sides.

**R Lecture Block**

```r
primate_ace_long<-data.frame(species=c(primate$species, #species names
                                       names(long_ace$ace)), #nodes names
                             longevity=c(primate$longevity, #longevity values
                                         long_ace$ace)) #estimated ancestral longevity

#Combine data and tree
primate_ace_long$node<-nodeid(primatetree,primate_ace_long$species)#Convert tips labels to numbers
View(primate_ace_long)#how does it look?
primate_ace_long_tree <- full_join(primatetree,primate_ace_long, by = 'node')#Combine with tree in ggtree object
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics import summarize_numeric_trait_readiness

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

readiness = summarize_numeric_trait_readiness(
    tree_path, traits_path, trait="longevity", taxon_column="species"
)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)

tip_rows = len(readiness.analysis_taxa)
internal_rows = sum(1 for row in ancestral.estimates if not row.is_tip)
assert tip_rows + internal_rows == len(ancestral.estimates)
```

**Comparison Evidence**

- R table rows: tips `75` + internal `74` = `149`
- bijux table rows: tips `75` + internal `74` = `149`

Raw comparison payload: [`block-payloads/ancestral-table-assembly.json`](./results/block-payloads/ancestral-table-assembly.json)

### ancestral-colored-tree-plot Ancestral colored tree plot

- status: `plot_only`
- script lines: `414-419`
- verdict: The ancestral-state branch-color figure is tracked as a visual rendering surface only.

**R Lecture Block**

```r
ggtree(primate_ace_long_tree)+
  geom_tree(aes(color=primate_ace_long_tree@extraInfo$longevity),size=1) + # changes colour of branches to indicate ancestral estimates
  geom_tiplab(aes(x=x+0.08),size=2) + #plot tiplabels
  xlim(0,1.6)+
  guides(color=guide_legend(title="Lifespan"),size=guide_legend(title="Lifespan"))+
  theme(legend.position = c(.1, .85))
```

**Python `bijux-phylogenetics` Block**

```python
# Plot-only block.
# The ancestral colored tree figure is tracked as a rendering surface.
```

**Comparison Evidence**

- reason: The ancestral-state branch-color figure is tracked as a visual rendering surface only.

Raw comparison payload: [`block-payloads/ancestral-colored-tree-plot.json`](./results/block-payloads/ancestral-colored-tree-plot.json)

### bonobo-gibbon-mrca-estimate Bonobo/Gibbon MRCA estimate

- status: `verified_with_tolerance`
- script lines: `421-429`
- verdict: The MRCA clade and the ancestral point estimate match.

**R Lecture Block**

```r
#Lets examine the estimated lifespans of Bonobos and Gibbons and their the most recent common ancestor
primate_ace_long %>% filter(species == "Pan_paniscus")
primate_ace_long %>% filter(species == "Hylobates_lar")

#Node number of most recent common ancestor
getMRCA(primatetree, tip=c("Pan_paniscus", "Hylobates_lar"))

#Estimate of node
primate_ace_long %>% filter(node == "107")
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

tree = load_tree(tree_path)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
estimate_by_node = {row.node: row.estimate for row in ancestral.estimates}

target = {"Pan_paniscus", "Hylobates_lar"}
mrca = min(
    (node for node in tree.iter_nodes() if target <= set(node_descendant_taxa(node))),
    key=lambda node: len(node_descendant_taxa(node)),
)
mrca_signature = node_signature(mrca)
mrca_estimate = estimate_by_node[mrca_signature]
```

**Comparison Evidence**

- MRCA node in R: `107`
- MRCA clade signature in bijux: `Hylobates_lar|Hylobates_pileatus|Nomascus_concolor|Pan_paniscus|Symphalangus_syndactylus`
- MRCA estimate in R: `520.0784230077601`
- MRCA estimate in bijux: `520.07842300776`

Raw comparison payload: [`block-payloads/bonobo-gibbon-mrca-estimate.json`](./results/block-payloads/bonobo-gibbon-mrca-estimate.json)

### lifespan-increase-counts How many times lifespan increased across primates

- status: `verified`
- script lines: `431-448`
- verdict: The branch-wise increase counts match the direct R reference path.

**R Lecture Block**

```r
#***********************************************************************************
#* How many times has lifespan increased across primates? ----
#***********************************************************************************

#Create a dataset of ancestors and descendants
tree_df<-tidytree::as_tibble(primatetree)

#Add in ancestral estimates matching by node id's
tree_df <- tree_df %>% mutate(ancestor_long=primate_ace_long$longevity[match(parent,primate_ace_long$node)],
                              descendant_long=primate_ace_long$longevity[match(node,primate_ace_long$node)])

tree_df

#How many increases are there?
tree_df <- tree_df %>% mutate(diff_long = descendant_long - ancestor_long) #Calculate the difference between ancestor and descendant
tree_df %>% filter(diff_long > 0) %>% summarise(n=n()) #Count the number of increases

tree_df %>% filter(diff_long > 12) %>% summarise(n=n()) #Count the number of increases
```

**Python `bijux-phylogenetics` Block**

```python
from pathlib import Path

from bijux_phylogenetics import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.io.trees import load_tree

tree_path = Path("PCM1_plots_signal/Lecture/R/Data/trimmed_primatetree.nex")
traits_path = Path("PCM1_plots_signal/Lecture/R/Data/primate.csv")

tree = load_tree(tree_path)
ancestral = reconstruct_continuous_ancestral_states(
    tree_path, traits_path, trait="longevity", taxon_column="species", model="brownian"
)
estimate_by_node = {row.node: row.estimate for row in ancestral.estimates}

increase_count = 0
increase_gt12_count = 0

def visit(node, parent=None):
    global increase_count, increase_gt12_count
    current = node_signature(node)
    if parent is not None:
        diff = estimate_by_node[current] - estimate_by_node[parent]
        if diff > 0:
            increase_count += 1
        if diff > 12:
            increase_gt12_count += 1
    for child in node.children:
        visit(child, current)

visit(tree.root)
```

**Comparison Evidence**

- increases > 0: R `72` vs bijux `72`
- increases > 12: R `55` vs bijux `55`

Raw comparison payload: [`block-payloads/lifespan-increase-counts.json`](./results/block-payloads/lifespan-increase-counts.json)

### final-workspace-artifact Save final analysis workspace

- status: `artifact_only`
- script lines: `450-451`
- verdict: The lecture script saves an `.RData` workspace; this report saves explicit machine-readable evidence artifacts instead.

**R Lecture Block**

```r
#You can save all your work in a RData object
save.image("./Results/primate_results.RData")
```

**Python `bijux-phylogenetics` Block**

```python
# Workflow/artifact block.
# The R lecture script saves an `.RData` workspace; this evidence pass saves
# explicit JSON artifacts under `evidence-book/studies/primate-longevity-signal/evidence-001/`.
```

**Comparison Evidence**

- R save target: `./Results/primate_results.RData`
- report artifact directory: `evidence-book/studies/primate-longevity-signal/evidence-001`

Raw comparison payload: [`block-payloads/final-workspace-artifact.json`](./results/block-payloads/final-workspace-artifact.json)

