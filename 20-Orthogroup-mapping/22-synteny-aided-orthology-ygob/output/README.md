`20240206-4sps-ygob-cgob-based-orthology-map.tsv`

 2024-02-06, synteny based 4 species orthology map

See `../script/20240203-extract-and-merge-ygob-cgob.R` for detailed steps.

**_Notes on dataset bias_**

Because _S. cerevisiae_ is the common "pillar" between the CGOB abd YGOB datasets, the resulting merged dataset is **"_S. cerevisiae_ centered" or "biased"**, in the sense that only rows containing at least one _S. cerevisiae_ gene are retained. In other words, the table contains many rows with only _S. cerevisiae_ orthologs (could be either gene gains in Sc or losses in all the other three, or failed mapping). In the meantime, genes that are lost in _S. cerevisiae_ will be missing. There are ways we can get around this, but I didn't attempt it here.

**_WGD and ohnologs_**

_S. cerevisiae_ and _C. glabrata_ belong to the post-WGD species, meaning that their ancestor had the whole genome duplicated, with two copies of every gene (this is a simplification, see PMID: 26252497. Most of the duplicated genes from the WGD were subsequently lost. The result is that the post-WGD species generally have ~2x the number of chromosomes but about the same number of genes (~5-6k ORFs) as the pre-WGD species, like _K. lactis_ and _C. albicans_. The YGOB mapping used sequence similarity and synteny to reconstruct the evolutionary histories of all genes. For those genes where the WGD copies were retained, such as _MSN4_ and _MSN2_, and many genes in the glycolysis pathway, there will be two genes mapped to each ancestral gene as well as the extent copy in the pre-WGD species. This is why the table has two separate "Scer" and "Cgla" columns. <p style="color:red">**Important**</p>: here the Scer1 and Scer2 do not correspond to the A and B copies on the YGOB websites. The order is random for our purposes. Also, during the processing steps, I've further altered the order such that they don't even correspond to the original `Pillars.tab` file.

**_How to use_**

The table is inversely sorted by the number of missing homologs (last column termed `n.NAs`). The simplest way to use the table is to retain only those rows that satisfy the following:

    There is one and only one ortholog in each of the four species

This essentially selects for the single gene orthogroup subset. For convenience, I generated this file, `20240206-4sps-ygob-cgob-one-gene-per-sps.tsv`.
