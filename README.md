# D_Cottons_USDA
Gossypium D Subgenomes - USDA, MSU, ISU Collaboration

## “Thirteen genome sequences representing the entire subgenus Houzingenia (_Gossypium_): insights into evolution of the New World diploid cottons”

Project Lead: Corrinne E. Grover - Iowa State University

Collaborators:
* Brian Scheffler - United States Department of Agriculture
* Jodi Scheffler - United States Department of Agriculture
* Amanda Hulse-Kemp - United States Department of Agriculture
* Mark A. Arick II - Mississippi State University
* Daniel G. Peterson - Mississippi State University
* Dinum Perera - Mississippi State University
* Jonathan F. Wendel - Iowa State University
* William S. Sanders - Mississippi State University / The Jackson Laboratory


**Summary**: write the summary for the paper here.

**Background literature** (crude summary): The American diploid D-genome cottons (subgenus *Houzingenia*) comprise a monophyletic subgenus which includes 13-14 species (see references in Ulloa 2013) distributed from Southwest Mexico to Arizona, with additional disjunct distributions in Peru and Galapagos. Among these species are *G. raimondii* (D5), the model diploid progenitor to allopolyploid cotton, and *G. harknessii* (D2-2), an important species for cytoplasmic male sterility in cotton (are there others? -- forgetting). Research using various markers have indicated several instances of hybridization, both within the subgenus and, in one remarkable case (*G. gossypioides*), between a member of *Houzingenia* and another, geographically isolated subgenus from Africa. (*G. gossypioides* is also multiply introgressed) 

Descriptions of the subgenus have largely centered around phylogenetics and/or describing the relationships among species. This includes cpDNA restriction sites, nuclear genes (few), ITS, SSR/EST-SSR markers, RAPDS. Many of these agree on certain relationships (e.g., D3d + D3k, D1 + D8, D2-1 + D2-2 + D10, D4 + D7 + D9 + D11), with some deviation in these, and conflict in the branching order among these groups. (Note: there are other studies where these conflict -- and the trees are totally crazy). TEs have been described a bit (Zhao 1998, some in Hawkins 2006), and aside from a crazy microsat/SSR study, not much is known about genetic distance. Some resources exist in Genbank (cpDNA, one mito genome).

Modern sequencing techniques make it easy to produce a substantial amount of genomic sequencing suitable for addressing basic questions in a more thorough manner. Here we use modest coverage Illumina sequencing to present an in-depth view of the subgenus *Houzingenia*. Genome and plastome assemblies are presented and available for use in genome evolutionary research, as well as cotton-specific questions. Comparative molecular evolutionary dynamics of the genic fraction are conducted, revealing the pace and pattern of evolution/substitution in this subgenus. In addition, the intergenic regions are evaluated for the first time to characterize the amount of divergence outside of genes, and due to indels or SNP.

To improve our understanding of the relationships among the D-genome species, as well as provide insight into the introgressed nature of some species, we conduct phylogenetic analyses on (100s, 1000s?) of nuclear genes for which orthology could be strictly identified.  Conflict among genes (indicative of possible hybrid past) is quantified. Traits are characterized and mapped to the phylogeny (TE amount, genome size, NUMTs others) to determine the evolution of these over time.

**Methods** (very rough)

All data was trimmed and filtered with Trimmomatic (v0.32) using the following criteria: remove all bases following a sequence adapter, remove leading and trailing bases with a quality score below 28, remove all bases after average quality over an eight base window falls below 28, remove all bases after a single base quality falls below 10, and filter reads with a length shorter than 85. 

The trimmed data was aligned to the CottonGen.org D5 JGI reference using BWA v0.7.10. For each species, the alignments were merged and processed with BamBam v1.3 to produce consensus sequences.

For each species, the trimmed data was assembled using ABySS v2.0.1, stepping through every 5th k-mer from 40 through 100. The best assemblies (highest N50) for each species was then annotated with MAKER v2.31.6. Evidence for the predicted genes came from the NCBI D5 EST database, the CottonGen.org D5 predicted proteins, and the ab initio gene predictors Genemark v4.30, SNAP v2013-11-29, and Augustus v3.0.3. SNAP and Augustus models will be trained using BUSCO v2.0.

RepeatExplorer was used to discover and compare transposable elements between all species. The trimmed reads cut at 85 bases and sampled at 1% of the genome size for each species. A cotton-specific repeat database was used to classify the clusters.

## Tasks

1. Genome assemblies
  1. reassembly underway (by Tony, 12/06/2016)
  1. gather assembly stats
  2. assembly methods summarized (see above)
  4. run MAKER on assemblies for *de novo* gene predictions/content evaluations
2. Chloroplast assemblies (on hold, not a priority)
  1. gather assembly stats
  2. summarize assembly methods
  3. annotate cpDNA (DOGMA, CpGAVAS, or ?)
3. Intergenic space alignments/synteny
  1. map contigs to D5 genome
  2. extract regions where all (>10? >8?) species are represented
  2. for aligned intergenic regions 
    1. determine nt diversity/SNP differences (A genome outgroup)
    2. characterize number/length of indels (A genome outgroup)
    3. summarize for all species 
3. Gene content
  1. summarize MAKER output
  2. create gene families via OrthoFinder
  3. summarize gene family differences
  4. GO annotation of gene families (?)
4. Phylogeny
  1. align reads against D5 reference (Tony has this underway, 12/06/2016)
  2. bam2consensus to get genes for phylogeny
  3. filter alignments (https://github.com/Wendellab/phylogenetics/blob/master/process_alignments or similar)
  4. Phylogenetics
    1. Ka/Ks
    1. Raxml (+ Astral coalescence?) for ML
    2. MrBayes (+ BUCKy concordance) for Bayesian
5. TE characterization
  1. Cluster via RepeatExplorer
  2. Annotate with newest cotton-enriched RepeatMasker database (created 09/2016)
  2. Hierarchical clustering of repeats to assess similarity among species
  3. Map character to tree and assess gain/loss (fastANC or similar) overall, and on a per character basis (however many subdivisions)
6. NUMT/NUPT characterization
  1. Use RepeatExplorer clusters
  2. Annotate against existing organellar sequences for cotton
  3. Determine extent of occupation, differences
  3. Map against phylogeny (if interesting)
  
