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

*Sequence generation and initial processing*. DNA was extracted from (leaves) using (what kit), and sent to (where) for library construction and sequencing.  Sequencing was completed on the Illumina (what machine) using (what type of sequencing). The data were trimmed and filtered with Trimmomatic v0.32 (citation) with the following options: (1) sequence adapter removal, (2) removal of leading and/or trailing bases when the quality score (Q) <28, (3) removal of bases after average Q <28 (8 nt window) or single base quality <10, and (4) removal of reads <85 nt. 

*Genome assembly and annotation.* The trimmed data was independently assembled for each species (combining accessions of the same species) via ABySS v2.0.1 (citation), using every 5th kmer value from 40 through 100. A single assembly with the highest E-size (an alternative statistic to N50; cite Salzberg 2011) was selected for each species and subsequently annotated with MAKER v2.31.6 (citation) using evidence from: (1) the NCBI G. raimondii EST database (citation), (2) *G. raimondii* reference genome predicted proteins, as hosted by CottonGen.org (citation, citation), and (3) three ab initio gene prediction programs, i.e. Genemark v4.30 (citation), SNAP v2013-11-29 (citation), and Augustus v3.0.3 (citation). Both the SNAP and Augustus models were trained using BUSCO v2.0 (citation). Gene orthology and family designations were determined via OrthoFinder (citation), and subsequently annotated via blast2go (citation).  

*Phylogenetics and molecular evolutionary analyses.* The trimmed reads were mapped against the *G. raimondii* reference sequence (cite Paterson 2012) using BWA v0.7.10 (citation), post-processed with samtools (which version) (citation), and individual genes were independently assembled for each species/accession via BamBam v 1.3 (citation) in conjunction with the G. raimondii reference annotation (cite Paterson 2012).  Alignments were pruned for genes and/or gene regions with insufficient coverage, i.e., too many ambiguous bases, using filter_alignments (https://github.com/Wendellab/diploidDphylogeny). Parameters were set to remove sequences with more than 30% ambiguous bases within species and to remove aligned positions with more than 20% ambiguity among species. Only those genes with (what species configuration) were retained for phylogenetic and molecular analyses. Those genes passing these filtering steps were also concatenated for phylogenetic analysis.

Phylogenetic relationships were inferred for individual genes and the concatenated alignment using both Maximum Likelihood and Bayesian analyses. Individual models were determined by jModelTest2 (citation) for each gene and binned according to their parameters. Maximum likelihood analyses were performed using RaxML v(which version) (cite Stamatakis 2014) using (the basic general time reversible model with gamma distribution (GTRGAMMA)), 1000 alternative runs on distinct starting trees, and rapid bootstrapping with consensus tree generation. Bayesian analyses were generated using MrBayes v (which version) (cite Ronquist and Huelsenbeck 2003) under (GTR gamma) with the following parameters: 4 runs with 4 chains for 10 million generations and using a burn-in fraction of 25%. Resulting trees were rooted with a member of subgenus Gossypium, G. arboreum. Concordance among individual gene trees was assessed via BUCKy (Ané et al., 2007; Larget et al., 2010) with 3 runs, each with 4 chains and 1 million iterations, and default parameters. Since BUCKy requires common taxa for all gene trees, a basic set of taxa was selected using a single representative accession for each species. Coalescence was assessed via ASTRALv (which version) (cite Mirarab 2014) under (which parameters).

Measures of molecular evolution were all calculated in R v(which version)(citation). Species divergence estimates were calculated via the {chronoMPL} package (citation), using the (which time estimates?) (citation). Trees were visualized using the {ape} package (citation). Synonymous and non-synonymous divergence estimates were calculated via the {seqinr} package (citation) and visualized using {ggplot2} (citation). Distance measures of aligned intergenic regions were estimated via {ape}, and indels were characterized by (???).

*Transposable element characterization and ancestral state reconstruction*. Transposable elements (TEs) and other repetitive DNAs were assessed via the RepeatExplorer pipeline (citation). The trimmed reads were uniformly cut to 85 nt and sampled to 1% genome size equivalents for each species using the estimates of Hendrix and Stewart (2005). The data were clustered according to standard parameters, substituting a cotton-repeat enriched database in the RepeatMasker (citation) step. This database is derived from all plant entries in RepBase v(which) (citation) combined with those cotton repeats previously identified (cite Hawkins, Grover, Paterson, etc). Custom (what language) scripts, available at some/github/repo, were used to summarize the data.

Total amounts of each category and subcategory of TEs were summed in R and used as input in to {fastAnc} (citation) for ancestral state reconstruction of each feature. The fitContinuous function of {Geiger} (citation) was also used to test different models of ancestral state reconstruction for each.


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
  
