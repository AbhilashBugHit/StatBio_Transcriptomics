# StatBio_Transcriptomics
Some functions I DIYed for Transcriptomic Analysis

There is a function to calculate a Benjamini-Hochberg corrected P-value which is used by the other function in the file
which is basically a simple Ontology enrichment tester. Given a certain table containing Gene ID's and their corresponding
ontology terms in two columns, you can pass the table to this function with a gene list and a background list. It will
test the proportion of ontology terms in the gene list against the proportion of the same terms in the background list.
If the proportions are significantly different, as determined by giving the option of the choice of test, Fisher Exact,
Binomial or Hypergeometric then we know that the gene list is enriched or depleted in genes from a particular ontology category.
The Binomial and Fisher Exact tests are set to two-sided by default so the determination of enrichment or depletion is done 
by visual inspection of stacked barplots that are generated for each ontology category that clears the statistical test and
multiple correction. Barplots where the gene list (red) have more genes than the background for that same cateogry are
considered to be enriched and if they are lesser in proportion to the green bar are considered depleted.
