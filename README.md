# StatBio_Transcriptomics
Some functions I DIYed for Transcriptomic Analysis

There is a function to calculate a Benjamini-Hochberg corrected P-value which is used by the other function in the file
which is basically a simple Ontology enrichment tester. Given a certain table containing Gene ID's and their corresponding
ontology terms in two columns, you can pass the table to this function with a gene list and a background list. It will
test the proportion of ontology terms in the gene list against the proportion of the same terms in the background list.
If the proportions are significantly different, as determined by giving the option of the choice of test, Fisher Exact,
Binomial or Hypergeometric then we know that the gene list is enriched in genes from a particular ontology category.
The Binomial and Fisher Exact tests are set to alternative="greater" by default so they will test for enrichment only. 

Visual inspection of stacked barplots should be done carefully as the axis is not linear and therefore the lengths of the stacked barplots are not geometrically proportional to the relative numbers of cateogry genes in the query list and the background. Refer to the scale on the x-axis to get an idea of the number of genes in the list. The logarithmic scale provides a visually appealing comparison between ontology categories that may vary orders of magnitude in the number of genes that they contain. 

The barplots are generated for each ontology category that clears the statistical test and multiple correction P-value correction under user-specified FDR. 
