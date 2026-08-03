# Methods — Interactive web application (CoPA)

To make the consensus co-expression modules and their analyses openly
explorable, we developed CoPA, an interactive web application built in R with
the Shiny framework (bslib UI, plotly interactive graphics) and deployed on
shinyapps.io. The application comprises four modules. The CoPA tab browses each
co-expression module, displaying its top genes by module membership (kME), the
coherence of co-expression across seven individual bulk prefrontal cortex
RNA-seq datasets, gene set enrichment analysis, and module projections (mean
expression across cortical cell types, as log UMI counts or REI) in four adult
human dorsolateral prefrontal cortex snRNA-seq datasets. The dCoPA tab
visualizes disease-associated module changes (AD and schizophrenia versus
control) resolved by cell type. The Core GBmap tab projects a user-supplied
module onto a reference atlas, and the Gene Projection tab reports per-gene
cell-type expression across reference datasets. All projections use pre-computed
cell-type means and variances for responsiveness.
