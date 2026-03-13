# Making flow chart
install.packages("igraph", repos="https://cloud.r-project.org", type="binary")
install.packages("DiagrammeR", repos="https://cloud.r-project.org", type="binary")
library(DiagrammeR)
grViz("
digraph lncRNA_pipeline {

  graph [layout = dot, rankdir = TB]

  node [shape = rectangle, style = filled, 
        fillcolor = LightBlue]

  A [label = 'All Transcripts:\n1,264,505 Transcripts']
  B [label = 'CPAT Filter:\n274,333 Transcripts']
  C [label = 'GffCompare Filter:\n154,635 Transcripts']
  D [label = 'Clustering Filter:\n73,195 Transcripts']
  E [label = 'Pfam Filter:\n44,281 Transcripts']
  F [label = 'Rfam Filter:\n43,998 Transcripts']
  G [label = 'Strain-specific lncRNAs:\n500 Transcripts']
  H [label = 'Haplotype-specific lncRNAs:\n112 Transcripts']

  A -> B -> C -> D -> E -> F -> G -> H 
}
")

# 
