# Identifying changes in RNA binding
A short workshop to demonstrate how to use linear modeling to identify changes in RNA binding from parallel `Total` and `RNA-bound` protein quantification.

Here, we assume the reader is aware of the basic principles of protein quantification using isobaric tags (in this case, TMT labeling). We further assume the reader is aware of the OOPS method used to obtain RNA-bound protein (). Finally, the reader will benefit from prior understanding of the `MsnSet` class. 

The workshop takes the form of R markdown notebooks (See [here](notebooks/)) with worked examples and questions and tasks to prompt the reader to consider important aspects of the analysis. When working through the workshop by oneself, the answers to the questions and tasks can be found [here](workshop_answers.txt).

The notebooks cover the following:

[1_simple_example.Rmd](notebooks/1_simple_example.Rmd)
Working up from an example of modeling the abundance of a single protein in a control vs treatment experiment to modeling the change in RNA binding for a single protein in a control vs treatment experiment. 

[2_identify_changes_in_RNA_binding.Rmd](notebooks/2_identify_changes_in_RNA_binding.Rmd)
Applying a linear model (using `lm`) to each protein in a real data set to identify those with a signficant change in RNA binding

[3_using_limma.Rmd](notebooks/3_using_limma.Rmd)
Applying a moderated linear model using `limma` to shrink the observed biological variance toward the expected values from proteins with a similar abundance. Comparing the results obtained with `lm` and `limma`.

[4_get_all_go_terms_h_sapiens.Rmd](notebooks/4_get_all_go_terms_h_sapiens.Rmd)
Expanding the set of annotated GO terms to the complete set of implictly annotated GO terms

[5_Identify_over_rep_GO_terms.Rmd](notebooks/5_Identify_over_rep_GO_terms.Rmd)
Using the above expanded GO term annotations and `goseq` to identify the GO terms over-represented in the proteins with a significant change in RNA binding. Using the `topTreat` function from `limma` to test an adjusted null hypothesis to additionally take into account the mimimum effect size which is deemed biologically relevant.
