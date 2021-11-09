# TCGA_PDAC
---
title: " Differntial gene expression of Pancreatic adenocarcnoma with limma package"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
```
#Get list of available projects of GDC

```{r}
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-PAAD")
query_TCGA = GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")
```

#To visualize the query results in a more readable way, we can use the command getResults.
```{r}
pc_res = getResults(query_TCGA) # make results as table
# head([pc_res) # data of the first 6 patients.
colnames(pc_res) # columns present in the table

head(pc_res$sample_type) # first 6 types of tissue.

summary(pc_res$sample_type) # summary of distinct tissues types present in this study
```
#Preprocessing step:
```{r}
query_TCGA = GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")) # Download solid and tumor samples

GDCdownload(query = query_TCGA) # Download the files in the query

tcga_data = GDCprepare(query_TCGA) #Finally, load the actual RNASeq data into R. 
```
```{r}
dim(tcga_data) #check the file size
```

#There are 3 functions that allow us to access to most important data present in this object, these are: colData(), rowData(), assays(). Use the command ?SummarizedExperiment to find more details. colData() allows us to access the clinical data associated with our samples. The functions colnames() and rownames() can be used to extract the column and rows names from a given table respectively.
```{r}
# In R (and other programming languages) chain functions allows us to save time and space.
colnames(colData(tcga_data))
```
#The table() function (in this context) produces a summary with the sum of each of the factors present in a given column.
```{r}
table(tcga_data@colData$vital_status)
```

```{r}
table(tcga_data@colData$tumor_stage)
```
```{r}
table(tcga_data@colData$definition)
```
```{r}
table(tcga_data@colData$tissue_or_organ_of_origin)
```
``{r}
table(tcga_data@colData$gender)
```
```{r}
table(tcga_data@colData$race)
```
#To distinguish normal vs tumor samples of PC datasets.
```{r}
dim(assay(tcga_data))     # gene expression matrices.
```
#We use assay and row functions
```{r}
head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples
```
```{r}
head(rowData(tcga_data))    # ensembl id and gene id of the first 6 genes.
```

#Finally, we can use some core functionality of R to save the TCGA_data as a .RDS file. This is faster than repeating the previous operations (GDC download) and useful to work in same datasets for several days.
```{r}
# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)
```

#The data is loaded in R with the following command
```{r}
tcga_data = readRDS(file = "tcga_data.RDS")
```

#Limma differential expression analysis in R
# RNA sequence normalization
#A typical task on RNA-Seq data is differential expression (DE) analysis, based on some clinical phenotypes. This, in turn, requires normalization of the data, as in its raw format it may have batch effects and other artifacts.

#A common approach to such complex tasks is to define a computational pipeline, performing several steps in sequence, allowing the user to select different parameters.
We will now define and run one such pipeline, through the use of an R function.
The function is called limma_pipeline(tcga_data, condition_variable, reference_group), where tcga_data is the data object we have gotten from TCGA and condition_variable is the interesting variable/condition by which you want to group your patient samples. You can also decide which one of the values of your conditional variable is going to be the reference group, with the reference_group parameter.This function returns a list with three different objects:

#A complex object, resulting from running voom, this contains the TMM+voom normalized data;
#A complex object, resulting from running eBayes, this contains the the fitted model plus a number of statistics related to each of the probes;
#A simple table, resulting from running topTable, which contains the top 100 differentially expressed genes sorted by p-value. This is how the code of this function:
```{r}
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){

  design_factor = colData(tcga_data)[, condition_variable, drop=T]

  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}

  design = model.matrix(~ group)

  dge = DGEList(counts=assay(tcga_data),
                 samples=colData(tcga_data),
                 genes=as.data.frame(rowData(tcga_data)))

  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)

  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)

  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)

  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")

  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}
```

#With the following command, we can obtain the DE analysis comparing Primary solid Tumor samples against Solid Tissue Normal. This will be used in the next section, on the classification task.

```{r}
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)
```
setwd("/github/home/mariam.zamani")
#Let’s save this object to file, like we did with tcga_data:
```{r}
# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)
```

#As an example, we perform limma_pipeline to perform DE analysis by grouping patients by gender instead of by tissue type.
```{r}
gender_limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="gender",
  reference_group="female"
)
```

#Visualization
#PCA Plot given the voom object created by the limma_pipeline function. 
```{r}
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
```

#By calling the function plot_PCA, we get a plot of the first two principal components: Condition variable = gender
```{r}
res_pca = plot_PCA(limma_res$voomObj, "definition")
```
#Condition variable = Solid Tissue Normal
```{r}
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
res_pca = plot_PCA(limma_res$voomObj, "definition")
```
##Limma differential expression analysis in R

#In our pipeline function, we use the package limma. We will select a particular clinical feature of the data to use as class for grouping the samples as either tumor vs normal tissue. This data is available under the column definition for tcga_data, but needs the use of the function colData to be accessed. In addition, limma requires this data to be a factor, so we convert it as such:
```{r}
clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)
```

#As seen before, we have 2 distinct groups of tissues defined in this column, Solid Tissue Normal (our control samples) and Primary Solid Tumor (our samples of interest). In this factor, we also want to define Solid Tissue Normal as being our base or reference level.
```{r}
group = relevel(group, ref="Solid Tissue Normal")
```

#Next, we need to create a design matrix, which will indicate the conditions to be compared by the DE analysis. The ~ symbol represents that we are constructing a formula.
```{r}
design = model.matrix(~group)
head(design)
```
#Before performing DE analysis, we remove genes, which have low amount of counts. We transform our tcga_data object as DGEList, which provides functions for filtering. By default genes with counts with less than 10 reads are removed.
```{r}
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore
```
#Before we fit a model to our data, we normalize the data to minimize batch effects and technical variation with the TMM (trimmed mean of M-values) normalization method. Moreover, to apply limma on RNA-seq, we need to convert the data to have a similar variance as arrays. This is done with the VOOM method.
```{r}
dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)
```
#Finally, using lmFit lets fit a series of linear models, one to each of the probes. These data will then be fed to eBayes to produce a complex object which holds a number of statistics that we can use to rank the differentially expressed genes.

```{r}
fit = lmFit(v, design)
fit = eBayes(fit)
```
#Using the function of topTable we can check the top10 genes classified as being differentially expressed. 

```{r}
topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)
```
#Now, we will explore the use of a machine learning method to classify an unseen sample as being a tumor or not. To achieve this goal we are going to build first a simple linear model (with feature selection), and then an Elastic Net model. For this, we need to split the data into two sets: a training set, which we will use to train our model, and a testing set. The test data serves as an independent dataset to validate our results. It is important that you do not use test data to optimize your results or this will include bias in the classifier.

#First extract the data that we are going to use to build our model. We want the expression data that has already been normalized and a clinical feature which divides our data into different groups, such as tumor vs. non-tumor or tumor stage. We can get the normalized expression values from limma_res$voomObj$E and the type of sample is determined by the definition column.
```{r}
# Transpose and make it into a matrix object
d_mat = as.matrix(t(limma_res$voomObj$E))

# As before, we want this to be a factor
d_resp = as.factor(limma_res$voomObj$targets$definition)
```
#With the data in the correct format we can now divide it into a train set and a test set. We will use the function createDataPartition which creates a vector of booleans (TRUE or FALSE) that we can then use to subset the matrix in this case leaving 75% of samples for training and 25% for testing.
```{r}
# Divide data into training and testing set

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(42)
train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]
```

#x_train and y_train are the data we will use to train our model, where x is the matrix with the normalized expression values and y is the vector with the response variable, Primary solid Tumor and Solid Tissue Normal.

#Following the same logic, x_test and y_test are the matrix with normalized expression values and the response variable respectively. Again, we will only use these (*_test) to perform a prediction and evaluate how good this prediction was after the training process has finished.


#Train Elastic Net model
Definition of Elastic model:
#We will train an Elastic Net model, which is a generalized linear model that combines the best of two other models: LASSO and Ridge Regression. Ridge Regression is often good at doing predictions but its results are not very interpretable. LASSO is good at picking up small signals from lots of noise but it tends to minimize redundancy so if there are two genes that are equally good predictors (features with high correlation) it will tend to select one. Elastic Net is a balance between both methods, it selects the genes or groups of genes (if they are correlated) that best predict each of the conditions and use these to build a model that will then be used for classification.

#We can then look at these genes individually to see if some interesting gene of biological relevance for the classification problem is selected. When using Elastic Net there are other parameters than we have to set, specifically alpha. This parameter will define if the Elastic Net will behave more like LASSO (alpha = 1) or like Ridge Regression (alpha = 0). For simplicity we will set it to 0.5 however in a real setting we would probably vary this value in order to find the best model (minimizing the error).
```{r}
# Train model on training dataset using cross-validation
res = cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "binomial" 
)
```
Definition of Confusion matrix:
#A confusion matrix is a simple table that compares the predictions from our model against their real values. Therefore, it shows us the true positives, true negatives, false positives and false negatives. We can use it to compute a number of accuracy metrics that we can then use to understand how good our model actually is.
```{r}
# Test/Make prediction on test dataset
y_pred = predict(res, newx=x_test, type="class", s="lambda.min")
View(y_pred)
```

```{r}
confusion_matrix = table(y_pred, y_test)
# Evaluation statistics
print(confusion_matrix)
```
```{r}
print(paste0("Sensitivity:", sensitivity(confusion_matrix))
print(paste0("Specificity: ",specificity(confusion_matrix)))
```
```{r}
# Getting genes that contribute for the prediction
res_coef = coef(res, s="lambda.min") # the "coef" function returns a sparse matrix
dim(res_coef)
```
#We can now take a look at the genes (coefficients), that Elastic Net selected to build its model.
```{r}
head(res_coef) # in a sparse matrix the "." represents the value of zero
```
#The number of coefficients is large (as there are many genes!). We only want to consider coefficients with non-zero values, as these represent variables (genes) selected by the Elastic Net.

```{r}
# get coefficients with non-zero values
res_coef = res_coef[res_coef[,1] != 0,]
# note how performing this operation changed the type of the variable
head(res_coef)
```
```{r}
# remove first coefficient as this is the intercept, a variable of the model itself
res_coef = res_coef[-1]

relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes
```

#Identify relevant genes
```{r}
head(relevant_genes) # few select genes
```
#Get the common gene name from limma_res$voomObj$genes as in this table we can find the ensembl to gene name correspondence.
```{r}
head(limma_res$voomObj$genes)
```
#Relevant_gene_names 
```{r}
relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"external_gene_name"]

head(relevant_gene_names) # few select genes (with readable names now)
```
#Question: Are there any genes of particular relevance?
#Hint's look up!
#Did limma and Elastic Net select some of the same genes? We can check the common genes between the two results by using the intersect function.
```{r}
print(intersect(limma_res$topGenes$ensembl_gene_id, relevant_genes)) # Shows relevant genes
```
#Hierarchical clustering
#Finally we can take a look at how our samples cluster together by running an hierarchical clustering algorithm. We will only be looking at the genes Elastic Net found and use these to cluster the samples. The genes highlighted in green are the ones that limma had also selected as we’ve seen before. The samples highlighted in red are Solid Tissue Normal, the samples highlighted in black are Primary solid Tumor.

```{r}
# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)

# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_genes %in% limma_res$topGenes$ensembl_gene_id),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)

# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)
```
#As you’ve seen, selected genes group into two classes: genes highly expressed in tumors vs. genes in control. Interestingly, genes also detected by DE analysis are only associated to high expression in the control group. One interesting question is, are selected genes (up in control or up in tumor) associated to any type of common biological problem? I Tred to do a GO analysis on them.

```{r}

# Extract the hierarchical cluster from heatmap to class "hclust"
hc = as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters = cutree(hc, k=2)
table(clusters)
```
```{r}
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "PITX1","CD5L","PLXNA2","FXYD3","MISP","BMF"), 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p
```
```{r}
# Group 2, up in control
library(gprofiler2)
get_base_url()
get_user_agent()
## Not run: 
gorth(c("PITX1","CD5L","PLXNA2","FXYD3","MISP","BMF","Nanog"), source_organism="mmusculus", target_organism="hsapiens")
query <- gprofiler(c("PITX1","CD5L","PLXNA2","FXYD3","MISP","BMF"),)
gprofiler(query, organism = "hsapiens", sort_by_structure = T,
          ordered_query = F, significant = T, exclude_iea = F,
          underrep = F, evcodes = F, region_query = F, max_p_value = 1,
          min_set_size = 0, max_set_size = 0, min_isect_size = 0,
          correction_method = "analytical", hier_filtering = "none",
          domain_size = "annotated", custom_bg = "", numeric_ns = "",
          png_fn = NULL, include_graph = F, src_filter = NULL)
## End(Not run)


```{r}
# retain only a small subset of the genes (see documentation for ?varFilter)
d_mat = varFilter(limma_res$voomObj$E, var.func=IQR, var.cutoff=0.95)

# transpose the matrix, so that it has the same shape as the d_mat we used at the beginning
d_mat = t(d_mat)

print(dim(d_mat))
```
Part 2:
##NB: set the working directory of RStudio correctly
```{r}
# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
setwd("/gpfs1/home/mariam.zamani")
tcga_data = readRDS(file = "tcga_data.RDS")
limma_res = readRDS(file = "limma_res.RDS")
```
#Survival analysis using gender variables

```{r}

# extract clinical data
clinical = tcga_data@colData

dim(clinical)
# we are only interested in the "Primary solid Tumor" cases for survival
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                    c("patient",
                      "vital_status",
                      "days_to_death",
                      "days_to_last_follow_up",
                      "gender",
                      "tumor_stage")]

```

#Now we have a new dataframe, clin_df, containing only the information that is relevant to survival analysis. In addition to gender, we have added vital_status (whether patient is alive or dead), tumor_stage (from stage 1 to 4) and two important variables: days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients), and days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

#We need to format the dataframe in a way that is readable to the methods from the survival package :
```{r}
# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                   clin_df$days_to_death,
                                   clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)
Surv(clin_df$overall_survival, clin_df$deceased)
```
#Kaplan-Meier plots
As a first step, we need to define a survival formula with the help of the Surv object.

#In R, formulas are special constructs of the form y ~ x, and in the context of linear models you can see x as the independent variable and y as the dependent variable.
This works also for multivariate models: age ~ gender + height is a formula that can be used to predict age from gender and height. You can refer to the documentation of formula for further examples and explanations, by typing ?formula in a R shell. In survival. We have a categorical variable, gender, that needs to be used to separate (or, more appropriately, stratify) the available death events. The survival package provides us with an object, Surv, to form a dependent variable out of the overall_survival and deceased information:


```{r}
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender
# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)

print(fit)
```
#This modifies our overall survival vector by adding censoring information (the + just after the time), which requires a small digression.

This data is right censored, meaning that for some patients we only have the time of the last follow up but we don’t know if they died at a later date or not.
These patients are kept in the early stages of the analysis (eg, they are part of the survival curve) but they are dropped (or as it is said, censored) when the time of their last follow up arrives. Now that the survival time has been tagged with the censoring, we can add the categorical independent variable gender, and effectively create a formula

```{r}
# we produce a Kaplan Meier plot
ggsurvplot(fit, data=clin_df)
```
```{r}
ggsurvplot(fit, data=clin_df, pval=T)
```
#The p-value is non-significant, so gender alone does not significantly sway prognosis in this dataset.

#Question = Can we see the number of patients dying (or being “censored”) as Time increases? Indeed we can, with what is called the “at risk table”.

```{r}
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, risk.table.col="strata")
```
#With the risk.table=T argument, we get the number of patients “at risk”, that is neither dead nor censored at a certain time point.

#The argument risk.table.col="strata" tells survminer to colour the table in the same way as the strata, or groups, are coloured.

#Question = how does tumor stage affect survival?

#The tumor_stage variable that TCGA provides for this tumor contains both stages and sub-stages, eg stage iiia or stage ivb. We want to join together the sub-stages, to increase the group size and reduce complexity (and thus increase the power of the logrank statistics).

```{r}
# remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[abc]$", "", clin_df$tumor_stage)

# we remove those with stage "not reported", since they are unknown
clin_df[which(clin_df$tumor_stage == "not reported"), "tumor_stage"] = NA

# finally, we also remove those with tumor stage 4, since they are too few
clin_df[which(clin_df$tumor_stage == "stage iv"), "tumor_stage"] = NA

table(clin_df$tumor_stage)
```
#We can now fit a new survival model with the tumor stage groups (one to four, plus the “not reported”):
```{r}
fit = survfit(Surv(overall_survival, deceased) ~ tumor_stage, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)
```

#Generate a Kaplan-Meier plot from the fitted model
```{r}
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T) # # we produce a Kaplan-Meier plot from the fitted model
```
#We get an overall p-value testing the null hypothesis that all the curves are similar at every time point. In this case, the p-value is small enough that we can reject the null hypothesis.

#What we saw here is an easy way of producing Kaplan-Meier plots to investigate survival, as well as evaluating whether the survival curves are significantly different or not. A more interesting application to this is using, for example, gene expression to divide the patients into groups, to see whether up or down regulation of genes affects survival. We’ll see this in the next section.

# Gene expression and survival:
#Ealier, we looked at the RNASeq data for this tumor. We found some genes that are differentially expressed between the healthy and disease samples, and we also trained an Elastic Net model and investigated those predictors that are important in discriminating healthy and disease tissue. So, taken one of the selected genes, we know that they are either up-regulated in the tumor tissue tissue and not in the normal tissue, or viceversa. But do these genes also provide an advantage, or disadvantage, regarding prognosis?
We already have the top differentially expressed genes, ordered by significance, in the limma_res$topGenes dataframe, so we just have to take the first one.

```{r}
# let's extract the table of differential expression we got earlier
expr_df = limma_res$topGenes

# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[1, ])
```
#This work was part of unpublished research of analysis of publicly available "TCGA pancreatic adenocarcinoma gene expression dataset" from Jansen lab of NDSU.

#Acknowledgemnt:https::costalab.ukaachen.de/
