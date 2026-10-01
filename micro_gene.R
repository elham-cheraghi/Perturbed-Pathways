library(KEGGdzPathwaysGEO)
library(Biobase)
library(limma)
library(AnnotationDbi)
library(BiocManager)

alpha <- 0.05
logfc_cutoff <- 1
dir.create("micro_gene_output", showWarnings = FALSE)

# Get all datasets in the package
datasets <- data(package = "KEGGdzPathwaysGEO")$results[, "Item"]

for (dataset in datasets) {

  cat("\n============================\n")
  cat("Processing:", dataset, "\n")

  tryCatch({

    # Load dataset
    data(list = dataset, package = "KEGGdzPathwaysGEO")
    expressionSet <- get(dataset)

    expressionMatrix <- exprs(expressionSet)

    # Detect platform
    platform <- annotation(expressionSet)
    cat("Platform:", platform, "\n")

    # Annotation package name
    pkg <- paste0(platform, ".db")

    # Install if necessary
    # if (!requireNamespace(pkg, quietly = TRUE)) {
    #   BiocManager::install(pkg, ask = FALSE)
    # }

    # Load package
    library(pkg, character.only = TRUE)

    # Get annotation database object
    anno_db <- get(pkg)

    # Probe -> Gene Symbol mapping
    geneSymbols <- mapIds(
      anno_db,
      keys = rownames(expressionMatrix),
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )

    # Remove probes without annotation
    keep <- !is.na(geneSymbols)
    expressionMatrix <- expressionMatrix[keep, ]
    geneSymbols <- geneSymbols[keep]

    # Average probes per gene
    expressionMatrix <- avereps(
      expressionMatrix,
      ID = geneSymbols
    )

    # Phenotype information
    phenotypeData <- pData(expressionSet)
    group <- factor(phenotypeData$Group)

    # Design matrix
    designMatrix <- model.matrix(~ group)

    # Differential expression
    fit <- lmFit(expressionMatrix, designMatrix)
    fit <- eBayes(fit)

    deg <- topTable(
      fit,
      coef = 2,
      number = Inf,
      adjust.method = "BH"
    )
    deg <- deg[order(rownames(deg)), , drop = FALSE]

    # sig <- subset(
    #   deg,
    #   adj.P.Val < alpha &
    #     abs(logFC) >= logfc_cutoff
    # )

    # Save results
    write.csv(
      deg,
      file.path("micro_gene_output", paste0(dataset, "_DEG_all_genelevel.csv"))
    )

    # write.csv(
    #   sig,
    #   file.path("micro_gene_output", paste0(dataset, "_DEG_significant_genelevel.csv"))
    # )

    cat("Genes tested:", nrow(deg), "\n")
    cat("Significant genes:", nrow(sig), "\n")

  }, error = function(e) {

    cat("ERROR in", dataset, "\n")
    cat(conditionMessage(e), "\n")

  })

}
print(datasets)
