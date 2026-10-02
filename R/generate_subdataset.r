library(KEGGdzPathwaysGEO)
library(Biobase)
library(limma)
library(AnnotationDbi)

# alpha <- 0.05
# logfc_cutoff <- 1
subdataset_sizes <- c(6, 12, 24, 48)
number_of_subdatasets <- 3

set.seed(42)
dir.create("subdataset_gene_output", showWarnings = FALSE)

datasets <- data(package = "KEGGdzPathwaysGEO")$results[, "Item"]

for (dataset in datasets) {

  cat("\nProcessing:", dataset, "\n")

  tryCatch({

    data(list = dataset, package = "KEGGdzPathwaysGEO")
    expressionSet <- get(dataset)
    expressionMatrix <- exprs(expressionSet)

    platform <- annotation(expressionSet)
    pkg <- paste0(platform, ".db")
    library(pkg, character.only = TRUE)
    anno_db <- get(pkg)

    geneSymbols <- mapIds(
      anno_db,
      keys = rownames(expressionMatrix),
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )

    keep <- !is.na(geneSymbols)
    expressionMatrix <- avereps(
      expressionMatrix[keep, ],
      ID = geneSymbols[keep]
    )

    phenotypeData <- pData(expressionSet)
    control_samples <- which(phenotypeData$Group == "c")
    disease_samples <- which(phenotypeData$Group == "d")

    for (subdataset_size in subdataset_sizes) {

      if (subdataset_size > 0.5 * ncol(expressionMatrix)) {
        cat("Skipping size", subdataset_size, "(more than 50% of samples)\n")
        next
      }

      minimum_controls <- max(2, subdataset_size - length(disease_samples))
      maximum_controls <- min(length(control_samples), subdataset_size - 2)

      if (minimum_controls > maximum_controls) {
        cat("Skipping size", subdataset_size, "(not enough samples per group)\n")
        next
      }

      controls_to_sample <- round(
        subdataset_size * length(control_samples) / ncol(expressionMatrix)
      )
      controls_to_sample <- min(
        maximum_controls,
        max(minimum_controls, controls_to_sample)
      )
      diseases_to_sample <- subdataset_size - controls_to_sample
      previous_draws <- character()

      for (subdataset_number in seq_len(number_of_subdatasets)) {

        repeat {
          selected_samples <- sort(c(
            sample(control_samples, controls_to_sample),
            sample(disease_samples, diseases_to_sample)
          ))
          draw <- paste(selected_samples, collapse = ",")
          if (!draw %in% previous_draws) break
        }
        previous_draws <- c(previous_draws, draw)

        subExpressionMatrix <- expressionMatrix[, selected_samples]
        group <- factor(
          phenotypeData$Group[selected_samples],
          levels = c("c", "d")
        )

        designMatrix <- model.matrix(~ group)
        fit <- eBayes(lmFit(subExpressionMatrix, designMatrix))

        deg <- topTable(
          fit,
          coef = 2,
          number = Inf,
          adjust.method = "BH"
        )
        deg <- deg[order(rownames(deg)), , drop = FALSE]

        # sig <- subset(
        #   deg,
        #   adj.P.Val < alpha & abs(logFC) >= logfc_cutoff
        # )

        output_name <- paste0(
          dataset,
          "_size", subdataset_size,
          "_repeat", subdataset_number
        )

        write.csv(
          deg,
          file.path(
            "subdataset_gene_output",
            paste0(output_name, "_DEG_all_genelevel.csv")
          )
        )

        # write.csv(
        #   sig,
        #   file.path(
        #     "subdataset_gene_output",
        #     paste0(output_name, "_DEG_significant_genelevel.csv")
        #   )
        # )

        cat(
          "Size:", subdataset_size,
          "Repeat:", subdataset_number, "\n"
        )
      }
    }

  }, error = function(e) {
    cat("ERROR in", dataset, "-", conditionMessage(e), "\n")
  })
}
