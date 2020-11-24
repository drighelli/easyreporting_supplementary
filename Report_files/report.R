# BiocManager::install("easyreporting")

library(easyreporting)
proj.path <- file.path(getwd(), "rnaseq_report_live")
bioEr <- easyreporting(filenamePath=proj.path, title="RNA-seq Analysis Report",
                       author=c("Dario Righelli"))

mkdTitle(bioEr, title="Loading Counts Data")

## first way CC creation
mkdCodeChunkSt(bioEr, sourceFilesList=system.file("script/importFunctions.R", package="easyreporting"), isComplete=TRUE)
mkdVariableAssignment(bioEr, "geneCounts", "as.matrix(importData(system.file('extdata/BMDC_counts_FeatureCounts.xlsx', package='easyreporting')))", show=FALSE)
mkdGeneralMsg(bioEr, "head(geneCounts, 20)")
mkdCodeChunkEnd(bioEr)

## one call CC creation
mkdCodeChunkComplete(object=bioEr, 
    message=paste("geneCounts <- as.matrix(importData(system.file('extdata/BMDC_counts_FeatureCounts.xlsx', package='easyreporting')))",
    "head(geneCounts, 20)", sep="\n"),
    sourceFilesList=system.file("script/importFunctions.R", 
                    package="easyreporting"), 
    optionList=makeOptionsList(evalFlag=FALSE))

mkdTitle(bioEr, title="Plot Boxplot on count data", level=2)
mkdCodeChunkComplete(bioEr, message=paste0("boxplot(log(geneCounts+1),",
                                           " col=c('red','red','orange','orange','purple','purple'),", 
                                           " main='Counts BoxPlot',las=1)"))

mkdTitle(bioEr, title="Filtering Low Abundant Features", level=1)
mkdCodeChunkComplete(object=bioEr, message=paste("fgeneCounts <- ",
              "NOISeq::filtered.data(dataset=geneCounts, ",
              "factor=c('D', 'D', 'E', 'E', 'C', 'C'), norm=FALSE, method=3, cv.cutoff=100, cpm=0.5)",
              "boxplot(log(fgeneCounts+1),col=c('red','red','orange','orange','purple','purple'), ",
              "main='Counts BoxPlot',las=1)", sep="\n"))


mkdTitle(bioEr, title="Normalizing Features Across Samples", level=1)
mkdCodeChunkComplete(object=bioEr, message="nfgeneCounts <- EDASeq::betweenLaneNormalization(fgeneCounts, which='upper')")


mkdTitle(bioEr, title="Plot PCA on count data", level=2)
mkdCodeChunkComplete(bioEr, message=paste("se <- SummarizedExperiment::SummarizedExperiment((log2(nfgeneCounts)+ 1), ",
                "colData=S4Vectors::DataFrame(rownames=colnames(nfgeneCounts), " ,
                "condition=c('DEC', 'DEC', 'E2', 'E2', 'CTRL', 'CTRL')))",
                "DESeq2::plotPCA(DESeq2::DESeqTransform(se))", sep="\n"))


## one call with title and CC, can also add comment
mkdCodeChunkTitledCommented(bioEr, title="Differential Expression Analysis",
                codeMsg=paste0("degList <- applyEdgeREx(counts=nfgeneCounts, ",
                "factors=c('DEC', 'DEC', 'E2', 'E2', 'UNTR', 'UNTR'), ",
                "contrasts=c('DEC - UNTR', 'E2 - UNTR'),",
                "p.threshold=1)"),
            commentMsg=paste0("As we saw from the PCA, the groups are well separated, ",
            "so we can perform a Differential Expression analysis with edgeR."),
            sourceFilesList=system.file("script/geneFunctions.R", package="easyreporting"))

mkdCodeChunkTitledCommented(bioEr, title="MA-Plot", level=2,
                                    codeMsg="MAedgeRMAPlotEx(degList=degList)",
                                    sourceFilesList=system.file("script/plotFunctions.R",
                                                        package="easyreporting"))

mkdCodeChunkTitledCommented(bioEr, title="Volcano-Plot DEC-UNTR", level=2,
                            codeMsg="VolcanoPlot(degs=degList[[1]])",
                            sourceFilesList=system.file("script/plotFunctions.R",
                                                        package="easyreporting"))


mkdCodeChunkTitledCommented(bioEr, title="Volcano-Plot E2-UNTR", level=2,
                            codeMsg="VolcanoPlot(degs=degList[[2]])",
                            sourceFilesList=system.file("script/plotFunctions.R",
                                                        package="easyreporting"))

compile(bioEr)



