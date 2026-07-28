###############################################################################
# Projeto:
# BreastCancer-LuminalB-Bioinformatics
#
# Aula 01
# Download dos dados TCGA-BRCA utilizando o pacote TCGAbiolinks
#
# Objetivos:
# 1. Conhecer o projeto TCGA
# 2. Aprender a utilizar o TCGAbiolinks
# 3. Baixar dados de RNA-seq do TCGA-BRCA
# 4. Explorar um objeto SummarizedExperiment
#
# Autor:
# Jean Resende
###############################################################################

#==============================================================================
# 1. Instalação dos pacotes (executar apenas uma vez)
#==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")

install.packages("tidyverse")

#==============================================================================
# 2. Carregar os pacotes
#==============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

#==============================================================================
# 3. Criar estrutura do projeto
#==============================================================================

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#==============================================================================
# 4. Construindo a consulta ao GDC
#==============================================================================

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Visualizar o objeto da consulta
query

#==============================================================================
# 5. Download dos dados
#==============================================================================

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 20
)

#==============================================================================
# 6. Preparar os dados
#==============================================================================

brca <- GDCprepare(query)

#==============================================================================
# 7. Explorando o objeto
#==============================================================================

class(brca)

brca

# Quantas linhas e colunas?

dim(brca)

#==============================================================================
# 8. Matriz de expressão
#==============================================================================

counts <- assay(brca)

dim(counts)

head(counts)

#==============================================================================
# 9. Informações dos genes
#==============================================================================

gene_info <- rowData(brca)

gene_info

head(gene_info)

#==============================================================================
# 10. Informações das amostras
#==============================================================================

sample_info <- colData(brca)

sample_info

head(sample_info)

#==============================================================================
# 11. Metadados do objeto
#==============================================================================

metadata(brca)

#==============================================================================
# 12. Nome das amostras
#==============================================================================

samples <- colnames(counts)

head(samples)

#==============================================================================
# 13. Nome dos genes
#==============================================================================

head(rownames(counts))

#==============================================================================
# 14. Removendo a versão do Ensembl
#==============================================================================

rownames(counts) <- sub("\\..*", "", rownames(counts))

head(rownames(counts))

#==============================================================================
# 15. Salvando os objetos
#==============================================================================

saveRDS(brca,
        "data/TCGA_BRCA_SummarizedExperiment.rds")

saveRDS(counts,
        "data/TCGA_BRCA_Counts.rds")

saveRDS(sample_info,
        "data/TCGA_BRCA_SampleInfo.rds")

saveRDS(gene_info,
        "data/TCGA_BRCA_GeneInfo.rds")

###############################################################################
# Fim da Aula 01
###############################################################################