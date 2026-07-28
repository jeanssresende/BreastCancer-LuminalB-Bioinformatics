###############################################################################
# Projeto:
# BreastCancer-LuminalB-Bioinformatics
#
# Aula 01B
# Download de matrizes de expressão do TCGA-BRCA via UCSC Xena
#
# Objetivos:
# 1. Baixar matrizes de expressão (Counts e TPM)
# 2. Baixar dados clínicos
# 3. Baixar informações de sobrevida
# 4. Baixar classificação molecular PAM50
#
# Fonte dos dados:
# https://xenabrowser.net/
#
# Autor:
# Jean Resende
###############################################################################

#==============================================================================
# 1. Criar diretórios
#==============================================================================

dir.create("data", showWarnings = FALSE)
dir.create("data/raw", showWarnings = FALSE)

#==============================================================================
# 2. URLs dos arquivos
#==============================================================================

counts_url <-
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_counts.tsv.gz"

tpm_url <-
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_tpm.tsv.gz"

clinical_url <-
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.clinical.tsv.gz"

survival_url <-
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.survival.tsv.gz"

pam50_url <-
  "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/brca%2Fpam50"

#==============================================================================
# 3. Download dos arquivos
#==============================================================================
options(timeout = 600)

download.file(
  url = counts_url,
  destfile = "data/raw/TCGA_BRCA_counts.tsv.gz",
  mode = "wb"
)

download.file(
  url = tpm_url,
  destfile = "data/raw/TCGA_BRCA_tpm.tsv.gz",
  mode = "wb"
)

download.file(
  url = clinical_url,
  destfile = "data/raw/TCGA_BRCA_clinical.tsv.gz",
  mode = "wb"
)

download.file(
  url = survival_url,
  destfile = "data/raw/TCGA_BRCA_survival.tsv.gz",
  mode = "wb"
)

download.file(
  url = pam50_url,
  destfile = "data/raw/TCGA_BRCA_PAM50.tsv",
  mode = "wb"
)

#==============================================================================
# 4. Importar os arquivos
#==============================================================================

counts <- read.delim(
  gzfile("data/raw/TCGA_BRCA_counts.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

tpm <- read.delim(
  gzfile("data/raw/TCGA_BRCA_tpm.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

clinical <- read.delim(
  gzfile("data/raw/TCGA_BRCA_clinical.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

survival <- read.delim(
  gzfile("data/raw/TCGA_BRCA_survival.tsv.gz"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

pam50 <- read.delim(
  "data/raw/TCGA_BRCA_PAM50.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

#==============================================================================
# 5. Explorando os dados
#==============================================================================

dim(counts)

dim(tpm)

dim(clinical)

dim(survival)

dim(pam50)

head(clinical)

head(survival)

head(pam50)

#==============================================================================
# 6. Informações importantes
#==============================================================================

cat("\n")
cat("=====================================\n")
cat("INFORMAÇÕES IMPORTANTES\n")
cat("=====================================\n\n")

cat("Counts: log2(count + 1)\n")

cat("TPM: log2(TPM + 1)\n")

cat("\n")

cat("Essas matrizes NÃO representam contagens brutas.\n")

cat("Os valores já estão transformados em log2(x + 1).\n")

cat("\n")

cat("Para análises de expressão diferencial (DESeq2/edgeR),\n")
cat("utilize contagens brutas obtidas pelo TCGAbiolinks.\n")

cat("\n")

cat("As matrizes do Xena são ideais para:\n")

cat("- Heatmaps\n")
cat("- Correlação\n")
cat("- PCA\n")
cat("- Machine Learning\n")
cat("- Visualização\n")
cat("- Sobrevida\n")

cat("\n")

###############################################################################
# Fim da Aula 01B
###############################################################################