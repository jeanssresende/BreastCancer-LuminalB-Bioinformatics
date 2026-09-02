###############################################################################
# Projeto:
# BreastCancer-LuminalB-Bioinformatics
#
# Aula 02
# Organização e seleção da coorte TCGA-BRCA
#
# Objetivos:
# 1. Carregar os dados obtidos pelo UCSC Xena
# 2. Explorar a estrutura dos arquivos
# 3. Entender os identificadores das amostras TCGA
# 4. Identificar tumores primários
# 5. Entender a relação entre amostra e paciente
# 6. Integrar expressão, PAM50 e sobrevida
# 7. Selecionar pacientes Luminal B
# 8. Preparar a coorte para as próximas análises
#
# Autor:
# Jean Resende
###############################################################################


#==============================================================================
# 1. CARREGAR OS DADOS
#==============================================================================

counts <- read.delim(
  "data/raw/TCGA_BRCA_counts.tsv.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

tpm <- read.delim(
  "data/raw/TCGA_BRCA_tpm.tsv.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

clinical <- read.delim(
  "data/raw/TCGA_BRCA_clinical.tsv.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

survival <- read.delim(
  "data/raw/TCGA_BRCA_survival.tsv.gz",
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
# 2. CONHECENDO OS OBJETOS
#==============================================================================

class(counts)

class(tpm)

class(clinical)

class(survival)

class(pam50)


#==============================================================================
# 3. DIMENSÕES DOS DADOS
#==============================================================================

dim(counts)

dim(tpm)

dim(clinical)

dim(survival)

dim(pam50)


#==============================================================================
# 4. CONHECENDO AS COLUNAS
#==============================================================================

names(clinical)

names(survival)

names(pam50)


#==============================================================================
# 5. VISUALIZANDO OS DADOS
#==============================================================================

head(clinical)

head(survival)

head(pam50)


#==============================================================================
# 6. ENTENDENDO A MATRIZ DE EXPRESSÃO
#==============================================================================

# A primeira coluna contém os identificadores dos genes.
# As demais colunas representam as amostras.

head(counts[, 1:5])


# Nome da primeira coluna

names(counts)[1]


# Nomes das primeiras amostras

head(names(counts)[-1])


#==============================================================================
# 7. CRIANDO A MATRIZ DE EXPRESSÃO
#==============================================================================

# A primeira coluna contém os identificadores dos genes
gene_id <- counts[, 1]


# Remover a primeira coluna
expression <- counts[, -1]


# Utilizar os identificadores dos genes como nomes das linhas
rownames(expression) <- gene_id


# Verificar
dim(expression)

head(rownames(expression))

head(colnames(expression))


#==============================================================================
# 8. VERIFICANDO OS IDENTIFICADORES ENSEMBL
#==============================================================================

# Os identificadores possuem uma versão após o ponto.
#
# Exemplo:
#
# ENSG00000000003.15
# ENSG00000000005.6
# ENSG00000000419.13
#
# Neste momento NÃO vamos remover a versão.
#
# O identificador original fornecido pelo Xena será mantido.


head(rownames(expression))


# Verificar se os identificadores originais são únicos

sum(duplicated(rownames(expression)))


#==============================================================================
# 9. CRIANDO UM ID ENSEMBL SEM VERSÃO
#==============================================================================

# Para futuras etapas de anotação, podemos criar uma versão
# do identificador sem a informação de versão.
#
# IMPORTANTE:
#
# Não vamos substituir os rownames da matriz.
#
# A matriz continuará utilizando os IDs originais do Xena.


gene_id_clean <- sub("\\..*", "", rownames(expression))


# Visualizar
head(gene_id_clean)


# Verificar se a remoção da versão gera duplicidades
sum(duplicated(gene_id_clean))


#==============================================================================
# 10. CONFERINDO A ESCALA DOS DADOS
#==============================================================================

# O Xena disponibiliza esta matriz como:
#
# log2(count + 1)
#
# Portanto, os valores NÃO são counts brutos.

summary(as.numeric(as.matrix(expression[, 1:10])))


# Valores mínimo e máximo

range(as.numeric(as.matrix(expression[, 1:10])), na.rm = TRUE)


#==============================================================================
# 11. CONHECENDO OS BARCODES DO TCGA
#==============================================================================

samples <- colnames(expression)


# Visualizar alguns barcodes

head(samples)


#==============================================================================
# 12. ENTENDENDO O BARCODE TCGA
#==============================================================================

# O barcode TCGA contém informações sobre o paciente
# e sobre a amostra.


# Identificador do paciente
patient_id <- substr(samples, 1, 12)


# Código do tipo de amostra
sample_type_code <- substr(samples, 14, 15)


# Criar tabela de informações das amostras
sample_info <- data.frame(
  sample_barcode = samples,
  patient_id = patient_id,
  sample_type_code = sample_type_code,
  stringsAsFactors = FALSE
)


# Visualizar
head(sample_info)

#==============================================================================
# 13. QUANTAS AMOSTRAS E PACIENTES TEMOS?
#==============================================================================

# Número de amostras
nrow(sample_info)


# Número de pacientes
length(unique(sample_info$patient_id))


# Número de amostras por paciente
table(table(sample_info$patient_id))


#==============================================================================
# 14. IDENTIFICANDO O TIPO DE AMOSTRA
#==============================================================================

table(sample_info$sample_type_code)


# Códigos utilizados pelo TCGA:
#
# 01 = Primary Solid Tumor
# 06 = Metastatic
# 11 = Solid Tissue Normal


sample_info$sample_type <- ifelse(
  sample_info$sample_type_code == "01",
  "Primary Tumor",
  ifelse(
    sample_info$sample_type_code == "11",
    "Normal",
    ifelse(
      sample_info$sample_type_code == "06",
      "Metastatic",
      "Other"
    )
  )
)


# Verificar a distribuição
table(sample_info$sample_type)


#==============================================================================
# 15. SELECIONANDO APENAS TUMORES PRIMÁRIOS
#==============================================================================

primary_samples <- sample_info[sample_info$sample_type == "Primary Tumor",]


# Número de amostras
nrow(primary_samples)

# Número de pacientes
length(unique(primary_samples$patient_id))

# Visualizar
head(primary_samples)


#==============================================================================
# 16. VERIFICANDO OS DADOS PAM50
#==============================================================================

names(pam50)

head(pam50)


# Distribuição dos subtipos
table(pam50$PAM50, useNA = "ifany")


#==============================================================================
# 17. PREPARANDO O PAM50
#==============================================================================

# O arquivo possui:
#
# sample
# PAM50
#
# Vamos criar o identificador do paciente a partir do sample.


pam50_info <- data.frame(
  patient_id = substr(pam50$sample, 1, 12),
  PAM50 = pam50$PAM50,
  stringsAsFactors = FALSE
)


# Visualizar
head(pam50_info)


#==============================================================================
# 18. VERIFICANDO DUPLICATAS NO PAM50
#==============================================================================

# Quantos pacientes aparecem mais de uma vez?
sum(duplicated(pam50_info$patient_id))


# Visualizar pacientes duplicados
pam50_info$patient_id[duplicated(pam50_info$patient_id)]


#==============================================================================
# 19. VERIFICANDO SE UM PACIENTE POSSUI MAIS DE UM PAM50
#==============================================================================

# Número de classificações diferentes por paciente
pam50_per_patient <- aggregate(
  PAM50 ~ patient_id,
  data = pam50_info,
  FUN = function(x) length(unique(x))
)

# Pacientes com mais de uma classificação
pam50_per_patient[pam50_per_patient$PAM50 > 1,]

#==============================================================================
# 20. VERIFICANDO A SOBREPOSIÇÃO ENTRE OS DATASETS
#==============================================================================

# Pacientes na expressão
expression_patients <- unique(sample_info$patient_id)

# Pacientes no PAM50
pam50_patients <- unique(pam50_info$patient_id)

# Número de pacientes presentes nos dois datasets

length(
  intersect(
    expression_patients,
    pam50_patients
  )
)


#==============================================================================
# 21. INTEGRANDO TUMORES PRIMÁRIOS COM PAM50
#==============================================================================

cohort <- merge(
  primary_samples,
  pam50_info,
  by = "patient_id"
)


# Verificar
dim(cohort)

head(cohort)


#==============================================================================
# 22. DISTRIBUIÇÃO DOS SUBTIPOS
#==============================================================================

table(cohort$PAM50, useNA = "ifany")


#==============================================================================
# 23. VERIFICANDO DUPLICAÇÃO DE PACIENTES APÓS A INTEGRAÇÃO
#==============================================================================

sum(duplicated(cohort$patient_id))


# Visualizar alguns pacientes duplicados

head(
  cohort$patient_id[
    duplicated(cohort$patient_id)
  ]
)


#==============================================================================
# 24. SELECIONANDO LUMINAL B
#==============================================================================

luminal_b <- cohort[cohort$PAM50 == "LumB",]


# Verificar
dim(luminal_b)

head(luminal_b)


# Número de amostras
nrow(luminal_b)

# Número de pacientes

length(unique(luminal_b$patient_id))


#==============================================================================
# 25. VERIFICAR PACIENTES REPETIDOS EM LUMINAL B
#==============================================================================

sum(duplicated(luminal_b$patient_id))


# Visualizar pacientes repetidos
luminal_b$patient_id[duplicated(luminal_b$patient_id)]


#==============================================================================
# 26. SELECIONAR UMA AMOSTRA POR PACIENTE
#==============================================================================

# O projeto trabalha no nível do paciente.
#
# Portanto, precisamos garantir que cada paciente
# seja representado uma única vez na coorte final.
#
# Neste momento estamos selecionando a primeira amostra
# disponível para cada paciente.
#
# A escolha definitiva da amostra será revisada caso
# existam múltiplas amostras primárias relevantes.


luminal_b <- luminal_b[!duplicated(luminal_b$patient_id),]


# Verificar novamente

nrow(luminal_b)

length(unique(luminal_b$patient_id))


#==============================================================================
# 27. PREPARANDO OS DADOS DE SOBREVIDA
#==============================================================================

# O arquivo survival possui:
#
# sample
# OS.time
# OS
# _PATIENT


head(survival)


# Criar tabela simplificada de sobrevida

survival_info <- data.frame(
  patient_id = survival$`_PATIENT`,
  OS.time = survival$OS.time,
  OS = survival$OS,
  stringsAsFactors = FALSE
)


# Visualizar
head(survival_info)


#==============================================================================
# 28. VERIFICANDO DUPLICATAS NA SOBREVIDA
#==============================================================================

sum(duplicated(survival_info$patient_id))

#==============================================================================
# 29. VERIFICANDO A SOBREPOSIÇÃO COM A SOBREVIDA
#==============================================================================

length(
  intersect(
    luminal_b$patient_id,
    survival_info$patient_id
  )
)


#==============================================================================
# 30. INTEGRANDO SOBREVIDA À COORTE LUMINAL B
#==============================================================================

luminal_b <- merge(
  luminal_b,
  survival_info,
  by = "patient_id",
  all.x = TRUE
)


# Verificar
dim(luminal_b)

head(luminal_b)


#==============================================================================
# 31. VERIFICANDO OS DADOS DE SOBREVIDA
#==============================================================================

# Pacientes com informação de OS
sum(!is.na(luminal_b$OS))


# Pacientes sem informação de OS

sum(is.na(luminal_b$OS))


# Distribuição do evento
table(luminal_b$OS, useNA = "ifany")


#==============================================================================
# 32. RECUPERANDO OS BARCODES DAS AMOSTRAS
#==============================================================================

head(luminal_b$sample_barcode)


#==============================================================================
# 33. CRIANDO VETOR DAS AMOSTRAS LUMINAL B
#==============================================================================

luminal_b_samples <- luminal_b$sample_barcode


# Número de amostras
length(luminal_b_samples)


#==============================================================================
# 34. SUBSET DA MATRIZ DE EXPRESSÃO
#==============================================================================

expression_luminal_b <- expression[, 
                                   colnames(expression) %in% luminal_b_samples,
                                   drop = FALSE]


# Verificar
dim(expression_luminal_b)


#==============================================================================
# 35. SUBSET DA MATRIZ TPM
#==============================================================================

# A matriz TPM possui a mesma estrutura:
#
# primeira coluna = gene
# demais colunas = amostras
#
# Os dados estão em:
#
# log2(TPM + 1)


gene_id_tpm <- tpm[, 1]

expression_tpm <- tpm[, -1]

rownames(expression_tpm) <- gene_id_tpm


# Manter os identificadores originais do Xena

head(rownames(expression_tpm))


# Verificar se os IDs originais são únicos

sum(duplicated(rownames(expression_tpm)))


# Criar versão dos IDs sem versão para futura anotação

gene_id_tpm_clean <- sub("\\..*", "", rownames(expression_tpm))


# Verificar possíveis duplicidades após remover a versão

sum(duplicated(gene_id_tpm_clean))


#==============================================================================
# 36. SELECIONAR SOMENTE AS AMOSTRAS LUMINAL B
#==============================================================================

expression_tpm_luminal_b <- expression_tpm[,
                                           colnames(expression_tpm) %in% 
                                             luminal_b_samples,
                                           drop = FALSE]


# Verificar
dim(expression_tpm_luminal_b)


#==============================================================================
# 37. CONFERIR AS AMOSTRAS DE EXPRESSÃO
#==============================================================================

all(colnames(expression_luminal_b) %in% luminal_b$sample_barcode)


# Conferir TPM

all(colnames(expression_tpm_luminal_b) %in% luminal_b$sample_barcode)


#==============================================================================
# 38. CONFERIR O NÚMERO DE AMOSTRAS
#==============================================================================

ncol(expression_luminal_b)

ncol(expression_tpm_luminal_b)

nrow(luminal_b)


#==============================================================================
# 39. RESUMO FINAL DA COORTE
#==============================================================================

cat("\n")
cat("============================================\n")
cat("RESUMO DA COORTE TCGA-BRCA LUMINAL B\n")
cat("============================================\n")

cat(
  "Amostras de expressão:",
  ncol(expression_luminal_b),
  "\n"
)

cat(
  "Pacientes:",
  length(unique(luminal_b$patient_id)),
  "\n"
)

cat(
  "Genes:",
  nrow(expression_luminal_b),
  "\n"
)

cat(
  "Pacientes com OS:",
  sum(!is.na(luminal_b$OS)),
  "\n"
)

cat(
  "Eventos de OS:",
  sum(luminal_b$OS == 1, na.rm = TRUE),
  "\n"
)

cat("\n")
cat("============================================\n")


#==============================================================================
# 40. SALVAR A COORTE
#==============================================================================

dir.create(
  "data/processed",
  showWarnings = FALSE
)


# Informações da coorte

saveRDS(
  luminal_b,
  "data/processed/TCGA_BRCA_LuminalB_samples.rds"
)


# Matriz de expressão

saveRDS(
  expression_luminal_b,
  "data/processed/TCGA_BRCA_LuminalB_expression.rds"
)


# Matriz TPM

saveRDS(
  expression_tpm_luminal_b,
  "data/processed/TCGA_BRCA_LuminalB_TPM.rds"
)


###############################################################################
# Fim da Aula 02
###############################################################################