setwd('C:/Users/jterr/Documents/USP trabalhos/FIgenomica/ucs_tcga_pan_can_atlas_2018')
dir('.')
install.packages("stringr")
install.packages("dplyr")
library(stringr)
library(dplyr)

clin <- read.delim('data_clinical_patient.txt', skip = 4)
clin[1:4, 1:6]

mut <- read.delim('data_mutations_extended.txt')
mut[1:10,1:3]

rna <- read.delim('data_RNA_Seq_v2_expression_median.txt')
rna[1:7, 1:4]

str_sub('12345', -2,5)
nchar('12345')
table(str_sub(colnames(rna), -2, nchar(colnames(rna)))[3:ncol(rna)])

isNotnas <- which(!rna$Hugo_Symbol == '')
isNotnas
rna <- rna[isNotnas, ]
dups <- which(duplicated(rna$Hugo_Symbol))
rna[dups, 1:4]
rna_f <- rna[-dups, ]
which(duplicated(rna_f$Hugo_Symbol))
rownames(rna_f) <- rna_f$Hugo_Symbol
rna_f[1:4,1:4]
rna_f <- rna_f[ , -c(1,2)]
class(rna_f)
rna_m <- as.matrix(rna_f)
rna_m[1:4,1:4]
is.numeric(rna_m)

colnames(rna_m) <- str_sub(colnames(rna_m),1,12)
head(clin$PATIENT_ID)
?gsub
colnames(rna_m) <- gsub('.', '-', colnames(rna_m), fixed = TRUE)
all(colnames(rna_m) %in% clin$PATIENT_ID)

rna_m <- rna_m[ , clin$PATIENT_ID]
all(colnames(rna_m) == clin$PATIENT_ID)

clin[1:4,1:4]

boxplot(log(rna_m['TP53', ] + 1))
boxplot(log(rna_m['RAD51', ] + 1))
plot(log(rna_m['TP53', ] + 1) ~ log(rna_m['RAD51', ] + 1), pch=19)
abline(lm(log(rna_m['TP53', ] + 1) ~ log(rna_m['RAD51', ] + 1)), col='red')
cor(log(rna_m['TP53', ] + 1), log(rna_m['RAD51', ] + 1))
head(clin)
table(clin$DSS_STATUS)
clin$DSS_STATUS[clin$DSS_STATUS == ''] <- NA
boxplot(log(rna_m['TP53', ] + 1) ~ clin$DSS_STATUS)
t.test((log(rna_m['TP53', ] + 1) ~ clin$DSS_STATUS))
boxplot(log(rna_m['RAD51', ] + 1) ~ clin$DSS_STATUS)
t.test((log(rna_m['RAD51', ] + 1) ~ clin$DSS_STATUS))

for (row in rownames(rna_m)[1:3]) {
  print(t.test(rna_m[row , ]  ~ clin$DSS_STATUS)$p.value)
}
myfun <- function(x)
  t.test(x  ~ clin$DSS_STATUS)$p.value
pvals <- apply(log(rna_m + 1), 1, myfun)
df_ex <- data.frame(gene=rownames(rna_m),
                    pval=pvals)
head(df_ex)
dim(df_ex)
df_de <- subset(df_ex, df_ex$pval < 0.05)
df_ex$pajust <- p.adjust(df_ex$pval, method = 'fdr')
df_adjust <- subset(df_ex, df_ex$pajust < 0.1)

#Na última aula nés carregamos 3 arquivos do banco de dados do TCGA de amostras de Carcinossarcoma Uterino. Um arquivo tem dados clínicos e da patologia, o outro tem dados sobre mutação e o outro é uma matriz de expressão de todos os genes (RNAseq). 

# 1. Crie uma função para mostrar os 10 genes mais frequentemente mutados nesses tumores.

fun10<-function()
  {freq<-as.data.frame(table(mut$Hugo_Symbol))
  mut10=filter(freq, row_number(desc(freq$Freq)) <= 10)
  return(mut10)}


# 2. Crie uma função que retorna qual a classificação das variantes (Variant_Classification) em cada um dos Top 10 genes.

class_mut <- function()
{
  top10 <- mut10
  names_variant <- top10$Var1
  tamanho <- length(names_variant)
  df <- data.frame()
  
for (i in 1:tamanho) 
  {
  frame<-subset.data.frame(mut, mut$Hugo_Symbol==muttop$Var1[i])
  df <- rbind(df, frame)
  }
  df <- df[ ,c(1,9)]
}
  
# 3. Veja que os genes tem diferentes proporções de variantes. Qual seria o significado biológico desta observação?


  
# 4.  Crie uma função que retorna quais os exons (Exon_Number) estão acometido em cada um dos Top 10 genes.

class_exon <- function()
{
  top10_exon <- mut10
  names_exon <- top10_exon$Var1
  tamanho <- length(names_exon)
  df_exon <- data.frame()

    for (i in 1:tamanho) 
     {
    frame<-subset.data.frame(mut, mut$Hugo_Symbol==muttop$Var1[i])
    df_exon <- rbind(df_exon, frame)
    }
  df_exon <- df_exon[ ,c(1,111)]
  df_exon
}
  
# 5. Você vai ver que em alguns genes as mutações ocorrem em múltiplos exons, em outros, eles ocorrem em um número limitado de exons. Você acha que isto tem alguma implicação biológica? Explique.

# 6. Classifique os tumores baseado na presença ou ausência de mutação em cada um dos 10 Top genes e coloque esta informação na tabela de dados clínicos. Ou seja, você terá que inserir 10 colunas na tabela e ela vai ter a informação se a amostra tem ou não a mutação no gene. Incluir somente se a mutação for dos tipos Frame_Shift_Del Frame_Shift_Ins, Missense_Mutation, Nonsense_Mutation e Splice_Site.

df_del <- df[ which(df$Variant_Classification=='Frame_Shift_Del'| df$Variant_Classification
          == 'Frame_Shift_Ins'|df$Variant_Classification=='Missense_Mutation'|
            df$Variant_Classification=='Nonsense_Mutation'|df$Variant_Classification=='Splice_Site'), ]

list <- df_del$Tumor_Sample_Barcode

df_del$id <- substr(list,1,nchar(list)-3)

# 7. Crie uma função para plotar a expressão dos genes baseado na presença ou não de mutação (boxplot) e para gerar uma análise estatística onde H0: a expressão do gene Xi é a mesma em amostras com e sem mutação e H1: a expressão do gene Xi é diferente entre as amostras com e sem mutação.


