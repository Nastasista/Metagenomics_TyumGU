---
title: "Tests_and_Hypotheses"
author: "Anastasia Poluzerova"
date: "2023-06-05"
output: 
  html_document: 
    keep_md: yes
---


```r
knitr::opts_chunk$set(fig.width = 10, fig.height = 6)
 
library('phyloseq')
library('tidyverse')
library('vegan')



set.seed(1609)
setwd('/home/nastasista/Metagenomics')
ps <- readRDS("ps.no.organells.RData")

# Select only samples with more than 8k reads per sample
ps <- prune_samples(sample_sums(ps) > 8000, ps)
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4460 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 4460 taxa by 6 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 4460 reference sequences ]
```

## Alpha-diversity


```r
# rarefy to minimal depth
ps.raref <- rarefy_even_depth(ps)
```

```r
plot_richness(ps.raref, x = "Source", measures=c("Observed", "Simpson"), color = "Site")
```

![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-1-1.png)<!-- -->
Создание таблицы альфа-разнообразия:


```r
alpha_div <- function(ps, measures){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  obs_sim <- estimate_richness(ps, measures = measures)
  Site <- ps@sam_data$Site
  Source <- ps@sam_data$Source
  alpha <- cbind(obs_sim, Site, Source)
  return(alpha)
}

alpha <- alpha_div(ps.raref, c("Observed", "Simpson"))
alpha
```

```
##                            Observed   Simpson Site             Source
## Self.growing.Dumps.B1.AY.1      353 0.9954960   B1 Self-growing Dumps
## Litostrat.B2.C.1                336 0.9949253   B2          Litostrat
## Litostrat.B2.C.2                307 0.9950977   B2          Litostrat
## Litostrat.B2.C.3                346 0.9946523   B2          Litostrat
## Litostrat.B2.C.4                312 0.9910318   B2          Litostrat
## Self.growing.Dumps.B1.AY.2      326 0.9954286   B1 Self-growing Dumps
## Coal.Mine.Terricon.B3.C.1       157 0.9770122   B3 Coal Mine Terricon
## Coal.Mine.Terricon.B3.C.2       117 0.9719444   B3 Coal Mine Terricon
## Coal.Mine.Terricon.B3.C.3       103 0.9706158   B3 Coal Mine Terricon
## Coal.Mine.Terricon.B3.C.4       119 0.9753717   B3 Coal Mine Terricon
## Self.growing.Dumps.B1.AY.3      477 0.9970006   B1 Self-growing Dumps
## Local.Reference.B4.AY.1         252 0.9941508   B4    Local Reference
## Local.Reference.B4.AY.2         255 0.9937443   B4    Local Reference
## Local.Reference.B4.AY.3         184 0.9917881   B4    Local Reference
## Self.growing.Dumps.B1.AY.4      335 0.9955142   B1 Self-growing Dumps
## Local.Reference.B4.AY.4         191 0.9922160   B4    Local Reference
## Embryo.Sand.B5.AY.1             403 0.9961900   B5        Embryo Sand
## Embryo.Sand.B5.AY.2             329 0.9953602   B5        Embryo Sand
## Embryo.Sand.B5.AY.3             344 0.9952160   B5        Embryo Sand
## Embryo.Sand.B5.AY.4             312 0.9950426   B5        Embryo Sand
## Regional.Reference.B6.AY.1      320 0.9954253   B6 Regional Reference
## Regional.Reference.B6.AY.2      290 0.9944651   B6 Regional Reference
## Regional.Reference.B6.AY.3      282 0.9941904   B6 Regional Reference
## Regional.Reference.B6.AY.4      252 0.9937457   B6 Regional Reference
```





```r
aov <- aov(Observed ~ Site + Source, data = alpha)
summary(aov)
```

```
##             Df Sum Sq Mean Sq F value   Pr(>F)    
## Site         5 172192   34438   21.41 5.31e-07 ***
## Residuals   18  28950    1608                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Результаты свидетельствуют о том, что переменная "Site" имеет статистически значимое влияние на переменную "Observed". 


## Beta-diversity


```r
ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))
ord.pcoa.bray <- ordinate(ps.prop, method='PCoA', distance='bray')
plot_ordination(ps.prop, ord.pcoa.bray, color = 'Source', shape = "Site") +
    geom_point(size=3, alpha=0.7) + 
    theme_light()
```

![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-4-1.png)<!-- -->



```r
# Calculate Bray-Curtis distance matrix
dist <- phyloseq::distance(ps, method = "bray")

sample_df <- data.frame(sample_data(ps))

# Perform PERMANOVA
permanova <- adonis2(dist ~ Site + Source, data = sample_df)
permanova
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = dist ~ Site + Source, data = sample_df)
##          Df SumOfSqs      R2      F Pr(>F)    
## Site      5   6.6519 0.64328 6.4919  0.001 ***
## Residual 18   3.6887 0.35672                  
## Total    23  10.3406 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Результаты говорят о том, что "Site" оказывает значимое влияние на структуру данных, поскольку его включение в модель дает значимый F-статистический тест и низкое p-значение (<0.001).

## CCA

Проверка связи между питательными факторами и микробными популяциями :


```r
# Collapse samples according to nutrition data. Add nutrition data to metadata section
ps@sam_data
```

```
##                                 Filename             Source Site Horizont
## Self-growing Dumps.B1.AY.1  Abacumov-B-1 Self-growing Dumps   B1       AY
## Litostrat.B2.C.1           Abacumov-B-13          Litostrat   B2        C
## Litostrat.B2.C.2           Abacumov-B-14          Litostrat   B2        C
## Litostrat.B2.C.3           Abacumov-B-15          Litostrat   B2        C
## Litostrat.B2.C.4           Abacumov-B-16          Litostrat   B2        C
## Self-growing Dumps.B1.AY.2  Abacumov-B-2 Self-growing Dumps   B1       AY
## Coal Mine Terricon.B3.C.1  Abacumov-B-25 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.2  Abacumov-B-26 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.3  Abacumov-B-27 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.4  Abacumov-B-28 Coal Mine Terricon   B3        C
## Self-growing Dumps.B1.AY.3  Abacumov-B-3 Self-growing Dumps   B1       AY
## Local Reference.B4.AY.1    Abacumov-B-37    Local Reference   B4       AY
## Local Reference.B4.AY.2    Abacumov-B-38    Local Reference   B4       AY
## Local Reference.B4.AY.3    Abacumov-B-39    Local Reference   B4       AY
## Self-growing Dumps.B1.AY.4  Abacumov-B-4 Self-growing Dumps   B1       AY
## Local Reference.B4.AY.4    Abacumov-B-40    Local Reference   B4       AY
## Embryo Sand.B5.AY.1        Abacumov-B-49        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.2        Abacumov-B-50        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.3        Abacumov-B-51        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.4        Abacumov-B-52        Embryo Sand   B5       AY
## Regional Reference.B6.AY.1 Abacumov-B-61 Regional Reference   B6       AY
## Regional Reference.B6.AY.2 Abacumov-B-62 Regional Reference   B6       AY
## Regional Reference.B6.AY.3 Abacumov-B-63 Regional Reference   B6       AY
## Regional Reference.B6.AY.4 Abacumov-B-64 Regional Reference   B6       AY
##                            Repeat                   SampleID
## Self-growing Dumps.B1.AY.1      1 Self-growing Dumps.B1.AY.1
## Litostrat.B2.C.1                1           Litostrat.B2.C.1
## Litostrat.B2.C.2                2           Litostrat.B2.C.2
## Litostrat.B2.C.3                3           Litostrat.B2.C.3
## Litostrat.B2.C.4                4           Litostrat.B2.C.4
## Self-growing Dumps.B1.AY.2      2 Self-growing Dumps.B1.AY.2
## Coal Mine Terricon.B3.C.1       1  Coal Mine Terricon.B3.C.1
## Coal Mine Terricon.B3.C.2       2  Coal Mine Terricon.B3.C.2
## Coal Mine Terricon.B3.C.3       3  Coal Mine Terricon.B3.C.3
## Coal Mine Terricon.B3.C.4       4  Coal Mine Terricon.B3.C.4
## Self-growing Dumps.B1.AY.3      3 Self-growing Dumps.B1.AY.3
## Local Reference.B4.AY.1         1    Local Reference.B4.AY.1
## Local Reference.B4.AY.2         2    Local Reference.B4.AY.2
## Local Reference.B4.AY.3         3    Local Reference.B4.AY.3
## Self-growing Dumps.B1.AY.4      4 Self-growing Dumps.B1.AY.4
## Local Reference.B4.AY.4         4    Local Reference.B4.AY.4
## Embryo Sand.B5.AY.1             1        Embryo Sand.B5.AY.1
## Embryo Sand.B5.AY.2             2        Embryo Sand.B5.AY.2
## Embryo Sand.B5.AY.3             3        Embryo Sand.B5.AY.3
## Embryo Sand.B5.AY.4             4        Embryo Sand.B5.AY.4
## Regional Reference.B6.AY.1      1 Regional Reference.B6.AY.1
## Regional Reference.B6.AY.2      2 Regional Reference.B6.AY.2
## Regional Reference.B6.AY.3      3 Regional Reference.B6.AY.3
## Regional Reference.B6.AY.4      4 Regional Reference.B6.AY.4
```

```r
ps <- subset_samples(ps, Site != "B6") #удаление образцов, для которых нет данных агрохимии


agro <- read.csv("agro.csv")
agro
```

```
##    Site Repeat  pH   P   K  NH4.
## 1    B1      1 7.4 178 265 26.99
## 2    B1      2 7.5 185 252 28.33
## 3    B1      3 7.4 258 219 19.92
## 4    B1      4 7.4 207 245 25.08
## 5    B2      1 4.3  16 223  3.23
## 6    B2      2 2.6  18 134  5.42
## 7    B2      3 3.0  12  92  8.29
## 8    B2      4 3.3  15 150  5.64
## 9    B3      1 5.2  80 370 46.54
## 10   B3      2 5.5  80 349 36.98
## 11   B3      3 5.4  52 332 27.29
## 12   B3      4 5.4  71 350 36.90
## 13   B4      1 5.2 170 126 16.14
## 14   B4      2 4.8 166 122  0.18
## 15   B4      3 5.1 229 126 13.16
## 16   B4      4 5.0 188 125  9.83
## 17   B5      1 4.3 104 303 57.87
## 18   B5      2 4.4 120 429 24.31
## 19   B5      3 4.3  83 235 19.74
## 20   B5      4 4.3 102 322 33.97
```

```r
ps@sam_data <- ps@sam_data[order(ps@sam_data$Site), ] #нужно упорядочить данные 
ps@sam_data
```

```
##                                 Filename             Source Site Horizont
## Self-growing Dumps.B1.AY.1  Abacumov-B-1 Self-growing Dumps   B1       AY
## Self-growing Dumps.B1.AY.2  Abacumov-B-2 Self-growing Dumps   B1       AY
## Self-growing Dumps.B1.AY.3  Abacumov-B-3 Self-growing Dumps   B1       AY
## Self-growing Dumps.B1.AY.4  Abacumov-B-4 Self-growing Dumps   B1       AY
## Litostrat.B2.C.1           Abacumov-B-13          Litostrat   B2        C
## Litostrat.B2.C.2           Abacumov-B-14          Litostrat   B2        C
## Litostrat.B2.C.3           Abacumov-B-15          Litostrat   B2        C
## Litostrat.B2.C.4           Abacumov-B-16          Litostrat   B2        C
## Coal Mine Terricon.B3.C.1  Abacumov-B-25 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.2  Abacumov-B-26 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.3  Abacumov-B-27 Coal Mine Terricon   B3        C
## Coal Mine Terricon.B3.C.4  Abacumov-B-28 Coal Mine Terricon   B3        C
## Local Reference.B4.AY.1    Abacumov-B-37    Local Reference   B4       AY
## Local Reference.B4.AY.2    Abacumov-B-38    Local Reference   B4       AY
## Local Reference.B4.AY.3    Abacumov-B-39    Local Reference   B4       AY
## Local Reference.B4.AY.4    Abacumov-B-40    Local Reference   B4       AY
## Embryo Sand.B5.AY.1        Abacumov-B-49        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.2        Abacumov-B-50        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.3        Abacumov-B-51        Embryo Sand   B5       AY
## Embryo Sand.B5.AY.4        Abacumov-B-52        Embryo Sand   B5       AY
##                            Repeat                   SampleID
## Self-growing Dumps.B1.AY.1      1 Self-growing Dumps.B1.AY.1
## Self-growing Dumps.B1.AY.2      2 Self-growing Dumps.B1.AY.2
## Self-growing Dumps.B1.AY.3      3 Self-growing Dumps.B1.AY.3
## Self-growing Dumps.B1.AY.4      4 Self-growing Dumps.B1.AY.4
## Litostrat.B2.C.1                1           Litostrat.B2.C.1
## Litostrat.B2.C.2                2           Litostrat.B2.C.2
## Litostrat.B2.C.3                3           Litostrat.B2.C.3
## Litostrat.B2.C.4                4           Litostrat.B2.C.4
## Coal Mine Terricon.B3.C.1       1  Coal Mine Terricon.B3.C.1
## Coal Mine Terricon.B3.C.2       2  Coal Mine Terricon.B3.C.2
## Coal Mine Terricon.B3.C.3       3  Coal Mine Terricon.B3.C.3
## Coal Mine Terricon.B3.C.4       4  Coal Mine Terricon.B3.C.4
## Local Reference.B4.AY.1         1    Local Reference.B4.AY.1
## Local Reference.B4.AY.2         2    Local Reference.B4.AY.2
## Local Reference.B4.AY.3         3    Local Reference.B4.AY.3
## Local Reference.B4.AY.4         4    Local Reference.B4.AY.4
## Embryo Sand.B5.AY.1             1        Embryo Sand.B5.AY.1
## Embryo Sand.B5.AY.2             2        Embryo Sand.B5.AY.2
## Embryo Sand.B5.AY.3             3        Embryo Sand.B5.AY.3
## Embryo Sand.B5.AY.4             4        Embryo Sand.B5.AY.4
```

```r
if (all(agro$Site == sample_data(ps)$Site)) {
  sample_data(ps) <- cbind(sample_data(ps)[, c("SampleID", "Site", "Source")], agro[, c("pH", "P", "K", "NH4.")])
  print("Replaced")
}
```

```
## [1] "Replaced"
```

```r
ps@sam_data
```

```
##                                              SampleID Site             Source
## Self-growing Dumps.B1.AY.1 Self-growing Dumps.B1.AY.1   B1 Self-growing Dumps
## Litostrat.B2.C.1                     Litostrat.B2.C.1   B2          Litostrat
## Litostrat.B2.C.2                     Litostrat.B2.C.2   B2          Litostrat
## Litostrat.B2.C.3                     Litostrat.B2.C.3   B2          Litostrat
## Litostrat.B2.C.4                     Litostrat.B2.C.4   B2          Litostrat
## Self-growing Dumps.B1.AY.2 Self-growing Dumps.B1.AY.2   B1 Self-growing Dumps
## Coal Mine Terricon.B3.C.1   Coal Mine Terricon.B3.C.1   B3 Coal Mine Terricon
## Coal Mine Terricon.B3.C.2   Coal Mine Terricon.B3.C.2   B3 Coal Mine Terricon
## Coal Mine Terricon.B3.C.3   Coal Mine Terricon.B3.C.3   B3 Coal Mine Terricon
## Coal Mine Terricon.B3.C.4   Coal Mine Terricon.B3.C.4   B3 Coal Mine Terricon
## Self-growing Dumps.B1.AY.3 Self-growing Dumps.B1.AY.3   B1 Self-growing Dumps
## Local Reference.B4.AY.1       Local Reference.B4.AY.1   B4    Local Reference
## Local Reference.B4.AY.2       Local Reference.B4.AY.2   B4    Local Reference
## Local Reference.B4.AY.3       Local Reference.B4.AY.3   B4    Local Reference
## Self-growing Dumps.B1.AY.4 Self-growing Dumps.B1.AY.4   B1 Self-growing Dumps
## Local Reference.B4.AY.4       Local Reference.B4.AY.4   B4    Local Reference
## Embryo Sand.B5.AY.1               Embryo Sand.B5.AY.1   B5        Embryo Sand
## Embryo Sand.B5.AY.2               Embryo Sand.B5.AY.2   B5        Embryo Sand
## Embryo Sand.B5.AY.3               Embryo Sand.B5.AY.3   B5        Embryo Sand
## Embryo Sand.B5.AY.4               Embryo Sand.B5.AY.4   B5        Embryo Sand
##                             pH   P   K  NH4.
## Self-growing Dumps.B1.AY.1 7.4 178 265 26.99
## Litostrat.B2.C.1           4.3  16 223  3.23
## Litostrat.B2.C.2           2.6  18 134  5.42
## Litostrat.B2.C.3           3.0  12  92  8.29
## Litostrat.B2.C.4           3.3  15 150  5.64
## Self-growing Dumps.B1.AY.2 7.5 185 252 28.33
## Coal Mine Terricon.B3.C.1  5.2  80 370 46.54
## Coal Mine Terricon.B3.C.2  5.5  80 349 36.98
## Coal Mine Terricon.B3.C.3  5.4  52 332 27.29
## Coal Mine Terricon.B3.C.4  5.4  71 350 36.90
## Self-growing Dumps.B1.AY.3 7.4 258 219 19.92
## Local Reference.B4.AY.1    5.2 170 126 16.14
## Local Reference.B4.AY.2    4.8 166 122  0.18
## Local Reference.B4.AY.3    5.1 229 126 13.16
## Self-growing Dumps.B1.AY.4 7.4 207 245 25.08
## Local Reference.B4.AY.4    5.0 188 125  9.83
## Embryo Sand.B5.AY.1        4.3 104 303 57.87
## Embryo Sand.B5.AY.2        4.4 120 429 24.31
## Embryo Sand.B5.AY.3        4.3  83 235 19.74
## Embryo Sand.B5.AY.4        4.3 102 322 33.97
```

Make CCA for the top 1000 most abundant ASVs


```r
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

ps.top1k <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:1000]) %>% 
  prune_taxa(ps)

otus.ps.vegan <- veganifyOTU(ps.top1k)
metadata <- as(sample_data(ps.top1k), "data.frame")




vare.cca <- vegan::cca(otus.ps.vegan ~  pH + P + K + NH4., data=metadata)
vare.cca
```

```
## Call: cca(formula = otus.ps.vegan ~ pH + P + K + NH4., data = metadata)
## 
##               Inertia Proportion Rank
## Total          7.3433     1.0000     
## Constrained    2.7181     0.3702    4
## Unconstrained  4.6252     0.6298   15
## Inertia is scaled Chi-square 
## 
## Eigenvalues for constrained axes:
##   CCA1   CCA2   CCA3   CCA4 
## 0.9306 0.8344 0.7233 0.2298 
## 
## Eigenvalues for unconstrained axes:
##    CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8    CA9   CA10   CA11 
## 0.9104 0.4947 0.4056 0.3838 0.3683 0.3528 0.3175 0.3096 0.2939 0.2704 0.2134 
##   CA12   CA13   CA14   CA15 
## 0.1755 0.0570 0.0489 0.0234
```

Check the CCA model


```r
anova(vare.cca)
```

```
## Permutation test for cca under reduced model
## Permutation: free
## Number of permutations: 999
## 
## Model: cca(formula = otus.ps.vegan ~ pH + P + K + NH4., data = metadata)
##          Df ChiSquare      F Pr(>F)    
## Model     4    2.7181 2.2038  0.001 ***
## Residual 15    4.6252                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(vare.cca, by="terms")
```

```
## Permutation test for cca under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## Model: cca(formula = otus.ps.vegan ~ pH + P + K + NH4., data = metadata)
##          Df ChiSquare      F Pr(>F)    
## pH        1    0.8899 2.8861  0.001 ***
## P         1    0.7785 2.5248  0.001 ***
## K         1    0.7738 2.5094  0.001 ***
## NH4.      1    0.2760 0.8950  0.688    
## Residual 15    4.6252                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
Основываясь на результате анализа, можно сказать, что pH, P и K имеют статистически значимое влияние, в то время как NH4. не является значимым фактором. 


```r
vif.cca(vare.cca) #проверка на мультиколинеарность
```

```
##       pH        P        K     NH4. 
## 3.567322 3.294638 2.630232 2.369083
```
На основе значений VIF можно сделать вывод, что факторы pH, P, K и NH4 не демонстрируют сильной мультиколлинеарности между собой.


```r
vare.cca <- vegan::cca(otus.ps.vegan ~  pH + P + K, data=metadata) #убираю NH4.
anova(vare.cca, by="terms")
```

```
## Permutation test for cca under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## Model: cca(formula = otus.ps.vegan ~ pH + P + K, data = metadata)
##          Df ChiSquare      F Pr(>F)    
## pH        1    0.8899 2.9052  0.001 ***
## P         1    0.7785 2.5415  0.001 ***
## K         1    0.7738 2.5260  0.001 ***
## Residual 16    4.9011                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



```r
species.data <- vare.cca$CCA$v %>% 
               data.frame() %>% 
               mutate(ASV = rownames(.)) %>% 
               inner_join(data.frame(ASV = names(taxa_sums(ps.top1k)),
                                     Total.abund = taxa_sums(ps.top1k),
                                     ps.top1k@tax_table[,2], # Phylum
                                     ps.top1k@tax_table[,3]), # Class
                          by = "ASV")
species.data %>% head(10)
```

```
##         CCA1       CCA2       CCA3   ASV Total.abund         Phylum
## 1  0.4904147  1.3717209  0.5132103  ASV1        2983    Chloroflexi
## 2  0.4876729  1.3848019  0.5472641  ASV2        2672    Chloroflexi
## 3  0.4884236  1.3865418  0.5367946  ASV3        1855    Chloroflexi
## 4  0.4897593  1.3871675  0.5352946  ASV4        1663    Chloroflexi
## 5  1.0009242 -0.9909102  0.7841744  ASV5        1643  Cyanobacteria
## 6  1.0046774 -1.0148785  0.7611772  ASV6        1483  Cyanobacteria
## 7  0.4829375  0.5822380 -1.3766719  ASV7         625 Proteobacteria
## 8  0.4916113  0.5818614 -1.4447804  ASV8         640 Proteobacteria
## 9  0.4846995  1.3949187  0.4786259  ASV9        1252 Actinobacteria
## 10 0.4929340  1.3687222  0.5538311 ASV10        1251    Chloroflexi
##                  Class
## 1                  AD3
## 2                  AD3
## 3                  AD3
## 4                  AD3
## 5     Oxyphotobacteria
## 6     Oxyphotobacteria
## 7  Alphaproteobacteria
## 8  Alphaproteobacteria
## 9       Acidimicrobiia
## 10                 AD3
```

```r
samples.data <- vare.cca$CCA$u %>% 
  data.frame() %>% 
  mutate(Names = rownames(.)) %>% 
  inner_join(ps@sam_data, by = c("Names" = "SampleID"))


samples.data
```

```
##          CCA1          CCA2        CCA3                      Names Site
## 1  -1.2730625  0.6506876558  1.17336076 Self-growing Dumps.B1.AY.1   B1
## 2   0.8636058  0.0009770234  1.25746384           Litostrat.B2.C.1   B2
## 3   1.1417695 -1.2791597668 -0.21324753           Litostrat.B2.C.2   B2
## 4   0.9345713 -1.5160815635  0.70182651           Litostrat.B2.C.3   B2
## 5   0.9835285 -0.9272212435  0.58693055           Litostrat.B2.C.4   B2
## 6  -1.3977156  0.5376135020  1.23404755 Self-growing Dumps.B1.AY.2   B1
## 7   0.5073659  1.3230785225 -0.09775418  Coal Mine Terricon.B3.C.1   B3
## 8   0.3480329  1.2218375692  0.41539007  Coal Mine Terricon.B3.C.2   B3
## 9   0.5467038  1.1468036687  1.02752001  Coal Mine Terricon.B3.C.3   B3
## 10  0.4529528  1.2358367536  0.48387303  Coal Mine Terricon.B3.C.4   B3
## 11 -2.0245673 -0.0555450217 -0.24860687 Self-growing Dumps.B1.AY.3   B1
## 12 -0.9012440 -1.1780107068 -0.32875631    Local Reference.B4.AY.1   B4
## 13 -0.7517014 -1.3121876350 -0.69765362    Local Reference.B4.AY.2   B4
## 14 -1.3234951 -1.4261007080 -1.74019631    Local Reference.B4.AY.3   B4
## 15 -1.5552230  0.3651627311  0.68156962 Self-growing Dumps.B1.AY.4   B1
## 16 -0.9777321 -1.3106965167 -0.95768981    Local Reference.B4.AY.4   B4
## 17  0.4195232  0.3834397632 -1.23649881        Embryo Sand.B5.AY.1   B5
## 18  0.6326644  1.4707365725 -2.36222796        Embryo Sand.B5.AY.2   B5
## 19  0.3822257 -0.1419891950 -0.29314352        Embryo Sand.B5.AY.3   B5
## 20  0.4905996  0.5595476362 -1.32803154        Embryo Sand.B5.AY.4   B5
##                Source  pH   P   K  NH4.
## 1  Self-growing Dumps 7.4 178 265 26.99
## 2           Litostrat 4.3  16 223  3.23
## 3           Litostrat 2.6  18 134  5.42
## 4           Litostrat 3.0  12  92  8.29
## 5           Litostrat 3.3  15 150  5.64
## 6  Self-growing Dumps 7.5 185 252 28.33
## 7  Coal Mine Terricon 5.2  80 370 46.54
## 8  Coal Mine Terricon 5.5  80 349 36.98
## 9  Coal Mine Terricon 5.4  52 332 27.29
## 10 Coal Mine Terricon 5.4  71 350 36.90
## 11 Self-growing Dumps 7.4 258 219 19.92
## 12    Local Reference 5.2 170 126 16.14
## 13    Local Reference 4.8 166 122  0.18
## 14    Local Reference 5.1 229 126 13.16
## 15 Self-growing Dumps 7.4 207 245 25.08
## 16    Local Reference 5.0 188 125  9.83
## 17        Embryo Sand 4.3 104 303 57.87
## 18        Embryo Sand 4.4 120 429 24.31
## 19        Embryo Sand 4.3  83 235 19.74
## 20        Embryo Sand 4.3 102 322 33.97
```



```r
# plot species
ggplot() +
  geom_point(data=species.data,
             aes(x=CCA1, y=CCA2, color=Phylum, size=Total.abund), alpha=0.9) +
  geom_segment(data = vare.cca$CCA$biplot %>% data.frame(), 
               aes(x = 0, xend = CCA1, y = 0, yend = CCA2), 
               alpha=0.8, color = "black",arrow = arrow(angle = 3)) +
  geom_text(data = vare.cca$CCA$biplot %>% 
                    data.frame() %>% 
                    mutate(Label = rownames(.)), 
            aes(x=CCA1, y=CCA2, label= Label,
                hjust = -0.5), size=4) +
  theme_light() +
  ggtitle("A: Species")
```

![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

Черные точки отражают размер(численность) филума



```r
major.phyla <- species.data %>% 
  group_by(Phylum) %>% 
  summarize(sum = sum(Total.abund)) %>% 
  arrange(desc(sum)) %>% 
  select(Phylum) %>% 
  head(10) %>% 
  as.vector()

# с помощью цикла рисуется график для каждой Phyla
for (i in major.phyla$Phylum) {
  p <- ggplot() +
    geom_point(data=species.data,
               aes(x=CCA1, y=CCA2, size=Total.abund), alpha=0.2, color="grey80") +
    geom_point(data=species.data %>% filter(Phylum == i),
               aes(x=CCA1, y=CCA2, color=Class, size=Total.abund), alpha=0.9) +
    geom_segment(data = vare.cca$CCA$biplot %>% data.frame(), 
                 aes(x = 0, xend = CCA1, y = 0, yend = CCA2), 
                 alpha=0.8, color = "black",arrow = arrow(angle = 3)) +
    geom_text(data = vare.cca$CCA$biplot %>% 
                      data.frame() %>% 
                      mutate(Label = rownames(.)), 
              aes(x=CCA1, y=CCA2, label= Label,
                  hjust = -0.5), size=4) +
    theme_light() +
    ggtitle(i)
  print(p)
}
```

![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-1.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-2.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-3.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-4.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-5.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-6.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-7.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-8.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-9.png)<!-- -->![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-14-10.png)<!-- -->




```r
# plot samples
ggplot() +
  geom_point(data=samples.data,
             aes(x=CCA1, y=CCA2, color=Source, shape = Site), size=3, alpha=0.7) +
  geom_segment(data = vare.cca$CCA$biplot %>% data.frame(), 
               aes(x = 0, xend = CCA1, y = 0, yend = CCA2), 
               alpha=0.8, color = "black",arrow = arrow(angle = 3)) +
  geom_text(data = vare.cca$CCA$biplot %>% 
                    data.frame() %>% 
                    mutate(Label = rownames(.)), 
            aes(x=CCA1, y=CCA2, label= Label,
                hjust = -0.5), size=4) +
  theme_light() +
  ggtitle("B. Samples")
```

![](Tests_and_Hypotheses_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

