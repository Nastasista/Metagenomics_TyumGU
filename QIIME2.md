## Полузёрова Анастасия 
## Магистратура ТюмГУ
## Практическое задание 3
# <div style="text-align:center;"><font color="blue">QIIME 2 </font></div>





```python
!conda --version
```

    conda 23.5.0



```python
!ls data_met/sequences
```

    Abacumov-B-13_S13_L001_R1_001.fastq.gz	Abacumov-B-39_S39_L001_R1_001.fastq.gz
    Abacumov-B-13_S13_L001_R2_001.fastq.gz	Abacumov-B-39_S39_L001_R2_001.fastq.gz
    Abacumov-B-14_S14_L001_R1_001.fastq.gz	Abacumov-B-3_S3_L001_R1_001.fastq.gz
    Abacumov-B-14_S14_L001_R2_001.fastq.gz	Abacumov-B-3_S3_L001_R2_001.fastq.gz
    Abacumov-B-15_S15_L001_R1_001.fastq.gz	Abacumov-B-40_S40_L001_R1_001.fastq.gz
    Abacumov-B-15_S15_L001_R2_001.fastq.gz	Abacumov-B-40_S40_L001_R2_001.fastq.gz
    Abacumov-B-16_S16_L001_R1_001.fastq.gz	Abacumov-B-49_S49_L001_R1_001.fastq.gz
    Abacumov-B-16_S16_L001_R2_001.fastq.gz	Abacumov-B-49_S49_L001_R2_001.fastq.gz
    Abacumov-B-1_S1_L001_R1_001.fastq.gz	Abacumov-B-4_S4_L001_R1_001.fastq.gz
    Abacumov-B-1_S1_L001_R2_001.fastq.gz	Abacumov-B-4_S4_L001_R2_001.fastq.gz
    Abacumov-B-25_S25_L001_R1_001.fastq.gz	Abacumov-B-50_S50_L001_R1_001.fastq.gz
    Abacumov-B-25_S25_L001_R2_001.fastq.gz	Abacumov-B-50_S50_L001_R2_001.fastq.gz
    Abacumov-B-26_S26_L001_R1_001.fastq.gz	Abacumov-B-51_S51_L001_R1_001.fastq.gz
    Abacumov-B-26_S26_L001_R2_001.fastq.gz	Abacumov-B-51_S51_L001_R2_001.fastq.gz
    Abacumov-B-27_S27_L001_R1_001.fastq.gz	Abacumov-B-52_S52_L001_R1_001.fastq.gz
    Abacumov-B-27_S27_L001_R2_001.fastq.gz	Abacumov-B-52_S52_L001_R2_001.fastq.gz
    Abacumov-B-28_S28_L001_R1_001.fastq.gz	Abacumov-B-61_S61_L001_R1_001.fastq.gz
    Abacumov-B-28_S28_L001_R2_001.fastq.gz	Abacumov-B-61_S61_L001_R2_001.fastq.gz
    Abacumov-B-2_S2_L001_R1_001.fastq.gz	Abacumov-B-62_S62_L001_R1_001.fastq.gz
    Abacumov-B-2_S2_L001_R2_001.fastq.gz	Abacumov-B-62_S62_L001_R2_001.fastq.gz
    Abacumov-B-37_S37_L001_R1_001.fastq.gz	Abacumov-B-63_S63_L001_R1_001.fastq.gz
    Abacumov-B-37_S37_L001_R2_001.fastq.gz	Abacumov-B-63_S63_L001_R2_001.fastq.gz
    Abacumov-B-38_S38_L001_R1_001.fastq.gz	Abacumov-B-64_S64_L001_R1_001.fastq.gz
    Abacumov-B-38_S38_L001_R2_001.fastq.gz	Abacumov-B-64_S64_L001_R2_001.fastq.gz



```python
!qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data_met/sequences \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza
```

    [32mImported data_met/sequences as CasavaOneEightSingleLanePerSampleDirFmt to demux-paired-end.qza[0m
    [0m

## Deblur


```python
#фильтрация последовательностей на основе оценки качества 
!qiime quality-filter q-score \
 --i-demux demux-paired-end.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza
```

    [32mSaved SampleData[SequencesWithQuality] to: demux-filtered.qza[0m
    [32mSaved QualityFilterStats to: demux-filter-stats.qza[0m
    [0m


```python
#с помощью deblur объединяются похожие последовательности и удаляются ошибки секвенирования
!qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 120 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza
```

    [32mSaved FeatureTable[Frequency] to: table-deblur.qza[0m
    [32mSaved FeatureData[Sequence] to: rep-seqs-deblur.qza[0m
    [32mSaved DeblurStats to: deblur-stats.qza[0m
    [0m


```python
# табулирование статистики фильтрации для визуализации. 
# Результаты сохраняются в файле demux-filter-stats.qzv.
!qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv
!qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv
```

    [32mSaved Visualization to: demux-filter-stats.qzv[0m
    [0m[32mSaved Visualization to: deblur-stats.qzv[0m
    [0m

Оценить эффективность фильтрации и качество данных можно с помощь. визуализации полученной таблицы:


![deblur-stats.qzv](Screenshot_20230618_145223.png)


## Add metadata


```python
!sed 's/,/\t/g' map.csv > map.tsv #преобразование файла метаданных map.csv в формат TSV 
```


```python
!sed -i '1s/Filename/#SampleID/' map.tsv #замена"Filename" на "SampleID"
```


```python
# Создание сводной информации о таблице признаков
!qiime feature-table summarize \
  --i-table table-deblur.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file map.tsv

# Создание таблицы с информацией о последовательностях (включая идентификаторы и количество вхождений)
!qiime feature-table tabulate-seqs \
  --i-data rep-seqs-deblur.qza \
  --o-visualization rep-seqs.qzv

```

    [32mSaved Visualization to: table.qzv[0m
    [0m[32mSaved Visualization to: rep-seqs.qzv[0m
    [0m

![sample-frequencies](sample-frequenc.png)
![feature-frequencies](feature-frequenc.png)
![rep-seqs](rep-seq.png)



## Add tree


```python
# QIIME выравнивает последовательности и строит филогенетическое дерево: 
!qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-deblur.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

    [32mSaved FeatureData[AlignedSequence] to: aligned-rep-seqs.qza[0m
    [32mSaved FeatureData[AlignedSequence] to: masked-aligned-rep-seqs.qza[0m
    [32mSaved Phylogeny[Unrooted] to: unrooted-tree.qza[0m
    [32mSaved Phylogeny[Rooted] to: rooted-tree.qza[0m
    [0m

# Diversity Analysis

Анализ разнообразия микробиома на уровне альфа- и бета-разнообразия 



```python
#QIIME рассчитывает различные метрики альфа- и бета-разнообразия:
!qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-deblur.qza \
  --p-sampling-depth 8000 \
  --m-metadata-file map.tsv \
  --output-dir core-metrics-results
```

    [32mSaved FeatureTable[Frequency] to: core-metrics-results/rarefied_table.qza[0m
    [32mSaved SampleData[AlphaDiversity] to: core-metrics-results/faith_pd_vector.qza[0m
    [32mSaved SampleData[AlphaDiversity] to: core-metrics-results/observed_features_vector.qza[0m
    [32mSaved SampleData[AlphaDiversity] to: core-metrics-results/shannon_vector.qza[0m
    [32mSaved SampleData[AlphaDiversity] to: core-metrics-results/evenness_vector.qza[0m
    [32mSaved DistanceMatrix to: core-metrics-results/unweighted_unifrac_distance_matrix.qza[0m
    [32mSaved DistanceMatrix to: core-metrics-results/weighted_unifrac_distance_matrix.qza[0m
    [32mSaved DistanceMatrix to: core-metrics-results/jaccard_distance_matrix.qza[0m
    [32mSaved DistanceMatrix to: core-metrics-results/bray_curtis_distance_matrix.qza[0m
    [32mSaved PCoAResults to: core-metrics-results/unweighted_unifrac_pcoa_results.qza[0m
    [32mSaved PCoAResults to: core-metrics-results/weighted_unifrac_pcoa_results.qza[0m
    [32mSaved PCoAResults to: core-metrics-results/jaccard_pcoa_results.qza[0m
    [32mSaved PCoAResults to: core-metrics-results/bray_curtis_pcoa_results.qza[0m
    [32mSaved Visualization to: core-metrics-results/unweighted_unifrac_emperor.qzv[0m
    [32mSaved Visualization to: core-metrics-results/weighted_unifrac_emperor.qzv[0m
    [32mSaved Visualization to: core-metrics-results/jaccard_emperor.qzv[0m
    [32mSaved Visualization to: core-metrics-results/bray_curtis_emperor.qzv[0m
    [0m

Visualise alpha-diversity


```python
!qiime diversity alpha-rarefaction \
  --i-table table-deblur.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 20000 \
  --m-metadata-file map.tsv \
  --o-visualization alpha-rarefaction.qzv
```

    [32mSaved Visualization to: alpha-rarefaction.qzv[0m
    [0m

![alpha-rarefaction](alpha-rarefaction.png)
![bray-curtis](bray.png)


```python
#статистическая оценка значимости альфа-разнообразия между группами (на основе векторов альфа-разнообразия):
!qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file map.tsv \
  --o-visualization core-metrics-results/observed_features-significance.qzv

!qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file map.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

```

    [32mSaved Visualization to: core-metrics-results/observed_features-significance.qzv[0m
    [0m[32mSaved Visualization to: core-metrics-results/evenness-group-significance.qzv[0m
    [0m

![alpha-diversity-boxplots](alpha-diversity-boxplots.png)
![kruskal-wallis](kruskal-wallis.png)

Тест Kruskal-Wallis:
Нулевая гипотеза в этом тесте заключается в том, что распределения значений во всех группах равны.

На основании полученных результатов можно сделать вывод, что существует статистически значимая разница между группами. Значение p-value меньше уровня значимости 0.05, что позволяет отвергнуть нулевую гипотезу. Это означает, что распределения значений в разных группах статистически отличаются, и существует вероятность, что эти различия не случайны.

Между группами "Coal Mine Terricon" и "Embryo Sand" наблюдается статистически значимая разница (p-value = 0.020921), с учетом поправки на множественные сравнения (q-value = 0.024140).
Между группами "Litostrat" и "Embryo Sand" также наблюдается статистически значимая разница (p-value = 0.083265), но после коррекции на множественные сравнения эта разница становится менее значимой (q-value = 0.089212).
Между другими парами групп также наблюдаются значения статистики H, p-value и q-value, которые указывают на различия между группами, но точность этих различий может быть несколько ниже из-за коррекции на множественные сравнения.


```python
# Beta-diversity significance
!qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file map.tsv \
  --m-metadata-column Source\
  --o-visualization core-metrics-results/bray_curtis-body-site-significance.qzv \
  --p-pairwise
```

    [32mSaved Visualization to: core-metrics-results/bray_curtis-body-site-significance.qzv[0m
    [0m

![bray_curtis_body_site_significance1](bray_curtis_body_site_significance1.png)
![bray_curtis_body_site_significance2](bray_curtis_body_site_significance2.png)
![permanova_results](permanova_results.png)

Анализ PERMANOVA подтверждает, что существует статистически значимая разница между группами на основе анализа многомерных данных.
