## 串联重复元件
    HipSTR: Genome-wide profiling of heritable and de novo STR variations, 2017, Nature Methods. (https://hipstr-tool.github.io/HipSTR)
citings: The impact of short tandem repeat variation on gene expression, 2019, Nature genetics 

## 转录组进化
    treeExp：http://hupi.fudan.edu.cn/people/guxun
## 基因融合 
### 1. star-fusion （在速度和准确率上均比较高，但需要stat比对结果）
#### 1. 构建reference lib （人类和鼠的可以直接下载：https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/； 其中plug-n是已经建立好的reference lib, 而source里面包含了所需的原始文件。从原始文件构建reference lib的命令如下）
    FusionFilter/prep_genome_lib.pl \
    --genome_fa ref_genome.fa \
    --gtf ref_annot.gtf \
    --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \
    --annot_filter_rule AnnotFilterRule.pm \
    --pfam_db PFAM.domtblout.dat.gz



    STAR-Fusion \
    --genome_lib_dir CTAT_resource_lib \
    --left_fq reads_1.fq \
    --right_fq reads_2.fq \
    --output_dir star_fusion_outdir


