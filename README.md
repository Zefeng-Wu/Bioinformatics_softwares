## 串联重复元件
    HipSTR: Genome-wide profiling of heritable and de novo STR variations, 2017, Nature Methods. (https://hipstr-tool.github.io/HipSTR)
citings: The impact of short tandem repeat variation on gene expression, 2019, Nature genetics 

## 转录组进化
    treeExp：http://hupi.fudan.edu.cn/people/guxun
## 基因融合 
### 1. star-fusion （在速度和准确率上均比较高，基于star比对结果）
#### 1. **构建reference lib**. 人类和鼠的可以从以下网址直接下载：https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/; 其中plug-n是已经建立好的reference lib, 而source里面包含了所需的原始文件。从原始文件构建reference lib的命令如下；默认会在当前目录生成一个名为    ctat_genome_lib_build_dir的目录

    FusionFilter/prep_genome_lib.pl \                   #集成于Star-Fucsion里的perl脚本
    --genome_fa ref_genome.fa \                         # 参考基因组         
    --gtf ref_annot.gtf \                               # 基因注释文件
    --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \     # fusion 注释数据库（不是必须）
    --annot_filter_rule AnnotFilterRule.pm \
    --pfam_db PFAM.domtblout.dat.gz
#### 2. **运行STAR-fusion：** STAR-fusion支持两种模式，第一种是直接从fastq开始，第二种是自己手动进行STAR比对，然后在运行STAR-fusion。第一种模式的用法如下：
##### 双端
    STAR-Fusion \
    --genome_lib_dir CTAT_resource_lib \
    --left_fq reads_1.fq \
    --right_fq reads_2.fq \
    --output_dir star_fusion_outdir
##### 单端
    STAR-Fusion \
    --genome_lib_dir CTAT_resource_lib \
    --left_fq reads_1.fq \
    --output_dir star_fusion_outdir


