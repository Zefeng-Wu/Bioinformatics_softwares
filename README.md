## 串联重复元件
    HipSTR: Genome-wide profiling of heritable and de novo STR variations, 2017, Nature Methods. (https://hipstr-tool.github.io/HipSTR)
citings: The impact of short tandem repeat variation on gene expression, 2019, Nature genetics 

## 转录组进化
    treeExp：http://hupi.fudan.edu.cn/people/guxun
## 基因融合 
### 1. star-fusion （在速度和准确率上均比较高，基于star比对结果）
#### 1. 构建reference lib. 人类和鼠的可以从以下网址直接下载：https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/; 其中plug-n是已经建立好的reference lib, 而source里面包含了所需的原始文件。从原始文件构建reference lib的命令如下；默认会在当前目录生成一个名为    ctat_genome_lib_build_dir的目录
    FusionFilter/prep_genome_lib.pl \                   #集成于Star-Fucsion里的perl脚本
    --genome_fa ref_genome.fa \                         # 参考基因组         
    --gtf ref_annot.gtf \                               # 基因注释文件
    --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \     # fusion 注释数据库（不是必须）
    --annot_filter_rule AnnotFilterRule.pm \
    --pfam_db PFAM.domtblout.dat.gz
#### 2. 运行STAR-fusion. STAR-fusion支持两种模式，第一种是直接从fastq开始，第二种是自己手动进行STAR比对，然后在运行STAR-fusion。
#### 第一种模式的用法如下：
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
#### 第二种模式
#### 1. Star比对
    STAR --genomeDir ${star_index_dir} \   
      --readFilesIn ${left_fq_filename} ${right_fq_filename} \ 
      --twopassMode Basic \      
      --outReadsUnmapped None \                                                                                                  
      --chimSegmentMin 12 \                                                                                                    
      --chimJunctionOverhangMin 12 \                                                                                           
      --alignSJDBoverhangMin 10 \                                                                                              
      --alignMatesGapMax 100000 \ 
      --alignIntronMax 100000 \   
      --chimSegmentReadGapMax 3 \                                                                                    
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --runThreadN ${THREAD_COUNT} \ 
      --outSAMstrandField intronMotif \
      --chimOutJunctionFormat 1
#### 2. Star-funsion
    STAR-Fusion \
      --genome_lib_dir CTAT_resource_lib \
      -J Chimeric.out.junction \
      --output_dir star_fusion_outdir

## BAM 文件统计和可视化
### 1. Deeptools （三大功能：1. BAM & Bigwig格式文件处理；2. QC检测 3. 热图和metaplot）
    # 处理器数目设定
     -p max/2
     
    # 针对指定区域进行处理
     --region chr2:10000-20000

    # ignoreDuplicates参数去除重复序列，针对匹配到同一方向同一起点的序列，只保留一个
     -- ignoreDuplicates

    # 匹配得分阈值设定
     --minMappingQuality

    # warning，deeptools是在scaling data做低质量数据去除和去重，所以如果数据质量较差及重复数据很多，尽量事先使用samtools进行提前处理
#### 功能一：BAM & bigwig file processing
    multiBamSummary
    multiBigwigSummary
    correctGCbias
    bamCoverage
    bamCompare
    bigwigCompare
    computeMatrix
##### 1. multiBamSummary：可以用来处理bam文件在基因组上覆盖情况，默认输出npz文件，衔接plotCorrelation和plotPCA进行作图
    # bin mode
    multiBamSummary bins --bamfiles file1.bam file2.bam -out results.npz
    
    # BED-file mode
    multiBamSummary BED-file --BED selection.bed --bamfiles file1.bam file2.bam -out results.npz

    deepTools2.0/bin/multiBamSummary bins \
      --bamfiles testFiles/*bam \ # using all BAM files in the folder
      --minMappingQuality 30 \
      --region 19 \ # limiting the binning of the genome to chromosome 19
      --labels H3K27me3 H3K4me1 H3K4me3 HeK9me3 input \
      -out readCounts.npz --outRawCounts readCounts.tab

     head readCounts.tab
     'chr'   'start' 'end'   'H3K27me3'      'H3K4me1'       'H3K4me3'       'HeK9me3'       'input'
     19 10000   20000   0.0     0.0     0.0     0.0     0.0
     19 20000   30000   0.0     0.0     0.0     0.0     0.0
     19 30000   40000   0.0     0.0     0.0     0.0     0.0
     19 40000   50000   0.0     0.0     0.0     0.0     0.0
     19 50000   60000   0.0     0.0     0.0     0.0     0.0
     19 60000   70000   1.0     1.0     0.0     0.0     1.0
     19 70000   80000   0.0     1.0     7.0     0.0     1.0
     19 80000   90000   15.0    0.0     0.0     6.0     4.0
     19 90000   100000  73.0    7.0     4.0     16.0    5.0

##### 2.bamCoverage:可以用来将bam file转换成bigwig file，同时可以设定binSize参数从而的获取不同的分辨率，在比较非一批数据的时候，还可以设定数据normalizeTo1X到某个值（一般是该物种基因组大小）从而方便进行比较。
    bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
      --binSize 10
      --normalizeTo1x 2150570000
      --ignoreForNormalization chrX
      --extendReads
##### 3.bamCompare：可以用来的处理treat组和control组的数据转换成bigwig文件，给出一个binsize内结合强度的比值（默认log2处理）。
    bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw --normalizeTo1x 2451960000
##### 4.computeMatrix：该功能可以计算每个基因区域的结合得分，生成中间文件用以给plotHeatmap和plotProfiles作图。computeMatrix有两种模式，scale-regions mode和reference-point mode

    computeMatrix scale-regions -p 10 \
      -R gene19.bed geneX.bed \
      -S test1.bw test2.bw \
      -b 3000 -a 3000 \
      --regionBodyLength 5000 \   
      --skipZeros \
      -o heatmap.gz 
reference-point mode则是给定一个bed file，以某个点为中心开始统计信号（TSS/TES/center）。但实际上我在尝试的时候regionBdoyLength参数也还是可以用的，所以估计和scale-regions区别也不是太大，主要是作图的一点区别。
    computeMatrix reference-point \ # choose the mode
       --referencePoint TSS \ # alternatives: TES, center
       -b 3000 -a 10000 \ # define the region you are interested in
       -R testFiles/genes.bed \
       -S testFiles/log2ratio_H3K4Me3_chr19.bw  \
       --skipZeros \
       -o matrix1_H3K4me3_l2r_TSS.gz \ # to be used with plotHeatmap and plotProfile
       --outFileSortedRegions regions1_H3K4me3_l2r_genes.bed
