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
    1. Deeptools （三大功能：1. BAM & Bigwig格式文件处理；2. QC检测 3. 热图和metaplot）
    # 处理器数目设定
     -p max/2
     
    # 针对指定区域进行处理
     --region chr2:10000-20000

    # ignoreDuplicates参数去除重复序列，针对匹配到同一方向同一起点的序列，只保留一个
     -- ignoreDuplicates

    # 匹配得分阈值设定
     --minMappingQuality

    # warning，deeptools是在scaling data做低质量数据去除和去重，所以如果数据质量较差及重复数据很多，尽量事先使用samtools进行提前处理
    
    功能一：BAM & bigwig file processing
        multiBamSummary
        multiBigwigSummary
        correctGCbias
        bamCoverage
        bamCompare
        bigwigCompare
        computeMatrix
    1. multiBamSummary：可以用来处理bam文件在基因组上覆盖情况，默认输出npz文件，衔接plotCorrelation和plotPCA进行作图
    bin mode
        multiBamSummary bins --bamfiles file1.bam file2.bam -out results.npz
    
    BED-file mode
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

    2.bamCoverage:可以用来将bam file转换成bigwig file，同时可以设定binSize参数从而的获取不同的分辨率，在比较非一批数据的时候，还可以设定数据normalizeTo1X到某个值（一般是该物种基因组大小）从而方便进行比较。
    bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2150570000
    --ignoreForNormalization chrX
    --extendReads
    
    3.bamCompare：可以用来的处理treat组和control组的数据转换成bigwig文件，给出一个binsize内结合强度的比值（默认log2处理）。
    bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw --normalizeTo1x 2451960000
    
    4.computeMatrix：该功能可以计算每个基因区域的结合得分，生成中间文件用以给plotHeatmap和plotProfiles作图。computeMatrix有两种模式，scale-regions mode和reference-point mode

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

    功能二：Tools for QC. 包括PCA作图，correlation作图等，都是运用multiBamSummary得到npz文件统计样本间的相关系数作图和PCA分析作图，没有需求故此处不做介绍。
        plotCorrelation
        plotPCA
        plotFingerprint
        bamPEFragmentSize
        computeGCBias
        plotCoverage

    功能三：Heatmaps and summary plots:主要用来画热图（并包含聚类功能);上游数据是computeMatrix得到的gz file
        plotHeatmap
        plotProfile
        plotEnrichment
        plotHeatmap
    
    1. plotHeatmap
    plotHeatmap -m matrix_two_groups.gz \ #输入gz文件
     -out ExampleHeatmap2.png \ 
     --colorMap RdBu \ #指定颜色
     --whatToShow 'heatmap and colorbar' \ #指定输出geatmap和colorbar
     --zMin -3 --zMax 3 \ #指定colorbar的范围
     --kmeans 4 #设定聚类个数

    2.plotProfile:主要用来画密度图,上游数据是computeMatrix得到的gz file;注意：默认针对单个bw文件作图或者把多个bw文件画在一个图里面（perGroup参数），同样也可以使用kmean或hclust聚类
      plotProfile -m matrix.mat.gz \
      --perGroup \
      --kmeans 2 \
      -out ExampleProfile3.png
    
    其他参数
        -z 给bed文件一个名称
        --samplesLabel  给bw文件一个名称
        --startLabel
        --endLabel
    
## 三代测序数据
    1reads过滤
    NanoFilt -q 9 -l 1000 > filter.fq (有问题) 
    filtlong --min_length 1000 --min_mean_q 9 SRR6924617.fastq >SRR6924617_filt_long_filter.fastq  #（正确）

## ATAC
    ATACseqQC Guide

## <font color="#006600">eQTL</font><br/>
    R MatrixEQTL
    Exploring regulation in tissues with eQTL networks
    
## Python 环境之anaconda. (尤其适用于在服务器上没有sduo权限的时候安装python包)
    1.下载并装miniconda
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod 777 Miniconda3-latest-Linux-x86_64.sh #给执行权限
    bash Miniconda3-latest-Linux-x86_64.sh #运行
    rm -rf ~/miniconda #卸载minicoda
    
    2.启动/退出
    cd miniconda3
    chmod 777 activate 
    . ./activate #这里的第一个点跟source是一样的效果
    conda list
    . ./deactivate #退出环境
    
    3.添加channel
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --get channels
    
    4.安装/卸载软件
    conda install gatk
    conda install gatk=3.7 #指定版本
    conda search gatk #搜索软件
    conda remove gatk #卸载软件
    conda update 软件名 # 更新软件
    
    5.创建conda环境
    conda env list #查当前环境
    conda create -n python2 python=2 #-n: 设置新的环境的名字, python=2 指定新环境的python的版本
    
## samtools
    
    samtools view  test.sorted.bam chr1:10000-20000 #提取特定区间的比对信息
    samtools view -c SAMPLE.bam ## get the total number of reads of a BAM file (may include unmapped and duplicated multi-aligned reads)
    samtools view -c -F 260 SAMPLE.bam #    ### counting only mapped (primary aligned) reads
    
    options
        -c  count reads and print the total number
        -f bitcode  output reads that fulfill the checked 'bitcode' criteria, see SAM bitcode fields
        -F bitcode  exclude reads that match one or more checked 'bitcode' criteria, see SAM bitcode fields
        -F 260  output primary aligned mapped reads
                       read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion
                       bit 3 + bit 9 = 4 + 256 = 260
                       
    samtools coverage test.sorted.bam  # 以染色体为单位统计覆盖度、测序深度、碱基数等
    samtools bedcov   test.sorted.bam  # 以bed文件为单位统计测序覆盖度（非reads数）
    samtools depth    test.sorted.bam  # 统计每个碱基的测序深度

  
 

    
## bedtools
    bedtools merge -i test.sorted.bam          # 合并overlappedd的reads,形成bed文件   
    bedtools coverage -a test.bed -b test.bam  #每个bed区间的reads数和总碱基数 
    
## 直系同源基因的鉴定
    https://orthovenn2.bioinfotoolkits.net/home


## 本地git新项目
    mkdir cv && cd
    git init 
    git config --global user.name "Zefeng-Wu"
    git config --global user.email "835102330@qq.com"
    
    ### download or self-make
    wget http://labfile.oss.aliyuncs.com/courses/624/cv-template.zip
    unzip cv-template
    mv cv-template/* .
    unzip cv-template
    mv cv-template/* 
    rm -rf cv-template* __MACOSX*（MACOSX前面是两根下划线）
    
    ## make a repository via net browser, member the git address
    git add .
    git commit -m 'commit my cv'
    git remote add origin https://github.com/Zefeng-Wu/cv.git
    git push -u origin master
    
## Cytoscape
    导入网络文件：file->import->network->file(net.txt)
    导入节点属性：file->import->table->file(node.txt)(此处为table而非network)
    

## Othofinder

    fasta_Modi_Uniq.R # modify the fasta header and isoforms
    for m in $(ls *.fa); do orthomclAdjustFasta $(basename ${m%.fa}) $m 1 ; done  # add species  in the header
    mv *.fasta > ../3adjust_fasta/

    conda create -n orthofinder
    conda activate orthofinder
    conda install -c bioconda orthofinder
    nohup orthofinder -f 3adjust_fasta/ -t 40 -M msa & 
