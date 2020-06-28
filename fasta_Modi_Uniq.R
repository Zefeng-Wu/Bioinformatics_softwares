library(Biostrings)
library(Biostrings)
library(tidyverse)

### 1.lotus data (ftp://ftp.kazusa.or.jp/pub/lotus/lotus_r3.0/)

pep_fa_file <-"1orignal_pep/Lj3.0_pep.fna"
pep_fa <- readAAStringSet(pep_fa_file)
names(pep_fa) <- str_split_fixed(names(pep_fa),pattern = " ",n =2 )[,1]

pep_iso_length <-data.frame(pep_name = names(pep_fa),length=width(pep_fa))
pep_iso_length$gene_name <-str_split_fixed(pep_iso_length$pep_name,pattern = "\\..",n =2 )[,1]
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))%>%filter(!grepl('chloro|mito',pep_name))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 39640

writeXStringSet(out_pep_fa,"2modified_pep/Lotus_japonicus.fa",format = "fasta")

### 2. Cicer arietinum ICC 4958 genome v2 (https://www.pulsedb.org/analysis/133)
pep_fa_file <-"1orignal_pep/Ca_Pep_v2.fasta"
pep_fa <- readAAStringSet(pep_fa_file)
names(pep_fa) <- str_split_fixed(names(pep_fa),pattern = "\\|",n =2 )[,1]

pep_iso_length <-data.frame(pep_name = names(pep_fa),length=width(pep_fa))
pep_iso_length$gene_name <-str_split_fixed(pep_iso_length$pep_name,pattern = "\\..",n =2 )[,1]
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 22368
writeXStringSet(out_pep_fa,"2modified_pep/Cicer_arietinum.fa",format = "fasta")

### 3. Phaseolus_vulgaris (ensembl plants)
pep_fa_file <-"1orignal_pep/Phaseolus_vulgaris.PhaVulg1_0.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- gene_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$gene_name] # 28134 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Phaseolus_vulgaris.fa",format = "fasta")

### 4. Medicago (ensembl plants)
pep_fa_file <-"1orignal_pep/Medicago_truncatula.MedtrA17_4.0.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- gene_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$gene_name] # 28134 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Medicago_truncatula.fa",format = "fasta")

### 5. pigeonpea (Cajanus cajan).http://gigadb.org/dataset/100028

pep_fa_file <-"1orignal_pep/PigeonPea_V5.0.gene.pep"  # 48680
pep_fa <- readAAStringSet(pep_fa_file)
names(pep_fa) <- str_split_fixed(names(pep_fa),pattern = " ",n =2 )[,1]
writeXStringSet(pep_fa,"2modified_pep/Cajanus_cajan.fa",format = "fasta")

### 6. Gmax (ensembl plants).
pep_fa_file <-"1orignal_pep/Glycine_max.Glycine_max_v2.1.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- gene_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$gene_name] # 55897 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Glycine_max.fa",format = "fasta")

### 7. Vigna_angularis (enmsebl plants)
pep_fa_file <-"1orignal_pep/Vigna_angularis.Vigan1.1.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:LR48_")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- gene_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$gene_name] # 33860same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Vigna_angularis.fa",format = "fasta")

### 8. Vigna_radiata (enmsebl plants)
pep_fa_file <-"1orignal_pep/Vigna_radiata.Vradiata_ver6.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 22368 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Vigna_radiata.fa",format = "fasta")

### 9. Vitis_vinifera.12X.pep.all.fa(enmsebl plants)
pep_fa_file <-"1orignal_pep/Vitis_vinifera.12X.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 29927 not same with ensembl states!
writeXStringSet(out_pep_fa,"2modified_pep/Vitis_vinifera.fa",format = "fasta")

### 10. Fragaria_vesca (rosaceae.org/species/fragaria_vesca/genome_v4.0.a1,)  GDR database)
pep_fa_file <-"1orignal_pep/Fragaria_vesca_v4.0.a2.proteins.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = "\\..",n = 10)[,1]
protein_name <- names(pep_fa)

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 34006 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Fragaria_vesca.fa",format = "fasta")

### 11. Pyrus communis Bartlett DH (V2) (https://www.rosaceae.org/analysis/340 GDR database)
pep_fa_file <-"1orignal_pep/PyrusCommunis_BartlettDHv2.0.pep.fasta"
pep_fa <- readAAStringSet(pep_fa_file)
names(pep_fa) <- str_split_fixed(names(pep_fa),pattern = "\t",n = 10)[,1] # # 37455
writeXStringSet(pep_fa,"2modified_pep/Pyrus_communis.fa",format = "fasta")

### 12. https://www.rosaceae.org/species/prunus_persica/genome_v2.0.a1  GDR database)
pep_fa_file <-"1orignal_pep/Prunus_persica_v2.0.a1.primaryTrs.pep.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,3]
gene_name <-str_remove_all(gene_name,pattern = "locus=")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 26873 
writeXStringSet(out_pep_fa,"2modified_pep/Prunus_persica.fa",format = "fasta")


### 13. http://cucurbitgenomics.org/organism/20 (Cucumis sativus)
pep_fa_file <-"1orignal_pep/ChineseLong_pep_v3.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = "\\..",n = 10)[,1]
protein_name <- names(pep_fa)
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))

out_pep_fa <- pep_fa[pep_iso_length$pep_name] # 34006 same with ensembl states
writeXStringSet(out_pep_fa,"2modified_pep/Cucumis_sativus.fa",format = "fasta")


### 14 TAIR (enmsebl plants)
pep_fa_file <-"1orignal_pep/Arabidopsis_thaliana.TAIR10.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))%>%filter(!grepl('MG|CG',pep_name))

out_pep_fa <- pep_fa[pep_iso_length$pep_name]# 27420 not same with ensembl states!
writeXStringSet(out_pep_fa,"2modified_pep/Arabidopsis_thaliana.fa",format = "fasta")

## 15 Juglans regia_V2 https://www.hardwoodgenomics.org/Genome-assembly/2539069?tripal_pane=group_downloads
pep_fa_file <-"1orignal_pep/Juglans regia_V2.fasta"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = "_p",n = 10)[,1]

pep_iso_length <-data.frame(pep_name = names(pep_fa),length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))%>%filter(!grepl('chloroplast',pep_name))

out_pep_fa <- pep_fa[pep_iso_length$pep_name]# 27420 
writeXStringSet(out_pep_fa,"2modified_pep/Juglans_regia.fa",format = "fasta")

## 16 Lupinus_angustifolius (ensembl plants)
pep_fa_file <-"1orignal_pep/Lupinus_angustifolius.LupAngTanjil_v1.0.pep.all.fa"
pep_fa <- readAAStringSet(pep_fa_file)

gene_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,4]
gene_name <-str_remove_all(gene_name,pattern = "gene:")
protein_name <- str_split_fixed(names(pep_fa),pattern = " ",n = 10)[,1]
names(pep_fa) <- protein_name

pep_iso_length <-data.frame(pep_name = protein_name,length=width(pep_fa))
pep_iso_length$gene_name <-gene_name
pep_iso_length <- pep_iso_length%>%group_by(gene_name)%>%slice(which.max(length))%>%filter(!grepl('MG|CG',pep_name))

out_pep_fa <- pep_fa[pep_iso_length$pep_name]# 33074 same with ensembl states
names(out_pep_fa)<-pep_iso_length$gene_name
writeXStringSet(out_pep_fa,"2modified_pep/Lupinus_angustifolius.fa",format = "fasta")
