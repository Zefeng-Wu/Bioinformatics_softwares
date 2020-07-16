from Bio import AlignIO

fa_align = AlignIO.read("test.fa",format="fasta")
print(fa_aling)
AlingIO.write(fa_align,"out.phy",format="phylip-relaxed")

