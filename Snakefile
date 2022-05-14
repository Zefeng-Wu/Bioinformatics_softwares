'''
This script is used to merge files
'''
rule all: 			# 如果流程里有2个或2个以上的rule，那么就一定要写rule all
	input: "ln.txt"  

rule merge:
	input:
		expand("{file}.txt",file=["hello","world"])
	output:
		"merged.txt"
	threads:
		2
	shell:
		"cat {input} > {output}"

rule countline:
	input:
		rules.merge.output #dependence
	output:
		"ln.txt"
	shell: 
		"wc -l {input} > {output}"


