from pandas import read_csv
from vcf import Reader as vcf_reader
from pyranges import PyRanges as pyrange
from prettytable import PrettyTable

genes = pyrange(
    read_csv(
        "refFlat.txt",
        sep="\t",
        names=[
            "geneName",
            "name",
            "Chromosome",
            "Strand",
            "Start",
            "End",
            "cdsStart",
            "cdsEnd",
            "exonCount",
            "exonStarts",
            "exonEnds",
        ],
    )
)["chr1"]

variations_in_genes_counter = dict()
for gen_name in genes.name:
    variations_in_genes_counter[gen_name] = 0

with open("coriell_chr1.vcf", "r") as f:
    variations = vcf_reader(f)
    for variation in variations:
        plus_strand_sel = (
            (genes.Strand == "+")
            & (genes.Start <= (variation.POS-1))
            & (genes.End >= (variation.POS-1))
        )
        minus_strand_sel = (
            (genes.Strand == "-")
            & (genes.Start >= (variation.POS-1))
            & (genes.End <= (variation.POS-1))
        )
        genes_containing_variation = genes[plus_strand_sel | minus_strand_sel]
        if  len(genes_containing_variation) > 0:
            for gen in genes_containing_variation.name.unique():
                variations_in_genes_counter[gen] += 1

results = PrettyTable(["Gene", "Number of Variations"])
for key, val in variations_in_genes_counter.items():
    results.add_row([key, val])

results_txt = results.get_string(sortby="Number of Variations", reversesort=True)
with open("results_tab.txt", "w") as f:
    f.write(results_txt)
print("Results (also saved into file results_tab.txt):\n{}".format(results_txt))