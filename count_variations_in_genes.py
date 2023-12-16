from pandas import read_csv
from vcf import Reader as vcf_reader
from pyranges import PyRanges as pyrange
from prettytable import PrettyTable

# wczytanie pliku ze współrzędnymi genów ...
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
)[
    "chr1"
]  # ... z jednoczesnym zawężeniem danych do analizowanej sekwencji

# inicjalizacja słownika liczników wariatnów dla poszczególnych genów
variations_in_genes_counter = dict()
for gen_name in genes.name:
    variations_in_genes_counter[gen_name] = 0

# sprawdzanie, w których genach występują kolejne warianty
with open("coriell_chr1.vcf", "r") as f:
    variations = vcf_reader(f)
    for variation in variations:
        plus_strand_sel = (
            (genes.Strand == "+")
            & (
                genes.Start <= (variation.POS - 1)
            )  # odejmujemy 1, bo indeksowanie w pliku vcf
            & (genes.End >= (variation.POS - 1))  # zaczyna się od 1, a w refFlat od 0
        )
        minus_strand_sel = (
            (genes.Strand == "-")
            & (genes.Start >= (variation.POS - 1))
            & (genes.End <= (variation.POS - 1))
        )
        genes_containing_variation = genes[plus_strand_sel | minus_strand_sel]
        if len(genes_containing_variation) > 0:  # zabezpieczenie przed "Empty PyRange"
            for gen in genes_containing_variation.name.unique():
                variations_in_genes_counter[gen] += 1

# utworzenie tabeli wyników
results = PrettyTable(["Gene", "Number of Variations"])
for key, val in variations_in_genes_counter.items():
    results.add_row([key, val])

# wyświetlenie wyników i zapis do pliku results_tab.txt
results_txt = results.get_string(sortby="Number of Variations", reversesort=True)
with open("results_tab.txt", "w") as f:
    f.write(results_txt)
print("Results (also saved into file results_tab.txt):\n{}".format(results_txt))
