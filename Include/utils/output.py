from sys import path

output_path = fr"{path[5]}\Include\data\output"


def all_gene_names(test: str, gene_names: str, manual: bool = False):
    finaloutput_path = output_path if not manual else f'{output_path}\manual'
    with open(fr'{finaloutput_path}\{test}_all_names.txt', 'w') as f:
        f.write(f'There are {len(gene_names)} genes\n')
        f.write(f"These are all {test}'s gene names:\n")
        for gene in gene_names:
            f.write(f'{gene}\n')


