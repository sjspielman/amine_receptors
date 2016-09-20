from Bio import SeqIO


def grabseq(id, refdict):
    for entry in refdict:
        if entry == id:
            sequence = str(refdict[entry].seq)
            break
    return sequence

mapfile = "sequence_descriptions.txt"
nucfile = "nucleotide_records.fasta"
protfile = "protein_records.fasta"


nucseqs = SeqIO.to_dict(SeqIO.parse(nucfile, "fasta"))
protseqs = SeqIO.to_dict(SeqIO.parse(protfile, "fasta"))

with open(mapfile, "r") as f:
    maplines = f.readlines()

chrm_nuc = {}
chrm_prot = {}

for line in maplines[1:]:
    parsedline = line.strip().split("\t")
    id = parsedline[0]
    gene = parsedline[1]
    if "CHRM" in gene:
        newid=id + "_" + gene
        print newid
        chrm_nuc[newid] = grabseq(id, nucseqs)
        chrm_prot[newid] = grabseq(id, protseqs)

with open("CHRM_nuc_raw.fasta", "w") as f:
    for entry in chrm_nuc:
        f.write(">" + entry + "\n" + chrm_nuc[entry] + "\n")
with open("CHRM_prot_raw.fasta", "w") as f:
    for entry in chrm_prot:
        f.write(">" + entry + "\n" + chrm_prot[entry] + "\n")
