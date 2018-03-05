import json, os, sys
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO, pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

Entrez.email = "antoniocolladopadilla12@gmail.com"

input_file = sys.argv[1]


# Input: Fasta file + Json with parameteres

def createinput(filename):
    dict = {}
    with open(filename) as f:
        for k, v in json.load(f).items():
            dict[k] = v
    print "Creating input dictionary."
    return dict


input = createinput(input_file)

# Blast search

results = []


def blastsearch():
    resultslist = []

    print "Performing Blast search"

    ids = open(input["protein_ids"]).read()
    handleresults = NCBIWWW.qblast(program='blastp', database='refseq_protein', sequence=ids,
                                   entrez_query=input["restriction_taxon"],
                                   expect=input["e_value_threshold"],
                                   hitlist_size=30)
    blast_records = NCBIXML.parse(handleresults)

    blast_records = list(blast_records)
    # each blast Lexa has several matches: compile

    count = input["e_value_threshold"]

    print "Choosing best query."

    for query in blast_records:
        scoring = 0
        for record in query.descriptions:
            scoring += record.e
        if scoring < count:
            count = scoring
            best_query = query

    print best_query.query

    for alignment in best_query.alignments:
        resultslist.append(alignment.title.split()[0].split('|')[3])

    return resultslist


# Entrez search



def entrezsearch():
    acc_list, position_list, strand_list, flanking_regions = [], [], [], []

    results = blastsearch()

    print "Retrieving Entrez queries."

    for acc in results:  # list with matches from the selected LexA reference

        handle = Entrez.efetch(db="protein", id=acc, rettype='ipg',
                               retmode='xml')  # for each match there are a lot of CDS
        records = Entrez.read(handle)  # list of the seqs of identical proteins of the match

        for record in records[0]['ProteinList']:

            if record.attributes['source'] == "RefSeq":

                try:

                    acc_list.append(record['CDSList'][0].attributes['accver'])

                    if record['CDSList'][0].attributes['strand'] == '+':

                        position_list.append(record['CDSList'][0].attributes['start'])
                        strand_list.append(1)

                    else:

                        position_list.append(record['CDSList'][0].attributes['stop'])
                        strand_list.append(2)

                except KeyError:

                    pass

    print "Retrieving flanking sequences."

    for i in range(len(acc_list)):

        if strand_list[i] == 1:
            flanking_handle = Entrez.efetch(db="nucleotide", id=acc_list[i], rettype="fasta", strand=strand_list[i],
                                            seq_start=int(position_list[i]) - int(input["downstream_region"]),
                                            seq_stop=int(position_list[i]) + int(input["upstream_region"]))

        else:
            flanking_handle = Entrez.efetch(db="nucleotide", id=acc_list[i], rettype="fasta", strand=strand_list[i],
                                            seq_start=int(position_list[i]) - int(input["upstream_region"]),
                                            seq_stop=int(position_list[i]) + int(input["downstream_region"]))

        flanking_parse = SeqIO.read(flanking_handle, "fasta")
        flanking_regions.append(flanking_parse.seq._data)

    return flanking_regions, acc_list


pass_region, pass_acc = [], []


def nucleotide_identity():
    flanking_regions, acc_list = entrezsearch()

    print "Removing similar sequences."

    for i in range(len(flanking_regions)):

        for k in range(i + 1, len(flanking_regions)):

            alignment = pairwise2.align.globalds(flanking_regions[i], flanking_regions[k], blosum62, -10, -1)
            seq1, seq2 = alignment[0][0], alignment[0][1]
            matches = ""

            for base1, base2 in zip(seq1, seq2):
                if base1 == base2 and base1 != '-':
                    matches += base1

            identity = float(len(matches)) / len(seq1) * 100

            if identity > input["identity_threshold"]:
                flanking_regions.pop(i)
                acc_list.pop(i)
                break

    output_fasta = 'filtered_flanking_regions.fasta'

    output_write = open(output_fasta, 'w')

    for i in range(len(flanking_regions)):
        output_write.write('>' + acc_list[i] + "\n" + flanking_regions[i] + "\n")

    output_write.close()

    return output_fasta


def motif_search():

    output_fasta = nucleotide_identity()

    cmd = "meme {input} -o {out} -dna -w {width} -mod {mode}".format(input=output_fasta, out=input["restriction_taxon"],
                                                                     width=input["motif_width"],
                                                                     mode=input["motif_distribution"])

    print "Performing motif search."

    os.system(cmd)


motif_search()
