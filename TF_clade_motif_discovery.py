from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez, AlignIO, Phylo
import json
from Bio.Align.Applications import ClustalwCommandline
from datetime import datetime


#Dictionary with the input parameters

def createinput(filename):

    dictionary = {}

    with open(filename) as f:

        for k, v in json.load(f).items():

            dictionary[k] = v

    print "Creating input dictionary."

    return dictionary


input = createinput('test_input.json')


#Create one dictionary for each type of sequences


def referencesequences(filename):

    dictionary = {}

    for sequence in SeqIO.parse(filename, 'fasta'):

        dictionary[sequence.id.split('|')[3]] = sequence.seq._data

    return dictionary


references = {'LexA' : referencesequences('ReferenceDataset_LexA.fasta'), 'UmuD' : referencesequences('ReferenceDataset_UmuD.fasta'),
              'CI' : referencesequences('ReferenceDataset_lambdaCI.fasta'), 'Similar' : referencesequences('ReferenceDataset_Similar.fasta')}


matrix = {}

def blastsearch(matrix, dictionary):

    for order in open(dictionary["bacteria_orders"], 'r'):

        startTime = datetime.now()

        print "Blast search query."

        print order

        taxon = "txid" + str(order.strip()) + "[orgn]"

        ids = open(dictionary["protein_ids"]).read()

        handleresults = NCBIWWW.qblast(program='blastp', database='refseq_protein',
                                       sequence=ids, entrez_query=taxon,
                                       expect=dictionary["e_value_threshold"], hitlist_size=2000)

        blast_records = list(NCBIXML.parse(handleresults))

        print datetime.now() - startTime

        startTime = datetime.now()

        for record in blast_records:

            record_id = record.query.split(' ')[0].split('|')[3]

            i = 2000

            if len(record.descriptions) == 0: continue

            while record.descriptions[-1].e < int(dictionary["e_value_threshold"]):

                print 'Repeat BLAST search with ' + record_id

                i += 500

                handlecorrection = NCBIWWW.qblast(program='blastp', database='refseq_protein',
                                               sequence=record_id, entrez_query=taxon,
                                               expect=int(dictionary["e_value_threshold"]), hitlist_size=i)

                blast_records_corrected = list(NCBIXML.parse(handlecorrection))

                record = blast_records_corrected[0]

            print 'Length is' + str(len(record.alignments))

            for hit in record.alignments:

                cov = float(hit.length) / record.query_length

                if cov > 0.7:

                    if hit.accession not in matrix.keys():

                        matrix[hit.accession] = {}

                    matrix[hit.accession][record_id] = {}

                    matrix[hit.accession][record_id]["evalue"], matrix[hit.accession][record_id]["coverage"], matrix[hit.accession][record_id]["order"] = hit.hsps[0].expect, cov, taxon


                    if record_id in references['LexA'].keys(): matrix[hit.accession][record_id]["source"] = "LexA"

                    elif record_id in references['UmuD'].keys(): matrix[hit.accession][record_id]["source"] = "UmuD"

                    elif record_id in references['CI'].keys(): matrix[hit.accession][record_id]["source"] = "repressor CI"

                    elif record_id in references['Similar'].keys() : matrix[hit.accession][record_id]["source"] = "Similar"

        print datetime.now() - startTime

    #Filter the matrix, excluding those values whose lowest E-value comes from a contamination query

    print len(matrix)

    matrix = {key: value for key, value in matrix.items() if sorted(value.values(), key=lambda k: k['evalue'])[0]['source'] == "LexA"}

    return matrix



matrix = blastsearch(matrix, input)



def entrezsearch(matrix):

    for hit in matrix.keys():

        protein_seq = SeqIO.read(Entrez.efetch(db="protein", id=hit, rettype="fasta"), 'fasta')

        matrix[hit]['sequence'] = protein_seq.seq._data


        records = Entrez.read(Entrez.efetch(db="protein", id=hit, rettype='ipg', retmode='xml'))


        for record in records[0]['ProteinList']:

            if record.attributes['source'] == "RefSeq":

                try:

                    accgenome = record['CDSList'][0].attributes['accver']

                    if record['CDSList'][0].attributes['strand'] == '+':

                        position = record['CDSList'][0].attributes['start']
                        strand = 1

                    else:

                        position = record['CDSList'][0].attributes['stop']
                        strand = 2

                except KeyError:

                    pass


        if strand == 1:

            flanking_handle = Entrez.efetch(db="nucleotide", id = accgenome, rettype = "fasta", strand = strand,
                                            seq_start = int(position) - int(input["upstream_region"]),
                                            seq_stop = int(position) + int(input["downstream_region"]))

        else:

            flanking_handle = Entrez.efetch(db="nucleotide", id=accgenome, rettype="fasta", strand=strand,
                                            seq_start=int(position) - int(input["downstream_region"]),
                                            seq_stop=int(position) + int(input["upstream_region"]))


        flanking_parse = SeqIO.read(flanking_handle, "fasta")

        matrix[hit]['flanking_region'] = flanking_parse.seq._data

    return matrix


matrix = entrezsearch(matrix)



def writefasta(matrix):

    lexa_fasta = open(input["output_name"] + '.fasta', 'w')

    #Add the references LexA sequences to the FASTA file

    for sequence in references['LexA']:

        lexa_fasta.write('>' + sequence + '\n' + references['LexA'][sequence] + '\n')

    #Add all the hits in the matrix to the FASTA file

    for hit in matrix.keys():

        lexa_fasta.write('>' + hit + '\n' + matrix[hit]['sequence'] + '\n')

    lexa_fasta.close()


#Run the function

writefasta(input, matrix)


def first_phylogeny(dictionary):

    clustal_exe = r"/Users/anto/Desktop/msc bioinformatics/module V - practicum/clustalw-2.1-macosx/clustalw2"

    clustal_cline = ClustalwCommandline(clustal_exe, infile = dictionary["output_name"] + ".fasta")

    stdout, stderr = clustal_cline()


first_phylogeny(input)
