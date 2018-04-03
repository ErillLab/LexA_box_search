from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez, AlignIO, Phylo
import json, os
from Bio.Align.Applications import ClustalwCommandline
from TF_clade_motif_discovery import writefasta


def createinput(filename):

    dictionary = {}

    with open(filename) as f:

        for k, v in json.load(f).items():

            dictionary[k] = v

    print "Creating input dictionary."

    return dictionary


input = createinput('test_input.json')


def lexa_phylogeny(dictionary, matrix):

    clustal_exe = r"/Users/anto/Desktop/msc bioinformatics/module V - practicum/clustalw-2.1-macosx/clustalw2"

    clustal_cline = ClustalwCommandline(clustal_exe, infile = dictionary["output_name"] + ".aln")

    stdout, stderr = clustal_cline()

    tree = Phylo.read(dictionary["output_name"] + ".dnd", "newick")

    external_nodes = tree.get_terminals()

    current_cluster, used_nodes, clusters, num_cluster = [], [], {}, 0

    for i in range(len(external_nodes)):

        if external_nodes[i] not in current_cluster and external_nodes[i] not in used_nodes:

            current_cluster = [external_nodes[i]]

            num_cluster += 1

            cluster_name = matrix[external_nodes[i].name]["order"]+ '_' + str(num_cluster)

            clusters[cluster_name] = current_cluster


        for k in range(i + 1, len(external_nodes)):

            if tree.distance(external_nodes[i], external_nodes[k]) < 0.5 and external_nodes[k] not in current_cluster:

                current_cluster.append(external_nodes[k])

                used_nodes.append(external_nodes[k])

                matrix[external_nodes[k]]["cluster"] = cluster_name

    return clusters



def clusterfiles(clusters, matrix):

    for cluster in clusters:

        cluster_file = open(cluster + '.fasta', 'w')

        for accession in clusters[cluster]:

            cluster_file.write('>' + accession + '\n' + matrix[accession]['sequence'] + '\n')



def motif_search(clusters):

    for cluster in clusters.keys():

        cmd = "meme {input} -o {out} -dna -w {width} -mod {mode}".format(input=cluster + '.fasta', out=input["restriction_taxon"],
                                                                         width=input["motif_width"],
                                                                         mode=input["motif_distribution"])

        print "Performing motif search."

        os.system(cmd)

