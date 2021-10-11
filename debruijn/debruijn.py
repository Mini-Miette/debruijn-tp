#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import itertools as it
import more_itertools as mit
import matplotlib.pyplot as plt
from operator import itemgetter
import pickle
import networkx as nx
import sys
import statistics as st
import time
import os
import argparse
import random
from random import randint
random.seed(9001)


__author__ = "Laura Xénard"
__copyright__ = "Universite de Paris"
__credits__ = ["Laura Xénard"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Laura Xénard"
__email__ = "laura.xenard@protonmail.com"
__status__ = "Developpement"


parse = False


def isfile(path):
    """Check if path is an existing file.


    Args:
        path (str): Path to the file.

    Raises:
        argparse: ArgumentTypeError

    Returns:
        path (str): Path to the file

    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieve the arguments of the program.

    Returns:
        Object: An object that contains the arguments

    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Read a fasta file.


    Args:
        fastq_file (str): Path of the fasta file

    Yields:
        str generator: DESCRIPTION.

    """

    with open(fastq_file, 'r') as f:
        for line in f:
            # We want to get 1 line out of 4.
            line = next(f)
            yield line.strip()
            # We need to strip to remove the trailing \n.
            line = next(f)
            line = next(f)


def cut_kmer(read, kmer_size):

    for i in range(len(read)-kmer_size+1):
        yield read[i:(i+kmer_size)]


def build_kmer_dict(fastq_file, kmer_size):
    """
    Retourne un dictionnaire ayant pour clé les k-mers et pour valeur le nombre
    d’occurrence de ce k-mer.
    """

    fastq_reader = read_fastq(fastq_file)
    kmer_dict = dict()
    for read in fastq_reader:
        kmer_reader = cut_kmer(read, kmer_size)
        for kmer in kmer_reader:
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict):

    graph = nx.DiGraph()
    for key, val in kmer_dict.items():
        n1 = key[:-1]  # Prefix node.
        n2 = key[1:]  # Suffix node.
        graph.add_edge(n1, n2, weight=val)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):

    for path in path_list:
        graph.remove_nodes_from(path[1:-1])
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
    return graph


def std(data):
    if len(data) == 1:
        data_std = 0
    else:
        data_std = st.stdev(data)
    return data_std


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):

    # Finding the best path
    weight_std = std(weight_avg_list)
    check_length = False
    if weight_std > 0:
        max_weight = max(weight_avg_list)
        paths_max_weight = [i for i, j in enumerate(weight_avg_list)
                            if j == max_weight]
        if len(paths_max_weight) == 1:
            best_path = paths_max_weight[0]
        else:
            check_length = True
    if weight_std == 0 or check_length:
        length_std = std(path_length)
        if length_std > 0:
            max_length = max(path_length)
            paths_max_length = [i for i, j in enumerate(path_length)
                                if j == max_length]
            if len(paths_max_length) == 1:
                best_path = paths_max_length[0]
            else:
                best_path = randint(0, len(path_length) - 1)
        else:
             best_path = randint(0, len(path_length) - 1)

    # Deleting the other paths
    path_list.pop(best_path)
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph, path):

    path_weight = 0
    for n1, n2 in mit.pairwise(path):
        path_weight += graph.edges[n1, n2]['weight']
    return path_weight / (len(path) - 1)


def solve_bubble(graph, ancestor_node, descendant_node):

    simple_path_reader = nx.all_simple_paths(graph, ancestor_node,
                                             descendant_node)
    path_list = list(simple_path_reader)
    if len(path_list) > 1:
        path_length = [len(path) for path in path_list]
        weight_avg_list = [path_average_weight(graph, path)
                           for path in path_list]
        graph = select_best_path(graph, path_list, path_length,
                                 weight_avg_list)
    return graph


def simplify_bubbles(graph):

    bubble = False

    for node in graph.nodes:
        if graph.in_degree(node) > 1:
            node_predecessors = list(graph.predecessors(node))
            for n1, n2 in it.combinations(node_predecessors, r=2):
                ancestor = nx.lowest_common_ancestor(graph, n1, n2)
                if ancestor is not None:
                    bubble = True
                    break
        if bubble:
            break
    if bubble:  # Recursive call.
        graph = simplify_bubbles(solve_bubble(graph, ancestor, node))

    return graph


def solve_entry_tips(graph, starting_nodes):

    tip = False

    for node in graph.nodes:

        # Searching for all the simple paths from a starting node to the
        # current one.
        all_paths = []
        for start in starting_nodes:
            path_list = list(nx.all_simple_paths(graph, start, node))
            if len(path_list) > 0:
                all_paths.append(*path_list)

        # If there's more than one such path, it means there's at least
        # an entry tip. We have to find the best path to keep.
        if len(all_paths) > 1:
            path_length = [len(path) for path in all_paths]
            weight_avg_list = [path_average_weight(graph, path)
                               for path in all_paths]
            graph = select_best_path(graph, all_paths, path_length,
                                     weight_avg_list, delete_entry_node=True)
            tip = True
            break

        if tip: # Recursive call.
            graph = solve_entry_tips(graph, starting_nodes)

    return graph


def solve_out_tips(graph, ending_nodes):

    tip = False

    for n1, n2 in it.combinations(ending_nodes, r=2):
        ancestor = nx.lowest_common_ancestor(graph, n1, n2)
        if ancestor is not None:

            # Searching for all the paths from ancestor to ending nodes.
            all_paths = [*nx.all_simple_paths(graph, ancestor, n1),
                         *nx.all_simple_paths(graph, ancestor, n2)]

            # If there's more than one such path, it means there's at least
            # an out tip. We have to find the best path to keep.
            if len(all_paths) > 1:
                path_length = [len(path) for path in all_paths]
                weight_avg_list = [path_average_weight(graph, path)
                                   for path in all_paths]
                # Selection of the best path.
                graph = select_best_path(graph, all_paths, path_length,
                                         weight_avg_list,
                                         delete_sink_node=True)
            tip = True
            break

        if tip: # Recursive call.
            graph = solve_out_tips(graph, ending_nodes)

    return graph


def get_starting_nodes(graph):
    """
    Retourne une liste de noeuds d’entrée i.e ceux qui n'ont pas de
    prédecesseurs.
    """

    return [node for node in graph.nodes if not graph.in_degree(node)]


def get_sink_nodes(graph):
    """
    Retourne une liste de noeuds de sortie i.e ceux qui n'ont pas de
    successeurs.
    """

    return [node for node in graph.nodes if not graph.out_degree(node)]


def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Retourne une liste de tuple(contig, longueur du contig).
    Un chemin entre un nœud d’entrée et un nœud de sortie
    correspond à une séquence contiguë (contig).
    """
    contig_list = []
    for start, end in it.product(starting_nodes, ending_nodes):
        if nx.has_path(graph, start, end):
            simple_path_reader = nx.all_simple_paths(graph, start, end)
            for path in simple_path_reader:
                contig = build_path(path)
                contig_list.append((contig, len(contig)))
    return contig_list


def build_path(path):
    """
    Build a string contig from a list of nodes.
    """

    path_str = ''
    for node in path[:-1]:
        path_str += node[:-1]
    # For the last node, we only need to add the last letter since the
    # rest of the kmer was added with the previous node.
    path_str += path[-1]
    return path_str


def save_contigs(contigs_list, output_file):

    with open(output_file, 'w') as f:
        for i, contig in enumerate(contigs_list):
            f.write(f'>contig_{i} len={contig[1]}\n')
            f.write(f'{fill(contig[0])}\n')


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    if parse:
        args = get_arguments()
        input_file = args.fastq_file
        kmer_size = args.kmer_size
        output_file = args.output_file
        graph_file = args.graphimg_file
    else:
        input_file = 'G:/RAID/Fac/M2_BI/Omiques/Assemblage - Génomique/debruijn-tp/data/eva71_hundred_reads.fq'
        input_file = 'G:/RAID/Fac/M2_BI/Omiques/Assemblage - Génomique/debruijn-tp/data/eva71_plus_perfect.fq'
        kmer_size = 22
        output_file = 'G:/RAID/Fac/M2_BI/Omiques/Assemblage - Génomique/debruijn-tp/output/eva71_plus_perfect_22-mer.fna'
        graph_file = 'graph.png'


    # Graph construction.
    kmer_dict = build_kmer_dict(input_file, kmer_size)
    graph = build_graph(kmer_dict)
    sources = get_starting_nodes(graph)
    sinks = get_sink_nodes(graph)
    contigs = get_contigs(graph, sources, sinks)

    # Bubble and tips resolution.
    print('Initial graph:')
    print(f'\tNb of nodes: {graph.number_of_nodes()}')
    print(f'\tNb of edges: {graph.number_of_edges()}')
    graph = simplify_bubbles(graph)
    print('Graph after bubbles simplification:')
    print(f'\tNb of nodes: {graph.number_of_nodes()}')
    print(f'\tNb of edges: {graph.number_of_edges()}')
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    print('Graph after entry tips removal:')
    print(f'\tNb of nodes: {graph.number_of_nodes()}')
    print(f'\tNb of edges: {graph.number_of_edges()}')
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    print('Graph after out tips removal:')
    print(f'\tNb of nodes: {graph.number_of_nodes()}')
    print(f'\tNb of edges: {graph.number_of_edges()}')


    # Writing the contig.
    # Since some nodes have been deleted, we have to reevaluate the source
    # and sink nodes.
    contigs = get_contigs(graph, get_starting_nodes(graph),
                          get_sink_nodes(graph))
    # We sort the contigs from longest to shortest.
    contigs.sort(key=itemgetter(1), reverse=True)
    print(f'Number of contigs found: {len(contigs)}')
    print(f'Length of longest contig: {contigs[0][1]}')
    save_contigs(contigs, output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':

    start_time = time.time()
    main()
    end_time = time.time() - start_time
    print(f'Assembly done in {end_time // 60:.0f} min {end_time % 60:.2f} s.')
