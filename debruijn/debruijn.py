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

import argparse
import os
import sys
import networkx as nx
import matplotlib 
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Aurelien GRAY"
__copyright__ = "CYtech"
__credits__ = ["Aurelien GRAY"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurelien GRAY"
__email__ = "grayaureli@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file,'rt') as file :
        i=0
        for line in file:
            line.strip('\n')
            if (i%4 ==1):
                
                yield line.replace('\n','')
            i+=1
    pass


def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size +1):
        yield read[i:i+kmer_size]
    pass


def build_kmer_dict(fastq_file, kmer_size):
    seq_dict = {}
    readfile = read_fastq(fastq_file)
    for j in readfile:    
        gen_seq = cut_kmer(j, kmer_size)
        for k in gen_seq:
            seq_dict[k] = seq_dict.get(k, 0) + 1
        
    return seq_dict
        
    pass


def build_graph(kmer_dict):
    G=nx.Graph()
    H=nx.DiGraph(G)
    
    for elt in kmer_dict:
        #print(kmer_dict[elt])
        n = len(elt)
        x,y = elt[0:n-1],elt[1:n]
        H.add_edge(x,y, weight=kmer_dict[elt])
    return H
    pass


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    if (delete_entry_node ==True and delete_sink_node== True):
        for i in range(len(path_list)):
            n = len(path_list[i])
            for k in range(n):    
                graph.remove_node(path_list[i][k])
        
    elif (delete_entry_node ==True and delete_sink_node== False):
        for i in range(len(path_list)):
            n = len(path_list[i])

            for k in range(0,n-1):    
                graph.remove_node(path_list[i][k])

        
    elif (delete_entry_node ==False and delete_sink_node== True):
        for i in range(len(path_list)):
            n = len(path_list[i])

            for k in range(1,n):    
                graph.remove_node(path_list[i][k])
        
    elif (delete_entry_node ==False and delete_sink_node== False):
        for i in range(len(path_list)):
            n = len(path_list[i])           
            for k in range(1,n-1):
                graph.remove_node(path_list[i][k])
    return graph
    pass

def std(data):
    return statistics.stdev(data)
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    
    max_weight = max(weight_avg_list)
    best_path_w =[]
    best_path_l=[]
    for i in range(len(path_list)):
        
        if(weight_avg_list[i]==max_weight):
            best_path_w.append(path_list[i])
        
        for j in range(len(best_path_w)):
            max_length = 0
            if (max_length < len(best_path_w[j])):
                max_length = len(best_path_w[j])
            if(len(best_path_w[j])==max_length):
                best_path_l.append(best_path_w[j])

    
    best_path = random.choice(best_path_l)
    print(random.choice(best_path_l))
    
    index = path_list.index(best_path)
    path_list.pop(index)
    
    
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph     
    pass


def path_average_weight(graph, path):
    all_w = []
    for i in range(len(path)-1):
        findNode = False
        for node, nbrs in graph.adj.items():
            if(findNode == True):
                break
            if(node == path[i]):
                for nbr, eattr in nbrs.items():
                    if(nbr == path[i+1]):
                        all_w.append(eattr['weight'])
                        findNode = True
                        break
            
    return statistics.mean(all_w)
    
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    p = nx.all_simple_paths(graph, source=ancestor_node,target=descendant_node)
    all_path=[]
    path_length=[]
    weight_list=[]
    for path in p:
        all_path.append(path)
        path_length.append(len(path))
        weight_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, all_path, path_length, weight_list)
    
    
    return graph
          
    pass

def simplify_bubbles(graph):
    ancestor_list =[]
    descendant_list=[]
    for node in graph.nodes():
        pre = graph.predecessors(node)
        suc = graph.successors(node)
        preList =[]
        sucList = []
        
        for preNode in pre:
            preList.append(preNode)
        for sucNode in suc:
            sucList.append(sucNode)
            
        if (len(sucList)>=2):
                ancestor_list.append(node)
                
        if (len(preList)>=2):
                descendant_list.append(node) 
    
    for i in range(len(ancestor_list)):
        for j in range(len(descendant_list)):

            if (ancestor_list[i] != descendant_list[j]):
                graph = solve_bubble(graph, ancestor_list[i], descendant_list[j])
                
    return graph
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    
    return [n for n, d in graph.in_degree() if d == 0]
    
    pass

def get_sink_nodes(graph):
    return [node for node, out_degree in graph.out_degree() if out_degree == 0]
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    liste = []
    n = len(starting_nodes[0])
    for i in starting_nodes:
        for j in ending_nodes:        
            p = nx.shortest_path(graph, source=i,target=j)
            res=p[0]
            for k in range(1,len(p)):
                res = res + p[k][n-1]
            length = len(res)
            liste.append((res,length))
            
    return liste
    pass

def save_contigs(contigs_list, output_file):
    with open(output_file,'w') as fichier:
        for i in range(len(contigs_list)):
            fichier.write('>contig_'+ str(i) + ' len=' + str(contigs_list[i][1]) + '\n' + fill(contigs_list[i][0]))
            fichier.write('\n')
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    #we get the file
    fastsq = args.fastq_file
    print('Making Graph')
    kmer_dict = build_kmer_dict(fastsq, 10)
    H = build_graph(kmer_dict)
    print('Graph Done')
    print('Simplify bubbles :')
    H = simplify_bubbles(H)
    print('Bubble Done')
    start = get_starting_nodes(H)
    sink = get_sink_nodes(H)
    print('Get contigs')
    test = get_contigs(H, start, sink)
    print('get contigs done')
    print('save contigs')
    save_contigs(test, 'file.fna')
    print('save contigs done')



if __name__ == '__main__':
    main()
