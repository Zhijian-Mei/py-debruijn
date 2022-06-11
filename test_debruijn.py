import copy
from collections import Counter
import sys
sys.setrecursionlimit(3000)
from tqdm import trange


class Node:
    """ Class Node to represent a vertex in the de bruijn graph """

    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0


class Edge:
    def __init__(self, lab):
        self.label = lab


def read_reads(fname):
    """ Read short reads in FASTA format. It is assumed that one line in the input file correspond to one read. """
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    reads = []

    for line in lines:
        if line[0] != '>':
            reads = reads + [line.rstrip()]
    return reads


def construct_graph(reads, k, threshold=3):
    """ Construct de bruijn graph from sets of short reads with k length word"""
    edges = dict()
    vertices = dict()
    pull_out_read = []

    for index in trange(len(reads)):
        read = reads[index]
        i = 0
        while i + k < len(read):
            v1 = read[i:i + k]
            v2 = read[i + 1:i + k + 1]
            if v1 in edges.keys():
                vertices[v1].outdegree += 1
                edges[v1] += [v2]
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edges[v1] = [v2]
            if v2 in edges.keys():
                vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1
                edges[v2] = []
            i += 1

    print(1)
    for edge in edges:
        previous = edges[edge]
        counter = Counter(previous)
        if len(counter) == 0:
            continue
        if len(counter) == 1:
            edges[edge] = list(counter)
        else:
            maxCountKmer = [counter.most_common(1)[0][0]]
            maxCount = counter.most_common(1)[0][1]
            for i in range(1, len(counter)):
                nextCount = counter.most_common()[i][1]
                if nextCount >= maxCount / threshold:
                    maxCountKmer += [counter.most_common()[i][0]]
            edges[edge] = maxCountKmer

    # current = edges['GPSVFP']
    # while current:
    #     if 'VTVSWN' in current:
    #         next = 'VTVSWN'
    #         current = edges[next]
    #         print(current)
    #         continue
    #     next = current[0]
    #     current = edges[next]
    #     print(current)
    #
    # quit()

    print(2)
    pull_out_kmer = []
    for edge in list(edges):
        if len(edges[edge]) > 1:
            for e in edges[edge]:
                pull_out_kmer.append([edge,e])
    print(3)
    if k != 6:
        for index in trange(len(reads)):
            read = reads[index]
            for kmer_pair in pull_out_kmer:
                if kmer_pair[0] in read and kmer_pair[1] in read:
                    if read not in pull_out_read:
                        pull_out_read.append(read)

    print(4)

    return (vertices, edges), pull_out_read


def DFS(current, E, vec, output, contig_copy):
    if current in vec:
        return
    vec.append(current)
    if len(E[current]) == 0:
        if vec not in output:
            result = vec[0]
            for i in range(1, len(vec)):
                result += vec[i][-1]
            # print(result,len(result))
            output.append(copy.deepcopy(vec))
            contig_copy.append(result)
        vec.pop()
        return
    for i in range(len(E[current])):
        DFS(E[current][i], E, vec, output, contig_copy)
    vec.pop()


def printPath(vec):
    # Print elements in the vector
    for ele in vec:
        print(ele, end=" ")
    print()


def output_contigs(g):
    """ Perform searching for Eulerian path in the graph to output genome assembly"""
    V = g[0]
    E = g[1]
    # Pick starting node (the vertex with zero in degree)
    starts = []
    for k in list(V.keys()):
        if V[k].indegree == 0:
            starts.append(k)
    print('Number of kmers have no income edges: ', len(starts))
    contig = []
    for i in trange(len(starts)):
        start = starts[i]
        # print('start=', start)
        current = start
        vec = []
        output = []
        contig_copy = []
        DFS(current, E, vec, output, contig_copy)
        contig.extend(contig_copy)
        # print('*' * 50)

    return contig


def print_graph(g):
    """ Print the information in the graph to be (somewhat) presentable """
    V = g[0]
    E = g[1]
    for k in V.keys():
        print("name: ", V[k].label, ". indegree: ", V[k].indegree, ". outdegree: ", V[k].outdegree)
        print("Edges: ")
        for e in E[k]:
            print(e)
        print()
