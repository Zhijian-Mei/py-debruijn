import copy
from collections import Counter


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


def construct_graph(reads, k):
    """ Construct de bruijn graph from sets of short reads with k length word"""
    edges = dict()
    vertices = dict()

    for read in reads:
        i = 0
        while i + k < len(read):
            v1 = read[i:i + k]
            v2 = read[i + 1:i + k + 1]
            if v1 in edges.keys():
                vertices[v1].outdegree += 1
                edges[v1] += [Edge(v2)]
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edges[v1] = [Edge(v2)]
            if v2 in edges.keys():
                vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1
                edges[v2] = []
            i += 1

    return (vertices, edges)


def DFS(current, E, vec, output,contig_copy):
    if current in vec:
        return
    vec.append(current)
    if len(E[current]) == 0:
        # print(vec)
        if vec not in output:
            result = vec[0]
            for i in range(1,len(vec)):
                result+=vec[i][-1]
            print(result,len(result))
            output.append(copy.deepcopy(vec))
            contig_copy.append(result)
        vec.pop()
        return
    for i in range(len(E[current])):
        DFS(E[current][i].label, E, vec, output,contig_copy)
    vec.pop()


def printPath(vec):
    # Print elements in the vector
    for ele in vec:
        print(ele, end=" ")
    print()


def output_contigs(g,step):
    """ Perform searching for Eulerian path in the graph to output genome assembly"""
    V = g[0]
    E = g[1]
    # Pick starting node (the vertex with zero in degree)
    all_starts = []
    for k in list(V.keys()):
        if V[k].indegree == 0:
            all_starts.append(k)
    # start = list(V.keys())[0]
    # for k in list(V.keys()):
    #     if V[k].indegree < V[start].indegree:
    #         start = k
    print('Number of kmers have no income edges: ', len(all_starts))
    contig = []
    for i in range(0,len(all_starts),step):
        try:
            starts = all_starts[i:i+step]
        except:
            starts = all_starts[i:]
        print('number of start kmer in this step: ',len(starts))
        for start in starts:
            print('start=', start)
            current = start
            vec = []
            output = []
            contig_copy = []
            DFS(current, E, vec, output,contig_copy)
            if len(contig) == 0:
                contig.extend(contig_copy)
            else:
                for j in range(len(contig)):
                    for k in range(len(contig_copy)):
                        if contig[j][-4:] == contig_copy[k][:4]:
                            contig[j] = contig[j]+contig_copy[k][4:]
                            contig_copy.remove(contig_copy[k])
                        elif contig[j][:4] == contig_copy[k][-4:]:
                            contig[j] = contig_copy[k] + contig[j][4:]
                            contig_copy.remove(contig_copy[k])
                contig.extend(contig_copy)
            print('*' * 50)
    for item in contig:
        print(item,len(item))
    print('Concatenation finished')
    quit()



def print_graph(g):
    """ Print the information in the graph to be (somewhat) presentable """
    V = g[0]
    E = g[1]
    for k in V.keys():
        print("name: ", V[k].label, ". indegree: ", V[k].indegree, ". outdegree: ", V[k].outdegree)
        print("Edges: ")
        for e in E[k]:
            print(e.label)
        print()
