# -*- coding: utf-8 -*-
"""
This module implements community detection.
"""
from __future__ import print_function

import array

import numbers
import warnings

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


class Status(object):
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}

    def __init__(self):
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])

    def __str__(self):
        return ("node2com : " + str(self.node2com) + " degrees : "
                + str(self.degrees) + " internals : " + str(self.internals)
                + " total_weight : " + str(self.total_weight))

    def copy(self):
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight

    def init(self, graph, weight, part=None):
        """Initialize the status of a graph with every node in one community"""
        count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight=weight)
        if part is None:
            for node in graph.nodes():
                self.node2com[node] = count
                deg = float(graph.degree(node, weight=weight))
                if deg < 0:
                    error = "Bad node degree ({})".format(deg)
                    raise ValueError(error)
                self.degrees[count] = deg
                self.gdegrees[node] = deg
                edge_data = graph.get_edge_data(node, node, default={weight: 0})
                self.loops[node] = float(edge_data.get(weight, 1))
                self.internals[count] = self.loops[node]
                count += 1
        else:
            for node in graph.nodes():
                com = part[node]
                self.node2com[node] = com
                deg = float(graph.degree(node, weight=weight))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg
                inc = 0.
                for neighbor, datas in graph[node].items():
                    edge_weight = datas.get(weight, 1)
                    if edge_weight <= 0:
                        error = "Bad graph type ({})".format(type(graph))
                        raise ValueError(error)
                    if part[neighbor] == com:
                        if neighbor == node:
                            inc += float(edge_weight)
                        else:
                            inc += float(edge_weight) / 2.
                self.internals[com] = self.internals.get(com, 0) + inc


__PASS_MAX = -1
__MIN = 0.0000001


def check_random_state(seed):
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError("%r cannot be used to seed a numpy.random.RandomState"
                     " instance" % seed)


def partition_at_level(dendrogram, level):
    partition = dendrogram[0].copy()
    for index in range(1, level + 1):
        for node, community in partition.items():
            partition[node] = dendrogram[index][community]
    return partition


def modularity(partition, graph, weight='weight'):
    if graph.is_directed():
        raise TypeError("Bad graph type, use only non directed graph")

    inc = dict([])
    deg = dict([])
    links = graph.size(weight=weight)
    if links == 0:
        raise ValueError("A graph without link has an undefined modularity")

    for node in graph:
        com = partition[node]
        deg[com] = deg.get(com, 0.) + graph.degree(node, weight=weight)
        for neighbor, datas in graph[node].items():
            edge_weight = datas.get(weight, 1)
            if partition[neighbor] == com:
                if neighbor == node:
                    inc[com] = inc.get(com, 0.) + float(edge_weight)
                else:
                    inc[com] = inc.get(com, 0.) + float(edge_weight) / 2.

    res = 0.
    for com in set(partition.values()):
        res += (inc.get(com, 0.) / links) - \
               (deg.get(com, 0.) / (2. * links)) ** 2
    return res


def best_partition(graph,
                   a_matrix,
                   partition=None,
                   weight='weight',
                   resolution=1.,
                   randomize=None,
                   random_state=None):
    dendo = generate_dendrogram(graph,
                                a_matrix,
                                partition,
                                weight,
                                resolution,
                                randomize,
                                random_state)
    return partition_at_level(dendo, len(dendo) - 1)


def generate_dendrogram(graph,
                        a_matrix,
                        part_init=None,
                        weight='weight',
                        resolution=1.,
                        randomize=None,
                        random_state=None):
    if graph.is_directed():
        raise TypeError("Bad graph type, use only non directed graph")

    # Properly handle random state, eventually remove old `randomize` parameter
    # NOTE: when `randomize` is removed, delete code up to random_state = ...
    if randomize is not None:
        warnings.warn("The `randomize` parameter will be deprecated in future "
                      "versions. Use `random_state` instead.", DeprecationWarning)
        # If shouldn't randomize, we set a fixed seed to get determinisitc results
        if randomize is False:
            random_state = 0

    # We don't know what to do if both `randomize` and `random_state` are defined
    if randomize and random_state is not None:
        raise ValueError("`randomize` and `random_state` cannot be used at the "
                         "same time")

    random_state = check_random_state(random_state)

    # special case, when there is no link
    # the best partition is everyone in its community
    if graph.number_of_edges() == 0:
        part = dict([])
        for i, node in enumerate(graph.nodes()):
            part[node] = i
        return [part]

    current_graph = graph.copy()
    status = Status()
    status.init(current_graph, weight, part_init)
    status_list = list()
    __one_level(current_graph, status, weight, resolution, random_state, a_matrix, status_list)
    new_mod = __modularity(status, a_matrix, status_list)
    partition = __renumber(status.node2com)
    status_list.append(partition)
    mod = new_mod
    current_graph = induced_graph(partition, current_graph, weight)
    status.init(current_graph, weight)

    while True:
        __one_level(current_graph, status, weight, resolution, random_state, a_matrix, status_list)
        new_mod = __modularity(status, a_matrix, status_list)
        if new_mod - mod < __MIN:
            break
        partition = __renumber(status.node2com)
        status_list.append(partition)
        mod = new_mod
        current_graph = induced_graph(partition, current_graph, weight)
        status.init(current_graph, weight)
    return status_list[:]


def induced_graph(partition, graph, weight="weight"):
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())

    for node1, node2, datas in graph.edges(data=True):
        edge_weight = datas.get(weight, 1)
        com1 = partition[node1]
        com2 = partition[node2]
        w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get(weight, 1)
        ret.add_edge(com1, com2, **{weight: w_prec + edge_weight})
    return ret


def __renumber(dictionary):
    count = 0
    ret = dictionary.copy()
    new_values = dict([])

    for key in dictionary.keys():
        value = dictionary[key]
        new_value = new_values.get(value, -1)
        if new_value == -1:
            new_values[value] = count
            new_value = count
            count += 1
        ret[key] = new_value

    return ret


def load_binary(data):
    data = open(data, "rb")

    reader = array.array("I")
    reader.fromfile(data, 1)
    num_nodes = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_nodes)
    cum_deg = reader.tolist()
    num_links = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_links)
    links = reader.tolist()
    graph = nx.Graph()
    graph.add_nodes_from(range(num_nodes))
    prec_deg = 0

    for index in range(num_nodes):
        last_deg = cum_deg[index]
        neighbors = links[prec_deg:last_deg]
        graph.add_edges_from([(index, int(neigh)) for neigh in neighbors])
        prec_deg = last_deg

    return graph


def __one_level(graph, status, weight_key, resolution, random_state, amatrix, status_list):
    modified = True
    nb_pass_done = 0
    cur_mod = __modularity(status, amatrix, status_list)
    new_mod = cur_mod

    while modified and nb_pass_done != __PASS_MAX:
        cur_mod = new_mod
        modified = False
        nb_pass_done += 1

        for node in random_state.permutation(list(graph.nodes())):
            com_node = status.node2com[node]
            degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight * 2.)  # NOQA
            neigh_communities = __neighcom(node, graph, status, weight_key)
            remove_cost = - resolution * neigh_communities.get(com_node, 0) + \
                          (status.degrees.get(com_node, 0.) - status.gdegrees.get(node, 0.)) * degc_totw
            remove_cost_entropy = __remove_cost_entropy(node, status, amatrix, status_list)

            __remove(node, com_node,
                     neigh_communities.get(com_node, 0.), status)
            best_com = com_node
            best_increase = 0
            for com, dnc in random_state.permutation(list(neigh_communities.items())):
                incr = remove_cost + resolution * dnc - \
                       status.degrees.get(com, 0.) * degc_totw
                incr_entropy = __inc_entropy(node, status, amatrix, com, status_list) + remove_cost_entropy

                if incr + incr_entropy > best_increase:
                    best_increase = incr + incr_entropy
                    best_com = com
            __insert(node, best_com,
                     neigh_communities.get(best_com, 0.), status)
            if best_com != com_node:
                modified = True
        new_mod = __modularity(status, amatrix)
        if new_mod - cur_mod < __MIN:
            break


def __neighcom(node, graph, status, weight_key):
    """
    Compute the communities in the neighborhood of node in the graph given
    with the decomposition node2com
    """
    weights = {}
    for neighbor, datas in graph[node].items():
        if neighbor != node:
            edge_weight = datas.get(weight_key, 1)
            neighborcom = status.node2com[neighbor]
            weights[neighborcom] = weights.get(neighborcom, 0) + edge_weight

    return weights


def __remove(node, com, weight, status):
    """ Remove node from community com and modify status"""
    status.degrees[com] = (status.degrees.get(com, 0.)
                           - status.gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) -
                                  weight - status.loops.get(node, 0.))
    status.node2com[node] = -1


def __insert(node, com, weight, status):
    """ Insert node into community and modify status"""
    status.node2com[node] = com
    status.degrees[com] = (status.degrees.get(com, 0.) +
                           status.gdegrees.get(node, 0.))
    status.internals[com] = float(status.internals.get(com, 0.) +
                                  weight + status.loops.get(node, 0.))


def __modularity(status, amatrix, status_list):
    """
    Fast compute the modularity of the partition of the graph using
    status precomputed
    """

    part = partition_at_level(status_list, len(status_list) - 1)

    links = float(status.total_weight)
    result = 0.
    result_a = 0.
    prob_a = dict([])
    prob_a[0] = np.sum(amatrix == 0, axis=0) * 1.0 / amatrix.shape[0]
    prob_a[1] = np.sum(amatrix == 1, axis=0) * 1.0 / amatrix.shape[0]
    ha = -(prob_a[0] * np.log2(prob_a[0]) + prob_a[1] * np.log2(prob_a[1]))

    for community in set(status.node2com.values()):
        in_degree = status.internals.get(community, 0.)
        degree = status.degrees.get(community, 0.)
        if links > 0:
            result += in_degree / links - ((degree / (2. * links)) ** 2)
        list_comm_nodes = np.array([item for item in status.node2com.keys() if status.node2com[item] == community])
        list_nodes = np.array([item for item in part.keys() if part[item] in list_comm_nodes])
        result_a += np.average(np.abs(ha - __modularity_community_attribute(list_nodes, amatrix)))

    return result + result_a


def __remove_cost_entropy(node, status, a_matrix, status_list):
    part = partition_at_level(status_list, len(status_list) - 1)
    node_com = status.node2com[node]
    list_comm_nodes = np.array([item for item in status.node2com.keys() if status.node2com[item] == node_com])
    list_nodes_origin = np.array([item for item in part.keys() if part[item] in list_comm_nodes])
    ha_ori = __modularity_community_attribute(list_nodes_origin, a_matrix)
    list_nodes_remove = np.array([item for item in part.keys() if part[item] in list_comm_nodes and part[item] != node])
    ha_remove = __modularity_community_attribute(list_nodes_remove, a_matrix)
    return ha_remove - ha_ori


def __inc_entropy(node, status, a_matrix, com, status_list):
    part = partition_at_level(status_list, len(status_list) - 1)
    list_comm_nodes = np.array([item for item in status.node2com.keys() if status.node2com[item] == com])
    list_nodes_origin = np.array([item for item in part.keys() if part[item] in list_comm_nodes])
    ha_ori = __modularity_community_attribute(list_nodes_origin, a_matrix)
    list_nodes_add = np.array([item for item in part.keys() if part[item] in list_comm_nodes and part[item] == node])
    ha_add = __modularity_community_attribute(list_nodes_add, a_matrix)
    return ha_add - ha_ori


def __modularity_community_attribute(list_nodes, amatrix):
    a_nodes = np.array([amatrix[int(item)] for item in list_nodes])
    prob_a = dict([])
    prob_a[0] = np.sum(a_nodes == 0, axis=0) * 1.0 / a_nodes.shape[0]
    prob_a[1] = np.sum(a_nodes == 1, axis=0) * 1.0 / a_nodes.shape[0]
    ha = -(prob_a[0] * np.log2(prob_a[0]) + prob_a[1] * np.log2(prob_a[1]))

    return ha


def readAttributes(graphfile):
    f = open(graphfile, "r")
    line = f.readline()

    while line == "\n":
        line = f.readline()

    n_nodes = int(line.split(" ")[0])
    n_attributes = int(line.split(" ")[2])

    a = np.zeros([n_nodes, n_attributes])

    for i in range(n_nodes):
        line = f.readline()
        start = line.index("[")
        end = line.index("]")

        attr = line[start + 1: end].split(", ")

        for j in range(n_attributes):
            a[i][j] = int(attr[j])
    return a


def plotGraph(filename):
    f = open(filename, "r")
    line = f.readline()

    while line == "\n":
        line = f.readline()

    n_nodes = int(line.split(" ")[0])
    n_edges = int(line.split(" ")[1])
    g = nx.Graph()

    for i in range(n_nodes):
        g.add_node("v" + 'i')
        f.readline()

    for j in range(n_edges):
        line = f.readline()

        source = line[1:-2].split(" : ")[0]
        target = line[1:-2].split(" : ")[1]
        """print(source + " " + target)"""
        g.add_edge(source, target)

    pos = nx.spring_layout(g)
    nx.draw_networkx_nodes(g, pos, node_size=300)
    nx.draw_networkx_edges(g, pos)

    nx.draw_networkx_labels(g, pos, fontsize=11)

    plt.show()
    return g


def loadGraph(filename):
    f = open(filename, "r")
    line = f.readline()

    while line == "\n":
        line = f.readline()

    n_nodes = int(line.split(" ")[0])
    n_edges = int(line.split(" ")[1])
    g = nx.Graph()

    for i in range(n_nodes):
        g.add_node("v" + 'i')
        f.readline()

    for j in range(n_edges):
        line = f.readline()

        source = line[1:-2].split(" : ")[0]
        target = line[1:-2].split(" : ")[1]
        """print(source + " " + target)"""
        g.add_edge(source, target)


if __name__ == '__main__':
    g = loadGraph('/home/xihui/eclipse-workspace/Graph2019/gnp20.txt')
    amatrix = readAttributes('/home/xihui/eclipse-workspace/Graph2019/gnp20.txt')
    print(amatrix)
