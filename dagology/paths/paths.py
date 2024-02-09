__author__ = "\n".join(["Nicolas Kozak (nicolas.kozak5@gmail.com)"])


def shortest_path(node, weights, distances, predecessors):
    for v, w in weights.items():
        new_distance = distances[node] + w
        if distances[v] > new_distance:
            distances[v], predecessors[v] = new_distance, node


def longest_path(node, weights, distances, predecessors):
    for v, w in weights.items():
        new_distance = distances[node] - w
        if distances[v] > new_distance:
            distances[v], predecessors[v] = new_distance, node

def greedy_shortest_path(node, weights, distances, predecessors):
    if weights:
        min_neighbor = min(weights, key=weights.get)
        new_distance = distances[node] + weights[min_neighbor]
        if distances[min_neighbor] > new_distance:
            distances[min_neighbor] = new_distance
            predecessors[min_neighbor] = node

def greedy_longest_path(node, weights, distances, predecessors):
    if weights:
        max_neighbor = max(weights, key=weights.get)
        new_distance = distances[node] - weights[max_neighbor]
        if distances[max_neighbor] > new_distance:
            distances[max_neighbor] = new_distance
            predecessors[max_neighbor] = node