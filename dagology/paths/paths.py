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