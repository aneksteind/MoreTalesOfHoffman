import os
import itertools
import click
import scipy.sparse

import numpy as np
import networkx as nx

from math import floor
from scipy.linalg import eigvalsh

import schrijver
import util


@click.group()
def cli():
    pass


def laplacian_pieces(G, nodelist=None, weight="weight"):
    """
        Get the diagonal and adjacency matrices of G
    """

    if nodelist is None:
        nodelist = G.nodes()

    A = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, weight=weight, format="csr")
    n, m = A.shape

    diags = A.sum(axis=1)
    D = scipy.sparse.spdiags(diags.flatten(), [0], m, n, format="csr")

    return D, A


def adjacency(G, nodelist=None, weight="weight"):
    """
        Returns the sparse adjacency matrix
        representation of the graph.
    """

    if nodelist is None:
        nodelist = G.nodes()

    A = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, weight=weight, format="csr")

    return A


def spectrum(G):
    """
        Gets the eigenvalues of the graph, represented
        with a sparse of dense adjacency matrix.
    """

    try:
        return eigvalsh(G.todense())
    except ValueError:
        return eigvalsh(G)


def lima(m, n, d):
    """
        Calculates the Lima bound of a graph given the
            - number of graph edges m
            - number of graph vertices n
            - largest eigenvalue of the signless laplacian
    """

    num = 2 * m
    den = num - n * d
    return 1 + (num / den)


def kolotina(mu_1, delta, theta):
    """
        Calculates the Kolotilina bound of a graph given the
            - largest eignvalue of the graph mu_1
    """

    return mu_1 / (mu_1 - delta + theta)


def ando_lin(eigen_values):
    """
        Calculates the Ando-Lin bound of a graph.
    """

    pos, neg, _ = cluster_eignvalues(eigen_values)

    sp = sumsq(pos)
    sn = sumsq(neg)

    try:
        return 1 + max([sp / sn, sn / sp])
    except:
        raise Exception(
            f"Division by zero in Ando-Lin bound construction: negative = {sn}, positive = {sp}"
        )


def hoffman(mu_1, mu_n):
    """
        Calculates the Hoffman bound of a graph given
        the largest and smallest eigenvalues of the adjacency matrix.
    """
    return 1 + (mu_1 / abs(mu_n))


def bound1(n, lambda_1, delta):
    return n * (lambda_1 - delta) / lambda_1 if lambda_1 else 0


def bound2(n, mu_n, delta):
    return n * abs(mu_n) / (delta + abs(mu_n)) if (delta + abs(mu_n)) > 0 else 0


def sumsq(values):
    """
        Calculates the sum of squares of a list of values.
    """

    return sum(map(lambda x: x ** 2, values))


def check(x, y):
    """
        A simple check determining whether x <= y,
        reports whether the two values were close.
        The bound check passes if x <= y.
    """

    close = is_close(x, y)
    passed = x <= y

    if not passed and not close:
        return None
    else:
        return (passed, close, x, y)


def cluster_eignvalues(vals):
    """
        Separates eignvalues into postive, negative, and zero groups
    """

    pos, neg, zer = [], [], []
    for val in vals:
        if is_close(val, 0.0):
            zer.append(val)
        elif val > 0:
            pos.append(val)
        elif val < 0:
            neg.append(val)
        else:
            raise Exception("The world is a weird place.")

    return pos, neg, zer


def inertia(mus):
    """
        Calculates the inertia of a graph.
    """
    pos, negs, zeros = cluster_eignvalues(mus)

    return len(zeros) + min(len(pos), len(negs))


def is_regular(G):
    """
        Checks if all nodes in the graph have the same degree.
    """

    deg = None

    for node in G:
        d = G.degree(node)
        if deg is None:
            deg = d
        if d != deg:
            return False

    return True


def prepare_write_status(graph_name, *args):
    return f'"{graph_name}",' + ",".join([str(arg) for arg in args]) + "\n"


def compare(G, name, tmp_file, verbose=False):
    strings = []

    if verbose:
        print(f"Testing {name}")

    try:
        x_vec = schrijver.thmbar(G)
    except:
        raise Exception(f"{name}: thmb 0 exception.")

    D, A = laplacian_pieces(G)

    # eigenvalues of A
    mus = spectrum(A)

    pos, negs, _ = cluster_eignvalues(mus)
    rank = len(pos) + len(negs)

    regularity = is_regular(G)

    al = ando_lin(mus)

    results = check(al, x_vec)
    if results is None:
        raise Exception(f"{name} failed the ando-lin bound!")

    strings += prepare_write_status(name, "ando-lin", *results, regularity)

    results = check(x_vec, rank)
    if results is None:
        raise Exception(f"{name} failed the rank bound!")

    strings += prepare_write_status(name, "sjbarrank", *results, regularity)

    with open(tmp_file, "a+") as logfile:
        for string in strings:
            logfile.write(string)

    return strings


def is_close(x, y, thresh=1e-8):
    """
        Tests if x is close to y as measured by some threshold.
    """

    diff = x - y
    return diff > (-thresh) and diff < thresh


def circulant_adj(js, num_vertices):
    """
        Returns the adjacency graph of a given
        circulant graph specification.
    """

    # make a fresh/empty adjacency matrix
    adj = np.zeros((num_vertices, num_vertices), dtype="uint8")

    # make a connect between each node and a node
    # that is j spaces to the left and to right of the node
    for i in range(num_vertices):
        for j in js:
            adj[i][(i - j) % num_vertices] = 1
            adj[i][(i + j) % num_vertices] = 1

    return adj


def name_circulant(num_vertices, j_value_set):
    """
        Gives a string representation of a given
        circulant graph.
    """

    return f"Cir [{num_vertices}] [{j_value_set}]"


def powerset(iterable):
    """
        Gets the powerset of a collection.
    """

    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(1, len(s) + 1)
    )


def from_dot(path):
    """
        Reads a dot file, returning a NetworkX graph
    """
    return nx.Graph(nx.drawing.nx_agraph.read_dot(path))


def circulant_gen(min_order, max_order):
    """
        A generator yielding all circulant graphs with
        number of vertices >= min_order and <= max_order.
    """

    for num_vertices in range(min_order, max_order + 1):
        all_j_values = [x for x in range(1, floor(num_vertices / 2.0))]
        j_values_iter = powerset(all_j_values)

        # for every possible offset combination
        for j_value_set in j_values_iter:
            # get the adjacency matrix of the circulant graph
            adj = circulant_adj(j_value_set, num_vertices)
            G = nx.from_numpy_matrix(adj)

            if G.size() > 0 and nx.is_connected(G):
                yield (G, name_circulant(num_vertices, j_value_set))


def named_graphs(min_order, max_order, dotdir, verbose=False):
    """
        Gets named graphs represented as dot files.
    """

    path = dotdir
    # graph the full list of graph file names
    files = os.listdir(path)

    for filename in files:
        G = from_dot(path + filename)
        if (
            G.size() > 0
            and len(G) >= min_order
            and len(G) <= max_order
            and nx.is_connected(G)
        ):
            if verbose:
                print(f"Reading {filename}")
            yield (G, filename)


@cli.command()
@click.option(
    "--max", "max_order", required=True, type=int, help="the maximum number of vertices"
)
@click.option(
    "--min", "min_order", type=int, default=1, help="the minimum number of vertices"
)
@click.option(
    "--graph-source",
    "source",
    required=True,
    type=click.Choice(["named", "circulant"]),
    help="named graphs from the mathematica database or generated circulant graphs",
)
@click.option(
    "--out", "outfile", required=True, help="the .csv file to write results to"
)
@click.option(
    "--verbose/--quiet", default=False, help="print the graph names during execution"
)
@click.option(
    "--wolfram/--custom",
    default=False,
    help="Flag specifying use of DOT files generated from Wolfram GraphData, names must match those in graphdata.csv",
)
@click.option(
    "--dot-dir",
    "dotdir",
    default="GraphData/",
    show_default=True,
    help="the directory containing DOT files for use in named graph bounds checking",
)
def test(min_order, max_order, source, outfile, verbose, wolfram, dotdir):
    # write headers of temp csv file
    with open(outfile, "w") as logfile:
        logfile.write(
            ",".join(
                [
                    "Name",
                    "Test Type",
                    "Passed Bound Check",
                    "Bound is Close",
                    "Lower",
                    "Upper",
                    "IsRegular",
                ]
            )
            + "\n"
        )

    # get the graphs to test
    if source == "named":
        graphs = named_graphs(min_order, max_order, dotdir, verbose=verbose)
    else:
        graphs = circulant_gen(min_order, max_order)

    # run and write the tests
    _ = list(map(lambda a: compare(*a, outfile, verbose=verbose), graphs))

    # write Wolfram graphs in a special way (meta data included)
    if wolfram:
        util.write_final(outfile, source == "named")


if __name__ == "__main__":
    cli()
