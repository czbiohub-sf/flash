import argparse
import glob
import itertools
import json
import os
from collections import defaultdict

import networkx

from common import (
    Component,
    Gene
)
from gurobipy import (
    GRB,
    Model,
    quicksum
)


def set_required_value_for_guides(guide_vars, m, guides, value):
    for guide in guides:
        if guide in guide_vars:
            m.addConstr(guide_vars[guide], GRB.EQUAL, value)


def optimize(genes, required_guides=[], excluded_guides=[], MIN_FRAGMENT_SIZE=200,
             verbose=False):
    print("optimizing %d genes starting with %s" % (len(genes), genes[0].name))

    m = Model()

    guide_vars = {}

    for g in genes:
        for t in g.targets:
            if t.guide not in guide_vars:
                guide_vars[t.guide] = m.addVar(vtype=GRB.BINARY, name=t.guide)

    # Set required and excluded guides.
    set_required_value_for_guides(guide_vars, m, required_guides, 1)
    set_required_value_for_guides(guide_vars, m, excluded_guides, 0)

    impossible_snps = []

    def targets_in_range(targets, r):
        return [t for t in targets if t.cut in r]

    # add constraints
    for i, g in enumerate(genes):
        if verbose and i % 100 == 0:
            print("Added constraints for {} genes.".format(i))

        if len(g.targets) < 1:
            continue

        longest_fragment_possible = g.targets[-1].cut - g.targets[0].cut

        if longest_fragment_possible > MIN_FRAGMENT_SIZE:
            # At least 2 cuts
            m.addConstr(quicksum(guide_vars[t.guide] for t in g.targets),
                        GRB.GREATER_EQUAL, 2)

            # No two cuts within 200 of each other
            for t in g.targets:
                for t2 in targets_in_range(g.targets, range(t.cut, t.cut + MIN_FRAGMENT_SIZE)):
                    if t != t2:
                        m.addConstr(quicksum([guide_vars[t.guide], guide_vars[t2.guide]]),
                                    GRB.LESS_EQUAL, 1)
        else:
            if longest_fragment_possible > 100:
                print(g.name, " is only ", longest_fragment_possible, " long.")
                print("Insisting on using longest possible cut.")
                # use both ends
                m.addConstr(quicksum([guide_vars[g.targets[0].guide],
                                      guide_vars[g.targets[-1].guide]]),
                            GRB.GREATER_EQUAL, 2)
                # exclude middle
                m.addConstr(quicksum(guide_vars[t.guide] for t in g.targets[1:-1]),
                            GRB.LESS_EQUAL, 0)
            else:
                print(g.name, " is only ", longest_fragment_possible, " long.")
                print("Insisting on using a single cut.")
                # Exactly one cut
                m.addConstr(quicksum(guide_vars[t.guide] for t in g.targets),
                            GRB.EQUAL, 1)

        # Cover SNPs
        for mutation, snp_range in g.get_mutation_ranges():
            ranges = [
                range(snp_range.start - 150, snp_range.stop + 150),
                range(snp_range.start - 400, snp_range.start),
                range(snp_range.stop, snp_range.stop + 400)
            ]
            for r in ranges:
                targets = targets_in_range(g.targets, r)
                if len(targets) > 0:
                    m.addConstr(quicksum(guide_vars[t.guide] for t in targets),
                                GRB.GREATER_EQUAL, 1)
                else:
                    print("Impossible SNP in " + g.name + " " + str(snp_range))
                    impossible_snps.append((g, mutation, snp_range))

                    # Force failure
                    m.addConstr(quicksum(guide_vars[t.guide] for t in targets),
                                GRB.GREATER_EQUAL, 1)

        # Exclude guides that overlap SNPs
        for t in g.targets:
            if g.target_overlaps_mutation(t):
                m.addConstr(guide_vars[t.guide], GRB.EQUAL, 0)

    guide_freq = defaultdict(int)
    for g in genes:
        for guide, cut in g.targets:
            guide_freq[guide] += 1

    m.setObjective(
        quicksum(
            guide_vars.values())-1.1*quicksum([v * guide_freq[g] for g, v in guide_vars.items()]),
        GRB.MINIMIZE
    )

    if not verbose:
        m.Params.LogToConsole = 0
    m.update()
    print('Done building model')
    m.optimize()

    library = []

    if m.Status in (GRB.Status.INF_OR_UNBD, GRB.Status.INFEASIBLE, GRB.Status.UNBOUNDED):
        print('The model cannot be solved because it is infeasible or unbounded')
    elif m.Status == GRB.Status.SUBOPTIMAL:
        print('The model may have a suboptimal solution.')
    else:
        for var in m.getVars():
            if var.X == 1.0:
                library.append(var.VarName)

    print("Library Size: ", len(library))
    return m, genes, library, impossible_snps


def get_guides_from_file_or_empty(f):
    if f:
        return [e.strip() for e in f.readlines()]
    else:
        return []

def find_components(genes):
    print("Finding components...")
    guide_to_genes = defaultdict(set)

    for g in genes:
        for target in g.targets:
            guide_to_genes[target.guide].add(g.name)

    graph = networkx.Graph()

    for g in genes:
        graph.add_node(g.name)

    def guide_to_node(guide):
        return 'guide: ' + guide
    def is_guide(node):
        if node[:6] == 'guide:':
            return True
        else:
            return False

    for guide in guide_to_genes:
        graph.add_node(guide_to_node(guide))
        for gene in guide_to_genes[guide]:
            graph.add_edge(guide_to_node(guide), gene)

    connected_components = networkx.connected_components(graph)
    connected_genes = [[node for node in c if not is_guide(node)] for c in connected_components]

    gene_dict = {g.name: g for g in genes}

    components = [Component([gene_dict[name] for name in sorted(c)]) for c in connected_genes]

    print("There are ", len(components), " components.")

    return components

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--include',
                        type=argparse.FileType('r'),
                        metavar="file")
    parser.add_argument('--exclude',
                        type=argparse.FileType('r'),
                        metavar="file")
    parser.add_argument('--output',
                        type=argparse.FileType('w'),
                        metavar="file",
                        default="library.txt")
    args = parser.parse_args()

    print("Loading genes...")
    gene_names = sorted([os.path.splitext(os.path.basename(f))[0] for f in glob.glob(
        'generated_files/under_version_control/genes/*.fasta')])

    genes = [Gene(name) for name in gene_names]

    existing_guides = get_guides_from_file_or_empty(args.include)
    excluded_guides = get_guides_from_file_or_empty(args.exclude)

    with open('inputs/additional/padding.json', 'r') as fp:
        padding_seqs = json.load(fp)

    for g in genes:
        # Targets are generated by the snakemake step
        g.load_targets("dna_good_5_9_18.txt")
        if g.name in padding_seqs:
            g.verify_padding(padding_seqs[g.name])
        else:
            assert g.padding is None

    genes_without_targets = [g for g in genes if g.targets is None]

    if len(genes_without_targets) > 0:
        print("The following genes did not have targets: %s" %
              ",".join(genes_without_targets))

    genes = [g for g in genes if g.targets is not None]

    components = find_components(genes)

    solved_guides = set()

    for comp in components:
        m, genes, library, impossible_snps = optimize(
            comp.genes,
            existing_guides,
            excluded_guides
        )
        if not library:
            m, genes, library, impossible_snps = optimize(
                comp.genes,
                existing_guides,
                excluded_guides,
                MIN_FRAGMENT_SIZE=178
            )
        if library:
            solved_guides = solved_guides.union(set(library))
        else:
            print("Component %s was not solved" % comp.name)

    for guide in sorted(solved_guides):
        args.output.write("{}\n".format(guide))


if __name__ == "__main__":
    main()
