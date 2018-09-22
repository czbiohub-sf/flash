#!/usr/bin/env python3

import json
import glob
import sys
from collections import defaultdict


def is_aro(value):
    return value.startswith("ARO:") and value.split("ARO:", 1)[1].isdigit()


def parse_aro_obo(filename):
    # these are ignored silently
    ignored_keys = set(["auto-generated-by", "date", "default-namespace",
        "format-version", "synonym", "xref"])
    # anything not in this list will be accepted but with a warning
    known_edge_types = set(["confers_resistance_to", "confers_resistance_to_drug",
        "derives_from", "has_part", "part_of", "participates_in", "regulates",
        "targeted_by", "targeted_by_drug"])
    # atom means not list
    known_atom_types = set(["auto-generated-by", "date", "def", "default-namespace",
        "format-version", "id", "name", "namespace", "ontology", "saved-by"])
    # silently ignore terms that start with these tokens
    ignored_tokens = set(["[Typedef]"])
    class Count(object):
        def __init__(self):
            self.lines = 0
            self.errors = 0
            self.warnings = 0
            self.ignored_edges = 0
            self.ignored_terms = 0
    count = Count()
    all_terms = {}
    with open(filename) as input_file:
        obo_raw_lines = input_file.read().split('\n')
    # skip header, this gets rid of spurious error trying to parse
    # header as term
    if len(obo_raw_lines) > 8:
        if obo_raw_lines[7].strip() == "[Term]":
            obo_raw_lines = obo_raw_lines[7:]
    # append empty line to make sure we emit last term
    obo_raw_lines.append('')
    term = {}
    def error(message="Term will be ignored", silent=False):
        if not silent:
            count.errors += 1
            print("Error:{}:{}: {}.".format(filename, count.lines, message))
        if "parse_error" not in term:
            if not silent:
                count.ignored_terms += 1
            term["parse_error"] = True
    def warn(message):
        count.warnings += 1
        print("WARNING:{}:{}: {}.".format(filename, count.lines, message))
    for raw_line in obo_raw_lines:
        count.lines += 1
        line = raw_line.strip()
        if not line or line == "[Term]":
            if term and "parse_error" not in term:
                if "id" not in term:
                    error("Missing term ID")
                    continue
                all_terms[term["id"]] = term
            term = {}
            continue
        if "parse_error" in term:
            continue
        parts = line.split(": ", 1)
        if len(parts) != 2:
            error("Unrecognized token '{}'".format(line),
                line in ignored_tokens)
            continue
        key, value = parts
        if key in ignored_keys:
            continue
        if key == "id":
            if not is_aro(value):
                error("Non-ARO ID")
                continue
            if value in all_terms:
                error("Duplicate ID")
                continue
        # part to the right of ' ! ' is a comment
        value = value.rsplit(" ! ", 1)[0].strip()
        # if relationship, it describes an edge in the ARO graph
        if key == "relationship":
            value = tuple(value.split(sep=None, maxsplit=1))
            if value[0] not in known_edge_types:
                warn("Unknown edge type '{}'".format(value[0]))
            if not is_aro(value[1]):
                count.ignored_edges += 1
                warn("Ignoring edge to non-ARO node '{}'".format(value[1]))
                continue
        if key in known_atom_types:
            if key in term:
                error("Duplicate key")
                continue
            term[key] = value
        else:
            if key not in term:
                term[key] = []
            # print(key, value)
            term[key].append(value)
    print("Parsed {} ARO terms.  Ignored {} terms and {} edges due to errors."
        .format(len(all_terms), count.ignored_terms, count.ignored_edges))
    return all_terms


def chase(all_terms, term, get_edges, level=0, visited=None):
    """Return the set of vertices reachable from term via get_edges()."""
    results = set()
    if term not in all_terms:
        print("WARNING: Referenced {} not found in ARO.OBO.".format(term))
        return results
    results.add((level, term))
    if visited == None:
        visited = set()
    visited.add(term)
    for parent in get_edges(all_terms[term]):
        if parent not in visited:
            results |= chase(all_terms, parent, get_edges, level + 1, visited)
    return results


def related_via(term_dict, *edge_labels):
    return [
        edge[1]
        for edge in term_dict.get("relationship", [])
        if edge[0] in edge_labels
    ]


def is_a(term_dict):
    return (
        term_dict.get("is_a", [])
    )


def is_sort_of(term_dict):
    return (
        term_dict.get("is_a", []) +
        related_via(term_dict, "part_of", "regulates", "participates_in")
    )


def conferred_resistances(term_dict):
    return related_via(term_dict,
        "confers_resistance_to", "confers_resistance_to_drug")


def antibiotic_fully_qualified_name(all_terms, term, level, add_aro=False, use_plural=False):
    # Get term superclasses sorted by level (in order of increasing generality)
    parents = sorted(chase(all_terms, term, is_a))
    assert parents[0][1] == term
    # The ultimate parent of every antibiotic is the term "antibiotic molecule"
    # whose parent is the term "process or components of antibiotic biology or chemistry"
    if (len(parents) < 3 or
        all_terms[parents[-1][1]]["name"] !=
            "process or component of antibiotic biology or chemistry" or
        all_terms[parents[-2][1]]["name"] !=
             "antibiotic molecule"):
        # This should never happen.
        result = "unsuppported__{}_{}".format(all_terms[term]["name"], term)
    else:
        result = all_terms[term]["name"]
        # Traverse superclasses up to but not incl. "process or component..."
        for p_level, p_term in parents[1:-2]:
            p_name = all_terms[p_term]["name"]
            # The following classifications are not useful in tags
            if p_name in ("antibiotic mixture", "antibiotic without defined classification"):
                continue
            result = p_name + "__" + result
    result += "__"
    result = result.replace(" ", "_")
    if use_plural:
        result = result.replace("_antibiotic__", "s__")
    result = result.strip("_")
    if add_aro:
        result += "_{}".format(term).replace(":", "_")
    return (result, level)


def fully_qualify(all_terms, resistances):
    return dict(antibiotic_fully_qualified_name(all_terms, r, l)
                for r, l in resistances.items())


def resistance_depths(all_terms, *genes):
    resistances = {}
    for g in genes:
        parents = sorted(chase(all_terms, g, is_sort_of))
        for p_level, p_term in parents:
            for r in conferred_resistances(all_terms[p_term]):
                if r not in resistances or resistances[r] > p_level:
                    resistances[r] = p_level
    return resistances


def fqn_prefixes(fully_qualified_name):
    parts = fully_qualified_name.split("__")
    for j in range(1, len(parts) + 1):
        yield "__".join(parts[:j])


def prefix_complete(fqn_resistances):
    infinity = 999
    if not fqn_resistances:
        return {}
    result = {}
    for fqn, depth  in fqn_resistances.items():
        result[fqn] = min(depth, result.get(fqn, depth))
        for prefix in fqn_prefixes(fqn):
            if prefix != fqn:
                result[prefix] = min(infinity, result.get(prefix, infinity))
    return result


def prefix_remove(fqn_resistances):
    # discard classes if there are specific instances of the class
    # unused
    resistances = set(fqn_resistances)
    rs = sorted(resistances)
    for i in range(1, len(resistances)):
        if rs[i].startswith(rs[i-1]):
            resistances.discard(rs[i-1])
    return resistances


def glob_aro_genes():
    files = glob.glob('generated_files/under_version_control/genes/*.fasta')
    all_genes = {}
    for f in files:
        gene = f.rsplit("/", 1)[1].rsplit(".", 1)[0]
        aro = gene.rsplit("__", 1)[1]
        if "ARO_" in aro:
            aro = aro.replace("ARO_", "ARO:")
            all_genes[aro] = gene
    print("Found {} gene aros.".format(len(all_genes)))
    return all_genes


def infer_resistances(all_terms, all_genes):
    results = {}
    for gene_aro in all_genes.keys():
        resistances = resistance_depths(all_terms, gene_aro)
        resistances = prefix_complete(fully_qualify(all_terms, resistances))
        results[gene_aro] = sorted(resistances.items())
    return results


def cross_reference_antibiotic_under_gene(all_terms, all_genes, antibiotics_for_gene):
    for gene_name, gene_aro in sorted((v, k) for k, v in all_genes.items()):
        print()
        print(gene_name)
        print(
            "    " +
            "\n    ".join(
                str(r[1]).rjust(3, " ") + "  " + str(r[0])
                for r in antibiotics_for_gene.get(gene_aro, [])
            )
        )


def inverse_index(antibiotics_for_gene):
    genes_for_antibiotic = defaultdict(list)
    for gene_aro in sorted(antibiotics_for_gene.keys()):
        for antibiotic_fqn, _ in antibiotics_for_gene[gene_aro]:
            genes_for_antibiotic[antibiotic_fqn].append(gene_aro)
    return genes_for_antibiotic


def  cross_reference_gene_under_antibiotic(all_terms, all_genes, genes_for_antibiotic):
    for antibiotic_fqn in sorted(genes_for_antibiotic.keys()):
        print()
        print(str(antibiotic_fqn).ljust(80),
              str(len(genes_for_antibiotic[antibiotic_fqn])).rjust(5))
        print(
            "    " +
            "\n    ".join(
                all_genes[gene_aro]
                for gene_aro in genes_for_antibiotic[antibiotic_fqn]
            )
        )


def list_antibiotics_for(all_terms, all_genes):
    print()
    print("FINDING ANTIBIOTICS FOR {} GENE(S)".format(len(all_genes)))
    print()
    print(
        "    " +
        "\n    ".join(gene for gene in sorted(all_genes.values(), key=str.lower))
    )
    resistances = resistance_depths(all_terms, *all_genes.keys())
    resistances = prefix_complete(fully_qualify(all_terms, resistances))
    print()
    print("FOUND {} ANTIBIOTICS".format(len(resistances)))
    print()
    print(
        "    " +
        "\n    ".join(
            str(resistances[fqn]).rjust(3, " ") + "  " + str(fqn)
            for fqn in sorted(resistances.keys())
        )
    )
    return 0


def list_genes_for(antibiotics, brief):
    all_genes = glob_aro_genes()
    all_terms = parse_aro_obo("inputs/card/aro.obo")
    antibiotics_for_gene = infer_resistances(all_terms, all_genes)
    hits = set()
    canonical = set()
    print()
    print("INPUT {} ANTIBIOTIC(S)".format(len(antibiotics)))
    print()
    print(
        "    " +
        "\n    ".join(gene for gene in sorted(antibiotics))
    )
    subset_genes = {}
    subset_resistances = {}
    for gene_aro, gene_antibiotics in antibiotics_for_gene.items():
        for needle in antibiotics:
            for hay in gene_antibiotics:
                # substring search!
                if needle in hay[0]:
                    hits.add(needle)
                    canonical.add(hay[0])
                    subset_resistances[gene_aro] = gene_antibiotics
                    subset_genes[gene_aro] = all_genes[gene_aro]
    print()
    print("MATCHING AGAINST {} CANONICALLY NAMED ANTIBIOTIC(S)".format(len(canonical)))
    print()
    print(
        "    " +
        "\n    ".join(ca for ca in sorted(canonical))
    )
    print()
    print("FOUND {} GENES".format(len(subset_genes)))
    if not brief:
        cross_reference_antibiotic_under_gene(all_terms, subset_genes, subset_resistances, brief)
    if brief:
        print("\n    " + "\n    ".join(sorted(subset_genes.values(), key=str.lower)))
    if hits != set(antibiotics):
        missing = set(antibiotics) - hits
        print()
        print("WARNING:  Found no genes covering {} antibiotics:".format(len(missing)))
        print("\n    " + "\n    ".join(m for m in missing))
        print()
        return -2
    return 0



def input_genes():
    all_genes = {}
    txt = sys.stdin.read()
    for gene in txt.split():
        if not gene:
            continue
        if "__" in gene:
            aro = gene.rsplit("__", 1)[1]
        else:
            aro = gene
        aro = aro.replace("ARO_", "ARO:")
        if is_aro(aro):
            all_genes[aro] = gene
        else:
            print("WARNING: Could not parse ARO from '{}'".format(gene))
    return all_genes



def input_antibiotics():
    antibiotics = set()
    txt = sys.stdin.read()
    for a in txt.split():
        if not a:
            continue
        a = a.lower().replace(" ", "_")
        if a in antibiotics:
            print("WARNING: Duplicate input '{}'".format(a))
        else:
            antibiotics.add(a)
    return antibiotics



def print_usage(me):
    print("""
         0.  List the genes that confer resistance to some clinically
         important drugs.

         python {me} list_genes_for < inputs/additional/clinically_important_drugs.txt

         1.  Cross reference antibiotics to genes.

         python3 {me} list_antibiotics_by_gene

         For every gene in generated_files/under_version_control/genes,
         list the antibiotics it confers resistance to, along with the
         ontology network path legnth connecting the gene to the
         antibiotic.  Path length 999 represents classes of antibiotics
         that are not connected to the gene in the ontology, but which
         have a constituent antibiotic connected to the gene.


         2.  Cross reference genes to antibiotics.

         python3 {me} list_genes_by_antibiotic

         For every antibiotic, list all genes that confer resistance to
         that antibiotic.  For every class of antibiotics, list all
         genes that confer resistance to any member of the class.


         3.  List all antibiotics covered by any of the given genes.

             python3 {me} list_antibiotics_for < file_with_gene_names.txt
         or
             echo gene_1 gene_2 ... | python3 {me} antibiotics_for

         Input list of gene names or ARO:XXXXXXs.
         Output list of antibiotics those genes confer resistance to.

         For every antibiotic or antibiotic class in the output, show
         the shortest path length in the ontology graph that connects
         that antibiotic or antibiotic class to any input gene.  See
         above for the meaning of path length 999.


         4.  List all genes that confer resistance to the input antibiotics.

             python3 {me} list_genes_for < file_with_antibiotic_names.txt
         or
             echo antibiotic_1 antibiotic_2 ... | python3 {me} genes_for

         Input list of fully qualified antibiotic or antibiotic class names.
         Output list of genes conferring resistance to any of the above.

         For every antibiotic or antibiotic class in the input, show
         the shortest path length in the ontology graph that connects
         that antibiotic or antibiotic class to an input gene.  See
         above for the meaning of path length 999.
    """.format(me=me))


def main():
    if len(sys.argv) == 1:
        print_usage(sys.argv[0])
        return 0
    if sys.argv[1] in ("list_genes_by_antibiotic", "list_antibiotics_by_gene"):
        all_genes = glob_aro_genes()
        all_terms = parse_aro_obo("inputs/card/aro.obo")
        antibiotics_for_gene = infer_resistances(all_terms, all_genes)
        if sys.argv[1] == "list_antibiotics_by_gene":
             print("\nINDEX OF ANTIBIOTICS BY GENE\n")
             cross_reference_antibiotic_under_gene(all_terms, all_genes, antibiotics_for_gene)
        if sys.argv[1] == "list_genes_by_antibiotic":
            print("\nINDEX OF GENES BY ANTIBIOTIC\n")
            genes_for_antibiotic = inverse_index(antibiotics_for_gene)
            cross_reference_gene_under_antibiotic(all_terms, all_genes, genes_for_antibiotic)
        return 0
    if sys.argv[1] == "list_antibiotics_for":
        all_terms = parse_aro_obo("inputs/card/aro.obo")
        all_genes = input_genes()
        return list_antibiotics_for(all_terms, all_genes)
    if sys.argv[1] == "list_genes_for":
        antibiotics = input_antibiotics()
        return list_genes_for(antibiotics, brief=True)
    print_usage(sys.argv[0])
    return -1


if __name__ == "__main__":
    retcode = main()
    sys.exit(retcode)
