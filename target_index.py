# targets module for flash library constructor

def offtag(ot):
    return "off_" + str(ot)


def read_all_targets(input_path):
    "Returns a list of the 20-mers in all_targets.txt."
    with open(input_path, "r") as f:
        targets = []
        for line in f:
            if line and line[0].isalpha():
                targets.append(line.strip())
        print("Read {} targets from {}.".format(len(targets), input_path))
        return targets


def read_all_targets_with_cut_sites(input_path):
    "Returns a map of 20-mer => gene => cut site"
    with open(input_path, "r") as f:
        targets = {}
        last_t = None
        for line in f:
            if line and line[0].isalpha():
                assert last_t == None
                last_t = line.strip()
                targets[last_t] = {}
            elif line.strip():
                gene, cut_site = line.strip().split()
                targets[last_t][gene] = int(cut_site)
            else:
                last_t = None
        print("Read {} targets from {}.".format(len(targets), input_path))
        return targets


def read_tagged_targets(input_path):
    "Returns a dict of target => [tag, ...] from parsing input_path."
    with open(input_path, "r") as f:
        tagged_targets = {}
        for line in f:
            if line and line[0].isalpha():
                try:
                    parts = line.strip().split()
                    target = parts[0]
                    tags = parts[1:]
                except:
                    print("'{}'".format(line))
                    raise
                assert len(target) == 20
                tagged_targets[target] = tags
        print("Read {} targets from {}.".format(len(tagged_targets), input_path))
        return tagged_targets
