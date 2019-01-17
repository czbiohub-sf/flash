#!/usr/bin/env python3

import sys, time
from multiprocessing import Pool, cpu_count
import flash, target_index, build


def filter_targets(all_targets, off_targets):
    result = {}
    for target in all_targets:
        poor_structure_reasons = flash.poor_structure(target, need_explanation=True)
        if poor_structure_reasons or off_targets.get(target):
            result[target] = (
                [target_index.offtag(ot) for ot in off_targets.get(target, [])] +
                poor_structure_reasons
            )
        else:
            result[target] = ["ok"]
    return result


def collated(subprocess_results):
    num_ok = 0
    filtered_targets = {}
    for ft_part in subprocess_results:
        for k, v in ft_part.items():
            if v == ["ok"]:
                num_ok += 1
            filtered_targets[k] = v
    num_all = len(filtered_targets)
    print("Only {} out of {} targets ({:3.1f}%) survived all filters."
          .format(num_ok, num_all, (100.0*num_ok/max(1, num_all))))
    return filtered_targets


def output(filtered_targets_output, filtered_targets, off_targets):
    with open(filtered_targets_output, "w") as f:
        for k in sorted(filtered_targets.keys()):
            v = filtered_targets[k]
            assert k not in off_targets or v != ["ok"]
            f.write(k + "    " + "    ".join(v) + "\n\n")


def main():
    t = time.time()
    all_targets = target_index.read_all_targets(build.all_targets_path)
    off_targets = target_index.read_tagged_targets(build.off_targets_path)
    assert set(off_targets) <= set(all_targets)
    num_subprocs = int(1.5 * cpu_count())
    with Pool(num_subprocs) as p:
        size = int((len(all_targets) + num_subprocs - 1) / num_subprocs)
        # Pool.starmap requires Python3.3+
        subprocess_results = p.starmap(filter_targets,
            [(all_targets[index : index + size], off_targets)
             for index in range(0, len(all_targets), size)])
    filtered_targets = collated(subprocess_results)
    assert set(filtered_targets.keys()) == set(all_targets)
    output(build.filtered_targets_path, filtered_targets, off_targets)
    print("Completed structure filtering in {:3.1f} seconds.".format(time.time() - t))


if __name__ == "__main__":
    retcode = main()
    sys.exit(retcode)
