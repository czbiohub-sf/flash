# Optimizer Hints

The optimizer attempts to generate an optimal library for a component given
certain constraints:

* Coverage of SNPs
* Minimum fragment sizes
* Minimum coverage of each gene

If those cannot be satisfied, then then the component will
not have a solution.

Such an unsolved component can be identified since its library is empty.

```
for c in components:
    if not c.library:
        print(c.name)
```

To investigate the cause of a component not being solved,
you can display the targets available to the solver.

```
for g in comp.genes:
  g.display_gene_targets()
```

If there are very few targets (eg because of E. Coli offtarget filtering) you may want to ignore the gene. If there is a SNP at the beginning or end of the gene that cannot be captured, you may want to add padding.

If, as in the case of one component of _dfr_ genes, the genes are simply too short, you can rerun the optimzier with a different MIN_FRAGMENT_SIZE.

```
unsolved_component = 'dfrB1__AY139601__ARO_3002864-comp-12'

fname = 'generated_files/untracked/components/%s.pickle' % unsolved_component
with open(fname,'rb') as f:
    dfr = pickle.load(f)
m, genes, library, impossible_snps = optimizer.optimize(dfr.genes, MIN_FRAGMENT_SIZE=178)
dfr.set_library(library)
with open(fname, 'wb') as fp:
    pickle.dump(dfr, fp)
```
