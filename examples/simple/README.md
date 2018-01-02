# Toy problem example

The commands in `sh commands.sh` illustrate WExT on a small toy problem.  If the dependencies for WExT are correctly installed, then the `commands.sh` script should run within several seconds on your computer.  These commands have three major parts.

1. Process mutation data [`process_mutations.py`](https://github.com/raphael-group/wext/wiki/Process-mutations)
2. Compute mutation weights/probabilities [`compute_mutation_probabilities.py`](https://github.com/raphael-group/wext/wiki/Compute-mutation-probabilities).
3. Find gene sets with significant patterns of mutually exclusive or co-occurring mutations [`find_sets.py`](https://github.com/raphael-group/wext/wiki/Find-sets).

This example applies the WR-exact test for mutual exclusivity and co-occurrence on pairs of genes in a small toy dataset using a saddlepoint approximation of each test.
The WExT [wiki](https://github.com/raphael-group/wext/wiki) describes the details of the above tests and the arguments for the above scripts in more detail.

## Input

The dataset has 8 genes, `a`, `b`, ..., `h`, and 14 patients, `1`, `2`, ..., `14`, where all mutations in genes `a` and `b` are mutually exclusive and all mutations in genes `c` and `d` are co-occurring.
```
1	a	c	d	e	f
2	b	c	d	g	h
3	a	c	d
4	b	c	d
5	a	c	d
6	b	c	d
7	a	c	d
8	b
9	a
10	b
11	a	e
12	b	g
13	c   d
14  c   d
```

## Output

The mutual exclusivity results (`exclusivity_results-k2.tsv`) should be very similar to the following table, which shows that only genes `a` and `b` have mutually exclusive mutations with the (default) FDR < 0.5.
```
#Gene set	WRE (Saddlepoint) P-value	WRE (Saddlepoint) FDR	WRE (Saddlepoint) Runtime	T	Z	t00	t01	t10	t11
a, b	0.000246182027855	0.0101775930969	0.00236916542053	12	0	2	6	6	0
```

The co-occurrence results (`any-co-occurrence_results-k2.tsv` and `all-co-occurrence_results-k2.tsv`) should be very similar to the following table, which shows that only genes `c` and `d` have co-occurring mutations with the (default) FDR < 0.5.
```
#Gene set	WRE (Saddlepoint) P-value	WRE (Saddlepoint) FDR	WRE (Saddlepoint) Runtime	T	Z	t00	t01	t10	t11
c, d	0.0037414655014	0.186225589558	0.00187397003174	0	9	5	0	0	9
```
