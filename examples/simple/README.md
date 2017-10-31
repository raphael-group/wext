# Minimal working example

This simple example applies the WR-exclusivity test on pairs of genes in a small toy dataset.  It uses the saddlepoint approximation.

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

The results (`example_results-k2.tsv`) should be very similar to the following table:
```
#Gene set	WRE (Saddlepoint) P-value	WRE (Saddlepoint) FDR	WRE (Saddlepoint) Runtime	T	Z	t00	t01	t10	t11
a, b	0.000127177896943	0.00473591769391	0.00224995613098	12	0	2	6	6	0
a, g	0.0590504141841	0.732983503396	0.00198483467102	8	0	6	6	2	0
b, e	0.0590504141841	0.732983503396	0.00172400474548	8	0	6	6	2	0
b, d	0.218284206842	1.0	0.00166797637939	7	3	4	3	4	3
b, c	0.218284206842	1.0	0.00179195404053	7	3	4	3	4	3
e, g	0.332644664728	1.0	0.00157308578491	4	0	10	2	2	0
c, e	0.344539115446	1.0	0.0016610622406	7	1	6	6	1	1
d, g	0.344539115446	1.0	0.0016450881958	7	1	6	6	1	1
c, g	0.344539115446	1.0	0.00165319442749	7	1	6	6	1	1
d, e	0.344539115446	1.0	0.00164103507996	7	1	6	6	1	1
a, d	0.629174565866	1.0	0.00207710266113	5	4	5	2	3	4
a, c	0.629174565866	1.0	0.00289821624756	5	4	5	2	3	4
```
