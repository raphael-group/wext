#!/usr/bin/env python
import numpy as np

# Add a y=x line to the given matplotlib axis
def add_y_equals_x(ax, c='k', line_style='--', alpha=0.75):
    # shamelessly stolen from http://goo.gl/9ZttXZ
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, line_style, c=c, alpha=alpha, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

def aligned_plaintext_table(table, sep='\t', spaces=2):
    """
    Create and return an aligned plaintext table.

    Arguments:
    table: DSV table; string
    sep: delimiter for DSV table; string; optional with tab as default
    spaces: minimum number of spaces between columns; integer; optional with 2 as default
    """
    # Separate entries.
    rows = [[entry.strip() for entry in row.split(sep)] for row in table.strip().split('\n')]

    # Find numbers of rows and columns.
    m = len(rows)
    lengths = map(len, rows)
    n = max(lengths)

    # Pad rows with a deficient number of columns.
    entries = [[rows[i][j] if j<lengths[i] else '' for j in range(n)] for i in range(m)]

    # Find column widths.
    sizes = [max(len(entries[i][j]) for i in range(m)) for j in range(n)]

    # Return results.
    return '\n'.join([''.join([entries[i][j].rjust(sizes[j]+spaces) for j in range(n)]).rstrip() for i in range(m)])

def rank(a, reverse=False, ties=2):
    """
    Find the ranks of the elements of a.

    Arguments:
    a: entries; may be np.ndarray, list, or typle
    reverse: reverse ranks; optional
    ties: options for breaking ties; optional; see below examples

    Example usage:
    In [1]: a = np.array([3.142, 2.718, 1.414, 1.618])
    In [2]: rank(a)
    Out[2]: array([3, 2, 0, 1])
    In [3]: rank(a, reverse=True)
    Out[3]: array([0, 1, 3, 2])
    In [4]: b = [7, 8, 9, 7, 8, 9]
    In [5]: rank(b, ties=0)
    Out[5]: [0, 2, 4, 1, 3, 5]
    In [6]: rank(b, ties=1)
    Out[6]: [0, 1, 2, 0, 1, 2]
    In [7]: rank(b, ties=2)
    Out[7]: [0, 2, 4, 0, 2, 4]
    """
    import numpy as np

    if type(a)==np.ndarray:
        x = a.flatten()
    else:
        x = a

    n = len(x)
    y = np.argsort(x)

    if ties==0:
        z = np.argsort(y)
        j = n-1
    elif ties==1 :
        z = np.zeros(n, dtype=y.dtype)
        j = 0
        for i in xrange(1, n):
            if x[y[i]]!=x[y[i-1]]:
                j += 1
            z[y[i]] = j
    elif ties==2:
        z = np.zeros(n, dtype=y.dtype)
        j = 0
        for i in xrange(1, n):
            if x[y[i]]!=x[y[i-1]]:
                j = i
            z[y[i]] = j
    else:
        raise NotImplementedError('Choice {} for ties has not been implemented.'.format(ties))

    if reverse:
        z = j-z

    if type(a)==np.ndarray:
        z = z.reshape(a.shape)
    elif type(a)==list:
        z = list(z)
    elif type(a)==tuple:
        z = tuple(z)

    return z
