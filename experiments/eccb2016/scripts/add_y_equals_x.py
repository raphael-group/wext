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

                                    
