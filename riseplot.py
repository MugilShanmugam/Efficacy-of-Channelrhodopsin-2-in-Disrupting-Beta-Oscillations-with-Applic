# -*- coding: utf-8 -*-
"""
riseplot.py by Marianne Bezaire August 2022
plot NEURON shapes using matplotlib.pyplot

USAGE:
import riseplot

riseplot.plot(soma)
riseplot.plot([cell1.soma, cell2.soma, cell3.soma])
riseplot.plot(soma, ptype)

NOTES:
Requires a reference to the soma of each neuron
to be included in the plot. If only one cell, the
soma section can be passed in as the first argument.
If more than one cell must be plotted, first create
a list of cell soma sections and pass that in.

Hint: use a list comprehension to build a list of
multiple cells' soma sections:
    
somalist = [c.soma for c in cells]
riseplot.plot(somalist)

The above works if you have a list of cells called *cells*
and the cell objects have a property called soma that
references the cell body section.

ptype optional argument:
-1: 3d-ish plot using matplotlib.pyplot via NEURON's functionality
 0: skeleton plot
 1: plot in XY plane respecting diameters
 2: 3d plot (not yet implemented, use -1 version instead)
"""

from neuron import h
from matplotlib import pyplot
import math
import matplotlib.patches as patches

xpts=[]
ypts=[]
xl = [0, 1]
yl = [0, 1]
colors = iter(['g','r','b','o','k','p']*100)

def main(): # example usage
    soma = h.Section(name='soma')
    dend = []
    dend.append(h.Section(name='dend[0]'))
    dend.append(h.Section(name='dend[1]'))
    
    soma.v = -65
    soma.L = 10
    soma.diam = 10
    
    for d in dend:
        d.v = -60
        d.L = 50
        d.diam = 1
        d.connect(soma(1))    
        
    plot(soma, ptype = 1)


def plot(sec, ptype = 1):
    global xl, yl
    h.define_shape()
    
    if ptype < 0:
        # plot using matplotlib, without diameters
        
        ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
        ps.plot(pyplot)
        pyplot.show()
        
    else:          
        # or take into account diameters using this custom code
        fig = pyplot.figure()
        ax = fig.add_subplot(111, aspect='equal')


        plot3d(sec, view = ptype, ax = ax)
        xl = [(xl[0]-1)*1.1, xl[1]*1.1]
        yl = [(yl[0]-1)*1.1, yl[1]*1.1]
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        pyplot.show()

def plot3d(sec = None, col = 'k', view = 1, ax = None):
    global xpts, ypts, xl, yl, colors

    if ax == None:
        fig = pyplot.figure()
        ax = fig.add_subplot(111, aspect='equal')
    
    if type(sec) == list:
        for s in sec:
            plot3d(s, ax = ax, view = view)
        return
    
    if sec == None:
        sec = h.cas()
    for x in range(sec.n3d() - 1):
        # print(soma.x3d(x), soma.y3d(x), soma.z3d(x), soma.diam3d(x))

        if view == 2:
        # parametric cylindrical 3d
            print("not implemented")
        
        elif view == 1:
            # proportional 2D projection in XY plane
            xseg = [sec.x3d(x), sec.x3d(x+1)]
            yseg = [sec.y3d(x), sec.y3d(x+1)]
            diams = [sec.diam3d(x), sec.diam3d(x+1)]
            
            dx = xseg[1] - xseg[0]
            dy = yseg[1] - yseg[0]
            
            norm1 = [-dy, dx]
            norm2 = [dy, -dx]
            
            scale = diams[0]/(2*math.sqrt(norm1[0]**2 + norm1[1]**2))
            
            norm1 = [x*scale for x in norm1]
            norm2 = [x*scale for x in norm2]
            
            xpts = [xseg[0] + norm1[0], xseg[1] + norm1[0], xseg[1] + norm2[0], xseg[0] + norm2[0]]
            ypts = [yseg[0] + norm1[1], yseg[1] + norm1[1], yseg[1] + norm2[1], yseg[0] + norm2[1]]
            
    
            xl = [min([xl[0], min(xpts)]), max([xl[1], max(xpts)])]
            yl = [min([yl[0], min(ypts)]), max([yl[1], max(ypts)])]

            ax.add_patch(patches.Polygon(xy=list(zip(xpts,ypts)), facecolor=col)) # , fill=False
    
        else: # view == 0
            # scaled skeleton plot
            pyplot.plot([sec.x3d(x), sec.x3d(x+1)],[sec.y3d(x), sec.y3d(x+1)], color=col, linewidth=sec.diam3d(x))
            xl = [min([xl[0], sec.x3d(x), sec.x3d(x+1)]), max([xl[1], sec.x3d(x), sec.x3d(x+1)])]
            yl = [min([yl[0], sec.y3d(x), sec.y3d(x+1)]), max([yl[1], sec.y3d(x), sec.y3d(x+1)])]

    tcol = next(colors)
    for s in h.SectionRef(sec).child:
        plot3d(sec = s, col = 'k', view = view, ax = ax)
        # for a more colorful plot, use col = tcol


#%%
if __name__ == "__main__":
    main() # run example code