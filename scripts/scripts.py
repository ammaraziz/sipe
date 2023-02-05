# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# from https://onestopdataanalysis.com/depth-plot/
 
from collections import defaultdict
import toyplot
import toyplot.pdf 
 
def parse_depth(depth_input):
    """Parse depth file.
 
    Args:
        depth_input (str): Path to depth file. 
    Returns:
        list: List with depth.
 
    """

    depths = defaultdict(dict)
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.rstrip('\n').split("\t")
            depths[genome_id].update({position : depth_count})
    return(depths)
 
def plot_depth(y, output_name, plot_title, normalize=False):
    
    '''
    Plot genome Depth across genome.
 
    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).
    '''

 
    y_label = "Normalized Depth" if normalize else "Depth"
    data = [xx / max(data) for xx in data] if normalize else data
 
    canvas = toyplot.Canvas(width=1200, height=900)
    axes = canvas.cartesian()
    mark = axes.plot(y)

    print("Done :)")