# !/usr/bin/env python3

# from https://onestopdataanalysis.com/depth-plot/
 
from collections import defaultdict
import argparse
import toyplot
import toyplot.pdf 


parser = argparse.ArgumentParser(description='plot depth')
parser.add_argument('--input',help='path to samtools depth file')
parser.add_argument('--output',help='path and name for output')
args = parser.parse_args()

def parse_depth(depth_input):
	"""Parse depth file.

	Args:
		depth_input (str): Path to depth file. 
	Returns:
		list: List with depth.
 
	"""
	# written to handle multiple references
	# depths = defaultdict(dict)
	# with open(depth_input) as depth_object:
	# 	for row in depth_object:
	# 		genome_id, position, depth_count = row.rstrip('\n').split("\t")
	# 		depths[genome_id].update({position : depth_count})
	# return(depths)
	depth = []
	references = set()

	with open(depth_input) as depth_object:
		for row in depth_object:
			genome_id, position, depth_count = row.split()

			references.add(genome_id)

			if len(references) > 1:
				raise Exception(' This script only handles one genome - contig.')
 
			depth.append(int(depth_count))

	return depth



def plot_depth(depth_data, plot_title):
	
	'''
	Plot genome Depth across genome.
	'''

	canvas = toyplot.Canvas(width=1200, height=900)
	axes = canvas.cartesian(label=plot_title, xlabel="Depth", ylabel="Genomic Position")
	mark = axes.plot(depth_data)
	return(canvas)


depth_data = parse_depth(args.input)
plot_title = args.output.split('.')[0]
canvas = plot_depth(depth_data, plot_title)
toyplot.pdf.render(canvas, args.output)
print("Plotting Done")