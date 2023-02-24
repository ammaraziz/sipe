from Bio.Blast import NCBIXML
from collections import defaultdict
from csv import DictWriter
import argparse

parser = argparse.ArgumentParser(description='Optional app description')
parser.add_argument('--input-xml',help='path to xml file')
parser.add_argument('--output-dir',help='output directory')
parser.add_argument('--sample',help='sample name')
args = parser.parse_args()

output_all = args.output_dir + "/" + args.sample + "_all.sorted.tsv"

result_handle = open(args.input_xml)
blast_records = NCBIXML.parse(result_handle)

reads_by_entero = defaultdict(list)
blast_results = []
for record in blast_records:
	if not record.alignments:
		continue
	read_name = record.query
	read_name_search = record.query.split(' ')[0]
	top_hit = record.alignments[0].hit_def
	serotype = top_hit.split("_")[0]

	query_start = record.alignments[0].hsps[0].query_start
	query_end = record.alignments[0].hsps[0].query_end
	subject_start = record.alignments[0].hsps[0].sbjct_start
	subject_end = record.alignments[0].hsps[0].sbjct_end
	score = record.alignments[0].hsps[0].score

	reads_by_entero[serotype].append(read_name_search)
	
	blast_results.append(
			{
			'read_name' : read_name,
			'read_name_search' : read_name_search,
			'top_hit' : top_hit,
			'query_start' : query_start,
			'query_end' : query_end,
			'subject_start' : subject_start,
			'subject_end' : subject_end,
			'score' : score
			}
			)

# write out all results
keys = ['read_name', 'read_name_search', 'top_hit', 'query_start', 'query_end', 'subject_start', 'subject_end', 'score']
with open(output_all, 'w', newline='') as outtsv:
	dict_writer = DictWriter(outtsv, fieldnames=keys, delimiter = "\t")
	dict_writer.writeheader()
	dict_writer.writerows(blast_results)


# write out read names into individual files
# only hits with > 500 reads
for serotype in reads_by_entero:
	if len(reads_by_entero[serotype]) > 500:
		with open(args.output_dir + f"/{args.sample}_{serotype}_reads.txt", 'w') as file:
			for item in reads_by_entero[serotype]:
				file.write(f'{item}\n')