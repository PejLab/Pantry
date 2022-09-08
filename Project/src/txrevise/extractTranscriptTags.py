import os
import argparse
import fileinput
import subprocess
import gzip

parser = argparse.ArgumentParser(description = "Extract transcript tags from GTF file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gtf", help = "Path to the GTF file.")
args = parser.parse_args()

#Construct header of the output file
valid_tags = ["CCDS", "basic", "cds_start_NF", "mRNA_start_NF", "cds_end_NF","mRNA_end_NF"]
header = "\t".join(["transcript_id"] + valid_tags)
print(header)

#Iterate over the GTF file
file_handle = gzip.open(args.gtf, "r")
for line in file_handle:
	line = line.decode("utf8").rstrip()
	#Remove headers
	if (line[0:2] == "#!"):
		pass
	else:
		fields = line.split("\t")
		if (fields[2] == "transcript"):
			tags = fields[8]
			tags = tags.rstrip('";')
			tag_fields = tags.split('"; ')

			#Construct tag pairs
			tag_pairs = list()
			for tag_field in tag_fields:
				tag_pairs.append(tag_field.split(' "'))

			#Iterate over pairs to extract tx_id and tags
			true_tags = dict()
			transcript_id = ""
			for tag_pair in tag_pairs:
				if(tag_pair[0] == "transcript_id"):
					transcript_id = tag_pair[1]
				if(tag_pair[0] == "tag"):
					true_tags[tag_pair[1]] = 1
			if len(true_tags) > 0:
    			#Convert tags to a binary vector
				tag_vector = list()
				for tag in valid_tags:
					if tag in true_tags:
						tag_vector.append("1")
					else:
						tag_vector.append("0")
				print("\t".join([transcript_id, "\t".join(tag_vector)]))

