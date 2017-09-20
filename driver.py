import os
from subprocess import call

def without_extension(filename):
	return os.path.splitext(filename)[0]

def get_hit_sites(path):
	TA = 0; Reads = 0; HitSites = 0

	input_file = open(path)
	for line in input_file:
	    split = line.split()
	    if len(split) > 3:
	        if split[3][0].isdigit() == True:
	            TA += 1
	            Reads += float(split[3])
	            if float(split[3]) > 0: HitSites += 1
	input_file.close()
	return HitSites

inputs = [path for path in os.listdir('.') if path.endswith('.tabular')]
igv_files = []
for data in inputs:
	name = without_extension(data)
	map_to_ta_plus = name + 'plus.igv'
	map_to_ta_minus = name + 'minus.igv'
	run_command = "python %s Staph_TA.txt windows_table_270_135.txt " + data
	os.system(run_command % "igv_staph_saouhsc_plus.py" + " > " + map_to_ta_plus)
	os.system(run_command % "igv_staph_saouhsc_minus.py" + " > " + map_to_ta_minus)
	igv_files.append(map_to_ta_plus)
	igv_files.append(map_to_ta_minus)


def get_promoter_analysis_name(filename):
	name = without_extension(filename)
	# names are usually like moe16_tufminus
	return name.split('_')[1]

# Normalize!
# We need to group the igv files by promoter analysis name so we can
# pass the control and experiment to normalize
group_labels = set(map(get_promoter_analysis_name, igv_files))
grouped_items = []
for group_label in group_labels:
	group = []
	for igv_file in igv_files:
		if group_label in igv_file: group.append((igv_file, get_hit_sites(igv_file)))
	grouped_items.append(group)

# We want to sort by the hit sit count. We already stored this count in the tuple:
# [[(ctrl_path1, hit_site), (expr_path1, hit_site)], [(ctrl_path2, hit_site), (expr_path2, hit_site)]]
normalized_args = [sorted(group, key=lambda item: item[1], reverse=True) for group in grouped_items]

# Extract only path from each tuple in the groups
paths = [(group[0][0], group[1][0]) for group in normalized_args]

# We will generate a list of the inputs for each pair to promoter_analysis
# will have to repair plus and minus later
# TODO make this not suck
promoter_analysis_inputs = []
for arg1, arg2 in paths:
	# output arg1 with norm before extension
	destination = without_extension(arg1) + ".norm.igv"
	os.system("python normalization.py %s %s > %s" % (arg1, arg2, destination))
	if arg1.startswith('37ctrl_'):
		promoter_analysis_inputs.append((destination, arg2))
	else:
		promoter_analysis_inputs.append((arg2, destination))

promoter_args = []
# We need to repair plus with minus before calling promoter analysis
plus_strands = [strand for strand in promoter_analysis_inputs if "plus" in strand[0]]
minus_strands = [strand for strand in promoter_analysis_inputs if "plus" not in strand[0]]
for plus_strand in plus_strands:
	promoter_name = str.replace(get_promoter_analysis_name(str.replace(plus_strand[0], '.norm', '')), 'plus', '')
	for minus_strand in minus_strands:
		if promoter_name in minus_strand[0]:
			promoter_args.append(plus_strand + minus_strand)

for arglist in promoter_args:
	output_name = str.replace(str.replace(str.replace(arglist[0].split('_')[1], 'norm.igv', ''), '.igv', ''), 'plus', '') + ".promoter.out"
	command = "python promoteranalysis.py %s %s %s %s > " + output_name 
	os.system(command % arglist)
