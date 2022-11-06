#!/usr/local/bin/python3

#import modules
import os, subprocess, shutil, sys, collections, argparse, glob
import pandas as pd
import proseqfunctions as psf

########## Parsing user arguments ##########
############################################

#use argparse to allow user to parse command line arguments for the protein family and taxonomic group 
#also generates help and usage messages with a summary of the programme and to allow user to see what arguments can be parsed
parser = argparse.ArgumentParser(description='Proseq.py performs protein sequence analysis on a user defined subset of sequences extracted from NCBI databases. This subset is defined by the protein family and the taxonomic group.')
parser.add_argument('--protein_family', '-pf', type=str, help='The protein family you would like to search for.')
parser.add_argument('--taxonomic_group', '-tg', type=str, help='The taxonomic group you would like to search for.')
#parse any arguments to the variables protein_family and taxonomic_group
args = parser.parse_args()
protein_family=args.protein_family
taxonomic_group=args.taxonomic_group

############################################
############################################


########## Welcome message ##########
print("Welcome to Proseq V1.0!\n\n")


########## Determine directory to store results ##########
##########################################################

#loop indefinitely
while True:
	#ask the user for a path to store the results and store the input as a variable
	results_path = input("Please specify the path to where you would like to store your results. \n\n\t>")
	#check the input is actually a directory
	if os.path.isdir(results_path):
		#if so, change the directory specified
		os.chdir(results_path)
		#break the loop
		break
	#else the input is invalid
	else:
		#inform the user and revert to the start of the loop
		print("Path not found. Please try again.")
		continue

#use the clear_directory function from proseqfunctions.py to clear existing results directories if they exist
psf.clear_directory('multiple_sequence_alignment')
psf.clear_directory('query_files')
psf.clear_directory('motif_searches')
psf.clear_directory('isoelectric_point_results')
psf.clear_directory('transmembrane_map')
psf.clear_directory('single_sequences')
psf.clear_directory('conservation_plot')
psf.clear_directory('hydropathy_plot')


##########################################################
##########################################################



########## Determine the protein family and taxonomic group ##########
######################################################################

#if no argument parsed to protein_family, ask the user to input a protein family and assign to protein_family
if protein_family is None:
	while True:
		protein_family = input("What protein family would you like to search for?\n\n\t> ")
		#if input is empty, inform the user and continue with the loop, asking for input again
		if not protein_family:
			print("No input detected.")
			continue
		#else break from the loop
		else:
			break

#if no argument parsed to taxonomic_group, ask the user to input a protein family and assign to taxonomic_group
if taxonomic_group is None:
	#loop indefinitely
	while True:
		taxonomic_group = input("What taxonomic group would you like to search for?\n\n\t> ")
		#if input is empty, inform the user and continue with the loop, asking for input again
		if not taxonomic_group:
			print("No input detected.")
			continue
		#else break from the loop
		else:
			break

######################################################################
######################################################################


########## Search for sequences ##########
##########################################

os.mkdir("query_files")
#tell the user that the specified sequences are being searched for
print("Searching for sequences...\n")
#store a bash command in a variable which uses esearch to search the protein database for sequences matching the users specified protein family and taxonomic group
#the command pipes search results to efetch to fetch the accessions of the sequences and stores them in a file
fetch_accessions_command = "esearch -db protein -query '{}[protein name] AND {}[organism] NOT partial' 2>/dev/null | efetch -format acc > query_files/query.acc".format(protein_family, taxonomic_group)
#use subprocess.call to run the command in the linux environment
subprocess.call(fetch_accessions_command, shell=True)

##########################################
##########################################

########## How many sequences were found? ##########
####################################################

#open the accessions file and count the number of accessions to determine the total number of sequences found
with open('query_files/query.acc', 'r') as accessions:
	sequence_count=len(accessions.readlines())
#if no sequences were found, tell the user and exit the programme
if sequence_count == 0:
	print("0 sequences were found using the protein family {} and the taxonomic group: {}.\nPlease try again using a different query.".format(protein_family, taxonomic_group))
	sys.exit("Exiting programme...")
#else if only one sequence was found, tell the user that this is not enough for later analyses and exit the programme
elif sequence_count == 1:
	print("Only 1 sequence was found using the protein family {} and the taxonomic group: {}.\nAt least 2 sequences are required for subsequent analyses.\nPlease try again using a different query.".format(protein_family, taxonomic_group))
	sys.exit("Exiting programme...")
#else tell the user how many sequences were found and continue
else:
	print("A total of {} sequences matching your query were found!\n".format(sequence_count))

####################################################
####################################################

########## User decision - user sequence set of over 1000 or limit sequences to 1000? (conditional) ##########
##############################################################################################################

#if the sequence count is over 1000
if sequence_count > 1000:
	#loop indefinitely
	while True:
		#warn the user that a large starting sequence set is not advised and ask whether they would like to fetch only the first 1000 sequences
		limit_sequences = input("Note: A starting sequence set of over 1000 is not advised. Would you like to fetch the first 1000 sequences only?(y/n)\n\n\t> ")
		#if the user replies with any form of yes
		if limit_sequences in ('Y', 'y', 'yes', 'YES', 'Yes'):
			#confrim the users decision and call the limited_fetch function
			print("Fetching first 1000 sequences out of {}...\n".format(sequence_count))
			psf.limited_fetch(protein_family, taxonomic_group)
			#break the loop
			break
		#else if the user replies with any form of no
		elif limit_sequences in ('N', 'n', 'no', 'NO', 'No'):
			#confirm the users decision and call the standard_fetch function
			print("Fetching all {} sequences...\n".format(sequence_count))
			psf.standard_fetch(protein_family, taxonomic_group)
			#break the loop
			break
		#else the wrong input has been given, inform the user and continue with the loop, asking for input again
		else: 
			print("Wrong input! Please answer either y or n...\n")
			continue
#if sequence count is less than or equal to 1000, tell the user that all sequences will be fetched and call the standard fetch function
else: 
	print("Fetching all {} sequences...\n".format(sequence_count))
	psf.standard_fetch(protein_family, taxonomic_group)

#################################################################################################
#################################################################################################

########## Extracting species and accessions ##########
#######################################################

#open the file with accessions and species
with open('query_files/accessions_species.txt', 'r') as spec_acc:
	#create a list to store accessions and a list to store species
	all_accessions = []
	all_species = []
	species_dict = {}
	#for each line in the file
	for line in spec_acc:
		#extract the species name from the line
		species = line.split("\t")[1].strip("\n")
		#extract the accession from the line
		accession = line.split("\t")[0]
		#add the species name to the species list
		all_species.append(species)
		#add the accession to the accession list
		all_accessions.append(accession)
		#add the accession as a key and the species asa value in the dictionary - to be used for dataframe
		species_dict[accession] = species

#######################################################
#######################################################

########## User decision - continue with dataset given the number of species represented in the sequences? ##########
#####################################################################################################################

#inform the user that the total number of species is being counted
print("Counting the total number of represented species...\n")
#use set to get all non-redundant species in the species list and count the length of the list = total species represented in sequences
total_species = len(set(all_species))
#inform the user of the total number of species
print("There are a total of",total_species,"species represented in the sequences.\n")

#loop indefinitely
while True:
	#ask the user whetehr they would still like to continue with this starting dataset given the number of species represented in the sequences and store the decision in a variable
	continue_decision = input("Would you still like to continue with this dataset?(y/n)\n\n\t> ")
	#if the variable is any form of yes, break from the loop
	if continue_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
		break
	#if the variable is any form of no, exit from the programme
	elif continue_decision in ('N', 'n', 'no', 'NO', 'No'):
		sys.exit("Exiting programme...")
	#else the input is invalid and the user is returned to the start of the loop
	else:
		print("Wrong input! Please answer either y or n.\n")
		continue

#####################################################################################################################
###################################################################################################################


########## Motif analysis ##########
####################################

##patmatmotifs only accepts fasta files with single sequences - need to separate sequences out into indivudal files
print("Starting motif analysis...\n")
#create a counter
count=0
#create a temporary directory to store single sequences
os.mkdir("motif_searches")
#open the fasta file containing all sequences
with open('query_files/query.fasta', 'r') as sequences:
		#for each line in the fasta file
		for line in sequences:
				#if the line starts with > indicating a fasta header
				if line.startswith('>'):
						#add 1 to the count
						count += 1
						#create a filename to store the fasta sequence using the count
						filename = "motif_searches/" + all_accessions[count-1] + ".fasta"
						#open the new file to write to it
						sequence_out = open(filename, "w")
						#write the header line to the file
						sequence_out.write(line)
						#close the file
						sequence_out.close()

				#else the line is a sequence line
				else:
						#create the filename
						filename = "motif_searches/" + all_accessions[count-1] + ".fasta"
						#open the file to append to it
						sequence_out = open(filename, "a")
						#write the sequence line to the end of the file
						sequence_out.write(line)
						#close the file
						sequence_out.close()

#loop indefinitely
while True:
	#ask the user if they would like to display any motif hits to the screen an store the decision in a variable
	motif_decision = input("Would you like to display any motif hits to the screen? (y/n)\n\n\t>")
	#if the variable is any form of yes, assign show_motifs variable to True and break the loop
	if motif_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
		show_motifs = True
		break
	#else if the variable is any form of no, assign show_motifs variable to False and break from the loop
	elif motif_decision in ('N', 'n', 'no', 'NO', 'No'):
		show_motifs = False
		break
	#else the input is invalid, infrom the user and revert to the start of the loop
	else:
		print("Wrong input! Please answer either y or n.\n")
		continue

print("Scanning protein sequences for known motifs...\n")
#call a bash command that loops through each of the individual sequence files and uses patmatmotifs to search for PROSITE motifs in the sequence storing the output in a report file
subprocess.call("for file in motif_searches/*; do patmatmotifs $file motif_searches/$(basename $file .fasta).patmatmotifs -auto ; done;", shell=True)

#create a variable with UNIX filename pattern for motif reports
motif_files = 'motif_searches/*.patmatmotifs'
#create an empty dictionary to store sequence accessions and motifs
motif_dict = {}
#for each filename found with the UNIX pattern
for filename in glob.glob(motif_files):
	#open the file and store as motif_results object
	with open(filename, 'r') as motif_results:
		#create a list to store an motifs found
		motif_list = []
		#for each line in the report file
		for line in motif_results:
			#if the line starts with a string indicating it is the line with the accession value
			if line.startswith("# Sequence:"):
				#split the line by space and retreive the accession storing it as a variable
				accession = (line.split(" "))[2]
			#else if the line starts with a string indicating it is a line containing a motif hit
			elif line.startswith("Motif"):
				#split the line to retreive the motif name and strip the newline character off the end storing it in a variable
				motif = (line.split("= "))[1].strip("\n")
				#if this motif is not already in the motif list, append it to the list 
				if motif not in motif_list:
					motif_list.append(motif)
		#if the show_motifs variable was assigned as True
		if show_motifs:
			#for each motif in motif_list, print the motif along with accession value of the sequence report
			for m in motif_list:
				print(m, "motif found in sequence", accession)
		#store the motif list in the motif_dict dictionary with the accession value of the sequence as a key
		motif_dict[accession] = motif_list
print("Done.\n")
print("Motif reports saved to motif_searches directory\n")


####################################
####################################


########## Calculating IEPs ##########
#####################################

print("Calculating Isoelectric point of sequences...\n")
#create a new directory to store the IEP analysis results
os.mkdir("isoelectric_point_results")
#call a bash comand that calculates the IEP of all sequences in the query.fasta file and creates a report file as an output
subprocess.call("iep query_files/query.fasta -outfile isoelectric_point_results/IEPreport.iep -auto", shell=True)
#create an empty dictionary to store the IEPs of all sequences
IEP_dict = {}
#open the IEP report and store it as a file object
with open("isoelectric_point_results/IEPreport.iep", 'r') as IEP_results:
	#for each line in the file object
	for line in IEP_results:
		#if the line starts with a string indicating it is the line that contains the corresponding sequence accession
		if line.startswith("IEP"):
			#split the line and extract the third field (the accession value) storing it as a string
			accession = (line.split(" "))[2]
		#else if the line starts with a string indicating it is the line containing the corresponding IEP value for this accession
		elif line.startswith("Isoelectric"):
			#split the line and extract the fourth field containing the IEP value and strip the newline character, storing the value as a float variable
			IEP = float((line.split(" "))[3].strip("\n"))
			#add the IEP value to the IEP_dict dictionary with the corresponding accession value as the key
			IEP_dict[accession] = IEP
print("Done.\n")
print("IEP report saved to isoelectric_point_results directory\n")

#####################################
#####################################

########## Summary datframe creation ##########
###############################################
#create indiviudal dataframes from the species_dict, IEP_dict and motif_dict
df_species = pd.DataFrame(species_dict.items(), columns=['Accession', 'Species'])
df_IEP = pd.DataFrame(IEP_dict.items(), columns=['Accession', 'Isoelectric Point'])
df_motif = pd.DataFrame(motif_dict.items(), columns=['Accession', 'Motif(s)'])
#merge the dataframes in steps of two on the Accession column to create a final dataframe
df_final = pd.merge(df_species, df_IEP, on = 'Accession')
df_final = pd.merge(df_final, df_motif, on = 'Accession')

#write the dataframe to a new file stored in the current directory
df_final.to_csv("all_sequences_species_IEP_motif_summary.txt", sep='\t')
print("Motif/IEP summary file saved to specified path.\n")

#######################################
#######################################


########## Selecting sequence subset for multiple sequence alignment and further anlayses ##########
####################################################################################################

#create a new directory to store sequences used for MSA and the alignment files
os.mkdir("multiple_sequence_alignment")
#loop indefinitely
while True:
	#ask user how they would like to select sequences for MSA and store input in variable
	sequences_decision = input("Subsequent analyses require a multiple sequence alignment. It is recommended that you select a sequence subset for this alignment to reduce processing time if your sequence set is large. How would you like to select sequences for this step? (please enter a number)\n\n\t(1) Select from most frequent species (max: 10)\n\n\t(2) Select from all species\n\n\t(3) Use all {} sequences (slow: not recommended for large sequence sets)\n\n\t>".format(sequence_count))
	#create a list to store any selected sequences
	selected_sequences = []
	#if the the input is 1
	if sequences_decision == '1':
		#find the most common species max: 10 and store in dictionary
		most_common = collections.Counter(all_species).most_common(10)
		#create a counter
		count = 0
		#for each species in most_common
		print("Most frequent species:\n")
		for i, species in enumerate(most_common):
			#add 1 to the count
			count += 1
			#display the count, the species name, and the number of sequences of that species
			print('({})'.format(count), most_common[i][0], most_common[i][1])
		#loop indefinitely
		while True:
			#ask the user which of these species they would like to analyse by inputing the corresponding count number(s) (if multiple given separate by space)
			species_decision = input("Please enter the corresponding number(s) of the species you would like to analyse.(please separate by space)\n\n\t>")
			#use the digit check function to make sure the input is either a single digit, or digits separated by spaces and, if true, use the number check function to ensure the numbers given are between 1 and the final count (one of the options)
			if psf.digit_check(species_decision) and psf.number_check(species_decision, count):
				#if both conditions are met, split the input into a list
				species_decision = species_decision.split()
				#create a counter to detemine how many sequences are represented by the species chosen
				selected_sequence_count = 0
				#for each number in the list
				for selection in species_decision:
					#take 1 from the number to get the correct index
					selection = int(selection)-1
					#add the number of sequences for species to the counter
					selected_sequence_count += most_common[selection][1]
					#locate the species name in the most common list using the index and append to the selected sequences list
					selected_sequences.append(most_common[selection][0])
				#if the total sequence count is only 1, inform user and revert to start of the loop
				if selected_sequence_count == 1:
					print("Only 1 sequence selected! Please ensure that you select enough species that the total number of sequences is >1.\n")
					#reset the selected_sequences list
					selected_sequences = []
					continue
				#for each species name in the selected sequences list
				for s in selected_sequences:
					#open a new file to write the selected accessions to
					with open("multiple_sequence_alignment/selected_sequence_accessions.txt", "a") as outfile:
						#use df_final dataframe to extract all entries with the species name 
						selected_entries = df_final[(df_final['Species'] == s)]
						#exract all accessions from these entries and create a list
						selected_accessions = list(selected_entries['Accession'])
						#write the selected accessions to the file
						outfile.write("\n".join(selected_accessions) + "\n")
				#store in a variable a bash command which uses pullseq to extract the selected sequence subset from the fasta file containing all the sequences using the accessions
				pullseq_command = "/localdisk/data/BPSM/Assignment2/pullseq -i query_files/query.fasta -n multiple_sequence_alignment/selected_sequence_accessions.txt > multiple_sequence_alignment/sequence_subset.fasta"
				subprocess.call(pullseq_command, shell=True)
				#immediately break from the loop 
				break
			#if the input does not pass both the conditions, inform the user and continue with the loop, asking for input again
			else:
				print("Wrong input! Please try again.\n")
				continue
		#immediately break from the loop having selected a sequence subset
		break
	#esle if the original input was 2
	elif sequences_decision == '2':
		#count all the species and store in a dictionary
		all_species_counted = collections.Counter(all_species).most_common()
		#create a counter
		count = 0
		#for each of the species in the dictionary
		print("All species and sequences numbers:\n")
		for i, species in enumerate(all_species_counted):
			#add 1 to the count
			count += 1
			#display the count, the species name, and the number of sequences of that species
			print('({})'.format(count), all_species_counted[i][0], all_species_counted[i][1])
		#loop indefinitely
		while True:
			#ask the user which of these species they would like to analyse by inputing the corresponding count number(s) (if multiple given separate by space)
			species_decision = input("Please enter the corresponding number(s) of the species you would like to analyse.(please separate by space)\n\n\t>")
			#use the digit check function to make sure the input is either a single digit, or digits separated by spaces and, if true, use the number check function to ensure the numbers given are between 1 and the final count (one of the options)
			if psf.digit_check(species_decision) and psf.number_check(species_decision, count):
				#if both conditions are met, split the input into a list
				species_decision = species_decision.split()
				#create a counter to detemine how many sequences are represented by the species chosen
				selected_sequence_count = 0
				#for each number in the list
				for selection in species_decision:
					#take 1 from the number to get the correct index
					selection = int(selection) - 1
					#add the number of sequences for species to the counter
					selected_sequence_count += all_species_counted[selection][1]
					#locate the species name in the dictionary using the index and append to the selected sequences list
					selected_sequences.append(all_species_counted[selection][0])
				#if the total sequence count is only 1, inform user and revert to start of the loop
				if selected_sequence_count == 1:
					print("Only 1 sequence selected! Please ensure that you select enough species that the total number of sequences is >1.\n")
					#reset the selected_sequences list
					selected_sequences = []
					continue
				#for each species name in the selected sequences list
				for s in selected_sequences:
					#open a new file to write the selected accessions to
					with open("multiple_sequence_alignment/selected_sequence_accessions.txt", "a") as outfile:
						#use df_final dataframe to extract all entries with the species name 
						selected_entries = df_final[(df_final['Species'] == s)]
						#exract all accessions from these entries and create a list
						selected_accessions = list(selected_entries['Accession'])
						#write the selected accessions to the file
						outfile.write("\n".join(selected_accessions) + "\n")
				#store in a variable a bash command which uses pullseq to extract the selected sequence subset from the fasta file containing all the sequences using the accessions
				pullseq_command = "/localdisk/data/BPSM/Assignment2/pullseq -i query_files/query.fasta -n multiple_sequence_alignment/selected_sequence_accessions.txt > multiple_sequence_alignment/sequence_subset.fasta"
				subprocess.call(pullseq_command, shell=True)
				#immediately break from the loop
				break
			#if the input does not pass both the conditions, inform the user and continue with the loop, asking for input again
			else:
				print("Wrong input! Please try again.\n")
				continue
		#immediately break from the loop having selected a sequence subset
		break
	#else if the original input was 3
	elif sequences_decision == '3':
		#stroe in a variable a bash command that copies the contents of the fasta file containing all sequences to a new file which will be used for conservation/transmembrane analysis
		rename_command = "cp query_files/query.fasta multiple_sequence_alignment/sequence_subset.fasta"
		subprocess.call(rename_command, shell=True)
		#immediately break from the loop having selected to use all the sequences for conservation/transmembrane analysis
		break
	#else the wrong input has been given, inform the user and continue with the loop, asking for the original input again
	else:
		print("Wrong input! Please choose from one of the options.\n")
		continue

#########################################################################
#########################################################################

########## Multiple sequence alignment ########## 
#################################################

print("Aligning selected sequences...\n")
#call a bash command that uses clustalo to peform a multiple sequence alignment on the chosen subset of sequences
subprocess.call("clustalo -i multiple_sequence_alignment/sequence_subset.fasta -o multiple_sequence_alignment/alignment.fasta", shell=True)
#inform the user when the proccess is complete
print("Sequence alignment complete!\n")

#################################################
#################################################


########## Conservation analysis ########## om
###########################################

#create a directory to store the conservation plot
os.mkdir("conservation_plot")
print("Starting conservation plotting...\n")

#obtain the desired window size using the proseqfunctions.py function window_size_input() and reccomend a window size of 4
window_size = psf.window_size_input(4)

#store in a variable a bash command that uses plotcon to plot the conservation of the sequences using the sequence alignment and the chosen window size, displaying the plot to the screen
plotcon_display_command = "plotcon multiple_sequence_alignment/alignment.fasta -auto -graph x11 -winsize {}".format(window_size)
#same command, but the plot is saved
plotcon_save_command = "plotcon multiple_sequence_alignment/alignment.fasta -auto -goutfile conservation_plot/plotcon_plot -graph svg -winsize {}".format(window_size)
#loop indefinitely
while True:
	#ask the user if they would like to display a conservation plot to the screen and store the users input in a variable
	display_decision = input("Would you like to display the conservation plot to the screen? (y/n) (the plot will automatically be saved.)\n\n\t>")
	#if the user replies with any form of yes
	if display_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
		#call the plotcon_display_command
		subprocess.call(plotcon_display_command, shell=True)
		#immediately break from the loop
		break
	#else if the user replies with any form of no
	elif display_decision in ('N', 'n', 'no', 'NO', 'No'):
		#immediately break from the loop
		break
	#else the wrong input has been given, inform the user and continue with the loop, asking for the original input again
	else:
		print("Wrong input! Please try again.\n")
		continue

#automatically call the plotcon_save_command, saving the plot to the conservation_plot directory
subprocess.call(plotcon_save_command, shell=True)
print("Conservation plot saved to conservation_plot directory\n")

###########################################
###########################################

########## Hydropathy analysis ##########
##########################################

#create a new directory to store hyropathy analysis results
os.mkdir("hydropathy_plot")
print("Starting hydropathy plotting...\n")

#obtain the desired window size using the proseqfunctions.py function window_size_input() and reccomend a window size of 19
window_size = psf.window_size_input(19)

pepwindowall_display_command = "pepwindowall multiple_sequence_alignment/alignment.fasta -auto -graph x11 -gxtitle='Residue' -gytitle='Hydropathy' -window {}".format(window_size)
pepwindowall_save_command = "pepwindowall multiple_sequence_alignment/alignment.fasta -auto -graph svg -gxtitle='Residue' -gytitle='Hydropathy' -goutfile2 hydropathy_plot/pepwindowall_plot -window {}".format(window_size)

while True:
	#ask the user if they would like to display a hydropathy plot to the screen and store the users input in a variable
	display_decision = input("Would you like to display the hydropathy plot to the screen? (y/n) (the plot will automatically be saved.)\n\n\t>")
	#if the user replies with any form of yes
	if display_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
		#call the pepwindowall_display_command
		subprocess.call(pepwindowall_display_command, shell=True)
		#immediately break from the loop
		break
	#else if the user replies with any form of no
	elif display_decision in ('N', 'n', 'no', 'NO', 'No'):
		#immediately break from the loop
		break
	#else the wrong input has been given, inform the user and continue with the loop, asking for the original input again
	else:
		print("Wrong input! Please try again.\n")
		continue


#automatically call the pepwindowall_save_command, saving the plot to the hydropathy_plot directory
subprocess.call(pepwindowall_save_command, shell=True)
print("Hydropathy plot saved to hydropathy_plot directory.\n")

##########################################
##########################################


########## Transmembrane analysis (optional wildcard) ##########
################################################################

#loop indefinitely
while True:
	#ask the user whether they would like to perform transmembrane analysis on the sequence subset and store the input as a string variable
	transmembrane_analysis_decision = input("Would you like to perform transmembrane analysis on this sequence subset? (y/n)\n\n\t>")
	#if the string variable is any form of yes 
	if transmembrane_analysis_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
		#create a new directory to store results of the analysis
		os.mkdir("transmembrane_map")
		#create a bash command that uses tmap tool on the alignment file and displays the output map to the screeen
		tmap_display_command = "tmap multiple_sequence_alignment/alignment.fasta -auto -outfile transmembrane_map/report.tmap -graph x11 -gsubtitle 'Plot of the propensities to form the middle (solid line) and the end (dashed line) of transmembrane regions' -gytitle='Propensity'"
		#create a bash command that uses the tmap tool on the alignment file and saves the output to a new file
		tmap_save_command = "tmap multiple_sequence_alignment/alignment.fasta -auto -outfile transmembrane_map/report.tmap -goutfile transmembrane_map/tmap_map -graph svg -gsubtitle 'Plot of the propensities to form the middle (solid line) and the end (dashed line) of transmembrane regions' -gytitle='Propensity'"
		#call the save command automatically
		subprocess.call(tmap_save_command, shell=True)
		#loop indefinitely
		while True:
			#ask the user if they would like to display the map to the screen and store the input as a variable
			display_decision = input("Would you like to display the transmembrane map to the screen? (y/n) (the map will automatically be saved.)\n\n\t>")
			#if the string variable is any form of yes
			if display_decision in ('Y', 'y', 'yes', 'YES', 'Yes'):
				#call the display command
				subprocess.call(tmap_display_command, shell=True)
				#break from the loop
				break
			#else if the string variable is any form of no
			elif display_decision in ('N', 'n', 'no', 'NO', 'No'):
				#immediately break from the loop
				break
			#else if the string variable is anything else the input is invalid
			else:
				#inform the user
				print("Wrong input! Please try again.")
				#revert to the start of the loop
				continue
		#immediately break from the loop
		break
	#else if the user replies with any form of no
	elif transmembrane_analysis_decision in ('N', 'n', 'no', 'NO', 'No'):
		#immediately break from the loop
		break
	#else the wrong input has been given, inform the user and continue with the loop, asking for the original input again
	else:
		print("Wrong input! Please try again.\n")
		continue
print("Transmembrane map and report saved to transmembrane_map directory.\n")

################################################################
################################################################

print("All analyses complete! All results stored in the specified path.")
