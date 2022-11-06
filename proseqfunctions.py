#!/usr/local/bin/python3

import os, subprocess, shutil, sys, collections, argparse

########## Function for clearing existing results directories ##########
########################################################################
#define a function that clears existing results directories if they exist
def clear_directory(directory_name):
	#try removing the directory
	try:
		shutil.rmtree(directory_name)
	#if a FileNotFoundError is returned, pass and suppress the error
	except FileNotFoundError:
		pass
########################################################################
########################################################################


########## Functions for fetching sequences and accessions/species names ##########
###################################################################################

#define a function that fetches the fasta sequences and accessions/species names of all found sequences
def standard_fetch(protein_family, taxonomic_group):
	#store in a variable a bash command that searches for the same query and fetches all the sequences in fasta format and stores them in an output file
	fetch_sequences_command = "esearch -db protein -query '{}[protein name] AND {}[organism] NOT partial' | efetch -format fasta > query_files/query.fasta ".format(protein_family, taxonomic_group)
	#use subprocess.call to run the command in linux environment
	subprocess.call(fetch_sequences_command, shell=True)
	#store in a variable a bash command that searches for the same query and fetches the accession values and species of all sequences and stores them in an output file
	species_names_command = "esearch -db protein -query '{}[protein name] AND {}[organism] NOT partial' | efetch -format gbc | xtract -pattern INSDSeq -element INSDSeq_accession-version -element INSDSeq_organism | cut -f1,2 > query_files/accessions_species.txt".format(protein_family, taxonomic_group)
	#use subprocess.call to run the command in linux environment
	subprocess.call(species_names_command, shell=True)

#define a function that fetches the fasta sequences of only 1000 sequences from the found sequences 
def limited_fetch(protein_family, taxonomic_group):
	#store in a variable a bash command that searches for the same query and fetches the first 1000 the sequences in fasta format and stores them in an output file
	fetch_sequences_command = "esearch -db protein -query '{}[protein name] AND {}[organism] NOT partial' | efetch -format fasta -stop 1000 > query_files/query.fasta ".format(protein_family, taxonomic_group)
	#use subprocess.call to run the command in linux environment
	subprocess.call(fetch_sequences_command, shell=True)
	#store in a variable a bash command that searches for the same query and fetches the accession values and species of the first 1000 sequences and stores them in an output file
	species_names_command = "esearch -db protein -query '{}[protein name] AND {}[organism] NOT partial' | efetch -format gbc -stop 1000 | xtract -pattern INSDSeq -element INSDSeq_accession-version -element INSDSeq_organism | cut -f1,2 > query_files/accessions_species.txt".format(protein_family, taxonomic_group)
	subprocess.call(species_names_command, shell=True)

###################################################################################
###################################################################################

########## Functions for checking inputs when selecting sequences for MSA/next steps ##########
###############################################################################################

def digit_check(user_input):
	#if the user input contains a space
	if ' ' in user_input:
		#if all inputs are space
		if user_input.isspace():
			#input is invalid
			valid = False
		#if not all inputs are space
		else:
			#for each character in the input
			for char in user_input:
				#if the character is a digit or a space, input remains valid
				if char.isdigit() or char.isspace():				
					valid = True
				#if the character is anything else, the input is no longer valid
				else:
					valid = False
					#break from loop immediately 
					break
	#if user input is all digits i.e only one number given, input is valid
	elif user_input.isdigit():
		valid = True
	#if user input is anything else, input is invalid
	else:
		valid = False
	#return validity
	return valid

#define function with user_input and the limit as arguments
def number_check(user_input, limit):
	#split user input into a list
	user_input=user_input.split()
	#for each value in the list
	for n in user_input:
		#make it an integer
		n = int(n)
		#if it is not between 1 and the limit
		if not 1 <= n <= limit:
			#input is invalid
			valid = False
			#break from loop immediately
			break
		#otherwise input remains valid
		valid = True
	#return validity
	return valid

###############################################################################################
###############################################################################################

########## Window size input function ##########
################################################
#define function with a reccommended window size as an argument 
def window_size_input(recommended_window):
	#loop indefinitely
	while True:
		#ask the user what window size they would like to use and store the input in a variable
		window_size = input("What window size would you like to use for the plot? (recommended: {})\n\n\t>".format(recommended_window))
		#if the input is only digits and the integer valye is greter than 0 break from the loop
		if window_size.isdigit() and int(window_size) > 0:
			break
		#if the input includes anything else but digits (including spaces), the wrong input has been given, inform the user and continue with the loop, asking for the original input again
		else:
			print("Wrong input! Please try again.\n")
			continue
	#return the window size as a variable
	return window_size

################################################
################################################
