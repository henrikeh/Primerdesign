# coding: utf-8

### Script to generate a graph of the free energy of a primer sequence
#	Input: script, primer name via command line
#	Output: Plot showing free energy of pentamers
#	Based on: Rychlik, Wojciech. "Selection of primers for polymerase chain reaction." In PCR Protocols, 31–40. Springer, 1993.http://link.springer.com/protocol/10.1385/0-89603-244-2:31.

from sys import argv
import matplotlib.pyplot as plt
import numpy as np

# Command line arguments: script name, primer name (sequence has to be stored as .txt file)
script, primer_name = argv

primer_sequence = open(primer_name + ".txt")
primer = primer_sequence.read() 
primer_sequence.close()
primer_length = len(primer) 

# Values given in paper by Rychlik stored as dictionary

delta_G = {'AA': -1.9, 'AC': -1.3, 'AG': -1.6, 'AT': -1.5, 'CA': -1.9, 'CC': -3.1, 'CG': -3.6, 'CT': -1.6, 'GA': -1.6, 'GC': -3.1, 'GG': -3.1, 'GT': -1.3, 'TA': -1.0, 'TC': -1.6, 'TG': -1.9, 'TT': -1.9}

# Slice the primer sequence into pairs
primer_slices = []

for x in range(0, (primer_length - 1)): # iteration from the second base of the primer to the end (= length-1)
	primer_slices.append(primer[x:x+2]) # each iteration adds the xth and x+1st bases as a string to the list
	
value_list = []

for element in primer_slices: # go through the elements in the sliced primer sequence 1 by 1
	value_list.append(delta_G[element])  # appends the current value to the values list based on the pairs given in the dictionary delta_G

pentamer_value = []

for x in range(0,(primer_length - 4)):
	pentamer_value.append(abs(sum(value_list[x:(x+4)])))

x_values = np.arange(1,primer_length + 1) # creates dummy xvalues for plotting

# fill y values so lengths of lists match
missing = primer_length - len(pentamer_value)
addition = [None] * missing # fill difference with None to suppress plotting

pentamer_value.extend(addition)
x_ticks = list(primer) # Slice primer sequence into single characters for plotting

# Generate & save plot

plt.plot(x_values, pentamer_value, '-o')
plt.xticks(x_values, x_ticks) 
plt.xlabel('Primer sequence')
plt.ylabel('Free energy [kcal/mol]')
plt.savefig('./Output/' + primer_name + '_free_energy')
