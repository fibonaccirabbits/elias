# import stuff
from Bio import SeqIO
from find_files import find_files as fifi
import sys

def fasta_parse_and_slice(file, slice_size, slize_overlap):
    """A functions that takes a full file path to a fasta file as input and returns a dictionary of that fasta file, a dictionary of all the slices from that file with a certain length and a shortened dictionary where all the duplicates have been removed."""

    fasta_dict = {}  # A dicitionary of all sequences imported.
    slice_dict = {}  # A dictionary of all slices from the sequences
    for record in SeqIO.parse(file, "fasta"):
        slice_start = 0
        slice_index = 0
        fasta_dict[record.id] = str(record.seq)
        while True:
            slice = record.seq[slice_start:slice_start + slice_size]
            slice_dict[record.id + '_s' + str(slice_index)] = str(slice)  # creates slices and names the dictionary keys by the fasta header plus "_sX" where X is the slice number
            slice_start += slize_overlap  # Jumps forward in the sequence by the slize overlap set in the function variables.
            if len(slice) < 50:  # Breaks the loop when the last slice is shorter than the desiered slice length, as an indicator of the complete sequenced having been parsed.
                break
            slice_index += 1

        reverse_slice = dict((v, k) for k, v in slice_dict.items())  # reverses the keys and values in the dictionary to remove duplicate sequences
        slice_dict_shortened = dict((v, k) for k, v in reverse_slice.items())  # Reverts the dictionary back to the original order with the removed duplicates

    return(fasta_dict, slice_dict, slice_dict_shortened)


def write_dict_to_fasta(dictionary, output_file):
    """A simple function that takes a dictionary structure and writes it as a fasta file with the dictionary key as a header and the corresponding value as a sequence"""

    with open(output_file, 'w') as f:
        for key, value in dictionary.items():
            f.write('>' + key + '\n' + value + '\n\n')
        f.close()
    return()


#f, s, sds = fasta_parse_and_slice("/Users/eliastjarnhage/Desktop/ProteinFastaResults.fasta", 50, 28)
#write_dict_to_fasta(sds, "/Users/eliastjarnhage/Desktop/Output_FASTA.fasta")

#print(f)
#print(len(f))
#print(s)
#print(len(s))
#print(sds)
#print(len(sds))


infiles = fifi('input', '.fasta')
for infile in infiles:
	f, s, sds = fasta_parse_and_slice(infile, 50, 28)
	infileparts = infile.split('/')
	infilename = infileparts[-1].split('.')
	outname = 'output/' + infilename[0] + '_output.fasta'
	print('Writing outfile: %s' % outname)
	write_dict_to_fasta(sds, outname)

# 15 slices per NA and 18 slices per HA
