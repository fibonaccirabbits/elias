#!/usr/local/bin/python3
from Bio import SeqIO
import logging
import sys

def fasta_parse_and_slice(file, slice_size, slize_overlap, monitor_level):
    """A functions that takes a full file path to a fasta file as input and returns a dictionary of that fasta file, a dictionary of all the slices from that file with a certain length and a shortened dictionary where all the duplicates have been removed."""

    logging.info('Fasta parser started')
    fasta_dict = {}  # A dicitionary of all sequences imported.
    slice_dict = {}  # A dictionary of all slices from the sequences
    monitor = 0
    for record in SeqIO.parse(file, "fasta"):

        slice_start = 0
        slice_index = 0
        fasta_dict[record.id] = str(record.seq)
        while True:
            slice = record.seq[slice_start:slice_start + slice_size]
            slice_dict[record.id + '_s' + str(slice_index)] = str(slice)  # creates slices and names the dictionary keys by the fasta header plus "_sX" where X is the slice number
            slice_start += slize_overlap  # Jumps forward in the sequence by the slize overlap set in the function variables.
            if len(slice) < slice_size:  # Breaks the loop when the last slice is shorter than the desiered slice length, as an indicator of the complete sequenced having been parsed.
                break
            slice_index += 1
        monitor += 1

        if monitor % monitor_level == 0:
            logging.info(str(monitor) + ' sequences have been sliced')

    reverse_slice = dict((v, k) for k, v in slice_dict.items())  # reverses the keys and values in the dictionary to remove duplicate sequences
    slice_dict_shortened = dict((v, k) for k, v in reverse_slice.items())  # Reverts the dictionary back to the original order without the removed duplicates

    print(len(list(slice_dict)))
    print(len(list(slice_dict_shortened)))

    return(fasta_dict, slice_dict, slice_dict_shortened)


def write_dict_to_fasta(dictionary, output_file):
    """A simple function that takes a dictionary structure and writes it as a fasta file with the dictionary key as a header and the corresponding value as a sequence"""

    logging.info('Fasta writer started')
    with open(output_file, 'w') as f:
        for key, value in dictionary.items():
            f.write('>' + key + '\n' + value + '\n\n')
        f.close()
        logging.info('Writing to file complete')
    return()


def main(file, slice_size, slize_overlap, monitor_level):
    logging.basicConfig(filename=str(file[0:-6]) + '_log.log', filemode='w', format='%(asctime)s: %(message)s', level=logging.INFO)
    logging.info('Script started')
    f, s, sds = fasta_parse_and_slice(file, slice_size, slize_overlap, monitor_level)
    outname = file.split('/')[-1].split('.')[0]
    outname = 'output/' + outname + '_output.fasta'
    write_dict_to_fasta(sds, outname)


main(file=sys.argv[1], slice_size=int(sys.argv[2]), slize_overlap=int(sys.argv[3]),monitor_level=int(sys.argv[4]))
