# Data loading functionalities
import polars as pl
from pathlib import Path
import logging as log
from .seq import Assembly
from .feature import Pileup
import sys
import os
from nanomotif.seq import DNAsequence

def load_fasta(path, trim_names=False, trim_character=" ") -> dict:
    """
    Reads a fasta file and returns a dictionary with the contig names as 
    keys and the sequences as values
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    data = {}
    active_sequence_name = "no_header"
    for line in lines:
        line = line.rstrip()
        if line[0] == '>':
            active_sequence_name = line[1:]
            if trim_names:
                active_sequence_name = active_sequence_name.split(trim_character)[0]
            if active_sequence_name not in data:
                data[active_sequence_name] = ''
        else:
            data[active_sequence_name] += line
    # Convert sequences to DNAsequence objects
    for key in data:
        data[key] = DNAsequence(data[key])
    return data



# Code copied from https://github.com/lh3/readfq/blob/master/readfq.py
class FastaReader:
    def readfq(self, fp): # this is a generator function
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break

def get_files_from_directory(directory, extension):
    """Get all files with specified extension from directory."""
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"Error: Directory '{directory}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not dir_path.is_dir():
        print(f"Error: '{directory}' is not a directory", file=sys.stderr)
        sys.exit(1)
    
    # Ensure extension starts with a dot
    if not extension.startswith('.'):
        extension = '.' + extension
    
    files = list(dir_path.glob(f"*{extension}"))
    if not files:
        print(f"Warning: No files with extension '{extension}' found in '{directory}'", file=sys.stderr)
    
    return [str(f) for f in files]

def process_fasta_file(filepath, reader):
    """Process a single FASTA file and return list of (contig, basename) tuples."""
    results = []
    basename = Path(filepath).stem  # Get filename without extension
    
    try:
        with open(filepath, 'r') as fp:
            for name, seq, qual in reader.readfq(fp):
                if name:  # Only process if we have a valid contig name
                    results.append((name, basename))
    except IOError as e:
        print(f"Error reading file '{filepath}': {e}", file=sys.stderr)
        return []
    
    return results


def generate_contig_bin(args):
    """
    Generate a DataFrame with contig names and their corresponding bins.
    If contig_bin is specified, it reads from the provided file.
    If files are specified, it processes those files to extract contig names.
    If a directory is specified, it reads all files with the given extension.
    """
    if args.contig_bin:
        log.debug(f"Reading contig-bin mapping from {args.contig_bin}")
        contig_bin_df = pl.read_csv(args.contig_bin, separator="\t", has_header=False, infer_schema_length=10000) \
            .rename({"column_1":"contig", "column_2":"bin"}) \
            .cast({'contig': pl.String, 'bin': pl.String})
        return contig_bin_df
    if args.files:
        log.debug(f"Processing specified files: {args.files}")
        files = args.files
        # Check if files exist
        for filepath in files:
            if not Path(filepath).exists():
                log.warning(f"Error: File '{filepath}' does not exist")
                sys.exit(1)
    else:
        log.debug(f"Processing files in directory {args.directory} with extension {args.extension}")
        files = get_files_from_directory(args.directory, args.extension)
    
    if not files:
        log.error("No files to process")
        sys.exit(1)
    
    # Process files
    reader = FastaReader()
    all_results = []
    
    for filepath in files:
        results = process_fasta_file(filepath, reader)
        all_results.extend(results)
    
    # Write output
    output_path = os.path.join(args.out, "temp", "contig_bin.tsv") if args.out else sys.stdout
    if not os.path.exists(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except OSError as e:
            log.error(f"Error creating output directory: {e}")
            sys.exit(1)
    if output_path:
        try:
            output_file = open(output_path, 'w')
        except IOError as e:
            log.error(f"Error opening output file '{args.out}': {e}")
            sys.exit(1)
    
    try:
        # Write data
        for contig, basename in all_results:
            output_file.write(f"{contig}\t{basename}\n")
        
    finally:
        if output_path:
            output_file.close()
    # Create DataFrame from results
    log.debug(f"Processed {len(files)} files, found {len(all_results)} contigs")
    return pl.DataFrame({
        "contig": [result[0] for result in all_results],
        "bin": [result[1] for result in all_results]
    })

