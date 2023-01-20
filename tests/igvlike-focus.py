#! /usr/bin/python3
import argparse
import itertools
import mmap
import sys

verbose = False

def main(arguments):
    global verbose
    if arguments.verbosity:
        verbose = len(arguments.verbosity)
    chromosome = "chr"+str(arguments.chr)
    start = arguments.index-arguments.context
    end   = arguments.index+arguments.context+1
    if arguments.before != None:
        start = arguments.index-arguments.before
    if arguments.after != None:
        end = arguments.index+arguments.after
    fasta_file_name = None
    mapping_file_names = []
    for f in arguments.files:
        extension = f.split(".")[-1]
        if extension == "fasta":
            if fasta_file_name == None:
                fasta_file_name = f
            else:
                print("Can't handle more than one fasta file")
                return 1
        elif extension == "map":
            mapping_file_names.append(f)
        else:
            print("unknown file extension :", extension)
            print("for :", f)
            print("please only use fasta and map files");

    insertions = [0]*(end-start)
    block = []

    if fasta_file_name != None:
        block.append(list(zip(range(start, end), get_genome_part(fasta_file_name, chromosome, start, end, arguments.index))))

    if verbose:
        print("parsing intersecting mappings")
    for mapping_start_address, mapping_string in get_intersecting_mappings(mapping_file_names, chromosome, start, end, arguments.max_read_size):
        if verbose>1:
            print("found mapping at :", mapping_start_address)
        this_line_insertions = [0]*(end-start+60)
        block.append([])
        address = max(start-1, mapping_start_address-1)
        i = address-mapping_start_address+1
        while address<end:
            if i<len(mapping_string):
                letter = mapping_string[i]
                if letter.isupper():
                    this_line_insertions[address-start] += 1
                else:
                    address += 1
                if start<=address<end:
                    block[-1].append((address, letter))
                i+=1
            else:
                break
        insertions = [max(ins1, ins2) for ins1, ins2 in zip(insertions, this_line_insertions)]
        if arguments.progressive:
            print_alignments(insertions, block, start, end)

    print_alignments(insertions, block, start, end)


def print_alignments(insertions, block, start, end):
    if verbose:
        print("generating output")
    for line in block:
        strings = ["" for _ in range(start, end)]
        for addr, c in line:
            strings[addr-start] += c
        strings = [s.ljust(1+insertions[i], " ") for i,s in enumerate(strings)]
        print("".join(strings))

def get_genome_part(file_name, chromosome, start, end, index=None):
    if verbose:
        print("getting reference genome extract")
        print("in chromosome :", chromosome)
        print("around index =", index)
    f = open(file_name, "r")
    address = 0
    current_chr = None
    for line in f:
        line_length = len(line.strip())
        if line[0] == ">":
            current_chr = line[1:].strip()
            continue
        endline_address = address+line_length
        if current_chr == chromosome and endline_address>start:
            for i in range(max(start-address, 0), min(end-address, line_length)):
                if address+i == index-1:
                    yield line[i].upper()
                else:
                    yield line[i].lower()
            if endline_address>end:
                return
        address = endline_address

def get_intersecting_mappings(file_names, chromosome, start, end, max_read_size):
    for file_name in file_names:
        with get_index(file_name) as index_file:
            with open(file_name, mode="r", encoding="utf-8") as map_file:
                with mmap.mmap(map_file.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_map_file:
                    last_address=0
                    last_chr = None
                    current_chr = None
                    address = 0
                    for line in index_file.readlines():
                        last_address = address
                        last_chr = current_chr
                        current_chr, index, address = line.split()
                        index = int(index)
                        address = int(address)
                        if current_chr == chromosome and index >= start or current_chr != chromosome and last_chr == chromosome:
                            if verbose>1:
                                print("starting reading ("+current_chr+":"+str(index)+") at "+str(address))
                            break
                    mmap_map_file.seek(int(last_address))
                    yielded_line = False
                    for line in mmap_readlines(mmap_map_file):
                        if line == b'':
                            continue
                        current_chr = line.split()[0].decode()
                        index = int(line.split()[1])
                        if not(yielded_line) and index+max_read_size<start-1:
                            continue
                        if current_chr==chromosome and index<end:
                            yield index, line.split()[2].decode()
                            yielded_line = True
                            continue
                        if current_chr==chromosome and index>end-1:
                            break
                        if current_chr!=chromosome and yielded_line:
                            break

def mmap_readlines(mmap_obj):
    while True:
        try:
            yield mmap_obj.readline()
        except:
            return

def get_index(file_name):
    try:
        return open(file_name+".index", mode="r")
    except IOError:
        create_index(file_name)
    return open(file_name+".index", mode="r")

def create_index(file_name):
    print("creating index for", file_name)
    with open(file_name+".index", mode="w") as index_file:
        with open(file_name, mode="r", encoding="utf-8") as map_file:
            with mmap.mmap(map_file.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_map_file:
                if verbose:
                    print("file openned for indexing")
                for i in itertools.count():
                    if verbose>1:
                        print("indexing :", i*10000000)
                    try:
                        mmap_map_file.seek(i*10000000)
                    except ValueError:
                        return
                    if i>0:
                        mmap_map_file.readline()
                    line_address = mmap_map_file.tell()
                    line = mmap_map_file.readline().split()
                    chromosome = line[0]
                    index = line[1]
                    index_file.write(chromosome.decode()+"\t"+index.decode()+"\t"+str(line_address)+"\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="show read mappings around a specific part of the genome.")
    parser.add_argument("files", nargs="+", type=str, help="Fasta and map files to read from")
    parser.add_argument("-c", "--chr", type=int, dest="chr", help="The chromosome in which to show the read mappings")
    parser.add_argument("-i", "--index", type=int, dest="index", help="The address around which to show read mappings")
    parser.add_argument("-A", "--after", default=None, type=int, dest='after', help="how many bases to show after address")
    parser.add_argument("-B", "--before", default=None, type=int, dest='before', help="how many bases to show before address")
    parser.add_argument("-C", "--context", default=40, type=int, dest='context', help="how many bases to show before and after address")
    parser.add_argument("-s", "--max-read-size", default=150, type=int, dest="max_read_size", help="the maximum length of a read; preferably exact")
    parser.add_argument("-v", action="append_const", dest="verbosity", const=True, help="add verbosity; flag may be set multiple times for more verbosity")
    parser.add_argument("-p", "--progressive", action='store_true', dest="progressive", help="print results everytime a new read is found")
    main(parser.parse_args())
