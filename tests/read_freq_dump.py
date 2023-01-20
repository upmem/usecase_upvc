#! /usr/bin/python3
import argparse
import itertools
import mmap
import sys
import struct

class freq_dump:
    def __init__(self, file_name):
        self.file_object = open(file_name, mode="rb")
        self.mmap_freq_dump = mmap.mmap(self.file_object.fileno(), length=0, access=mmap.ACCESS_READ)
        self.n_seq = self.mmap_freq_dump.read_byte()
        self.seq_starts = [int.from_bytes(self.mmap_freq_dump.read(8), byteorder='little') for _ in range(self.n_seq)]

    def read_address(self, sequence_id, index):
        self.mmap_freq_dump.seek(self.seq_starts[sequence_id] + index*30)
        freqs = [struct.unpack('f', self.mmap_freq_dump.read(4))[0] for _ in range(5)]
        scores = [int.from_bytes(self.mmap_freq_dump.read(2), byteorder='little') for _ in range(5)]
        return freqs, scores

    def get_string_from_address(self, seq_n, index):
        freqs, scores = self.read_address(seq_n-1, index-1)
        return " ".join(str(v) for v in [f'{seq_n: >2}', f'{index: >8}', "", " ".join(f'{v:7.2f}' for v in freqs), "\t", " ".join(f'{v: >4}' for v in scores)])

    def __del__(self):
        self.close()

    def close(self):
        self.mmap_freq_dump.close()
        self.file_object.close()


def main(args):
    context = args.context
    dump_file_name = args.dump_file
    freq_dump_obj = freq_dump(dump_file_name)
    for address_file_name in args.address_files:
        with open(address_file_name, "r") as address_file:
            focuses = []
            last_dump = (0,0)
            for _ in range(context):
                line = address_file.readline()
                if line=="":
                    break
                split_line = line.split(" ")
                seq_n = int(split_line[0])
                index = int(split_line[1])
                focuses.append((seq_n, index))
            while focuses!=[]:
                line = address_file.readline()
                if line!="":
                    split_line = line.split(" ")
                    seq_n = int(split_line[0])
                    index = int(split_line[1])
                    focuses.append((seq_n, index))
                seq_n, index = focuses[0]
                start = index-context
                end = index+context
                if last_dump[0] == seq_n and last_dump[1]+1 >= start:
                    start = max(start, last_dump[1]+1)
                else:
                    print()
                for i in range(start, end+1):
                    if (seq_n, i) in focuses:
                        print("* "+freq_dump_obj.get_string_from_address(seq_n, i))
                    else:
                        print("  "+freq_dump_obj.get_string_from_address(seq_n, i))
                last_dump = (seq_n, end)
                focuses.pop(0)
    freq_dump_obj.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="show frequency table dump in human readable format around specific indices in the genome")
    parser.add_argument("dump_file", type=str, help="Frequency table bin file to read from")
    parser.add_argument("address_files", nargs="+", type=str, help="files containing addresses around which to show frequency table")
    parser.add_argument("-C", "--context", default=2, type=int, dest="context", help="how many addresses to show before and after the focused addresses")
    try:
        main(parser.parse_args())
    except BrokenPipeError:
        pass
