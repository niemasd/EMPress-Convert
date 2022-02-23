#!/usr/bin/env python3
from bp import parse_newick
from gzip import open as gopen
from itertools import zip_longest
from os.path import isfile
from struct import pack
from sys import stdin
import argparse
DEFAULT_BUFSIZE = 1048576

# https://github.com/biocore/empress/blob/bf6f755c5187e543ff13c2c624b153fb97918b23/empress/tools.py#L342-L407
def shifting(bitlist, size=51):
    """Takes a list of 0-1s, splits in size and converts it to a list of int
    Parameters
    ----------
    bitlist: list of int
        The input list of 0-1
    size: int
        The size of the buffer
    Returns
    -------
    list of int
        Representation of the 0-1s as a list of int
    Raises
    ------
    ValueError
        If any of the list values is different than 0 or 1
    References
    ----------
    Borrowed from https://stackoverflow.com/a/12461400
    Example
    -------
    shifting([1, 0, 0, 0, 0, 1], size=3) => [4, 1]
    """
    if not all(x in [0, 1] for x in bitlist):
        raise ValueError('Your list has values other than 0-1s')

    values = [iter(bitlist)] * size
    ints = []
    for num in zip_longest(*values):
        out = 0
        init_zeros = []
        seen_one = False
        for bit in num:
            if bit is None:
                continue

            if not seen_one:
                if bit == 0:
                    init_zeros.append(0)
                else:
                    seen_one = True

            out = (out << 1) | bit

        # if out == 0, everything was zeros so we can simply add init_zeros
        if out == 0:
            ints.extend(init_zeros)
        else:
            ints.append(out)

    # we need to check init_zeros for the last loop in case the last value
    # had padded zeros
    if init_zeros and out != 0:
        # rm last value
        ints = ints[:-1]
        # add zeros
        ints.extend(init_zeros)
        # readd last value
        ints.append(out)

    return ints

if __name__ == "__main__":
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-on', '--output_lengths', required=True, type=str, help="Output Lengths File (.gz)")
    parser.add_argument('-ol', '--output_labels', required=True, type=str, help="Output Labels File (.txt.gz)")
    parser.add_argument('-ot', '--output_tree', required=True, type=str, help="Output Tree File (.gz)")
    args = parser.parse_args()
    if args.input != 'stdin' and not isfile(args.input):
        raise ValueError("File not found: %s" % args.input)
    if isfile(args.output_lengths):
        raise ValueError("File exists: %s" % args.output_lengths)
    if not args.output_lengths.lower().endswith('.gz'):
        raise ValueError("Should be .gz: %s" % args.output_lengths)
    if isfile(args.output_labels):
        raise ValueError("File exists: %s" % args.output_labels)
    if not args.output_labels.lower().endswith('.txt.gz'):
        raise ValueError("Should be .txt.gz: %s" % args.output_labels)
    if isfile(args.output_tree):
        raise ValueError("File exists: %s" % args.output_tree)
    if not args.output_tree.lower().endswith('.gz'):
        raise ValueError("Should be .gz: %s" % args.output_tree)

    # load input tree
    if args.input == 'stdin':
        nwk = stdin.read().strip()
    elif args.input.lower().endswith('.gz'):
        f = gopen(args.input, 'r'); nwk = f.read().decode().strip(); f.close()
    else:
        f = open(args.input, 'r', buffering=DEFAULT_BUFSIZE); nwk = f.read().strip(); f.close()

    tree = parse_newick(nwk)

    # create output files
    out_lengths = gopen(args.output_lengths, 'w', 9)
    out_labels = gopen(args.output_labels, 'w', 9)
    out_tree = gopen(args.output_tree, 'w', 9)

    for i in range(1, len(tree) + 1):
        node = tree.postorderselect(i)
        name = tree.name(node)
        length = tree.length(node)

        # write current length
        if length is None:
            length = 0
        tmp_num = pack('>d', length)
        out_lengths.write(tmp_num)

        # write current label
        if name is None:
            tmp_label = ""
        else:
            tmp_label = name

        out_labels.write(tmp_label.encode())

    # convert tree structure
    for num in shifting(tree.B):
        tmp_num = pack('>Q', num)
        out_tree.write(tmp_num)

    # close output files
    out_lengths.close()
    out_labels.close()
    out_tree.close()
