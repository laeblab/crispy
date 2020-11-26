#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import bz2
import gzip
import sys

from os import fspath

from pathlib import Path


def open_ro(filename, mode="rt"):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    filename = fspath(filename)
    if mode not in ("rt", "rb", "r"):
        raise ValueError(mode)
    elif mode == "r":
        # Ensure uniform behavior between open/gzip.open/bz2.open
        mode = "rt"

    with open(filename, "rb") as handle:
        header = handle.read(2)

    if header == b"\x1f\x8b":
        return gzip.open(filename, mode)
    elif header == b"BZ":
        return bz2.open(filename, mode)
    else:
        return open(filename, mode)


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type=Path, help="GFF file to be filtered")
    parser.add_argument(
        "genes",
        type=Path,
        help="Text file containing gene names, one per line",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    id_fields = ("gene", "id")

    input_found = set()
    input_whitelist = set()
    with open_ro(args.genes) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                input_whitelist.add(line.lower())

    in_header = True
    warned_about = set([None])

    whitelist = set(input_whitelist)
    with open_ro(args.gff) as handle:
        for line in handle:
            if line.startswith("#"):
                if in_header:
                    print(line, end="")
            else:
                in_header = False
                fields = line.rstrip("\r\n").split("\t")

                props = {}
                for attr in fields[8].lower().split(";"):
                    key, value = attr.split("=", 1)
                    props[key] = value

                if any(props.get(key, "") in whitelist for key in id_fields):
                    whitelist.add(props["id"])
                    input_found.update(props[key] for key in id_fields if key in props)
                    print(line, end="")

                    parent_id = props.get("parent")
                    if parent_id not in whitelist and parent_id not in warned_about:
                        sys.stderr.write(
                            "WARNING: Parent {!r} of {!r} not whitelisted\n".format(
                                parent_id, props["id"]
                            )
                        )
                        warned_about.add(parent_id)

    for missing in (input_whitelist - input_found):
        sys.stderr.write(f"WARNING: {missing} not found in GFF!\n")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
