from __future__ import absolute_import
import sys

from .. import argparse, io, readers, check_metadata_problems

parser = argparse.ArgumentParser(
    description="Merge multiple FPS files into a single FPS file",
    )
parser.add_argument(
    "-o", "--output", metavar="FILENAME",
    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument("filenames", nargs="+", help="list of FPS files", default=[])

def main(args=None):
    args = parser.parse_args(args)

    if not args.filenames:
        return

    metadata = None
    sources = []
    for filename in args.filenames:
        fps_reader = readers.open_fps(filename)
        if metadata is None:
            metadata = fps_reader.metadata
        else:
            problems = check_metadata_problems(metadata, fps_reader.metadata)
            if problems:
                for (severity, error, msg_template) in problems:
                    msg = msg_template.format(metadata1 = repr(args.filenames[0]),
                                              metadata2 = repr(filename))
                    if severity == "warning":
                        sys.stderr.write("WARNING: " + msg + "\n")
                    elif severity == "error":
                        sys.stderr.write("ERRORR: " + msg + "\n")
                        raise SystemExit(1)
                    elif severity == "info":
                        sys.stderr.write("INFO: " + msg + "\n")
                    else:
                        raise AssertionError(severity)
        sources.extend(fps_reader.metadata.sources)

    metadata.sources = sources
    metadata.date = io.utcnow()

    with io._closing_output(args.output) as outfile:
        with io.ignore_pipe_errors:
            io.write_fps1_magic(outfile)
            io.write_fps1_header(outfile, metadata)
            
            for filename in args.filenames:
                fps_reader = readers.open_fps(filename)
                for block in fps_reader.iter_blocks():
                    outfile.write(block)

if __name__ == "__main__":
    main()
