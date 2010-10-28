"""Support module containing some common functionality

This should only be used by other modules in the chemfp package.
"""

import sys
from datetime import datetime
import errno

FPSv1_HEADER = (
    "#FPS1\n"
    "#num_bits={num_bits}\n"
    "{software_line}"
    "{type_line}"
    "{source_line}"
    "#date={date}\n")

def now_in_isoformat():
    "Current date and time in ISO format"
    # Exclude fractional seconds
    return datetime.now().isoformat().split(".", 1)[0]

def source_to_source_line(source):
    if source is None:
        return ""
    # This is likely a filename, which may be a unicode string.
    # Encode in utf-8 then remove any newlines.
    if isinstance(source, unicode):
        source = source.encode("utf8")
    source = "".join(source.splitlines())
    return "#source=" + source + "\n"

def format_fpsv1_header(num_bits, software=None, type=None, source=None, date=None):
    if date is None:
        date = now_in_isoformat()

    source_line = source_to_source_line(source)
    if software is None:
        software_line = ""
    else:
        software_line = "#software="+software+"\n"
    if type is None:
        type_line = ""
    else:
        type_line = "#type="+type+"\n"

    return FPSv1_HEADER.format(
        num_bits = num_bits,
        software_line = software_line,
        type_line = type_line,
        source_line = source_line,
        date = date,
        )

def convert_structures(reader, fingerprinter, outfile):
    for title, mol in reader:
        fp = fingerprinter(mol)
        outfile.write("%s %s\n" % (fp.encode("hex"), title))

def open_output(output_filename):
    if output_filename is None:
        return sys.stdout
    elif isinstance(output_filename, basestring):
        # Here's where I could check for a .gz extension and compress the output
        return open(output_filename, "w")
    return output_filename

def write_to_pipe(outfile, text):
    try:
        outfile.write(text)
    except IOError, err:
        # Catch any pipe errors, like when piping the output to "head -10"
        if err.errno != errno.EPIPE:
            raise
        try:
            outfile.close()
        except IOError:
            pass
        raise SystemExit()
        

def generate_fpsv1_output(header_kwargs,
                          reader,
                          fingerprinter,
                          output_filename):

    outfile = open_output(output_filename)
    try:
        outfile.write(format_fpsv1_header(**header_kwargs))
        convert_structures(reader, fingerprinter, outfile)
    except IOError, err:
        # Catch any pipe errors, like when piping the output to "head -10"
        if err.errno != errno.EPIPE:
            raise
        try:
            outfile.close()
        except IOError:
            pass
        
