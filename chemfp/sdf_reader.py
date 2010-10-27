"""sdf_reader - iterate through records in an SD file"""

# This is used by the command-line SDF reader, used when the source
# fingerprints are already encoded in fields of an SD file.

# It's also used by the RDKit parser - see the comments there.


from __future__ import absolute_import
from . import decompressors

__all__ = ["open_sdfile", "iter_sdf_records", "iter_two_tags", "iter_title_and_single_tag"]

import re
import chemfp


# Do a quick check that the SD record is in the correct format
_sdf_check_pat = re.compile(r"""
.*\n     # line 1
.*\n     # line 2
.*\n     # line 3
         # next are the number of atoms and bonds. This pattern allows
         #  '  0', ' 00', '000', '00 ', '0  ' and ' 0 ', which is the
         # same as checking for field.strip().isdigit()
         # The escape '\040' is for the space character
(\040\040\d|\040\d\d|\d\d\d|\d\d\040|\d\040|\040\d\040)  # number of bonds
(\040\040\d|\040\d\d|\d\d\d|\d\d\040|\d\040|\040\d\040)  # number of bonds
        # Only space and digits are allowed before the required V2000 or V3000
[\0400-9]{28}V(2000|3000)
""", re.X)

class FileLocation(object):
    """A mutable instance used to track record title and position information

    You may pass one of these to the readers to that information. I'm
    a bit unsure about doing it this way but I liked the options of
    passing back a 2-ple of (record, position) or records with an
    attribute like record.position even less.

    WARNING: the attributes '.lineno', and '._record' be set-able.
    Passing in a non-FileLocation object may cause interoperability
    problems in the future.
    """
    def __init__(self, filename=None):
        self.filename = filename
        self.lineno = 1
        self._record = None  # internal variable; it only valid enough to get the title
    @property
    def title(self):
        # The title isn't needed for most cases so don't extract it unless needed
        if self._record is None:
            return None
        return self._record[:self._record.find("\n")].strip()
    
    def where(self):
        s = "line {self.lineno}".format(self=self)
        if self.filename is not None:
            s += " of {self.filename!r}".format(self=self)
        title = self.title
        if title is not None:
            s += " ({title})".format(title=title)
        return s

    def info(self):
        return dict(filename=self.filename, lineno=self.lineno, title=self.title)

def default_reader_error(msg, loc):
    "This is the default error handler. It raises an informative TypeError"
    raise TypeError(msg + " at " + loc.where())

def open_sdfile(source=None, decompressor=decompressors.AutoDetectDecompression,
                loc=None, reader_error=default_reader_error):
    """Open an SD file and return an iterator over the SD records, as blocks of text

    If 'source' is None, use sys.stdin. If 'source' is a string, then
    open it for reading in universal mode. Otherwise, 'source' must be
    an input stream.

    'decompressor' can be 'none', 'gzip', or 'bzip2' to indicate a specific
    compression type. If 'auto' and the source is a filename then use the
    extension to determine the compression type.

    'loc' is for the experimental FileLocation tracker

    'reader_error' is a callback used to report failures, some of
    which may be ignored. This feature is also experimental.
    """
    infile = decompressors.open_and_decompress_universal(source, decompressor)
    return iter_sdf_records(infile, loc, reader_error)

# My original implementation used a slow line-oriented parser.  That
# was decently fast, but this version, which reads a block at a time
# and works directly on those blocks, is over 3 times as fast.

def iter_sdf_records(infile, loc=None, reader_error=default_reader_error):
    """Iterate over records in an SD file, returning records as strings

    'infile' is an input stream. If the 'name' attribute exists then it's
    used in the file location.

    'loc' is for the experimental FileLocation tracker

    'reader_error' is a callback used to report failures, some of
    which may be ignored. This feature is also experimental.
    """
    if loc is None:
        loc = FileLocation()
    if loc.filename is None:
        loc.filename = getattr(infile, "name", None)
    pushback_buffer = ''
    records = None
    while 1:
        if not records:
            read_data = infile.read(32768)
            if not read_data:
                # No more data from the file. If there is something in the
                # pushback buffer then it's an incomplete record
                if pushback_buffer:
                    if pushback_buffer.endswith("\n$$$$"):
                        # The file is missing the terminal newline. Compensate.
                        # This will cause an extra read after the known end-of-file.
                        # (Is that a problem? It's not supposed to be one.)
                        pushback_buffer += "\n"
                    else:
                        loc._record = None
                        if loc.lineno == 1:
                            # No records read. Wrong format.
                            reader_error("Could not find a valid SD record", loc)
                        else:
                            reader_error(
   "unexpected content at the end of the file (perhaps the last record is truncated?)",
                                loc)
                        break
                else:
                    # We're done!
                    break
            # Join the two blocks of text. This should be enough to
            # have lots of records, so split, and the last term is
            # either a partial record or the empty string. Keep track
            # of that for use in the next go-around.
            records = (pushback_buffer + read_data).split("\n$$$$\n")
            pushback_buffer = records.pop()  # either '' or a partial record

            # It is possible though unlikely that the merged blocks of
            # text contains only an incompelte record, so this might
            # loop again. However, the joining and searching is an
            # O(n**2) operation, so I don't want to do that too often.
            # While it's possible to fix this, there should be no
            # reason to support huge records - they don't exist unless
            # you are really stretching and doing things like storing
            # images or other large data in the SD tags.
            # To prevent timing problems, don't allow huge records.
            if len(pushback_buffer) > 200000:
                loc._record = None
                reader_error("record is too large for this reader", loc)
                return
        else:
            # We have a set of records, one string per record. Pass them back.
            for record in records:
                # A simple, quick check that it looks about right
                if not _sdf_check_pat.match(record):
                    loc._record = record
                    reader_error("incorrectly formatted record", loc)
                    # If the reader_error callback returns then just skip the record
                else:
                    record += "\n$$$$\n"  # restore the split text
                    loc._record = record
                    yield record
                loc.lineno += record.count("\n")
            records = None

### Is something like this needed? Perhaps when used as a co-process?
# def iter_sdf_records_line_buffered(infile, loc=None, reader_error=default_reader_error):
    

# I tried implementing this search with a regular expression but it
# was about 30% slower than this more direct (and complicated) search.

# Note: tag_substr must contain the "<" and ">"
def _find_tag_data(rec, tag_substr):
    "Return the first data line for the given tag substring, or return None"
    startpos = 0
    while 1:
        tag_start = rec.find(tag_substr, startpos)
        if tag_start == -1:
            return None

        # rfind cannot return -1 because _sdf_check_pat verified there
        # are at least 3 newlines. The +1 is safe because there's at
        # least the "<" and ">" from the tag.
        tag_line_start = rec.rfind("\n", None, tag_start) + 1
        if rec[tag_line_start] != ">":
            # This tag is not on a data tag line. It might be the value for
            # some of the text field.
            startpos = tag_start + 1
            continue

        # This is an actual tag line. Find the start of the next line.
        # The record must end with "\n$$$$\n" so find will never return -1
        # and never return the last character position.
        next_line_start = rec.find("\n", tag_start) + 1

        # These might occur if there is no data content
        if rec[next_line_start]==">" or rec[next_line_start:next_line_start+4]=="$$$$":
            return ""

        # Otherwise, get up to the end of the line
        return rec[next_line_start:rec.find("\n", next_line_start)]

# These are not legal tag characters (while others may be against the
# SD file spec, these will break the parser)
_bad_char = re.compile(r"[<>\n\r\t]")

def iter_two_tags(sdf_iter, tag1, tag2):
    """Iterate over SD records to get the first data line of the two named tags
    
    The values are returned as a (tag1_value, tag2_value) 2-ple,
    where the value comes from the first line after the given tag (if
    present). If there is no line then return the empty string. If the
    tag isn't present, return None. If the tag exists mutliple times
    then only the first value is returned.

    infile - an input stream (its 'name' attribute should be the source filename)
    tag1 - the name of the first tag
    tag2 - the name of the second tag
    """
    m = _bad_char.search(tag1)
    if m:
        raise TypeError("tag1 must not contain the character %r" % (m.group(0),))
    m = _bad_char.search(tag2)
    if m:
        raise TypeError("tag2 must not contain the character %r" % (m.group(0),))
        
    tag1_substr = "<" + tag1 + ">"
    tag2_substr = "<" + tag2 + ">"
    for rec in sdf_iter:
        yield _find_tag_data(rec, tag1_substr), _find_tag_data(rec, tag2_substr)

def iter_title_and_single_tag(sdf_iter, tag):
    """Iterate over SD records to get the title line and data line for the specified tag

    (NOTE: the title line is the first line of the SD file, and not an SD tag.)

    The values are returned as a (title, tag_value) 2-ple, where the
    value comes from the first line after the specified tag (if
    present). If there is no line then return the empty string. If the
    tag isn't present, return None. If the tag exists multiple times
    then only the first value is returned.

    infile - an input stream (its 'name' attribute should be the source filename)
    tag - the name of the tag value to return
    """
    m = _bad_char.search(tag)
    if m:
        raise TypeError("tag must not contain the character %r" % (m.group(0),))
    
    tag_substr = "<" + tag + ">"
    for rec in sdf_iter:
        yield rec[:rec.find("\n")].strip(), _find_tag_data(rec, tag_substr)
