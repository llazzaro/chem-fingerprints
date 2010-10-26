"""chemfp.decompressors - support APIs for decompressing different compression schemes

All decompressors implement the following API:

  open_filename_binary(filename) -- open filename for reading, in binary mode.
  open_filename_universal(filename) -- open filename for reading, in universal newline mode.
  open_fileobj(fileobj) -- take input from fileobj, converting as appropriate
  strip_extension(filename) -- remove the compression extension, if present


"""

import os
import gzip
import subprocess

# Decompressors must implement the following API:
#   open_filename_binary(filename) -- open the given file 

class Uncompressed(object):
    @staticmethod
    def open_filename_binary(filename):
        return open(filename, "rb")
    
    @staticmethod
    def open_filename_universal(filename):
        return open(filename, "rU")
    
    @staticmethod
    def open_fileobj(fileobj):
        return fileobj
    
    @staticmethod
    def strip_extension(filename):
        return filename
    
class GzipDecompressor(object):
    @staticmethod
    def open_filename_binary(filename):
        return gzip.GzipFile(filename)

    # gzip does not implement universal mode. Cross your fingers and hope this works.
    open_filename_universal = open_filename_binary
    
    @staticmethod
    def open_fileobj(fileobj):
        return gzip.GzipFile(fileobj=fileobj)
    
    @staticmethod
    def strip_extension(filename):
        if filename[-3:].lower() == ".gz":
            return filename[:-3]
        return filename

# TODO: Detect if the 3rd party bzip module is installed and use that
# instead of the subprocess call.

class Bzip2Decompressor(object):
    @staticmethod
    def open_filename_binary(filename):
        p = subprocess.Popen(["bzcat", filename],
                             stdout = subprocess.PIPE)
        return p.stdout

    @staticmethod
    def open_filename_universal(filename):
        p = subprocess.Popen(["bzcat", filename],
                             stdout = subprocess.PIPE,
                             universal_newlines = True)
        return p.stdout

    @staticmethod
    def open_fileobj(filename):
        raise NotImplementedError("bzip2 decompression does not yet support input streams")

    @staticmethod
    def strip_extension(filename):
        if filename[-4:].lower() == ".bz2":
            return filename[:-4]
        return filename
    

class AutoDetectDecompression(object):
    """Choose the right decompressor based on the filename extension"""
    @staticmethod
    def open_filename_binary(filename):
        decompressor = detect_decompressor(filename)
        return decompressor.open_filename_binary(filename)

    @staticmethod
    def open_filename_universal(filename):
        decompressor = detect_decompressor(filename)
        return decompressor.open_filename_universal(filename)

    @staticmethod
    def open_fileobj(fileobj):
        # No support for auto-sniffing this; assume uncompressed.  (To
        # support sniffing I think I would have to read the first few
        # bytes then implement some sort of file-like object which
        # merges the already read text with the rest of the partially
        # read file object. That should be possible with this API.)
        return Uncompressed.open_fileobj(fileobj)

    @staticmethod
    def strip_extension(filename):
        return detect_decompressor(filename).strip_extension(filename)

def detect_decompressor(filename):
    """Use the filename extension to pick the right decompressor.

    If the extension is unknown, assume it's uncompressed.
    """
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext == ".gz":
        return GzipDecompressor
    if ext == ".bz2":
        return Bzip2Decompressor
    return Uncompressed

decompressor_names = {
    "gzip": GzipDecompressor,
    "bzip2": Bzip2Decompressor,
    "uncompressed": Uncompressed,
    "autodetect": AutoDetectDecompression,
    }

def get_named_decompressor(name_or_obj):
    """if a string, return the named decompressor, otherwise return the input itself

    This is used to simplify decompressor selection
    """
    if isinstance(name_or_obj, basestring):
        return decompressor_names[name_or_obj]
    return name_or_obj
