import os

class UnsupportedPatternError(KeyError):
    def __init__(self, pattern, reason=None):
        KeyError.__init__(self, pattern)
        self.pattern = pattern
        if reason is None:
            reason = "Cannot interpret pattern definition"
        self.reason = reason
        self.filename = None
        self.lineno = None
    def __str__(self):
        msg = self.reason  + " " + repr(self.pattern)
        if self.lineno is not None:
            msg += " at line %d" % (self.lineno,)
        if self.filename is not None:
            msg += " in file %r" % (self.filename,)
        return msg

class PatternFile(object):
    def __init__(self, filename, max_bit, bit_definitions):
        assert max_bit >= 0, max_bit
        self.filename = filename
        self.max_bit = max_bit
        self.bit_definitions = bit_definitions
        self._bit_to_bit_definition = dict((bitdef.bit, bitdef) for bitdef in bit_definitions)

    def __getitem__(self, bit):
        return self._bit_to_bit_definition[bit]

    def __iter__(self):
        return iter(self._bit_to_bit_definition)

    
class BitDefinition(object):
    __slots__ = ("bit", "count", "pattern", "description", "lineno")
    def __init__(self, bit, count, pattern, description, lineno):
        self.bit = bit
        self.count = count
        self.pattern = pattern
        self.description = description
        self.lineno = lineno
    

def load_patterns(infile):
    if isinstance(infile, basestring):
        infile = open(infile, "rU")
    filename = getattr(infile, "name", "<unknown>")
    bit_definitions = list(read_patterns(infile))
    max_bit = max(bitdef.bit for bitdef in bit_definitions)
    
    return PatternFile(filename, max_bit, bit_definitions)

def read_patterns(infile):
    seen_bits = {}

    for lineno, line in enumerate(infile):
        lineno += 1

        # Leading and trailing whitespace is ignored
        line = line.strip()
        
        # Ignore blank lines or those with a leading "#"
        if not line or line.startswith("#"):
            continue

        # The first three columns, plus everything else for the description
        fields = line.split(None, 3)
        if len(fields) != 4:
            raise TypeError("Not enough fields on line %d: %r" % (lineno, line))

        # Normalize whitespace for the description
        fields[3] = " ".join(fields[3].split())

        # Do some type checking and error reporting
        bit, count, pattern, description = fields
        if not bit.isdigit():
            raise TypeError(
                "First field of line %d must be a non-negative bit position, not %r" %
                            (lineno, bit))
        bit = int(bit)

        if not count.isdigit() or int(count) == 0:
            raise TypeError(
                "Second field of line %d must be a positive minimum match count, not %r" %
                            (lineno, bit))
        count = int(count)
        
        if bit in seen_bits:
            raise TypeError("Line %d redefines bit %d, already set by line %d" %
                            (lineno, bit, seen_bits[bit]))
        seen_bits[bit] = lineno

        yield BitDefinition(bit, count, pattern, description, lineno)

class CountInfo(object):
    __slots__ = ("count", "bit", "byteno", "bitmask")
    def __init__(self, count, bit):
        self.count = count  # minimum count needed to enable this bit
        self.bit = bit  # used to set not_implemented, and useful for debugging
        
        # These simplify the fingerprint generation code
        self.byteno = bit//8
        self.bitmask = 1<<(bit%8)


def _bit_definition_to_pattern_definition(bit_definitions):
    "Helper function to organize the bit defintions based on pattern instead of bit"
    
    # A pattern definition is of the form:
    #  (pattern string, count_info_list)
    #     where the count_info list elements are sorted by count

    # I want to preserve the pattern order so that patterns which
    # are defined first are evaluated first
    ordered_patterns = []
    pattern_info = {}

    # Find all of the bit definitions for a given pattern
    for bitdef in bit_definitions:
        if bitdef.pattern not in pattern_info:
            pattern_info[bitdef.pattern] = []
            ordered_patterns.append(bitdef.pattern)
        pattern_info[bitdef.pattern].append( CountInfo(bitdef.count, bitdef.bit) )

    # Put them into a slighly more useful form
    #  - sorted now makes it easier to test when done
    #  - knowing the max match count lets some matchers optmize how to match
    for pattern in ordered_patterns:
        count_info_list = pattern_info[pattern]
        count_info_list.sort(key=lambda count_info: count_info.count)
        yield (pattern,
               count_info_list[-1].count,  # the largest count
               tuple(count_info_list)
               )
        
def _build_matchers(patterns, pattern_definitions, compile_pattern):
    not_implemented = set()
    matcher_definitions = []
    for (pattern, largest_count, count_info_tuple) in pattern_definitions:
        if pattern == "<0>":
            # Special case support for setting (or rather, ignoring) the 0 bit
            continue

        matcher = compile_pattern(pattern, largest_count)

        if matcher is NotImplemented:
            for count_info in count_info_tuple:
                not_implemented.add(count_info.bit)
            continue
        if matcher is NotImplementedError:
            raise AssertionError("Use 'NotImplemented' not 'NotImplementedError'")

        if matcher is None:
            raise UnsupportedPatternError(pattern)
        
        matcher_definitions.append( (matcher, largest_count, count_info_tuple) )

    return not_implemented, tuple(matcher_definitions)

def make_matchers(patterns, compile_pattern):
    pattern_definitions = _bit_definition_to_pattern_definition(patterns.bit_definitions)
    try:
        return _build_matchers(patterns, pattern_definitions, compile_pattern)
    except UnsupportedPatternError, err:
        err.filename = patterns.filename
        
        pattern = err.args[0]
        for bitdef in patterns.bit_definitions:
            if bitdef.pattern == pattern:
                err.lineno = bitdef.lineno
                raise
        raise
        

class PatternFingerprinter(object):
    def __init__(self, patterns, compile_pattern):
        self.patterns = patterns

        self.num_bytes = (patterns.max_bit // 8) + 1
        self.not_implemented, self.matcher_definitions = (
            make_matchers(patterns, compile_pattern)   )

    def describe(self, bit):
        description = self.patterns[bit].description
        if bit in self.not_implemented:
             description + " (NOT IMPLEMENTED)"
        return description

    def fingerprint(self, mol):
        raise NotImplemented("Must be implemented by a derived class")


def _load_named_patterns(name):
    filename = os.path.join(os.path.dirname(__file__), name + ".patterns")
    return load_patterns(filename)
