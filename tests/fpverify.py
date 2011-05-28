# Verify
import re
import itertools
import string

from chemfp.decoders import from_cactvs
from chemfp import types


# Parse terms which look like "34","82-88", and "0-4,9,23-880"

_range_pat = re.compile("(\d+)-(\d+)")
_bit_pat = re.compile("\d+")

def parse_bitset_definition(s):
    if s.endswith("=="):
        num_bits, decoded = from_cactvs(s)
        return set(bitno for (bitno, bit) in iter_bits(decoded) if bit)

    if s.startswith("0x"):
        decoded = s[2:].decode("hex")
        return set(bitno for (bitno, bit) in iter_bits(decoded) if bit)

    bits = set()
    if s == "X":
        return bits

    start = 0
    while 1:
        m = _range_pat.match(s, start)
        if m:
            # Ranges are inclusive, so "3-5" means {3, 4, 5}
            start = int(m.group(1))
            end = int(m.group(2))
            assert start <= end, (start, end)
            for bit in range(start, end+1):
                bits.add(bit)
            start = m.end()
        else:
            m = _bit_pat.match(s, start)
            if m:
                # a single bit
                bits.add(int(m.group()))
                start = m.end()
            else:
                raise ValueError("Unknown term %r position %d" % (s, start+1))
        if start == len(s):
            break
        if s[start] != ",":
            raise ValueError("Must have a ',' in %r at position %d" % (s, start+1))
        start += 1

    return bits

def parse_true(true_term, num_bits):
    if true_term.endswith("?"):
        assert "!" not in true_term
        true_bits = parse_bitset_definition(true_term[:-1])
        ignore_bits = set(range(num_bits)) - true_bits
        return true_bits, ignore_bits

    if "!" not in true_term:
        return parse_bitset_definition(true_term), set()

    left, right = true_term.split("!")
    return parse_bitset_definition(left), parse_bitset_definition(right)
    

def parse_false(false_term, num_bits, true_bits):
    if false_term == "*":
        return set(range(num_bits)) - true_bits
    return parse_bitset_definition(false_term)


_bit_offset_table = {}
for i in range(256):
    _bit_offset_table[chr(i)] = tuple(offset for offset in range(8) if i & (1<<offset))

def iter_bits(fp):
    for byteno, byte in enumerate(fp):
        base = byteno*8
        word = ord(byte)
        for j in (0,1,2,3,4,5,6,7):
            yield base+j, bool(word&(1<<j))


def get_fingerprinter(line, toolkit):
    fields = line.split()
    family = fields[0].format(TOOLKIT=toolkit)
    fields[0] = family
    return types.parse_type(" ".join(fields))

def evaluate_test_cases(fingerprinter, test_cases):
    f = open("tmp.smi", "w")
    for i, test_case in enumerate(test_cases):
        f.write("%s %d\n" % (test_case[0], i))
    f.close()

    errors = []
    tested_true = set()
    tested_false = set()

    reader = fingerprinter.read_structure_fingerprints("tmp.smi", "smi")
    for i, (test_case, fp_result) in enumerate(
                           itertools.izip_longest(test_cases, reader) ):
        assert test_case is not None, "too many results?"
        smiles, true_bits, false_bits, ignore_bits = test_case
        print "SMILES", smiles
        assert fp_result is not None, "not enough results?"
        fp, title = fp_result
        assert str(i) == title, (i, title)


        for bitno, val in iter_bits(fp):
            if bitno in ignore_bits:
                continue
            if bitno in true_bits:
                if not val:
                    errors.append( (bitno, val, True) )
                else:
                    tested_true.add(bitno)
            elif bitno in false_bits:
                if val:
                    errors.append( (bitno, val, False) )
                else:
                    tested_false.add(bitno)

        if errors:
            print "ERROR: Fingerprint failure:", smiles
            print " bit#  got  expected"
            print " ----  ---  --------"
            for (bitno, got, expected) in errors:
                print " %3d    %d      %d" % (bitno, got, expected)
            print "bit pattern", ",".join(str(bitno) for (bitno, val) in iter_bits(fp) if val)
            print "hex pattern", fp.encode("hex")
            raise SystemExit()

    all = set(range(fingerprinter.num_bits))
    print "Missing true bits:", " ".join(str(b) for b in sorted(all - tested_true))
    print "Missing false bits:", " ".join(str(b) for b in sorted(all - tested_true))
    
        

def main():
    import sys
    if len(sys.argv) > 1:
        toolkit = sys.argv[1]
        main2(toolkit)
    else:
        for toolkit in ("OpenEye", "RDKit", "OpenBabel"):
            main2(toolkit)

def main2(toolkit):
    fingerprinter = None

    num_bits = None
    toolkit_version = None
    
    test_cases = []
    skip = False

    def toolkit_selected(options):
        print toolkit, toolkit_version, options
        return toolkit in options or toolkit_version in options
    
    for line in open("rdmaccs.fpverify", "U"):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith("#type="):
            assert fingerprinter is None
            fingerprinter = get_fingerprinter(line[6:], toolkit)
            num_bits = fingerprinter.num_bits
            toolkit_version = fingerprinter.software.split()[0]
            continue

        if line.startswith("#"):
            continue

        fields = line.split()
        
        if fields[0] == "=skip":
            if toolkit_selected(fields[1:]):
                print "skip", fields[1:]
                skip = True
            continue

        if fields[0] == "=only":
            if not toolkit_selected(fields[1:]):
                skip = True
            continue

        if fields[0].startswith("="):
            if fields[0][1:2] in string.ascii_uppercase:
                # Some sort of toolkit-specific directive. Ignore.
                continue
            raise AssertionError(line)
        
        assert fingerprinter is not None # haven't defined the fingerprint type

        if skip:
            skip = False
            continue


        if len(fields) == 1:
            raise AssertionError(line)
        if len(fields) == 2:
            smiles, true_term = fields
            false_term = "*"
        else:
            smiles, true_term, false_term = fields # must have 2 or 3 terms

        true_bits, ignore_bits = parse_true(true_term, num_bits)
        false_bits = parse_false(false_term, num_bits, true_bits)

        true_bits = true_bits - ignore_bits
        false_bits = false_bits - ignore_bits

        assert not (true_bits & false_bits), (true_bits & false_bits)

        test_cases.append( (smiles, true_bits, false_bits, ignore_bits) )

    evaluate_test_cases(fingerprinter, test_cases)

if __name__ == "__main__":
    main()
