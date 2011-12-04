import chemfp
from chemfp import substruct
fps = chemfp.open_fps(None)

for fp, id in fps:
    print "==", id, fp.encode("hex")
    for i, c in enumerate(fp):
        if c == "\0":
            continue
        w = ord(c)
        for j in range(8):
            if w & (1<<j):
                bit = i*8+j
                print bit, substruct.describe_bit(bit)
