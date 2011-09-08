import itertools

def mutual_exclusion(parser, args, default, groups):
    true_groups = []
    for g in groups:
        if getattr(args, g):
            true_groups.append(g)

    if not true_groups:
        setattr(args, default, True)
    elif len(true_groups) == 1:
        pass
    else:
        parser.error("Cannot specify both --%s and --%s" % (true_groups[0], true_groups[1]))

def sys_exit_opener(opener, source, format, id_tag, aromaticity):
    try:
        return opener.read_structure_fingerprints(source, format, id_tag, aromaticity)
    except (IOError, oe.UnknownFormat, TypeError), err:
        raise SystemExit("Cannot read structure fingerprints: %s\n" % err)

def iter_all_sources(opener, filenames, format, id_tag, aromaticity):
    for filename in filenames:
        reader = sys_exit_opener(opener, filename, format, id_tag, aromaticity)
        for x in reader:
            yield x

def read_multifile_structure_fingerprints(opener, filenames, format, id_tag, aromaticity):
    if not filenames:
        reader = sys_exit_opener(opener, None, format, id_tag, aromaticity)
        return reader.metadata, reader

    reader = sys_exit_opener(opener, filenames[0], format, id_tag, aromaticity)
    if len(filenames) == 1:
        return reader.metadata, reader

    reader = sys_exit_opener(opener, filenames[0], format, id_tag, aromaticity)
    reader.metadata.sources = filenames
    return reader.metadata, itertools.chain(reader, iter_all_sources(opener, filenames[1:], format, id_tag, aromaticity))

