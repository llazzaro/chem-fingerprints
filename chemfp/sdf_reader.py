class Location(object):
    def __init__(self, lineno=0, title=None, filename=None):
        self.lineno = lineno
        self.title = title
        self.filename = filename
    def where(self):
        if self.filename is None:
            return "line {self.lineno}".format(self=self)
        else:
            return "line {self.lineno} of {self.filename!r}".format(self=self)
    def message(self):
        return "SD record ({title!r}) at {where}".format(
            title=self.title, where=self.where())

class LinesWithLocation(list):
    def __init__(self, filename):
        list.__init__(self)
        self.location = Location(filename=filename)
            

def read_sdf_records(infile, filename=None):
    """read_sdf_records(infile) -> iterator over SD record lines

    Iterate through each record in the SD file, returning the list of
    lines for each record.
    """
    if filename is None:
        filename = getattr(infile, "name")
    infile_reader = enumerate(infile)
    lines = LinesWithLocation(filename)
    while 1:
        loc = lines.location
        # Do a simple check that we are reading an SD file.
        # Get the 4th line. It must contain atom and bond counts
        # and the string "V2000" or "V3000"
        for i, (lineno, line) in zip("1234", infile_reader):
            if i == "1":
                loc.lineno = lineno+1
                loc.title = line.strip()
            lines.append(line)
        if not lines:
            # No lines read. End of file. Stop.
            break

        if len(lines) != 4:
            raise TypeError(
      "Unexpected end of file shortly after starting SD record at {where}".format(
                    where=loc.where()))

        count_line = lines[-1]
        if (("V2000" not in count_line and "V3000" not in count_line) or
            not count_line[0:3].strip().isdigit() or
            not count_line[3:6].strip().isdigit()):
            loc.lineno = lineno+1
            raise TypeError("Record at {where} is not a valid SD counts line".format(
                    where=loc.where()))

        # Likely good - read the rest of this record.
        for lineno, line in infile_reader:
            lines.append(line)
            if line.startswith("$$$$"):
                yield lines
                lines = LinesWithLocation(filename)
                break
        else:
            if lines:
                raise TypeError(
       "Unexpected end of file while parsing SD record {message}".format(loc.message()))

def _find_tag_value(record_lines, data_pattern):
    for lineno, line in enumerate(record_lines):
        if line.startswith(">") and data_pattern in line:
            data_value = record_lines[lineno+1]
            if data_value[:1] == ">" or data_value[:4] == "$$$$":
                # Empty tag value?
                return None
            return data_value.rstrip()

    # No matching tag found
    return None


def read_title_and_data_tag(infile, data_tag):
    """read_sdf_title_and_tag(infile) -> yield record titles and a tag value"""
    data_pattern = "<" + data_tag + ">"
    for record_lines in read_sdf_records(infile):
        title = record_lines[0].strip()
        if not title:
            title = None
        data_tag_value = _find_tag_value(record_lines, data_pattern)
        yield record_lines.location, title, data_tag_value

def _find_two_tag_values(record_lines, title_pattern, data_pattern):
    # Use -1 as a flag to mean "has not been found yet"
    title_value = data_value = -1
    for lineno, line in enumerate(record_lines):
        if line.startswith(">"):
            if title_pattern in line:
                value = record_lines[lineno+1]
                if value[:1] == ">" or value[:4] == "$$$$":
                    title_value = None
                else:
                    title_value = value.rstrip()
                if data_value is not -1:
                    # Have seen both fields at least once
                    return title_value, data_value

            elif data_pattern in line:
                value = record_lines[lineno+1]
                if value[:1] == ">" or value[:4] == "$$$$":
                    data_value = None
                else:
                    data_value = value.rstrip()
                if title_value is not -1:
                    # Have seen both fields at least once
                    return title_value, data_value

    # Didn't find both
    if title_value is -1:
        title_value = None
    if data_value is -1:
        data_value = None
    return title_value, data_value

def read_title_tag_and_data_tag(infile, title_tag, data_tag):
    title_tag_pattern = "<" + title_tag + ">"
    data_tag_pattern = "<" + data_tag + ">"
    for record_lines in read_sdf_records(infile):
        title_value, data_value = _find_two_tag_values(
            record_lines, title_tag_pattern, data_tag_pattern)
        yield record_lines.location, title_value, data_value
