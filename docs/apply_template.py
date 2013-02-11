import sys
import os
from os.path import join, exists, getmtime
import inspect

import chemfp
import chemfp.arena
import chemfp.bitops

name = sys.argv[1]
template_name = name + ".txt"
rst_name = name + ".rst"

from jinja2 import Template, Environment, BaseLoader, TemplateNotFound

class LocalDirLoader(BaseLoader):
    def __init__(self):
        self.path = "."
    
    def get_source(self, environment, template):
        path = join(self.path, template)
        if not exists(path):
            raise TemplateNotFound(template)
        mtime = getmtime(path)
        with file(path) as f:
            source = f.read().decode('utf-8')
        return source, path, lambda: mtime == getmtime(path)
        

def docstring(f):
    doc = inspect.getdoc(f)
    if doc is None:
        raise AssertionError("Missing docstring for " + f.__name__)
    doc = inspect.getdoc(f).rstrip("\n")
    doc = "    " + doc.replace("\n", "\n    ")
    return doc

env = Environment(loader = LocalDirLoader())
env.filters['docstring'] = docstring
template = env.get_template(template_name)

page = template.render({"chemfp": chemfp})

with open(rst_name, "w") as f:
    f.write(page)
