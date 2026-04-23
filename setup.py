# -*- coding: utf-8 -*-

#from distutils.core import setup
from setuptools import setup

import ast
import re
import sys
import shutil


def read_package_metadata():
    with open('hapi/hapi.py', 'r', encoding='utf-8') as f:
        module_ast = ast.parse(f.read(), filename='hapi/hapi.py')

    metadata = {}
    for node in module_ast.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name) and target.id in {'HAPI_VERSION', 'HAPI_HISTORY'}:
            metadata[target.id] = ast.literal_eval(node.value)

    return metadata['HAPI_VERSION'], metadata['HAPI_HISTORY']


HAPI_VERSION, HAPI_HISTORY = read_package_metadata()

# Update the README file
shutil.copy('README.md','README.bak')

with open('README.md') as f:
    README = f.read()
    
main_body = re.search('(## Introduction.+)$',README,re.MULTILINE|re.DOTALL).groups(0)[0].rstrip()
    
HISTORY_LIST = ''
for i,item in enumerate(HAPI_HISTORY):
    HISTORY_LIST += '  %d) %s\n'%(i+1,item)
    
README = """
# HITRAN Application Programming Interface (HAPI)
===============================================

Current version: {hapi_version}

## Version history

{history}

{body}

""".format(hapi_version=HAPI_VERSION,body=main_body,history=HISTORY_LIST.rstrip())

# Python 2 and 3 encoding support
if sys.version_info[0]<3:
    with open('README.md','w') as f:
        f.write(README.lstrip())
else:
    with open('README.md','wb') as f:
        f.write(README.lstrip().encode('utf-8'))
    
# Install the package
setup(
    name='hitran-api',
    version=HAPI_VERSION,
    packages=['hapi',],
    install_requires=['numpy'],
    license='MIT',
)