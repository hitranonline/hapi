# -*- coding: utf-8 -*-

#from distutils.core import setup
from setuptools import setup
from hapi.hapi import HAPI_VERSION, HAPI_HISTORY

import re
import sys
import shutil

# Update the README file
shutil.copy('README.md','README.bak')

with open('README.md') as f:
    README = f.read()
    
main_body = re.search('(## Introduction.+)$',README,re.MULTILINE|re.DOTALL).groups(0)[0].rstrip()
    
README = """
# HITRAN Application Programming Interface (HAPI)
===============================================

Current version: {hapi_version}

## Version history

{history}

{body}

""".format(hapi_version=HAPI_VERSION,body=main_body,history=HAPI_HISTORY.strip())

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
    license='MIT',
)