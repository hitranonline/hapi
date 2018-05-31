# -*- coding: utf-8 -*-

#from distutils.core import setup
from setuptools import setup
from hapi.hapi import HAPI_VERSION

setup(
    name='hitran-api',
    version=HAPI_VERSION,
    packages=['hapi',],
    license='MIT',
)