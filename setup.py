from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

import secmet

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='secmet',
    version=secmet.__version__,
    url='http://github.com/kblin/python-secmet/',
    license='Apache Software License',
    author='Kai Blin',
    tests_require=['pytest'],
    install_requires=['biopython>=1.65'],
    cmdclass={'test': PyTest},
    author_email='kblin@biosustain.dtu.dk',
    keywords = "bioinformatics",
    description='A library to handle secondary metabolite annotations',
    long_description=long_description,
    packages=['secmet'],
    include_package_data=True,
    platforms='any',
    test_suite='tests',
    classifiers = [
        'Programming Language :: Python',
        "Development Status :: 3 - Alpha",
        'Natural Language :: English',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
