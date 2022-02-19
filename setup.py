import os
from setuptools import setup, find_packages, Command
import subprocess
import yaml
from modules._version import __version__
from setuptools.command.install import install
import pip
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

class CustomInstall(install):
    '''
    '''

    def __init__(self, dist):
        super(install, self).__init__(dist)

    # user_options = []
    user_options = install.user_options + [
        ('custom-option=', None, 'Path to something')
    ]
    # This method must be implemented
    def initialize_options(self):
        pass

    # This method must be implemented
    def finalize_options(self):
        pass

    def run(self):
        '''
        '''
        install.run(self)

setup(
    name = "grapes2",
	spython_requires='>=3.6',
    version = __version__,
    author = "Bernat del Olmo",
    author_email = "bioinformatix@gencardio.com, bdelolmo@gencardio.com",
    description = ("GRAPES2 CNV/SV detection on gene panel sequencing"),
    keywords = "example documentation tutorial",
    url = "https://github.com/GENCARDIO/GRAPES2",
    # packages=['an_example_pypi_project', 'tests'],
    install_requires=[
        'pybedtools>=0.8.2',
        'pysam>=0.16.0.1',
        'matplotlib>=3.3.3',
        'numpy>=1.19.5',
        'scipy>=1.5.4',
        'seaborn>=0.11.1',
        'pandas>=1.1.5',
        'hmmlearn==0.2.7',
        'docx>=0.2.4',
        'PyYAML',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    # cmdclass={'install': CustomInstall},
)
