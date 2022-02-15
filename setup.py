#!/usr/bin/env python
# -*- encoding: utf-8 -*-

# Dynamic metadata (setup.py): possibly non-deterministic. Any items that are
# dynamic or determined at install-time, as well as extension modules or
# extensions to setuptools, need to go into setup.py.

import io
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='gepard',
    version='0.9.9',
    license='AGPL-3.0',
    description='Tool for studying the 3D quark and gluon distributions in the nucleon',
    long_description='{}'.format(
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.md'))
    ),
    author='Krešimir Kumerički',
    author_email='kkumer@phy.hr',
    url='https://gepard.phy.hr',
    project_urls={
        "Sources": "https://github.com/kkumer/gepard",
        "Bug Tracker": "https://github.com/kkumer/gepard/issues",
        'Changelog': 'https://github.com/kkumer/gepard/blob/master/CHANGELOG.rst',
    },
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    # include_package_data=True,
    package_data = {
        '': ['*.dat'],
    },
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: Implementation :: CPython',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: PyPy',
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=[
        'Physics', 'Particle Physics',
    ],
    python_requires='>=3.7',
    install_requires=['importlib-resources'
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'pytest-runner',
    ],
    # entry_points={
        # 'console_scripts': [
            # 'nameless = nameless.cli:main',
        # ]
    # },
)
