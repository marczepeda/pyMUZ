from setuptools import find_packages, setup
import os
work_dir = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(work_dir,'README.md'),'r') as f:
    long_description = f.read()

setup(
    name='pyMUZ',
    version='0.0.1',
    packages=find_packages(),
    description='Python Tools for MUZ',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/marczepeda/python',
    author='Marcanthony Zepeda',
    author_email='marczepeda17@gmail.com',
    license='GPLv3+',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3.11',
        'Operating System :: OS Independent',
    ],
    install_requires=['numpy',
                      'matplotlib',
                      'seaborn',
                      'pandas',
                      'biopython',
                      'requests',
                      'selenium'],
    python_requires='>=3.11'
)