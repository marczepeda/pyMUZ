from setuptools import find_packages, setup

with open('README.md','r') as f:
    long_description = f.read()

setup(
    name='MUZ',
    version='0.0.0',
    description='README.md',
    package_dir={'general':'MUZ/general',
                 'biology':'MUZ/biology'
                 },
    packages=find_packages(where='MUZ',
                           exclude=['old','revise'],
                           include=['general','biology']),
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/marczepeda/python',
    author='Marcanthony Zepeda',
    author_email='marczepeda17@gmail.com',
    license='GPLv3+',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3.10',
        'Operating System :: OS Independent',
    ],
    install_requires=["numpy",
                      "matplotlib",
                      "seaborn",
                      "pandas",
                      "biopython"],
    python_requires='>=3.10'
)