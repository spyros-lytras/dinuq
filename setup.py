from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()
        
        
setup(
    name='dinuq',
    version='1.0.0',
    description='The dinucleotide quantification python package',
    #long_description='long_description',
    url='',
    classifiers=[
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python",
    ],
    keywords='bioinformatics dinucleotides viruses',
    project_urls={
    'CVR Bioinformatics': 'https://bioinformatics.cvr.ac.uk/',
    'Funding': 'https://mrc.ukri.org/',
    'MRC - University of Glasgow Centre for Virus Research': 'https://www.gla.ac.uk/researchinstitutes/iii/cvr/',
    },
    author='Spyros Lytras',
    author_email='s.lytras.1@research.gla.ac.uk',
    license='',
    packages=['dinuq'],
    install_requires=['biopython'],
    python_requires="~=3.5")