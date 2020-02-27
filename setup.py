from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()
        
        
setup(
    name='dinuq',
    version='1.0.0',
    description='The Dinucleotide Quantification Python package',
    #long_description='long_description',
    url='https://github.com/spyros-lytras/dinuq',
    download_url = 'https://github.com/spyros-lytras/dinuq/archive/v1.0.0.tar.gz',
    classifiers=[
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python",
    ],
    keywords='bioinformatics dinucleotides viruses',
    project_urls={
    'CVR Bioinformatics': 'https://bioinformatics.cvr.ac.uk/',
    },
    author='Spyros Lytras',
    author_email='s.lytras.1@research.gla.ac.uk',
    license='',
    packages=['dinuq'],
    install_requires=['biopython'],
    python_requires="~=3.5")