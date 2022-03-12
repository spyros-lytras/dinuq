from setuptools import setup, Extension


def readme():
    with open('README.md') as f:
        return f.read()
        
        
setup(
    name='dinuq',
    version='1.2.0',
    description='The Dinucleotide Quantification Python package',
    long_description= readme(),
    long_description_content_type='text/markdown',
    url='https://spyros-lytras.github.io/dinuq/',
    download_url = 'https://github.com/spyros-lytras/dinuq/dist/dinuq-1.2.0.tar',
    classifiers=[
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python",
    ],
    keywords='bioinformatics dinucleotides viruses',
    project_urls={
    'Documentation': 'https://spyros-lytras.github.io/dinuq/',
    'CVR Bioinformatics': 'https://bioinformatics.cvr.ac.uk/',
    },
    author='Spyros Lytras',
    author_email='s.lytras.1@research.gla.ac.uk',
    license='MIT',
    packages=['dinuq'],
    install_requires=['biopython'],
    python_requires="~=3.6")