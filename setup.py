"""
automated installation instructions
"""


from setuptools import find_packages, setup


with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()

setup(
    name='automated-interpretation-pipeline',
    description='CPG Variant Prioritisation',
    long_description=readme,
    version='0.4.0',
    author='Matthew Welland, CPG',
    author_email=(
        'matthew.welland@populationgenomics.org.au, '
        'cas.simons@populationgenomics.org.au'
    ),
    package_data={'reanalysis': ['templates/*.jinja']},
    url='https://github.com/populationgenomics/automated-interpretation-pipeline',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=['peddy==0.4.8', 'cyvcf2==0.30.18'],
    extras_require={
        'full': [
            'click',
            'cloudpathlib[all]>=0.9.0',
            'cpg-utils>=4.10.3',
            'cpg_workflows>=1.10.0',
            'dill>=0.3.5.1',
            'hail>=0.2.105',
            'Jinja2>=3.0.3',
            'networkx>=2.8.3',
            'obonet>=0.3.1',
            'pandas>=1.4.3',
            'requests>=2.25.1',
            'sample-metadata>=5.2.1',
            'seqr-loader>=1.2.5',
            'tabulate>=0.9.0',
            'toml>=0.10.2',
        ],
        'test': [
            'pytest>=7.2.1',
            'pytest-cov>=3.0.0',
            'pytest-xdist>=2.5.0',
            'requests_mock>=1.10.0',
        ],
    },
)
