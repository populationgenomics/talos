"""
automated installation instructions
moving away from requirements.txt
"""
from setuptools import find_packages, setup

setup(
    name='automated-interpretation-pipeline',
    version='0.0.1',
    author='Matthew Welland, CPG',
    author_email=(
        'matthew.welland@populationgenomics.org.au, '
        'cas.simons@populationgenomics.org.au'
    ),
    url='https://github.com/populationgenomics/automated-interpretation-pipeline',
    description='Variant annotation and prioritisation process',
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
    install_requires=[
        'click==8.0.4',
        'cloudpathlib[all]==0.7.0',
        'cpg-utils==4.4.0',
        'cyvcf2==0.30.14',
        'dill==0.3.4',
        'hail==0.2.93',
        'Jinja2==3.0.3',
        'pandas==1.3.5',
        'peddy==0.4.8',
        'requests==2.25.1',
        'sample-metadata==4.10.0',
        'seqr-loader==1.2.5',
    ],
    extras_require={
        'test': [
            'pytest>=7.0.0',
            'pytest-cov==3.0.0',
            'tabulate==0.8.9',
            'pytest-xdist==2.5.0',
            'requests_mock==1.9.3',
        ]
    },
)
