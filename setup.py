"""
setup.py for the talos package
"""

from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list


    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [line.strip() for line in filehandler if line.strip() and not line.startswith('#')]


setup(
    name='talos',
    description='Centre for Population Genomics Variant Prioritisation',
    long_description=readme,
    version='4.0.1',
    author='Matthew Welland, CPG',
    author_email='matthew.welland@populationgenomics.org.au, cas.simons@populationgenomics.org.au',
    package_data={'talos': ['templates/*.jinja', 'reanalysis_global.toml']},
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
    install_requires=read_reqs('requirements.txt'),
    extras_require={
        'test': read_reqs('requirements-dev.txt'),
        'cpg': read_reqs('requirements-cpg.txt'),
    },
    entry_points={
        'console_scripts': [
            'generate_pedigree = talos.cpg_generate_pheno_ped:main',
            'vcf_to_mt = talos.vep_vcf_to_mt:main',
            'report_hunter = helpers.report_hunter:run_both',
            'hail_label = talos.hail_filter_and_label:cli_main',
            'hail_label_sv = talos.hail_filter_sv:cli_main',
            'hpo_panel_match = talos.hpo_panel_match:cli_main',
            'query_panelapp = talos.query_panelapp:cli_main',
            'validate_categories = talos.validate_categories:cli_main',
            'build_html = talos.html_builder:main',
        ]
    },
)
