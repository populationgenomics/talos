FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows

COPY requirements.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY reanalysis reanalysis/
RUN pip install .
