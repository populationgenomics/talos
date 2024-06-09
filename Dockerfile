FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:1.24.5

COPY requirements.txt .
COPY requirements-dev.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY talos talos/
RUN pip install .
