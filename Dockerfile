FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:511771e83afe64547a32eacb61ae4a216e4bcd85

COPY requirements.txt .
COPY requirements-dev.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY reanalysis reanalysis/
RUN pip install .
