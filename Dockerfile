FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:07f8422e67880a684e593e5e22cdf7ec5b566ae0

COPY requirements.txt .
COPY requirements-dev.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY reanalysis reanalysis/
RUN pip install .
