FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

COPY requirements.txt .
COPY requirements-dev.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY reanalysis reanalysis/
RUN pip install .
