FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:2d3b87a27ffa686fafb1e7216ca99a4f029d5585

COPY requirements.txt .
COPY requirements-dev.txt .
RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY helpers helpers/
COPY reanalysis reanalysis/
RUN pip install .
