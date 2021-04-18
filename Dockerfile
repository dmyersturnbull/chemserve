# https://stackoverflow.com/questions/53835198/integrating-python-poetry-with-docker/54763270#54763270

# Add a build argument that should either be 'true' or 'false'
# If true, installs without dev dependencies (poetry --no-dev)
ARG NO_DEV

FROM continuumio/miniconda3


# --------------------------------------
# ------------- Set labels -------------

# See https://github.com/opencontainers/image-spec/blob/master/annotations.md
LABEL name="chemserve"
LABEL version="0.2.0"
LABEL vendor="dmyersturnbull"
LABEL org.opencontainers.image.title="chemserve"
LABEL org.opencontainers.image.version="0.2.0"
LABEL org.opencontainers.image.url="https://github.com/dmyersturnbull/chemserve"
LABEL org.opencontainers.image.documentation="https://github.com/dmyersturnbull/chemserve"


# --------------------------------------
# ---------- Install from environment ----------

# ENV no longer adds a layer in new Docker versions,
# so we don't need to chain these in a single line
ENV NO_DEV=${NO_DEV}
ENV PYTHONFAULTHANDLER=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONHASHSEED=random
ENV PIP_NO_CACHE_DIR=off
ENV PIP_DISABLE_PIP_VERSION_CHECK=on
ENV PIP_DEFAULT_TIMEOUT=120

# Install system deps
RUN pip install "poetry>=1.1.6,<2.0"

COPY environment.yml .
RUN conda env create -f server-environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]


# --------------------------------------
# ---------- Copy and install ----------

# Copy only requirements to cache them in docker layer
WORKDIR /code
COPY poetry.lock pyproject.toml /code/

# Copy to workdir
COPY . /code

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "chemserve", "chemserve"]
