# Development and testing Docker image for PubChem Finder
FROM continuumio/miniconda3:latest

LABEL maintainer="pubchem-finder"
LABEL description="Development environment for PubChem Finder with RDKit"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install RDKit via conda (easiest method)
RUN conda install -c conda-forge rdkit python=3.11 -y \
    && conda clean -a -y

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:${PATH}"

# Set working directory
WORKDIR /app

# Copy project files
COPY pyproject.toml ./
COPY README.md ./

# Install project dependencies
# Note: We use --no-root to install deps without the package itself yet
RUN poetry config virtualenvs.create false \
    && poetry install --no-interaction --no-ansi --no-root

# Copy source code
COPY src/ ./src/
COPY tests/ ./tests/

# Install the package
RUN poetry install --no-interaction --no-ansi

# Default command: run tests
CMD ["pytest", "--cov", "-v"]
