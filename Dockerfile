FROM ubuntu:25.04

# Set environment variables to avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    python3-full \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libhts-dev \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python packages in the virtual environment
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pysam pandas

# Set up the app directory
WORKDIR /app

# Copy your application files
COPY . /app

# Default command (using the Python from the virtual environment)
# Using the full path to the Python binary in the virtual environment
CMD ["/opt/venv/bin/python", "LOD_edit.py"]