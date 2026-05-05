# platelet-wcm runtime image
# Python + all pip dependencies; no project source code.
#
# Build locally from the platelet-wcm/ project root:
#
#     docker build -f docker/runtime/Dockerfile -t platelet-wcm-runtime .

ARG from=python:3.11.5-slim
FROM ${from}

RUN echo "alias ls='ls --color=auto'" >> ~/.bashrc \
    && echo "alias ll='ls -l'" >> ~/.bashrc \
    && cp ~/.bashrc /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       gcc g++ gfortran make curl \
    && rm -rf /var/lib/apt/lists/*

# Consistent OpenBLAS results
ENV OPENBLAS_NUM_THREADS=1

COPY requirements.txt /
RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
    && pip install --no-cache-dir -r /requirements.txt

CMD ["/bin/bash"]
