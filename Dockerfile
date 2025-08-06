FROM ghcr.io/osgeo/gdal:alpine-normal-latest
WORKDIR /app

# Install build dependencies and pip
RUN apk update && apk add --no-cache \
    build-base \
    geos-dev \
    python3-dev \
    py3-pip \
    musl-dev \
    linux-headers

COPY requirements.txt .

RUN python3 -m venv /venv && \
    . /venv/bin/activate && \
    pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["/venv/bin/python", "combine_capa_with_streets.py"]