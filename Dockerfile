# syntax=docker/dockerfile:1

ARG PYTHON_VERSION=3.11.9

FROM python:${PYTHON_VERSION}-slim

LABEL fly_launch_runtime="flask"

WORKDIR /app

COPY requirements.txt 
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8080

CMD ["gunicorn", "-b", "0.0.0.0:8080", "app:app", "--workers", "1", "--timeout", "120"]
