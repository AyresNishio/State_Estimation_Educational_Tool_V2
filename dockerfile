From python:3.9-slim

# set working directory
WORKDIR /app

# install dependencies
COPY ./requirements.txt /app
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install build-essential -y
RUN pip install -r requirements.txt
# RUN pip install --upgrade numpy

# copy python app to the folder
COPY . /app

# start the server
CMD ["gunicorn", "-b", "0.0.0.0:8050", "main:server"]
