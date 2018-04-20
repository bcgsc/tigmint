FROM linuxbrew/linuxbrew
LABEL maintainer="Shaun Jackman <sjackman@gmail.com>"

WORKDIR /home/linuxbrew/tigmint
ADD . .
RUN sudo chown -R linuxbrew: .
RUN brew bundle && pip3 install -r requirements.txt
ENV PATH="/home/linuxbrew/tigmint/bin:$PATH"
