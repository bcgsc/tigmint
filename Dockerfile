FROM linuxbrew/linuxbrew
LABEL maintainer="Shaun Jackman <sjackman@gmail.com>"

WORKDIR /home/linuxbrew/tigmint
COPY . .
RUN sudo chown -R linuxbrew: . \
	&& brew bundle \
	&& rm -rf /home/linuxbrew/.cache \
	&& pip3 install cython git+https://github.com/daler/pybedtools.git \
	&& pip3 install -r requirements.txt
ENV PATH="/home/linuxbrew/tigmint/bin:$PATH"
