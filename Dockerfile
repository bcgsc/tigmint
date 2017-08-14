FROM linuxbrew/linuxbrew
LABEL maintainer="Shaun Jackman <sjackman@gmail.com>"

RUN brew tap homebrew/science \
	&& brew install bedtools bwa gawk gnu-sed miller pigz r samtools seqtk
RUN Rscript -e 'install.packages(c("ggplot2", "rmarkdown", "tidyverse", "uniqtag"), repos = c(CRAN = "https://cran.rstudio.com"))'
WORKDIR /home/linuxbrew/tigmint
ADD . .
RUN sudo chown -R linuxbrew: .
ENV PATH="/home/linuxbrew/tigmint:$PATH"
