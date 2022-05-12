---
title: Q&A
nav_exclude: false
nav_order: 5
---

<!-- prettier-ignore-start -->
# Q&A
{: .fs-9 }
<!-- prettier-ignore-end -->

## How to set up Singularity env?

{: .note }
- SingularityCE must be installed using root privilage, but you do not need root privilage to run Singularity images. For HPC user, please contact your administrator to install the `singularityCE` first.
- More detailed instruction and follow:
  https://sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps

There are 3 steps to install singularity runtime.

### 1. install dependencies

On Debian-based systems, including Ubuntu:

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup
```

### 2. install golang

```bash
# adjust this as necessary
export VERSION=1.17.2 OS=linux ARCH=amd64
```

- downlaod golang release

```bash
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \ # Downloads the required Go package
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \ # Extracts the archive
    rm go$VERSION.$OS-$ARCH.tar.gz    # Deletes the ``tar`` file
```

- set ENV

```bash
echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc && \
    source ~/.bashrc
```

### 3. download and install SingularityCE from a release

- select version

```bash
# adjust this as necessary
# By default Singularity will be installed in the /usr/local directory hierarchy. You can specify a custom directory with the --prefix option
export VERSION=3.9.5 INSTALL_PREFIX=/usr/local
```

```bash
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION} &&
    ./mconfig --prefix=${INSTALL_PREFIX} && \
    make -C builddir && \
    sudo make -C builddir install
```

## How to use singularity?

The help command gives an overview of SingularityCE options and subcommands as follows:

```bash
singularity help
```

You can use the pull and build commands to download pre-built images from an external resource like the Container Library or Docker Hub.

```bash
singularity pull library://lolcow
```

You can also use pull with the `docker://` uri to reference Docker images served from a registry. In this case pull does not just download an image file.

```bash
singularity pull docker://sylabsio/lolcow
```

## Is it possible to run on HPC envirentment?

You don't need root private to run the pipeline, and dependency will be configured automatically, which make it easy for HPC user. But you need to ask the system adminstraotr to install the SingularityCE.
