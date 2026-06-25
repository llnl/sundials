FROM ubuntu:latest AS bootstrap

ARG spack_version
ARG spack_packages_ref=develop
ARG spack_packages_commit
ARG spack_yaml

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    SPACK_ROOT=/opt/spack

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        bash \
        build-essential \
        bzip2 \
        ca-certificates \
        coreutils \
        csh \
        curl \
        file \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg2 \
        gzip \
        iproute2 \
        locales \
        lua-posix \
        make \
        patch \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-venv \
        python3-wheel \
        tar \
        unzip \
        xz-utils \
        zstd \
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*

RUN test -n "${spack_version}" \
 && git clone --depth 1 --branch "${spack_version}" \
        https://github.com/spack/spack.git "${SPACK_ROOT}" \
 && git clone --depth 1 --branch "${spack_packages_ref}" \
        https://github.com/spack/spack-packages.git /opt/spack-packages \
 && if [ -n "${spack_packages_commit}" ]; then \
        git -C /opt/spack-packages checkout "${spack_packages_commit}"; \
    fi \
 && mkdir -p "${SPACK_ROOT}/etc/spack" \
 && printf "repos:\n  builtin: /opt/spack-packages/repos/spack_repo/builtin\n" \
        > "${SPACK_ROOT}/etc/spack/repos.yaml"

ENV PATH="/opt/view/bin:${SPACK_ROOT}/bin:${PATH}"

RUN spack --version \
 && spack repo list \
 && spack compiler find \
 && spack compiler list

FROM bootstrap AS builder

RUN python3 -m venv /opt/python-venv/sundocs \
 && . /opt/python-venv/sundocs/bin/activate \
 && pip install wheel \
 && pip install "sphinx>=4.0.0" \
 && pip install sphinx-fortran \
 && pip install sphinx_rtd_theme \
 && pip install sphinxcontrib.bibtex \
 && pip install sphinx-copybutton

RUN mkdir /opt/spack-environment
COPY ${spack_yaml} /opt/spack-environment/spack.yaml

RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    cd /opt/spack-environment                                           \
 && spack env activate .                                                \
 && spack buildcache keys -it                                           \
 && status=0                                                            \
 && spack install --reuse --no-check-signature                          \
    || spack -d -v install --reuse --no-check-signature                 \
    || spack -d -v install --reuse --no-check-signature                 \
    || status=$?                                                        \
 && mkdir -p /var/cache/spack-mirror/build_cache                        \
 && spack buildcache update-index -d /var/cache/spack-mirror            \
 && spack buildcache list --allarch --very-long                         \
    | sed '/^$/d;/^--/d;s/@.\+//;s/\([a-z0-9]*\) \(.*\)/\2\/\1/'        \
    | sort > tmp.buildcache.txt                                         \
 && spack find --format {name}/{hash} | sort                            \
    | comm -23 - tmp.buildcache.txt                                     \
    | xargs --no-run-if-empty                                           \
      spack buildcache create --allow-root --only package --unsigned    \
                              --directory /var/cache/spack-mirror       \
                              --rebuild-index                           \
 && spack clean -a                                                      \
 && exit $status

RUN cd /opt/spack-environment \
 && spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    HWLOC_HIDE_ERRORS=2

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]
