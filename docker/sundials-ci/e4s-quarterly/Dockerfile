FROM ghcr.io/llnl/sundials-ci-e4s-base:e4s-22.05 as bootstrap

ARG spack_yaml

# Build stage with Spack pre-installed and ready to be used
FROM bootstrap as builder

# Install python environment for docs
RUN python3 -m venv /opt/python-venv/sundocs && \
    . /opt/python-venv/sundocs/bin/activate && \
    pip install wheel && \
    pip install sphinx >=4.0.0 && \
    pip install sphinx-fortran && \
    pip install sphinx_rtd_theme && \
    pip install sphinxcontrib.bibtex && \
    pip install sphinx-copybutton

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment
COPY $spack_yaml /opt/spack-environment/spack.yaml

# Install the spack environment
RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    cd /opt/spack-environment                                           \
 && spack env activate .                                                \
 && spack buildcache keys -it                                           \
 && cd /opt/spack-environment                                           \
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

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Add spack installed binaries to path
ENV PATH="/opt/view/bin:${PATH}"

# Allow mpirun as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
# Suppress HWLOC errors
ENV HWLOC_HIDE_ERRORS=2

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]
