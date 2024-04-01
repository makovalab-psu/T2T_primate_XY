## StitchOrthologs.py

This README describes the usage for `StitchOrthologs.py` which was used to identify single-copy ortholog segments in a multispecies, chromosome-scale nucleotide alignment and output them as a `fasta`-formatted alignment. 

**Usage:** `StitchOrthologs.py <in.maf> <out.fa> <out.csv>` 

To generate the input `MAF` alignment for `StitchOrthologs.py`, `hal2maf` (included in the `cactus` [github repositiory](https://github.com/ComparativeGenomicsToolkit/cactus)) was executed on a `cactus` alignment graph as the following: 

`hal2maf --noAncesetors --onlyOrthologs --refGenome <reference genome> <cactus.hal> <out.maf>`

An example of the `MAF` output from `hal2maf` can be found below, where `--refGenome` was set to `chm13` (first alignment block truncated for visualization):
```
##maf version=1 scoring=N/A
# hal (((Gorilla_gorilla:1,((Pan_troglodytes:1,Pan_paniscus:1)Anc5:1,(chm13:1,GRCh38:1)Anc6:1)Anc4:1)Anc2:1,(Pongo_abelii:1,Pongo_pygmaeus:1)Anc3:1)Anc1:1,Symphalangus_syndactylus:1)Anc0;

a
s       chm13.chrX      2394410 376     +       154259566       tcacttatgtctttagataaatgcacacacatatatccacatagcttggaaggtatataagctctggaaaactataattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       GRCh38.chrX     2781479 376     +       156040895       tcacttatgtctttagataaaTGCACACACATATCTCCACATAGCTTGGAAGGTATATAAGCTCTGGAAAACtataattttgagttagtctggtgataatttccaggccttctccctgtaaCA...
s       Gorilla_gorilla.chrX    12048775        376     +       177553137       tcacttatgtctttagataaatgcaCACACATATATCCACATAGCTTGGAAGGTATATAAGCTCTGGAAAACtataattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       Pan_paniscus.chrX       2527836 376     +       160249396       tcacttatgtctttagataaatgcacacacatatatccacataGCTTGGAAGGTATATAAGCTCTGGAAAACtAtaattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       Pan_troglodytes.chrX    3166326 376     +       153892427       tcagttatgtctttagataaatgcacacacatatatccacataGCTTGGAAGGTATATAAGCTCTGGAAAACTAtaattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       Pongo_abelii.chrX       2382228 376     +       162586350       tcacttatgtctttagataaatgcacacacatatatccacataGCTTGGAAGGTATATAAGCTCTGGAAAAAtataattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       Pongo_pygmaeus.chrX     2356735 375     +       161100531       tcacttatgtctttagataaatgcacacacatatatccacataGCTTGGAAGGTATATAAGCTCTGGAAAAAtataattttgagttagtctggtgataatttccaggccttctccctgtaaca...
s       Symphalangus_syndactylus.chrX   16744543        376     +       165588543       tcatttatgtctttaggtaaatgcacacACATATATCCACATAGCTTGGAAGGTGTATAAGCTCTGGAAAACtgtaattttgagttagtctggtgataatttccaggccttctccctgtaaca...

a
s       chm13.chrX      2394786 2       +       154259566       cc
s       GRCh38.chrX     2781855 2       +       156040895       CC
s       Gorilla_gorilla.chrX    12049151        2       +       177553137       CC
s       Pan_paniscus.chrX       2528212 2       +       160249396       CC
s       Pan_troglodytes.chrX    3166702 2       +       153892427       CC
s       Pongo_abelii.chrX       2382605 2       +       162586350       CT
s       Pongo_pygmaeus.chrX     2357111 2       +       161100531       CT
s       Symphalangus_syndactylus.chrX   16744919        2       +       165588543       CC


```

`StitchOrthologs.py` will use the `hal` guide-tree field in the `MAF` header to identify the species in the alignment. For each alignment block in the `MAF`, `StitchOrthologs.py` will classify the block as a *single-copy ortholog segment:* if and only if the block contains a single sequence from each species encoded as a tip in the guide-tree.

`StitchOrthologs.py` yeilds two outputs from the input `MAF`:

1. A single `fasta`-formatted alignment containing a all spliced single-copy orthologs per species
2. A `CSV` containing the local coordinate information for each single-copy ortholog segment in the input `MAF`

## Python environment versioning:

```
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                  2_kmp_llvm    conda-forge
addict                    2.4.0                    pypi_0    pypi
amqp                      5.1.0                    pypi_0    pypi
antlr4-python3-runtime    4.8                      pypi_0    pypi
apache-libcloud           2.8.3                    pypi_0    pypi
argcomplete               2.0.0                    pypi_0    pypi
aria2                     1.36.0               h319415d_0    conda-forge
art                       6.0                      pypi_0    pypi
attrs                     21.4.0             pyhd8ed1ab_0    conda-forge
bagit                     1.8.1                    pypi_0    pypi
bcbio-gff                 0.7.0              pyh7cba7a3_0    bioconda
bcftools                  1.9                  ha228f0b_4    bioconda
billiard                  3.6.4.0                  pypi_0    pypi
binutils_impl_linux-64    2.33.1               he6710b0_7  
binutils_linux-64         2.33.1              h9595d00_15  
biopython                 1.79                     pypi_0    pypi
blas                      1.0                         mkl  
blast                     2.14.1          pl5321h6f7f691_0    bioconda
bleach                    4.1.0                    pypi_0    pypi
blessed                   1.19.1                   pypi_0    pypi
boltons                   21.0.0                   pypi_0    pypi
boto                      2.49.0                   pypi_0    pypi
boto3                     1.21.23                  pypi_0    pypi
boto3-stubs               1.21.23                  pypi_0    pypi
botocore                  1.24.23                  pypi_0    pypi
botocore-stubs            1.24.23                  pypi_0    pypi
bottleneck                1.3.4            py39hce1f21e_0  
brotli                    1.0.9                he6710b0_2  
bx-python                 0.8.13           py39h6471ffd_1    bioconda
bzip2                     1.0.8                h7b6447c_0  
c-ares                    1.18.1               h7f8727e_0  
ca-certificates           2023.7.22            hbcca054_0    conda-forge
cachecontrol              0.12.10                  pypi_0    pypi
cachetools                4.2.4                    pypi_0    pypi
cactus                    1.0                      pypi_0    pypi
celery                    5.2.3                    pypi_0    pypi
certifi                   2023.7.22          pyhd8ed1ab_0    conda-forge
cffi                      1.15.0                   pypi_0    pypi
charset-normalizer        2.0.12                   pypi_0    pypi
cigar                     0.1.3                    pypi_0    pypi
click                     8.0.4                    pypi_0    pypi
click-didyoumean          0.3.0                    pypi_0    pypi
click-plugins             1.1.1                    pypi_0    pypi
click-repl                0.2.0                    pypi_0    pypi
clickclick                20.10.2                  pypi_0    pypi
cmake                     3.18.2               ha30ef3c_0    anaconda
coloredlogs               15.0.1                   pypi_0    pypi
connexion                 2.12.0                   pypi_0    pypi
curl                      7.80.0               h7f8727e_0  
cwltool                   3.1.20211107152837          pypi_0    pypi
cycler                    0.11.0             pyhd3eb1b0_0  
cython                    3.0.2                    pypi_0    pypi
dbus                      1.13.18              hb2f20db_0  
decorator                 5.1.1                    pypi_0    pypi
demes                     0.2.2              pyhd8ed1ab_0    conda-forge
dendropy                  4.6.1              pyhdfd78af_0    bioconda
diamond                   2.1.8                h43eeafb_0    bioconda
dill                      0.3.4                    pypi_0    pypi
docker                    5.0.3                    pypi_0    pypi
docutils                  0.18.1                   pypi_0    pypi
enlighten                 1.10.2                   pypi_0    pypi
entrez-direct             16.2                 he881be0_1    bioconda
expat                     2.4.4                h295c915_0  
fastme                    2.1.6.1              h031d066_3    bioconda
fasttree                  2.1.11               h031d066_2    bioconda
flask                     2.0.3                    pypi_0    pypi
flask-cors                3.0.10                   pypi_0    pypi
fontconfig                2.13.1               h6c09931_0  
fonttools                 4.25.0             pyhd3eb1b0_0  
freetype                  2.11.0               h70c0345_0  
future                    0.18.2                   pypi_0    pypi
futures                   3.1.1                    pypi_0    pypi
galaxy-containers         21.9.0                   pypi_0    pypi
galaxy-tool-util          21.9.2                   pypi_0    pypi
galaxy-util               21.9.0                   pypi_0    pypi
gawk                      5.1.0                h7f98852_0    conda-forge
gcc_impl_linux-64         7.3.0                habb00fd_1  
gcc_linux-64              7.3.0               h553295d_15  
gettext                   0.21.1               h27087fc_0    conda-forge
gffread                   0.12.7               hdcf5f25_3    bioconda
giflib                    5.2.1                h7b6447c_0  
glib                      2.69.1               h4ff587b_1  
gmp                       6.2.1                h58526e2_0    conda-forge
google-api-core           0.1.4                    pypi_0    pypi
google-auth               1.35.0                   pypi_0    pypi
google-cloud-core         0.28.1                   pypi_0    pypi
google-cloud-storage      1.6.0                    pypi_0    pypi
google-crc32c             1.3.0                    pypi_0    pypi
google-resumable-media    2.3.2                    pypi_0    pypi
googleapis-common-protos  1.56.0                   pypi_0    pypi
gsl                       2.7.1                hd82f3ee_0  
gst-plugins-base          1.14.0               h8213a91_2  
gstreamer                 1.14.0               h28cd5cc_2  
gunicorn                  20.1.0                   pypi_0    pypi
gxx_impl_linux-64         7.3.0                hdf63c60_1  
gxx_linux-64              7.3.0               h553295d_15  
h5py                      3.6.0            py39ha0f2276_0  
hdf5                      1.10.6               hb1b8bf9_0  
hmmer                     3.3.2                h1b792b2_1    bioconda
htop                      3.2.2                h8228510_0    conda-forge
http-parser               0.9.0                    pypi_0    pypi
humanfriendly             10.0                     pypi_0    pypi
icu                       58.2                 he6710b0_3  
idna                      3.3                      pypi_0    pypi
importlib-metadata        4.11.4           py39hf3d152e_0    conda-forge
importlib_resources       5.7.1              pyhd8ed1ab_1    conda-forge
inflection                0.5.1                    pypi_0    pypi
iniconfig                 1.1.1                    pypi_0    pypi
intel-openmp              2021.4.0          h06a4308_3561  
iqtree                    2.0.3                h176a8bc_1    bioconda
isodate                   0.6.1                    pypi_0    pypi
itsdangerous              2.1.1                    pypi_0    pypi
jinja2                    3.0.3                    pypi_0    pypi
jmespath                  1.0.0                    pypi_0    pypi
jpeg                      9d                   h7f8727e_0  
jsonschema                4.4.0                    pypi_0    pypi
kazoo                     2.8.0                    pypi_0    pypi
kiwisolver                1.3.2            py39h295c915_0  
kombu                     5.2.4                    pypi_0    pypi
krb5                      1.19.2               hac12032_0  
kubernetes                21.7.0                   pypi_0    pypi
lcms2                     2.12                 h3be6417_0  
ld_impl_linux-64          2.33.1               h53a641e_7  
libcurl                   7.80.0               h0b77cf5_0  
libdeflate                1.0                  h14c3975_1    bioconda
libedit                   3.1.20210910         h7f8727e_0  
libev                     4.33                 h7f8727e_1  
libffi                    3.3                  he6710b0_2  
libgcc-ng                 13.1.0               he5830b7_0    conda-forge
libgfortran-ng            7.5.0               ha8ba4b0_17  
libgfortran4              7.5.0               ha8ba4b0_17  
libidn2                   2.3.4                h166bdaf_0    conda-forge
libnghttp2                1.46.0               hce63b2e_0  
libnl                     3.7.0                h166bdaf_0    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libpng                    1.6.37               hbc83047_0  
libssh2                   1.9.0                h1ba5d50_1  
libstdcxx-ng              13.2.0               h7e041cc_1    conda-forge
libtiff                   4.2.0                h85742a9_0  
libunistring              0.9.10               h7f98852_0    conda-forge
libuuid                   1.0.3                h7f8727e_2  
libuv                     1.40.0               h7b6447c_0  
libwebp                   1.2.2                h55f646e_0  
libwebp-base              1.2.2                h7f8727e_0  
libxcb                    1.14                 h7b6447c_0  
libxml2                   2.9.12               h03d6c58_0  
libzlib                   1.2.13               hd590300_5    conda-forge
llvm-openmp               14.0.6               h9e868ea_0  
lockfile                  0.12.2                   pypi_0    pypi
lxml                      4.8.0                    pypi_0    pypi
lz4-c                     1.9.3                h295c915_1  
lzo                       2.10              h516909a_1000    conda-forge
macse                     2.07                 hdfd78af_0    bioconda
mafft                     7.310                h1b792b2_4    bioconda
markupsafe                2.1.1                    pypi_0    pypi
matplotlib                3.5.1            py39h06a4308_1  
matplotlib-base           3.5.1            py39ha18d171_1  
mcl                       22.282          pl5321h031d066_0    bioconda
mistune                   0.8.4                    pypi_0    pypi
mkl                       2021.4.0           h06a4308_640  
mkl-service               2.4.0            py39h7f8727e_0  
mkl_fft                   1.3.1            py39hd3c417c_0  
mkl_random                1.2.2            py39h51133e4_0  
mmseqs2                   14.7e284        pl5321h6a68c12_2    bioconda
mpi                       1.0                     openmpi    conda-forge
msgpack                   1.0.3                    pypi_0    pypi
msprime                   1.2.0            py39h78aa3b4_0    conda-forge
munkres                   1.1.4                      py_0  
muscle                    5.1                  h4ac6f70_3    bioconda
mypy-boto3-s3             1.21.7                   pypi_0    pypi
mypy-extensions           0.4.3                    pypi_0    pypi
ncbi-datasets-cli         15.19.1              ha770c72_0    conda-forge
ncbi-vdb                  3.0.8                hdbdd923_0    bioconda
ncurses                   6.3                  h7f8727e_2  
networkx                  2.7.1                    pypi_0    pypi
numexpr                   2.8.1            py39h6abb31d_0  
numpy                     1.21.2           py39h20f2e39_0  
numpy-base                1.21.2           py39h79a1101_0  
oauthlib                  3.2.0                    pypi_0    pypi
openjdk                   8.0.382              hd590300_0    conda-forge
openmpi                   4.1.5                external_1    conda-forge
openssl                   1.1.1w               hd590300_0    conda-forge
orthofinder               2.5.5                hdfd78af_1    bioconda
ossuuid                   1.6.2             hf484d3e_1000    conda-forge
packaging                 21.3               pyhd3eb1b0_0  
pandas                    1.4.1            py39h295c915_1  
pcre                      8.45                 h295c915_0  
perl                      5.32.1          4_hd590300_perl5    conda-forge
perl-alien-build          2.48            pl5321hec16e2b_0    bioconda
perl-alien-libxml2        0.17            pl5321hec16e2b_0    bioconda
perl-archive-tar          2.40            pl5321hdfd78af_0    bioconda
perl-business-isbn        3.007           pl5321hd8ed1ab_0    conda-forge
perl-business-isbn-data   20210112.006    pl5321hd8ed1ab_0    conda-forge
perl-capture-tiny         0.48            pl5321ha770c72_1    conda-forge
perl-carp                 1.50            pl5321hd8ed1ab_0    conda-forge
perl-common-sense         3.75            pl5321hd8ed1ab_0    conda-forge
perl-compress-raw-bzip2   2.201           pl5321h166bdaf_0    conda-forge
perl-compress-raw-zlib    2.202           pl5321h166bdaf_0    conda-forge
perl-constant             1.33            pl5321hd8ed1ab_0    conda-forge
perl-encode               3.19            pl5321h166bdaf_0    conda-forge
perl-exporter             5.74            pl5321hd8ed1ab_0    conda-forge
perl-exporter-tiny        1.002002        pl5321hd8ed1ab_0    conda-forge
perl-extutils-makemaker   7.70            pl5321hd8ed1ab_0    conda-forge
perl-ffi-checklib         0.28            pl5321hdfd78af_0    bioconda
perl-file-chdir           0.1011          pl5321hd8ed1ab_0    conda-forge
perl-file-path            2.18            pl5321hd8ed1ab_0    conda-forge
perl-file-temp            0.2304          pl5321hd8ed1ab_0    conda-forge
perl-file-which           1.24            pl5321hd8ed1ab_0    conda-forge
perl-importer             0.026           pl5321hd8ed1ab_0    conda-forge
perl-io-compress          2.201           pl5321hdbdd923_2    bioconda
perl-io-zlib              1.14            pl5321hdfd78af_0    bioconda
perl-json                 4.10            pl5321hdfd78af_0    bioconda
perl-json-xs              2.34            pl5321h4ac6f70_6    bioconda
perl-list-moreutils       0.430           pl5321hdfd78af_0    bioconda
perl-list-moreutils-xs    0.430           pl5321h031d066_2    bioconda
perl-parent               0.241           pl5321hd8ed1ab_0    conda-forge
perl-path-tiny            0.124           pl5321hd8ed1ab_0    conda-forge
perl-pathtools            3.75            pl5321h166bdaf_0    conda-forge
perl-scalar-list-utils    1.63            pl5321h166bdaf_0    conda-forge
perl-scope-guard          0.21            pl5321hd8ed1ab_0    conda-forge
perl-storable             3.15            pl5321h166bdaf_0    conda-forge
perl-sub-info             0.002           pl5321hd8ed1ab_0    conda-forge
perl-term-table           0.016           pl5321hdfd78af_0    bioconda
perl-test-fatal           0.016           pl5321ha770c72_0    conda-forge
perl-test-warnings        0.031           pl5321ha770c72_0    conda-forge
perl-test2-suite          0.000145        pl5321hdfd78af_0    bioconda
perl-try-tiny             0.31            pl5321ha770c72_0    conda-forge
perl-types-serialiser     1.01            pl5321hdfd78af_0    bioconda
perl-uri                  5.17            pl5321ha770c72_0    conda-forge
perl-xml-libxml           2.0207          pl5321h661654b_0    bioconda
perl-xml-namespacesupport 1.12            pl5321hd8ed1ab_0    conda-forge
perl-xml-sax              1.02            pl5321hd8ed1ab_0    conda-forge
perl-xml-sax-base         1.09            pl5321hd8ed1ab_0    conda-forge
pillow                    9.0.1            py39h22f2fdc_0  
pip                       21.3.1                   pypi_0    pypi
plink                     1.90b6.21            h031d066_5    bioconda
pluggy                    1.0.0                    pypi_0    pypi
prefixed                  0.3.2                    pypi_0    pypi
prompt-toolkit            3.0.28                   pypi_0    pypi
protobuf                  3.19.4                   pypi_0    pypi
prov                      1.5.1                    pypi_0    pypi
psutil                    5.9.0                    pypi_0    pypi
py                        1.11.0                   pypi_0    pypi
py-tes                    0.4.2                    pypi_0    pypi
pyasn1                    0.4.8                    pypi_0    pypi
pyasn1-modules            0.2.8                    pypi_0    pypi
pycparser                 2.21                     pypi_0    pypi
pycryptodome              3.14.1                   pypi_0    pypi
pydantic                  1.9.0                    pypi_0    pypi
pydot                     1.4.2                    pypi_0    pypi
pymesos                   0.3.15                   pypi_0    pypi
pynacl                    1.5.0                    pypi_0    pypi
pyparsing                 3.0.4              pyhd3eb1b0_0  
pypubsub                  4.0.3                    pypi_0    pypi
pyqt                      5.9.2            py39h2531618_6  
pyrsistent                0.18.1           py39hb9d737c_1    conda-forge
pysam                     0.21.0                   pypi_0    pypi
pytest                    7.1.1                    pypi_0    pypi
python                    3.9.11               h12debd9_1  
python-dateutil           2.8.2              pyhd3eb1b0_0  
python-lzo                1.15             py39heaa0706_0    conda-forge
python_abi                3.9                      2_cp39    conda-forge
pytz                      2022.1                   pypi_0    pypi
pyyaml                    5.4.1                    pypi_0    pypi
qt                        5.9.7                h5867ecd_1  
raxml                     8.2.12               h031d066_6    bioconda
raxml-ng                  1.2.0                h6d1f11b_1    bioconda
rdflib                    6.0.2                    pypi_0    pypi
readline                  8.1.2                h7f8727e_1  
repoze-lru                0.7                      pypi_0    pypi
requests                  2.27.1                   pypi_0    pypi
requests-oauthlib         1.3.1                    pypi_0    pypi
rhash                     1.4.1                h3c74f83_1  
routes                    2.5.1                    pypi_0    pypi
rsa                       4.8                      pypi_0    pypi
ruamel-yaml               0.17.16                  pypi_0    pypi
ruamel.yaml               0.17.21          py39hb9d737c_1    conda-forge
ruamel.yaml.clib          0.2.6            py39hb9d737c_1    conda-forge
s3transfer                0.5.2                    pypi_0    pypi
schema-salad              8.2.20220204150214          pypi_0    pypi
scipy                     1.7.3            py39hc147768_0  
seqkit                    2.5.0                h9ee0642_0    bioconda
setuptools                59.6.0                   pypi_0    pypi
shellescape               3.8.1                    pypi_0    pypi
sip                       4.19.13          py39h295c915_0  
six                       1.16.0             pyhd3eb1b0_1  
snp-sites                 2.5.1                h5bf99c6_1    bioconda
sonlib                    2.0                      pypi_0    pypi
sortedcontainers          2.4.0                    pypi_0    pypi
sqlite                    3.38.0               hc218d9a_0  
svgwrite                  1.4.2              pyhd8ed1ab_0    conda-forge
swagger-ui-bundle         0.0.9                    pypi_0    pypi
tk                        8.6.11               h1ccaba5_0  
toil                      5.11.0a1                 pypi_0    pypi
tomli                     2.0.1                    pypi_0    pypi
tornado                   6.1              py39h27cfd23_0  
tskit                     0.4.1            py39hce5d2b2_0    conda-forge
typing-extensions         4.1.1                    pypi_0    pypi
tzdata                    2021e                hda174b7_0  
urllib3                   1.26.9                   pypi_0    pypi
vcftools                  0.1.16               h9a82719_5    bioconda
vine                      5.0.0                    pypi_0    pypi
wcwidth                   0.2.5                    pypi_0    pypi
wdlparse                  0.1.0                    pypi_0    pypi
webencodings              0.5.1                    pypi_0    pypi
websocket-client          1.3.1                    pypi_0    pypi
werkzeug                  2.0.3                    pypi_0    pypi
wget                      1.20.3               ha56f1ee_1    conda-forge
wheel                     0.37.1             pyhd3eb1b0_0  
xz                        5.2.5                h7b6447c_0  
zipp                      3.8.0              pyhd8ed1ab_0    conda-forge
zipstream-new             1.1.8                    pypi_0    pypi
zlib                      1.2.13               hd590300_5    conda-forge
zstd                      1.4.5                h9ceee32_0  
```


