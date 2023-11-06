FROM python:3.10

RUN pip install --no-cache matplotlib-inline python-dateutil pyzmq celltypist leidenalg llvmlite loompy numba numpy \
    numpy-groupies numpyro pytorch-lightning scikit-image scipy scrublet sctk scvi-tools tensorstore \
    torch torchaudio torchmetrics torchvision tqdm pandas scikit-learn anndata matplotlib sympy umap-learn \
    statsmodels seaborn

RUN pip install --no-cache scanpy==1.9.3

RUN pip install --no-cache docopt

RUN pip install --no-cache annoy==1.15.1

RUN pip install --no-cache torch==1.13.1+cu116 torchvision==0.14.1+cu116 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu116

RUN pip install --no-cache scikit-misc

RUN pip install --no-cache pymde
