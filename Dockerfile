FROM python:3.11

RUN pip install --no-cache matplotlib-inline python-dateutil pyzmq celltypist leidenalg llvmlite loompy numba numpy \
    numpy-groupies numpyro pytorch-lightning scikit-image scipy scrublet sctk scvi-tools tensorstore \
    torch torchaudio torchmetrics torchvision tqdm pandas scikit-learn anndata matplotlib sympy umap-learn \
    statsmodels seaborn

RUN pip install --no-cache scanpy==1.9.3

RUN pip install --no-cache docopt

RUN pip install --no-cache annoy==1.15.1

RUN pip install --no-cache torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

RUN pip install --no-cache scikit-misc

RUN pip install --no-cache pymde

RUN pip install --no-cache rapids-singlecell