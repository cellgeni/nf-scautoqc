FROM python:3.11

RUN pip install --no-cache matplotlib-inline matplotlib-venn python-dateutil pyzmq celltypist leidenalg llvmlite loompy numba numpy \
    numpy-groupies numpyro pytorch-lightning scikit-image scipy scrublet scvi-tools tensorstore \
    torch torchaudio torchmetrics torchvision tqdm pandas scikit-learn anndata matplotlib sympy umap-learn \
    statsmodels seaborn docopt

RUN pip install --no-cache scanpy==1.9.3

RUN pip install --no-cache annoy==1.15.1

RUN pip install --no-cache torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

RUN pip install --no-cache scikit-misc pymde rapids-singlecell

RUN pip install --no-cache git+https://github.com/Teichlab/sctk

RUN pip install --no-cache scikit-learn==1.5.0 numpy==2.2.6 --force-reinstall