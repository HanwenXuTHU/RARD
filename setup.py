"""
Created on Wen Apr 8 10:18 2020

@author: Hanwen Xu

E-mail: hw-xu16@mails.tsinghua.edu.cn

New update: process data with unknown content
"""
from setuptools import setup, find_packages

setup(
    name="DeconCCN",
    version="0.1.5",
    keywords=("pip", "DeconCCN"),
    description="A cell subtype deconvolution algorithm based on component-wise condition number",
    long_description="a cell subtype deconvolution algorithm based on component-wise condition number",
    license="Hanwen",

    url="https://github.com/HanwenXuTHU/DeconCCN/",
    author="Hanwen Xu",
    author_email="xuhw20@mails.tsinghua.edu.cn",

    packages=find_packages(),
    include_package_data=True,
    platforms="any",
    install_requires=['numpy', 'tqdm', 'scipy', 'sklearn', 'statsmodels', 'pandas']
)