#!/usr/bin/env python

from zipfile import ZipFile
from ipywidgets import interactive, Layout, FloatSlider
import ipywidgets as widgets
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


with ZipFile('./PBE0_impurities_analysis/ALL_SR_MOLOPT/ADMM_UZH_HF27.zip', 'r') as UZH_HF27zip:
    UZH_HF27zip.extractall(path='./PBE0_impurities_analysis/ALL_SR_MOLOPT')


UZH = FloatSlider(min=-0.09, max=0.09, step=0.005, description="UZH",
                          layout=Layout(margin='-5px 0px 0px 225px'))

ui = widgets.HBox([UZH], layout=Layout(margin='0px 5000px 0px 80px'))


def wfn_plot(UZH):

    img_fileUZH = str(
        "./PBE0_impurities_analysis/All_SR_MOLOPT/ADMM_UZH_HF27/Al_Si/hex/neutral/wfn_pictures/wfn{}.tga".format(UZH))

    fig, ax = plt.subplots(figsize =(12,8))
    imgUZH = mpimg.imread(img_fileUZH)
    ax.imshow(imgUZH)
    ax.axis('off')
    plt.show()


out = widgets.interactive_output(wfn_plot, {'UZH': UZH})

display(out, ui)



# !rm -r. /PBE0_impurities_analysis/ALL_SR_MOLOPT/ADMM_UZH_HF27