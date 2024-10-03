#!/usr/bin/env python

# May 2022 - updated Sept 2024

import numpy as np
import tkinter as tk
from tkinter import *
from tkinter.ttk import *
from tkinter import ttk
from ipywidgets import interactive, FloatSlider
import ipywidgets as widgets
from PIL import ImageTk, Image
import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from mpl_toolkits import mplot3d
import ase.io  # Atomic structure tools library (Atomic Simulation Environment)
from ase.io.cube import read_cube_data
import csv


## Property classes
class DropdownMenuStatus:
    def __init__(self):
        self._blank_shown, self._A_notebook_shown, self._B_notebook_shown = True, False, False
        self._label_shown, self._previous_tab = False, ""

    @property
    def blank_shown(self):
        return self._blank_shown

    @blank_shown.setter
    def blank_shown(self, bool):
        self._blank_shown = bool

    @property
    def A_notebook(self):
        return self._A_notebook_shown

    @A_notebook.setter
    def A_notebook(self, bool):
        self._A_notebook_shown = bool

    @property
    def B_notebook(self):
        return self._B_notebook_shown

    @B_notebook.setter
    def B_notebook(self, bool):
        self._B_notebook_shown = bool

    @property
    def label(self):
        return self._label_shown

    @label.setter
    def label(self, bool):
        self._label_shown = bool

    @property
    def previous_tab_type(self):
        return self._previous_tab

    @previous_tab_type.setter
    def previous_tab_type(self, str):
        self._previous_tab = str


class Frames:
    def __init__(self):
        self._blank_pic_frame, self._A_pic_frame, self._B_pic_frame = None, None, None
        self._blank_pic_label, self._A_pic_label, self._B_pic_label = None, None, None
        self._aA, self._aB, self._bA, self._bB, self._cA, self._cB = None, None, None, None, None, None
        self._dA, self._dB, self._eA, self._eB = None, None, None, None

        self._left_aA, self._right_aA, self._left_bA, self._right_bA = None, None, None, None
        self._left_cA, self._right_cA = None, None

    @property
    def blank_pic_frame(self):
        return self._blank_pic_frame

    @blank_pic_frame.setter
    def blank_pic_frame(self, tk_frame):
        self._blank_pic_frame = tk_frame

    @property
    def A_pic_frame(self):
        return self._A_pic_frame

    @A_pic_frame.setter
    def A_pic_frame(self, tk_frame):
        self._A_pic_frame = tk_frame

    @property
    def B_pic_frame(self):
        return self._B_pic_frame

    @B_pic_frame.setter
    def B_pic_frame(self, tk_frame):
        self._B_pic_frame = tk_frame

    @property
    def blank_pic_label(self):
        return self._blank_pic_label

    @blank_pic_label.setter
    def blank_pic_label(self, tk_label):
        self._blank_pic_label = tk_label

    @property
    def A_pic_label(self):
        return self._A_pic_label

    @A_pic_label.setter
    def A_pic_label(self, tk_label):
        self._A_pic_label = tk_label

    @property
    def B_pic_label(self):
        return self._B_pic_label

    @B_pic_label.setter
    def B_pic_label(self, tk_label):
        self._B_pic_label = tk_label

    @property
    def frame_aA(self):
        return self._aA

    @frame_aA.setter
    def frame_aA(self, tk_frame):
        self._aA = tk_frame

    @property
    def frame_aB(self):
        return self._aB

    @frame_aB.setter
    def frame_aB(self, tk_frame):
        self._aB = tk_frame

    @property
    def frame_bA(self):
        return self._bA

    @frame_bA.setter
    def frame_bA(self, tk_frame):
        self._bA = tk_frame

    @property
    def frame_bB(self):
        return self._bB

    @frame_bB.setter
    def frame_bB(self, tk_frame):
        self._bB = tk_frame

    @property
    def frame_cA(self):
        return self._cA

    @frame_cA.setter
    def frame_cA(self, tk_frame):
        self._cA = tk_frame

    @property
    def frame_cB(self):
        return self._cB

    @frame_cB.setter
    def frame_cB(self, tk_frame):
        self._cB = tk_frame

    @property
    def frame_dA(self):
        return self._dA

    @frame_dA.setter
    def frame_dA(self, tk_frame):
        self._dA = tk_frame

    @property
    def frame_dB(self):
        return self._dB

    @frame_dB.setter
    def frame_dB(self, tk_frame):
        self._dB = tk_frame

    @property
    def frame_eA(self):
        return self._eA

    @frame_eA.setter
    def frame_eA(self, tk_frame):
        self._eA = tk_frame

    @property
    def frame_eB(self):
        return self._eB

    @frame_eB.setter
    def frame_eB(self, tk_frame):
        self._eB = tk_frame


class HandlingPictures:
    def __init__(self):
        self._blank, self._A, self._B = None, None, None
        self._blank_deleted, self._A_deleted, self._B_deleted = False, False, False

    @property
    def blank(self):
        return self._blank

    @blank.setter
    def blank(self, pyimage):
        self._blank = pyimage

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, pyimage):
        self._A = pyimage

    @property
    def B(self):
        return self._B

    @B.setter
    def B(self, pyimage):
        self._B = pyimage

    @property
    def blank_deleted(self):
        return self._blank_deleted

    @blank_deleted.setter
    def blank_deleted(self, bool):
        self._blank_deleted = bool

    @property
    def A_deleted(self):
        return self._A_deleted

    @A_deleted.setter
    def A_deleted(self, bool):
        self._A_deleted = bool

    @property
    def B_deleted(self):
        return self._B_deleted

    @B_deleted.setter
    def B_deleted(self, bool):
        self._B_deleted = bool


class DoublePageWidgets:
    def __init__(self):
        self._left_frame, self._left_canvas, self._right_frame, self._right_canvas = None, None, None, None
        self._left_txt, self._left_table_frame, self._left_table, self._right_txt = None, None, None, None
        self._right_table_frame, self._right_table = None, None
        self._left_p_fig, self._left_d_fig, self._left_extended = None, None, None
        self._right_p_fig, self._right_d_fig, self._right_extended = None, None, None
        self._left_L, self._left_H, self._right_L, self._right_H = None, None, None, None
        self._left_LUMO, self._left_HOMO, self._right_LUMO, self._right_HOMO = None, None, None, None
        self._left_slider, self._right_slider = None, None
        self._cur_l_canvs, self._cur_r_canvs = [], []

    @property
    def left_frame(self):
        return self._left_frame

    @left_frame.setter
    def left_frame(self, tk_frame):
        self._left_frame = tk_frame

    @property
    def left_canvas(self):
        return self._left_canvas

    @left_canvas.setter
    def left_canvas(self, tk_canvas):
        self._left_canvas = tk_canvas

    @property
    def right_frame(self):
        return self._right_frame

    @right_frame.setter
    def right_frame(self, tk_frame):
        self._right_frame = tk_frame

    @property
    def right_canvas(self):
        return self._right_canvas

    @right_canvas.setter
    def right_canvas(self, tk_canvas):
        self._right_canvas = tk_canvas

    @property
    def left_text(self):
        return self._left_txt

    @left_text.setter
    def left_text(self, tk_frame):
        self._left_txt = tk_frame

    @property
    def left_table_frame(self):
        return self._left_table_frame

    @left_table_frame.setter
    def left_table_frame(self, tk_frame):
        self._left_table_frame = tk_frame

    @property
    def left_table(self):
        return self._left_table

    @left_table.setter
    def left_table(self, t_):
        self._left_table = t_

    @property
    def right_text(self):
        return self._right_txt

    @right_text.setter
    def right_text(self, tk_frame):
        self._right_txt = tk_frame

    @property
    def right_table_frame(self):
        return self._right_table_frame

    @right_table_frame.setter
    def right_table_frame(self, tk_frame):
        self._right_table_frame = tk_frame

    @property
    def right_table(self):
        return self._right_table

    @right_table.setter
    def right_table(self, t_):
        self._right_table = t_

    @property
    def left_p_fig(self):
        return self._left_p_fig

    @left_p_fig.setter
    def left_p_fig(self, tk_frame):
        self._left_p_fig = tk_frame

    @property
    def left_d_fig(self):
        return self._left_d_fig

    @left_d_fig.setter
    def left_d_fig(self, tk_frame):
        self._left_d_fig = tk_frame

    @property
    def left_extended(self):
        return self._left_extended

    @left_extended.setter
    def left_extended(self, tk_frame):
        self._left_extended = tk_frame

    @property
    def right_p_fig(self):
        return self._right_p_fig

    @right_p_fig.setter
    def right_p_fig(self, tk_frame):
        self._right_p_fig = tk_frame

    @property
    def right_d_fig(self):
        return self._right_d_fig

    @right_d_fig.setter
    def right_d_fig(self, tk_frame):
        self._right_d_fig = tk_frame

    @property
    def right_extended(self):
        return self._right_extended

    @right_extended.setter
    def right_extended(self, tk_frame):
        self._right_extended = tk_frame

    @property
    def left_L(self):
        return self._left_L

    @left_L.setter
    def left_L(self, tk_frame):
        self._left_L = tk_frame

    @property
    def left_H(self):
        return self._left_H

    @left_H.setter
    def left_H(self, tk_frame):
        self._left_H = tk_frame

    @property
    def right_L(self):
        return self._right_L

    @right_L.setter
    def right_L(self, tk_frame):
        self._right_L = tk_frame

    @property
    def right_H(self):
        return self._right_H

    @right_H.setter
    def right_H(self, tk_frame):
        self._right_H = tk_frame

    @property
    def left_LUMO(self):
        return self._left_LUMO

    @left_LUMO.setter
    def left_LUMO(self, tk_frame):
        self._left_LUMO = tk_frame

    @property
    def left_HOMO(self):
        return self._left_HOMO

    @left_HOMO.setter
    def left_HOMO(self, tk_frame):
        self._left_HOMO = tk_frame

    @property
    def right_LUMO(self):
        return self._right_LUMO

    @right_LUMO.setter
    def right_LUMO(self, tk_frame):
        self._right_LUMO = tk_frame

    @property
    def right_HOMO(self):
        return self._right_HOMO

    @right_HOMO.setter
    def right_HOMO(self, tk_frame):
        self._right_HOMO = tk_frame

    @property
    def left_slider(self):
        return self._left_slider

    @left_slider.setter
    def left_slider(self, tk_scale):
        self._left_slider = tk_scale

    @property
    def right_slider(self):
        return self._right_slider

    @right_slider.setter
    def right_slider(self, tk_scale):
        self._right_slider = tk_scale

    @property
    def cur_l_canvas(self):
        return self._cur_l_canvs

    @property
    def cur_r_canvas(self):
        return self._cur_r_canvs


class DoublePageNotebooks:
    def __init__(self):
        self._dict = {}

    @property
    def dictionary(self):
        return self._dict

    @dictionary.setter
    def dictionary(self, dict):
        self._dict = dict


## classes
class Notebook_Handler:
    def __init__(self, notebook_name, letter):
        self.notebook, self.letter = notebook_name, letter

    def add_frames(self):
        self.notebook.grid(column=0, row=19, columnspan=90, rowspan=34)
        for let in 'a', 'b', 'c', 'd', 'e':
            frame = eval("Frames().frame_{}{}".format(let, self.letter))
            self.notebook.add(frame, text=let)
        if self.letter == 'A':
            # populate double notebook pages with widgets.
            for i in 'a', 'b', 'c':
                Double_Notebook_page(i)


class NotebookTab:
    def __init__(self, notebook_name):
        self.notebook = notebook_name
        self.tab_num = notebook_name.index('current')

    def get_tab_name(self):
        for indx, name in enumerate(['a', 'b', 'c', 'd', 'e']):
            if self.tab_num == indx:
                tab_name = name
        return tab_name


class SliderPictureUpdater:
    def __init__(self, fr, slider):
        self.notebook_name = A_notebook if str(dropdown.get()).find('A') != -1 else B_notebook

        self.tab, self.isovalue = NotebookTab(self.notebook_name).get_tab_name(), 0

        self.dct = DoublePageNotebooks().dictionary[self.tab]
        self.slider = self.dct[str("_{}_slider".format(slider))]

        self.H_fr, self.L_fr = None, None

        self.H_fr = [self.dct[str(f)] for f in fr if f.find('H') != -1][0] if type(fr) == list else \
            [self.dct[str(fr)] if fr.find('H') != -1 else self.H_fr][0]

        self.L_fr = [self.dct[str(f)] for f in fr if f.find('L') != -1][0] if type(fr) == list else \
            [self.dct[str(fr)] if fr.find('L') != -1 else self.L_fr][0]

    def set_slider(self, value=None):
        if value:
            self.slider.set(value)
        else:
            self.slider.set(0.01)

    def get_isovalue(self):
        self.isovalue = float(self.slider.get())

    def forget_canvas(self, fr):
        if self.dct[str("_cur_{}_canvs".format(fr))] != []:
            for canvas in self.dct[str("_cur_{}_canvs".format(fr))]:
                canvas.get_tk_widget().grid_remove()

            self.dct[str("_cur_{}_canvs".format(fr))] = []

    def plotter(self, file):
        figure = plt.Figure(figsize=(4.51, 2.78), dpi=100)
        ax = figure.add_subplot(111)
        img_wfn = mpimg.imread(file)
        ax.imshow(img_wfn, interpolation='nearest')  # , aspect='equal'
        ax.axis('off')
        figure.tight_layout()

        return figure

    def canvas(self, figure, frame):
        canvas = FigureCanvasTkAgg(figure, master=frame)
        canvas.get_tk_widget().grid(sticky="nsew")

        return canvas


class old_isovalues(SliderPictureUpdater):
    def __init__(self, event):
        super().__init__(["_left_H", "_left_L"], "left")
        super().get_isovalue()
        super().forget_canvas('l')

        H_file = str(f"./A_notebook/{self.tab}/old/Hold{self.isovalue}.jpg")
        L_file = str(f"./A_notebook/{self.tab}/old/Lold{self.isovalue}.jpg")
        for W in 'H', 'L':
            fig = super().plotter(eval("{}_file".format(W)))
            exec(f'{W}_fig_old{self.tab}A = fig')
            canvas = super().canvas(eval("{}_fig_old{}A".format(W, self.tab)), eval("self.{}_fr".format(W)))
            self.dct['_cur_l_canvs'].append(canvas)
        self.set_slider(self.isovalue)


class new_isovalues(SliderPictureUpdater):
    def __init__(self, event):
        super().__init__(["_right_H", "_right_L"], "right")
        super().get_isovalue()
        super().forget_canvas('r')

        H_file = str(f"./A_notebook/{self.tab}/new/Hnew{self.isovalue}.jpg")
        L_file = str(f"./A_notebook/{self.tab}/new/Lnew{self.isovalue}.jpg")
        for W in 'H', 'L':
            name = str(f"{W}_{self.isovalue}")

            fig = super().plotter(eval("{}_file".format(W)))
            exec(f'{W}_fig_new{self.tab}A = fig')
            canvas = super().canvas(eval("{}_fig_new{}A".format(W, self.tab)), eval("self.{}_fr".format(W)))
            self.dct['_cur_r_canvs'].append(canvas)


class Tables:
    def __init__(self, frame, wdth, lst):
        total_rows = len(lst)
        total_columns = len(lst[0])
        for i in range(total_rows):
            for j in range(total_columns):
                self.e = tk.Entry(frame, width=wdth[j], fg='white', font=('Raleway', 11))
                self.e.grid(row=i, column=j)
                self.e.insert(END, lst[i][j])


class Double_Notebook_page:
    def __init__(self, lett):
        self.lett = lett
        self.props = DoublePageWidgets()
        self.setup_right_left_frames()
        self.make_widget_frames()
        self.position_widgets()
        self.populate_table()
        self.populate_extended()

        if DoublePageNotebooks().dictionary == {}:
            DoublePageNotebooks.dictionary = {self.lett: self.props.__dict__}
        else:
            DoublePageNotebooks().dictionary[self.lett] = self.props.__dict__

        self.setup_first_pictures()

    def setup_right_left_frames(self):
        for side, col in zip(['left', 'right'], [3, 52]):
            exec(f'self.props.{side}_frame = ttk.Frame(Frames().frame_{self.lett}A, width=675, height=540)')
            exec(f'self.props.{side}_frame.grid(column={col}, columnspan=45, row=1, rowspan=57)')
            exec(f'self.props.{side}_canvas = tk.Canvas(self.props.{side}_frame, width=675, height=540, bg="#303030")')
            exec(f'self.props.{side}_canvas.grid(columnspan=43, rowspan=57)')

    def make_widget_frames(self):
        for side, txt in zip(["left", "right"], ["older version 1", "newer version 2"]):
            exec(f'self.props.{side}_text =  tk.Label(self.props.{side}_frame, text=txt)')
            exec(f'self.props.{side}_table_frame = ttk.Frame(self.props.{side}_frame)')
            exec(f'self.props.{side}_p_fig = ttk.Frame(self.props.{side}_frame, width=270, height=240)')
            exec(f'self.props.{side}_d_fig = ttk.Frame(self.props.{side}_frame, width=270, height=240)')
            exec(f'self.props.{side}_extended = ttk.Frame(self.props.{side}_frame, width=675, height=230)')

    def position_widgets(self):
        for side in 'left', 'right':
            exec(f'self.props.{side}_text.grid(column=2, row=0, rowspan=2)')
            exec(f'self.props.{side}_table_frame.grid(column=8, row=0, rowspan=3, columnspan=29)')
            exec(f'self.props.{side}_p_fig.grid(column=2, row=4, columnspan=14, rowspan=16)')
            exec(f'self.props.{side}_d_fig.grid(column=19, row=4, columnspan=20, rowspan=16)')
            exec(f'self.props.{side}_extended.grid(column = 0, row = 38, columnspan=45, rowspan = 22)')

    def populate_table(self):
        for side, age in zip(['left', 'right'], ['old', 'new']):
            lst = [('', "Energy(eV)", "HOMO(eV)", "LUMO(eV)", "Bandgap(eV)"),
                   (f"{self.lett}A defect {age}", " ", " ", " ", " ")]
            WDTs = [12, 10, 8, 8, 12]
            exec(f'self.props.{side}_table = Tables(self.props.{side}_table_frame, WDTs, lst)')

    def populate_extended(self):
        for side, age in zip(['left', 'right'], ['old', 'new']):
            exec(
                f'self.props.{side}_slider = tk.Scale(self.props.{side}_frame, from_=0.01, to=0.09, digits=2, resolution=0.005, orient=tk.VERTICAL, length=200, command={age}_isovalues)')
            exec(f'self.props.{side}_slider.grid(column=35, row=40, columnspan=9, rowspan=20)')
            for let, word, col1, col2 in zip(['H', 'L'], ['HOMO', 'LUMO'], [0, 15], [1, 16]):
                exec(f'self.props.{side}_{let} = ttk.Frame(self.props.{side}_extended, width=325, height=200)')
                exec(f'self.props.{side}_{let}.grid(column = col1, row = 3, columnspan=12, rowspan = 17)')

                exec(f'self.props.{side}_{word} = tk.Label(self.props.{side}_extended, text=word)')
                exec(f'self.props.{side}_{word}.grid(column=col2, columnspan = 2, row = 0, rowspan = 2)')

    def setup_first_pictures(self):
        for side, let, age in zip(['left', 'right'], ['l', 'r'], ['old', 'new']):
            setup = SliderPictureUpdater([f"_{side}_H", f"_{side}_L"], side)
            H_file = str(f"./A_notebook/{self.lett}/{age}/H{age}0.01.jpg")
            L_file = str(f"./A_notebook/{self.lett}/{age}/L{age}0.01.jpg")
            for W in 'H', 'L':
                fig = setup.plotter(eval("{}_file".format(W)))
                exec(f'{W}_fig_old{self.lett}A = fig')
                canvas = setup.canvas(eval("{}_fig_old{}A".format(W, self.lett)),
                                      eval("self.props.{}_{}".format(side, W)))
                DoublePageNotebooks().dictionary[self.lett][f"_cur_{let}_canvs"].append(canvas)
            setup.set_slider()


## functions
def setheaderpictures():
    blank_path, A_path, B_path = str("./blank.jpeg"), str("./A.jpeg"), str("./B.jpeg")
    blank_img, A_img, B_img = Image.open(blank_path), Image.open(A_path), Image.open(B_path)
    blank_resized = blank_img.resize((190, 205), Image.LANCZOS)
    A_resized = A_img.resize((190, 205), Image.LANCZOS)
    B_resized = B_img.resize((190, 205), Image.LANCZOS)
    HandlingPictures.blank = ImageTk.PhotoImage(blank_resized)
    HandlingPictures.A = ImageTk.PhotoImage(A_resized)
    HandlingPictures.B = ImageTk.PhotoImage(B_resized)


def heading_pictures(current, replace=None, new=None):
    if current == 'start':
        frame, label = Frames().blank_pic_frame, Frames().blank_pic_label
        Frames().blank_pic_frame.grid(column=28, row=0, columnspan=10, rowspan=8, sticky="nsew")
        Frames().blank_pic_label.pack(expand=1, fill=BOTH)
        HandlingPictures.blank_deleted = None

    elif replace:
        old_frame, new_frame = eval("Frames().{}_pic_frame".format(current)), eval("Frames().{}_pic_frame".format(new))
        old_frame.grid_remove()
        exec(f'HandlingPictures.{current}_deleted = True')
        if eval("HandlingPictures().{}_deleted".format(new)) is True:
            new_frame.grid(sticky="nsew")
            exec(f'HandlingPictures.{new}_deleted = None')
        elif eval("HandlingPictures().{}_deleted".format(new)) is False:
            frame, label = eval("Frames().{}_pic_frame".format(new)), eval("Frames().{}_pic_label".format(new))
            frame.grid(column=28, row=0, columnspan=10, rowspan=8, sticky="nsew")
            label.pack(expand=1, fill=BOTH)
            exec(f'HandlingPictures.{new}_deleted = None')


def dropdown_menu_command(event):
    opt = ['A' if str(dropdown.get()).find('A') != -1 else 'B'][0] if str(dropdown.get()) != '--' else 'blank'
    if DropdownMenuStatus().blank_shown is True:
        DropdownMenuStatus.blank_shown = False
        heading_pictures('blank', 'yes', opt)
        to_start.grid_remove()

    elif DropdownMenuStatus().A_notebook is True:
        DropdownMenuStatus.A_notebook = False
        heading_pictures('A', 'yes', opt)
        A_notebook.grid_remove()

    elif DropdownMenuStatus().label is True:
        heading_pictures(DropdownMenuStatus().previous_tab_type, 'yes', opt)
        DropdownMenuStatus.label = False
        C_label.grid_remove()

    elif DropdownMenuStatus().B_notebook is True:
        DropdownMenuStatus.B_notebook = False
        heading_pictures('B', 'yes', opt)
        B_notebook.grid_remove()

    if str(dropdown.get()).find('2') != -1:
        C_label.grid(column=15, columnspan=10, row=7, rowspan=2)
        DropdownMenuStatus.label = True
        exec(f'DropdownMenuStatus.{opt}_notebook = True')
    else:
        if DropdownMenuStatus().label is True:
            DropdownMenuStatus.label = False
            C_label.grid_remove()
        if str(dropdown.get()) == '--':
            to_start.grid(column=9, row=5, columnspan=19)
            DropdownMenuStatus.blank_shown = True
        else:
            Notebook_Handler(eval("{}_notebook".format(opt)), opt).add_frames()
            exec(f'DropdownMenuStatus.{opt}_notebook = True')
    DropdownMenuStatus.previous_tab_type = opt


def comparison_button_pressed():
    pop_up = tk.Toplevel()
    pop_up.geometry("440x200")
    newframe = tk.Frame(pop_up)
    newframe.pack()
    newlabel = tk.Label(pop_up, text="select analysis option to compare")
    newlabel.pack()
    listbox = tk.Listbox(pop_up)
    listbox.insert(1, "Total energies")
    listbox.insert(2, "Geometry and displacements")
    listbox.insert(3, "defect VBM, HOMO, and LUMO wfns")
    listbox.insert(4, "Pdos, HOMO/LUMO eigenvalues and defect level")
    listbox.insert(5, "Eigenvalue defect levels, incorperation and ionization energies")
    listbox.pack(expand=True, fill=BOTH)


def _quit():
    root.quit()
    root.destroy()


if __name__ == '__main__':
    root = tk.Tk()
    root.minsize(1515, 850)
    root.maxsize(1515, 850)

    canvas = tk.Canvas(root, width=1500, height=850)
    canvas.grid(columnspan=83, rowspan=47)

    button_quit = tk.Button(root, text="quit", command=_quit)
    button_quit.grid(column=1, row=0, rowspan=2)

    # Header
    A_frame = tk.Frame(root)
    A_frame.grid(column=1, row=2, columnspan=20, rowspan=3)
    lst = [('', "Energy(eV)", "HOMO(eV)", "LUMO(eV)", "Bandgap(eV)"),
           ("Bulk perfect", " ", " ", " ", " ")]
    WDTs = [12, 10, 8, 8, 12]
    t_A = Tables(A_frame, WDTs, lst)

    setheaderpictures()

    to_start = tk.Label(root, text="<- To start, select an option from the dropdown menu")
    to_start.grid(column=9, row=5, columnspan=19)

    for pic in 'blank', 'A', 'B':
        exec(f'Frames.{pic}_pic_frame = tk.Frame(root)')
        exec(
            f'Frames.{pic}_pic_label = tk.Label(Frames().{pic}_pic_frame, image=HandlingPictures().{pic}, compound=CENTER)')

    heading_pictures('start')

    corrections_label = tk.Label(root, text="Here are your corrections:")
    corrections_frame = tk.Frame(root)
    corrections_label.grid(column=73, row=0, columnspan=10)
    corrections_frame.grid(column=70, row=1, columnspan=25, rowspan=6)
    lst = [("Charge", "Point Charge \eV", "Lany-Zunger \eV", "Results  V_M^{scr}"),
           ("1", " ", " ", " "),
           ("2", " ", " ", " "),
           ("3", " ", " ", " "),
           ("4", " ", " ", " ")]
    WDTs = [6, 13, 13, 15]
    t_corrections = Tables(corrections_frame, WDTs, lst)

    # dropdown menu
    options = ["--", "A1", "A2", "B1", "B2"]
    dropdown = ttk.Combobox(root, values=options)
    dropdown.current(0)
    dropdown.bind("<<ComboboxSelected>>", dropdown_menu_command)
    dropdown.grid(column=1, row=5, columnspan=8)

    comparison_text = tk.StringVar()
    comparison_btn = tk.Button(root, textvariable=comparison_text, command=comparison_button_pressed,
                               font="Raleway", fg="black", bg="white")
    comparison_text.set("Comparison")
    comparison_btn.grid(column=80, row=7)

    C_label = tk.Label(root, text="testing, notebook not displayed for option yet")

    A_notebook = ttk.Notebook(root)
    B_notebook = ttk.Notebook(root)

    # set up notebook pages
    for i in 'a', 'b', 'c', 'd', 'e':
        exec(f'Frames.frame_{i}A = tk.Frame(A_notebook, width=1450, height=570)')
        exec(f'Frames.frame_{i}B = tk.Frame(B_notebook, width=1450, height=570)')
        for l in 'A', 'B':
            exec(f'Frames().frame_{i}{l}.grid(columnspan=100, rowspan=60)')
            exec(f'Frames().frame_{i}{l}.grid(sticky="nsew")')
            exec(f'canvas{i}{l} = tk.Canvas(Frames().frame_{i}{l}, width=1450, height=570)')
            exec(f'canvas{i}{l}.grid(columnspan=100, rowspan=60)')