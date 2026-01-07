#!/usr/bin/env python3
"""
Quick check of VPIC output
"""
import errno
import math
import os
import sys

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtWidgets

from mainwindow import Ui_MainWindow

mpl.use('Qt5Agg')


def get_vpic_info():
    """Get information of the VPIC simulation
    """
    with open('./info') as f:
        content = f.readlines()
    f.close()
    vpic_info = {}
    for line in content[1:]:
        if "=" in line:
            line_splits = line.split("=")
        elif ":" in line:
            line_splits = line.split(":")

        tail = line_splits[1].split("\n")
        vpic_info[line_splits[0].strip()] = float(tail[0])
    return vpic_info


vpic_info = get_vpic_info()
hdf5_fields = True  # whether data is in HDF5 format
smoothed_data = False  # whether data is smoothed
if smoothed_data:
    smooth_factor = 24  # smooth factor along each direction
else:
    smooth_factor = 1
dir_smooth_data = "data_smooth"
momentum_field = True  # whether momentum and kinetic energy data are dumped
time_averaged_field = False  # whether it is time-averaged field
turbulence_mixing = False  # whether it has turbulence mixing diagnostics
tmin, tmax = 0, 78
if time_averaged_field:
    tmin = 1
animation_tinterval = 100  # in msec
nt = tmax - tmin + 1


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def div0(a, b):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0]

    From: http://stackoverflow.com/a/35696047/2561161

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self,
                 parent=None,
                 width=5,
                 height=4,
                 dpi=100,
                 axes_h="X",
                 plot_type="Contour"):
        self.fig = Figure(figsize=(width, height),
                          constrained_layout=True,
                          dpi=dpi)
        self.create_axes(axes_h, plot_type)
        super(MplCanvas, self).__init__(self.fig)

    def update_axes(self, axes_h="X", plot_type="Contour"):
        self.create_axes(axes_h, plot_type)

    def create_axes(self, axes_h="X", plot_type="Contour"):
        if plot_type == "Contour":
            widths = [4.8, 0.2]
            spec = self.fig.add_gridspec(nrows=1, ncols=2, width_ratios=widths)
            self.ax_main = self.fig.add_subplot(spec[0, 0])
            self.ax_cbar = self.fig.add_subplot(spec[0, 1])
        elif plot_type in [
                "Contour+" + axes_h + "-Average",
                "Contour+" + axes_h + "-Slice"
        ]:
            widths = [0.2, 4.8, 1.2]
            spec = self.fig.add_gridspec(nrows=1,
                                         ncols=3,
                                         width_ratios=widths,
                                         wspace=0.0,
                                         hspace=0.0)
            self.ax_main = self.fig.add_subplot(spec[0, 1])
            self.ax1d = self.fig.add_subplot(spec[0, 2], sharey=self.ax_main)
            self.ax_cbar = self.fig.add_subplot(spec[0, 0])
        else:
            widths = [4.8, 0.2]
            heights = [4.8, 1.2]
            spec = self.fig.add_gridspec(nrows=2,
                                         ncols=2,
                                         width_ratios=widths,
                                         height_ratios=heights,
                                         wspace=0.0,
                                         hspace=0.0)
            self.ax_main = self.fig.add_subplot(spec[0, 0])
            self.ax_main.tick_params(axis='x', labelbottom=False)
            self.ax1d = self.fig.add_subplot(spec[1, 0], sharex=self.ax_main)
            self.ax_cbar = self.fig.add_subplot(spec[0, 1])


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, *args, obj=None, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)

        # simulation domain
        self.get_domain()

        # check if the simulation is 2D
        self.coords = ["x", "y", "z"]
        self.normal = "y"  # normal direction
        self.is_2d = False  # whether is a 2D simulation
        for c in self.coords:
            if self.vpic_domain["n" + c] == 1:
                self.normal = c
                self.is_2d = True
        self.hv = [c for c in self.coords if c != self.normal]  # in-plane

        # the file type (HDF5 or gda)
        if hdf5_fields:
            self.filetype_comboBox.setCurrentText("HDF5")
        else:
            self.filetype_comboBox.setCurrentText("gda")
            if smoothed_data:
                self.gda_path = dir_smooth_data
            else:
                self.gda_path = "data/"
            self.tframe_loaded = -1
            self.var_loaded = "random_var"

        # whether to automatically update the plot
        self.autoplot_checkBox.setChecked(False)
        self.auto_update = False
        self.autoplot_checkBox.stateChanged.connect(
            self.autoplot_checkBox_change)

        # plot type
        self.plottype_comboBox.currentTextChanged.connect(
            self.plottype_comboBox_change)

        # Raw plot variables
        self.raw_plot_variables()
        self.rawplot_comboBox.currentTextChanged.connect(
            self.rawplot_comboBox_vchange)

        # Diagnostics plot variables
        self.diag_plot = False  # whether to plot diagnostics
        self.diag_var_name = ""
        self.diag_plot_variables()
        self.diagplot_comboBox.currentTextChanged.connect(
            self.diagplot_comboBox_vchange)

        # Create toolbar and canvas
        self.margin = 30
        self.top = 320
        self.middle = self.width() // 2
        self.width_max = self.width() - 2 * self.margin
        self.height_max = self.height() - 2 * self.margin - self.top
        self.plot_vLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.plot_vLayoutWidget.setGeometry(
            QtCore.QRect(self.margin, self.top, self.width_max,
                         self.height_max))
        self.plot_vLayoutWidget.setObjectName("plot_vLayoutWidget")
        self.plot_verticalLayout = QtWidgets.QVBoxLayout(
            self.plot_vLayoutWidget)
        self.plot_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.plot_verticalLayout.setObjectName("plot_verticalLayout")
        self.canvas = MplCanvas(self,
                                width=8,
                                height=8,
                                dpi=100,
                                axes_h="X",
                                plot_type=self.plottype_comboBox.currentText())
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.plot_verticalLayout.addWidget(self.toolbar)
        self.plot_verticalLayout.addWidget(self.canvas)

        # Time slice
        self.tframe_hSlider.setMinimum(tmin)
        self.tframe_hSlider.setMaximum(tmax)
        self.tframe_SpinBox.setMinimum(self.tframe_hSlider.minimum())
        self.tframe_SpinBox.setMaximum(self.tframe_hSlider.maximum())
        self.tframe_hSlider.valueChanged.connect(self.tframe_hSlider_vchange)
        self.tframe_hSlider.sliderPressed.connect(
            self.tframe_hSlider_disconnect)
        self.tframe_hSlider.sliderReleased.connect(
            self.tframe_hSlider_reconnect)
        self.tframe_SpinBox.setKeyboardTracking(False)
        self.tframe_SpinBox.valueChanged.connect(self.tframe_SpinBox_vchange)
        self.var_name = self.rawplot_comboBox.currentText()
        self.plot_type = self.plottype_comboBox.currentText()
        self.tframe = self.tframe_hSlider.value()
        if "fields_interval" not in vpic_info:
            vpic_info["fields_interval"], _ = QtWidgets.QInputDialog.getInt(
                self, "Get fields interval", "Fields interval:", 100, 0,
                10000000, 1)
        self.tindex = self.tframe * int(vpic_info["fields_interval"])

        # x-range
        self.xmin_hScrollBar.valueChanged.connect(self.xmin_hScrollBar_vchange)
        self.xmin_SpinBox.valueChanged.connect(self.xmin_SpinBox_vchange)
        self.xmax_hScrollBar.valueChanged.connect(self.xmax_hScrollBar_vchange)
        self.xmax_SpinBox.valueChanged.connect(self.xmax_SpinBox_vchange)
        self.xmin_hScrollBar.setMinimum(int(self.vpic_domain["xmin"]))
        self.xmin_hScrollBar.setMaximum(int(self.vpic_domain["xmax"]))
        self.xmin_SpinBox.setMinimum(self.vpic_domain["xmin"])
        self.xmin_SpinBox.setMaximum(self.vpic_domain["xmax"])
        self.xmax_hScrollBar.setMinimum(self.xmin_hScrollBar.minimum())
        self.xmax_hScrollBar.setMaximum(self.xmin_hScrollBar.maximum())
        self.xmax_SpinBox.setMinimum(self.xmin_SpinBox.minimum())
        self.xmax_SpinBox.setMaximum(self.xmin_SpinBox.maximum())
        self.xmin_hScrollBar.setSliderPosition(self.xmin_hScrollBar.minimum())
        self.xmax_hScrollBar.setSliderPosition(self.xmax_hScrollBar.maximum())

        # y-range
        self.ymin_hScrollBar.valueChanged.connect(self.ymin_hScrollBar_vchange)
        self.ymin_SpinBox.valueChanged.connect(self.ymin_SpinBox_vchange)
        self.ymax_hScrollBar.valueChanged.connect(self.ymax_hScrollBar_vchange)
        self.ymax_SpinBox.valueChanged.connect(self.ymax_SpinBox_vchange)
        self.ymin_hScrollBar.setMinimum(int(self.vpic_domain["ymin"]))
        self.ymin_hScrollBar.setMaximum(int(self.vpic_domain["ymax"]))
        self.ymin_SpinBox.setMinimum(self.vpic_domain["ymin"])
        self.ymin_SpinBox.setMaximum(self.vpic_domain["ymax"])
        self.ymax_hScrollBar.setMinimum(self.ymin_hScrollBar.minimum())
        self.ymax_hScrollBar.setMaximum(self.ymin_hScrollBar.maximum())
        self.ymax_SpinBox.setMinimum(self.ymin_SpinBox.minimum())
        self.ymax_SpinBox.setMaximum(self.ymin_SpinBox.maximum())
        self.ymin_hScrollBar.setSliderPosition(self.ymin_hScrollBar.minimum())
        self.ymax_hScrollBar.setSliderPosition(self.ymax_hScrollBar.maximum())

        # z-range
        self.zmin_hScrollBar.valueChanged.connect(self.zmin_hScrollBar_vchange)
        self.zmin_SpinBox.valueChanged.connect(self.zmin_SpinBox_vchange)
        self.zmax_hScrollBar.valueChanged.connect(self.zmax_hScrollBar_vchange)
        self.zmax_SpinBox.valueChanged.connect(self.zmax_SpinBox_vchange)
        self.zmin_hScrollBar.setMinimum(int(self.vpic_domain["zmin"]))
        self.zmin_hScrollBar.setMaximum(int(self.vpic_domain["zmax"]))
        self.zmin_SpinBox.setMinimum(self.vpic_domain["zmin"])
        self.zmin_SpinBox.setMaximum(self.vpic_domain["zmax"])
        self.zmax_hScrollBar.setMinimum(self.zmin_hScrollBar.minimum())
        self.zmax_hScrollBar.setMaximum(self.zmin_hScrollBar.maximum())
        self.zmax_SpinBox.setMinimum(self.zmin_SpinBox.minimum())
        self.zmax_SpinBox.setMaximum(self.zmin_SpinBox.maximum())
        self.zmin_hScrollBar.setSliderPosition(self.zmin_hScrollBar.minimum())
        self.zmax_hScrollBar.setSliderPosition(self.zmax_hScrollBar.maximum())

        # create a dictionary for the ranges
        self.range_dist = {
            "xmin_bar": self.xmin_hScrollBar,
            "xmin_box": self.xmin_SpinBox,
            "ymin_bar": self.ymin_hScrollBar,
            "ymin_box": self.ymin_SpinBox,
            "zmin_bar": self.zmin_hScrollBar,
            "zmin_box": self.zmin_SpinBox,
            "xmax_bar": self.xmax_hScrollBar,
            "xmax_box": self.xmax_SpinBox,
            "ymax_bar": self.ymax_hScrollBar,
            "ymax_box": self.ymax_SpinBox,
            "zmax_bar": self.zmax_hScrollBar,
            "zmax_box": self.zmax_SpinBox
        }

        # 1D slices
        self.xslice_hScrollBar.valueChanged.connect(
            self.xslice_hScrollBar_vchange)
        self.xslice_SpinBox.valueChanged.connect(self.xslice_SpinBox_vchange)
        self.yslice_hScrollBar.valueChanged.connect(
            self.yslice_hScrollBar_vchange)
        self.yslice_SpinBox.valueChanged.connect(self.yslice_SpinBox_vchange)
        self.zslice_hScrollBar.valueChanged.connect(
            self.zslice_hScrollBar_vchange)
        self.zslice_SpinBox.valueChanged.connect(self.zslice_SpinBox_vchange)
        self.xslice_hScrollBar.setMinimum(int(self.vpic_domain["xmin"]))
        self.xslice_hScrollBar.setMaximum(int(self.vpic_domain["xmax"]))
        self.xslice_SpinBox.setMinimum(self.vpic_domain["xmin"])
        self.xslice_SpinBox.setMaximum(self.vpic_domain["xmax"])
        xmid = 0.5 * (self.vpic_domain["xmin"] + self.vpic_domain["xmax"])
        self.xslice_hScrollBar.setSliderPosition(int(xmid))
        self.yslice_hScrollBar.setMinimum(int(self.vpic_domain["ymin"]))
        self.yslice_hScrollBar.setMaximum(int(self.vpic_domain["ymax"]))
        self.yslice_SpinBox.setMinimum(self.vpic_domain["ymin"])
        self.yslice_SpinBox.setMaximum(self.vpic_domain["ymax"])
        ymid = 0.5 * (self.vpic_domain["ymin"] + self.vpic_domain["ymax"])
        self.yslice_hScrollBar.setSliderPosition(int(ymid))
        self.zslice_hScrollBar.setMinimum(int(self.vpic_domain["zmin"]))
        self.zslice_hScrollBar.setMaximum(int(self.vpic_domain["zmax"]))
        self.zslice_SpinBox.setMinimum(self.vpic_domain["zmin"])
        self.zslice_SpinBox.setMaximum(self.vpic_domain["zmax"])
        zmid = 0.5 * (self.vpic_domain["zmin"] + self.vpic_domain["zmax"])
        self.zslice_hScrollBar.setSliderPosition(int(zmid))

        # create a dictionary for the slices
        self.slice_dist = {
            "xslice_bar": self.xslice_hScrollBar,
            "xslice_box": self.xslice_SpinBox,
            "yslice_bar": self.yslice_hScrollBar,
            "yslice_box": self.yslice_SpinBox,
            "zslice_bar": self.zslice_hScrollBar,
            "zslice_box": self.zslice_SpinBox
        }

        # 2D plane
        if self.is_2d:
            self.plane_comboBox.setDisabled(True)
        self.plane_comboBox.currentTextChanged.connect(
            self.plane_comboBox_vchange)
        self.plane_hScrollBar.sliderPressed.connect(
            self.plane_hScrollBar_disconnect)
        self.plane_hScrollBar.sliderReleased.connect(
            self.plane_hScrollBar_reconnect)
        self.plane_hScrollBar.valueChanged.connect(
            self.plane_hScrollBar_vchange)
        self.plane_SpinBox.setKeyboardTracking(False)
        self.plane_SpinBox.valueChanged.connect(self.plane_SpinBox_vchange)
        self.set_normal_plane()
        self.set_plane_index()

        # read the field data
        self.read_data(self.var_name, self.tindex)

        # plot button
        self.plotButton.clicked.connect(self.update_plot)

        # checkbox for saving JPEGs during animation
        self.savejpg_checkBox.setChecked(False)
        self.save_jpegs = False
        self.savejpg_checkBox.stateChanged.connect(
            self.savejpg_checkBox_change)

        # animation buttons
        self.start_animateButton.clicked.connect(self.start_animation)
        self.stop_animateButton.clicked.connect(self.stop_animation)
        self.continue_animateButton.clicked.connect(self.continue_animation)
        self.stop_animateButton.setDisabled(True)
        self.continue_animateButton.setDisabled(True)
        self.is_animation = False

        # checkbox for fixing colormap
        self.cmap_checkBox.setChecked(False)
        self.fix_cmap = False
        self.cmap_checkBox.stateChanged.connect(self.cmap_checkBox_change)
        self.cmap_LineEdit.setText("coolwarm")

        # checkbox for fixing the colorbar range
        self.cbar_checkBox.setChecked(False)
        self.fix_cbar_range = False
        self.cbar_checkBox.stateChanged.connect(self.cbar_checkBox_change)

    def plottype_comboBox_change(self, value):
        self.plot_type = value
        self.canvas.fig.clf()
        self.canvas.update_axes(self.hv[0].upper(), value)
        self.update_plot()

    def rawplot_comboBox_vchange(self, value):
        if not self.diag_plot:
            self.var_name = self.rawplot_comboBox.currentText()
            self.read_data(self.var_name, self.tindex)
            self.update_plot()

    def diagplot_comboBox_vchange(self, value):
        self.diag_var_name = self.diagplot_comboBox.currentText()
        if self.diag_var_name != "":
            self.diag_plot = True
            self.var_name = self.diag_var_name
            self.read_data(self.var_name, self.tindex)
            self.update_plot()
        else:
            self.diag_plot = False
            self.var_name = self.rawplot_comboBox.currentText()
            self.read_data(self.var_name, self.tindex)
            self.update_plot()

    def savejpg_checkBox_change(self):
        self.save_jpegs = not self.save_jpegs

    def cmap_checkBox_change(self):
        self.fix_cmap = not self.fix_cmap

    def cbar_checkBox_change(self):
        self.fix_cbar_range = not self.fix_cbar_range
        if self.fix_cbar_range:
            self.cbar_min_LineEdit.setText(str(np.min(self.field_2d)))
            self.cbar_max_LineEdit.setText(str(np.max(self.field_2d)))

    def autoplot_checkBox_change(self):
        self.auto_update = not self.auto_update

    def tframe_hSlider_vchange(self, value):
        self.tframe_SpinBox.setValue(value)
        self.tframe = value
        self.tindex = self.tframe * int(vpic_info["fields_interval"])
        self.read_data(self.var_name, self.tindex)
        if self.auto_update:
            self.update_plot()

    def tframe_hSlider_disconnect(self):
        self.sender().valueChanged.disconnect()

    def tframe_hSlider_reconnect(self):
        self.sender().valueChanged.connect(self.tframe_hSlider_vchange)
        self.sender().valueChanged.emit(self.sender().value())

    def tframe_SpinBox_vchange(self, value):
        self.tframe_hSlider.setValue(value)

    def xmin_hScrollBar_vchange(self, value):
        self.xmin_SpinBox.setValue(float(value))

    def xmin_SpinBox_vchange(self, value):
        self.xmin_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def ymin_hScrollBar_vchange(self, value):
        self.ymin_SpinBox.setValue(float(value))

    def ymin_SpinBox_vchange(self, value):
        self.ymin_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def zmin_hScrollBar_vchange(self, value):
        self.zmin_SpinBox.setValue(float(value))

    def zmin_SpinBox_vchange(self, value):
        self.zmin_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def xmax_hScrollBar_vchange(self, value):
        self.xmax_SpinBox.setValue(float(value))

    def xmax_SpinBox_vchange(self, value):
        self.xmax_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def ymax_hScrollBar_vchange(self, value):
        self.ymax_SpinBox.setValue(float(value))

    def ymax_SpinBox_vchange(self, value):
        self.ymax_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def zmax_hScrollBar_vchange(self, value):
        self.zmax_SpinBox.setValue(float(value))

    def zmax_SpinBox_vchange(self, value):
        self.zmax_hScrollBar.setValue(int(value))
        if self.auto_update:
            self.update_plot()

    def xslice_hScrollBar_vchange(self, value):
        self.xslice_SpinBox.setValue(float(value))

    def xslice_SpinBox_vchange(self, value):
        self.xslice_hScrollBar.setValue(int(value))
        if self.plot_type == "Contour+X-Slice":
            if self.auto_update:
                self.update_plot()

    def yslice_hScrollBar_vchange(self, value):
        self.yslice_SpinBox.setValue(float(value))

    def yslice_SpinBox_vchange(self, value):
        self.yslice_hScrollBar.setValue(int(value))
        if self.plot_type == "Contour+Y-Slice":
            if self.auto_update:
                self.update_plot()

    def zslice_hScrollBar_vchange(self, value):
        self.zslice_SpinBox.setValue(float(value))

    def zslice_SpinBox_vchange(self, value):
        self.zslice_hScrollBar.setValue(int(value))
        if self.plot_type == "Contour+Z-Slice":
            if self.auto_update:
                self.update_plot()

    def plane_comboBox_vchange(self, value):
        current_plane = self.plane_comboBox.currentText()
        self.normal = current_plane[0].lower()
        self.set_normal_plane()

    def set_normal_plane(self):
        self.hv = [c for c in self.coords if c != self.normal]  # in-plane
        self.plane_comboBox.setCurrentText(self.normal.upper() + "-Plane")
        self.plane_hScrollBar.setMinimum(
            int(self.vpic_domain[self.normal + "min"]))
        self.plane_hScrollBar.setMaximum(
            int(self.vpic_domain[self.normal + "max"]))
        old_plane = self.plane_hScrollBar.value()
        self.plane_SpinBox.setMinimum(self.vpic_domain[self.normal + "min"])
        self.plane_SpinBox.setMaximum(self.vpic_domain[self.normal + "max"])
        mid = 0.5 * (self.vpic_domain[self.normal + "min"] +
                     self.vpic_domain[self.normal + "max"])
        self.plane_hScrollBar.setSliderPosition(int(mid))
        if int(mid) == old_plane:  # force update
            self.set_plane_index()
            self.read_data(self.var_name, self.tindex)
            if self.auto_update:
                self.update_plot()
        self.plottype_comboBox.setItemText(
            1, "Contour+" + self.hv[0].upper() + "-Average")
        self.plottype_comboBox.setItemText(
            2, "Contour+" + self.hv[0].upper() + "-Slice")
        self.plottype_comboBox.setItemText(
            3, "Contour+" + self.hv[1].upper() + "-Slice")

    def set_plane_index(self):
        plane_coord = self.plane_SpinBox.value()
        cmin = self.vpic_domain[self.normal + "min"]
        self.plane_index = int(
            (plane_coord - cmin) / self.vpic_domain["d" + self.normal])
        if self.plane_index == self.vpic_domain["n" + self.normal]:
            self.plane_index = self.vpic_domain["n" + self.normal] - 1

    def plane_hScrollBar_disconnect(self):
        self.sender().valueChanged.disconnect()

    def plane_hScrollBar_reconnect(self):
        self.sender().valueChanged.connect(self.plane_hScrollBar_vchange)
        self.sender().valueChanged.emit(self.sender().value())

    def plane_hScrollBar_vchange(self, value):
        self.plane_SpinBox.setValue(float(value))

    def plane_SpinBox_vchange(self, value):
        self.plane_hScrollBar.setValue(int(value))
        self.set_plane_index()
        self.read_data(self.var_name, self.tindex)
        if self.auto_update:
            self.update_plot()

    def raw_plot_variables(self):
        if hdf5_fields:
            self.fields_list = ["cbx", "cby", "cbz", "absb", "ex", "ey", "ez"]
            self.jlist = ["jx", "jy", "jz", "absj"]
            self.ehydro_list = [
                "ne", "vex", "vey", "vez", "pexx", "pexy", "pexz", "peyx",
                "peyy", "peyz", "pezx", "pezy", "pezz"
            ]
            self.Hhydro_list = [
                "ni", "vix", "viy", "viz", "pixx", "pixy", "pixz", "piyx",
                "piyy", "piyz", "pizx", "pizy", "pizz"
            ]
            if momentum_field:
                self.ehydro_list += [
                    "uex",
                    "uey",
                    "uez",
                ]
                self.Hhydro_list += [
                    "uix",
                    "uiy",
                    "uiz",
                ]
            self.hydro_list = self.jlist + self.ehydro_list + self.Hhydro_list
            self.var_list = self.fields_list + self.hydro_list
            self.var_list = sorted(set(self.var_list))
        else:
            flist = [
                _ for _ in os.listdir(self.gda_path) if _.endswith(".gda")
            ]
            if len([name for name in flist if name[:2] == "bx"]) == 1:
                # all frames are save in the same file
                self.var_list = sorted([f[:-4] for f in flist])
                self.single_gda = True  # all time frames are in the same file
            else:
                var_list = []
                for f in flist:
                    var_list.append(f.split("_")[0])
                self.var_list = sorted(set(var_list))
                self.single_gda = False
        _translate = QtCore.QCoreApplication.translate
        for ivar, var in enumerate(self.var_list):
            self.rawplot_comboBox.addItem("")
            self.rawplot_comboBox.setItemText(ivar,
                                              _translate("MainWindow", var))

    def diag_plot_variables(self):
        self.diag_var_list = ["", "jdotE"]
        if turbulence_mixing:
            self.diag_var_list.append("emix")
        _translate = QtCore.QCoreApplication.translate
        for ivar, var in enumerate(self.diag_var_list):
            self.diagplot_comboBox.addItem("")
            self.diagplot_comboBox.setItemText(ivar,
                                               _translate("MainWindow", var))

    def get_domain(self):
        """Get VPIC simulation domain
        """
        self.vpic_domain = {}
        self.vpic_domain["xmin"] = 0.0
        self.vpic_domain["xmax"] = vpic_info["Lx/de"]
        self.vpic_domain["ymin"] = -0.5 * vpic_info["Ly/de"]
        self.vpic_domain["ymax"] = 0.5 * vpic_info["Ly/de"]
        self.vpic_domain["zmin"] = -0.5 * vpic_info["Lz/de"]
        self.vpic_domain["zmax"] = 0.5 * vpic_info["Lz/de"]
        self.vpic_domain["nx"] = int(vpic_info["nx"]) // smooth_factor
        self.vpic_domain["ny"] = int(vpic_info["ny"]) // smooth_factor
        self.vpic_domain["nz"] = int(vpic_info["nz"]) // smooth_factor
        self.vpic_domain["dx"] = vpic_info["dx/de"] * smooth_factor
        self.vpic_domain["dy"] = vpic_info["dy/de"] * smooth_factor
        self.vpic_domain["dz"] = vpic_info["dz/de"] * smooth_factor
        for i in ["x", "y", "z"]:
            n = "n" + i
            d = "d" + i
            if vpic_info[n] < smooth_factor:
                self.vpic_domain[n] = 1
                self.vpic_domain[d] = vpic_info[d + "/de"] * vpic_info[n]
        hdx = 0.5 * self.vpic_domain["dx"]
        hdy = 0.5 * self.vpic_domain["dy"]
        hdz = 0.5 * self.vpic_domain["dz"]
        self.vpic_domain["xgrid"] = np.linspace(self.vpic_domain["xmin"] + hdx,
                                                self.vpic_domain["xmax"] - hdx,
                                                self.vpic_domain["nx"])
        self.vpic_domain["ygrid"] = np.linspace(self.vpic_domain["ymin"] + hdy,
                                                self.vpic_domain["ymax"] - hdy,
                                                self.vpic_domain["ny"])
        self.vpic_domain["zgrid"] = np.linspace(self.vpic_domain["zmin"] + hdz,
                                                self.vpic_domain["zmax"] - hdz,
                                                self.vpic_domain["nz"])

    def read_gda_file(self, vname, tindex):
        """read fields or hydro in gda format

        Args:
            vname (string): variable name
            tindex (int): time index
        """
        if self.is_2d:
            nh = self.vpic_domain["n" + self.hv[0]]
            nv = self.vpic_domain["n" + self.hv[1]]
            if self.single_gda:
                fname = self.gda_path + vname + ".gda"
                fdata = np.fromfile(fname,
                                    dtype=np.float32,
                                    count=nh * nv,
                                    offset=nh * nv * self.tframe * 4)
            else:
                fname = self.gda_path + vname + "_" + str(self.tindex) + ".gda"
                fdata = np.fromfile(fname, dtype=np.float32, count=-1)
            field_2d = fdata.reshape([nv, nh]).T
            return field_2d, field_2d
        else:
            nx = self.vpic_domain["nx"]
            ny = self.vpic_domain["ny"]
            nz = self.vpic_domain["nz"]
            ntot = nx * ny * nz
            # Since slicing of gda file is incontinent, we read all the 3D
            # data cube and take slices in the memory.
            if self.tframe_loaded != self.tframe or self.var_loaded != vname:
                if self.single_gda:
                    fname = self.gda_path + vname + ".gda"
                    field_3d = np.fromfile(fname,
                                           dtype=np.float32,
                                           count=ntot,
                                           offset=ntot * self.tframe *
                                           4).reshape([nz, ny, nx])
                else:
                    fname = self.gda_path + vname + "_" + str(
                        self.tindex) + ".gda"
                    field_3d = np.fromfile(fname, dtype=np.float32,
                                           count=-1).reshape([nz, ny, nx])
            self.tframe_loaded = self.tframe
            self.var_loaded = vname
            if self.normal == 'x':
                field_2d = field_3d[:, :, self.plane_index].T
            elif self.normal == 'y':
                field_2d = field_3d[:, self.plane_index, :].T
            else:
                field_2d = field_3d[self.plane_index, :, :].T
            return field_2d, field_3d

    def read_fields(self, vname, tindex):
        """read electric and magnetic fields in HDF5 format

        Args:
            vname (string): variable name
            tindex (int): time index
        """
        if smoothed_data:
            fname = ("./" + dir_smooth_data + "/fields_" + str(tindex) + ".h5")
        else:
            if time_averaged_field:
                fdir = "./fields-avg-hdf5/T." + str(tindex) + "/"
            else:
                fdir = "./field_hdf5/T." + str(tindex) + "/"
            fname = fdir + "fields_" + str(tindex) + ".h5"
        with h5py.File(fname, 'r') as fh:
            group = fh["Timestep_" + str(tindex)]
            if vname == "absb":
                bvec = {}
                for var in ["cbx", "cby", "cbz"]:
                    dset = group[var]
                    if self.normal == 'x':
                        bvec[var] = dset[self.plane_index, :, :]
                    elif self.normal == 'y':
                        bvec[var] = dset[:, self.plane_index, :]
                    else:
                        bvec[var] = dset[:, :, self.plane_index]
                field_2d = np.sqrt(bvec["cbx"]**2 + bvec["cby"]**2 +
                                   bvec["cbz"]**2)
            else:
                dset = group[vname]
                if self.normal == 'x':
                    field_2d = dset[self.plane_index, :, :]
                elif self.normal == 'y':
                    field_2d = dset[:, self.plane_index, :]
                else:
                    field_2d = dset[:, :, self.plane_index]
        return field_2d

    def read_current_density_species(self, vname, tindex, species):
        """read current density associated with one species

        Args:
            vname (string): variable name
            tindex (int): time index
            species (string): particle species
        """
        if smoothed_data:
            fname = ("./" + dir_smooth_data + "/hydro_" + species + "_" +
                     str(tindex) + ".h5")
        else:
            if time_averaged_field:
                fdir = "./hydro-avg-hdf5/T." + str(tindex) + "/"
            else:
                fdir = "./hydro_hdf5/T." + str(tindex) + "/"
            fname = fdir + "hydro_" + species + "_" + str(tindex) + ".h5"
        j2d = {}
        with h5py.File(fname, 'r') as fh:
            group = fh["Timestep_" + str(tindex)]
            if vname == "absj":
                for var in ["jx", "jy", "jz"]:
                    dset = group[var]
                    if self.normal == 'x':
                        j2d[var] = dset[self.plane_index, :, :]
                    elif self.normal == 'y':
                        j2d[var] = dset[:, self.plane_index, :]
                    else:
                        j2d[var] = dset[:, :, self.plane_index]
            else:
                dset = group[vname]
                if self.normal == 'x':
                    j2d["jdir"] = dset[self.plane_index, :, :]
                elif self.normal == 'y':
                    j2d["jdir"] = dset[:, self.plane_index, :]
                else:
                    j2d["jdir"] = dset[:, :, self.plane_index]
        return j2d

    def read_current_density(self, vname, tindex):
        """read current density

        Args:
            vname (string): variable name
            tindex (int): time index
        """
        # Electron
        if turbulence_mixing:
            j2d = self.read_current_density_species(vname, tindex,
                                                    "electronTop")
            jtmp = self.read_current_density_species(vname, tindex,
                                                     "electronBot")
            for var in j2d:
                j2d[var] += jtmp[var]
        else:
            j2d = self.read_current_density_species(vname, tindex, "electron")

        # Ion
        if turbulence_mixing:
            jtmp = self.read_current_density_species(vname, tindex, "ionTop")
            for var in j2d:
                j2d[var] += jtmp[var]
            jtmp = self.read_current_density_species(vname, tindex, "ionBot")
            for var in j2d:
                j2d[var] += jtmp[var]
        else:
            jtmp = self.read_current_density_species(vname, tindex, "ion")
            for var in j2d:
                j2d[var] += jtmp[var]

        if vname == "absj":
            field_2d = np.sqrt(j2d["jx"]**2 + j2d["jy"]**2 + j2d["jz"]**2)
        else:
            field_2d = j2d["jdir"]

        return field_2d

    def read_hydro_species(self, vname, tindex, species):
        """Read the hydro data of one species

        Args:
            vname (string): variable name
            tindex (int): time index
            species (string): particle species
        """
        if vname in self.ehydro_list:
            if smoothed_data:
                fname = ("./" + dir_smooth_data + "/hydro_" + species + "_" +
                         str(tindex) + ".h5")
            else:
                if time_averaged_field:
                    fdir = "./hydro-avg-hdf5/T." + str(tindex) + "/"
                else:
                    fdir = "./hydro_hdf5/T." + str(tindex) + "/"
                fname = fdir + "hydro_" + species + "_" + str(tindex) + ".h5"
        else:
            if smoothed_data:
                fname = ("./" + dir_smooth_data + "/hydro_" + species + "_" +
                         str(tindex) + ".h5")
            else:
                if time_averaged_field:
                    fdir = "./hydro-avg-hdf5/T." + str(tindex) + "/"
                else:
                    fdir = "./hydro_hdf5/T." + str(tindex) + "/"
                fname = fdir + "hydro_" + species + "_" + str(tindex) + ".h5"

        hydro = {}
        keys = []
        with h5py.File(fname, 'r') as fh:
            group = fh["Timestep_" + str(tindex)]
            if vname[0] == "n":
                keys.append("rho")
            elif vname[0] == 'v':  # for velocity
                keys.append("rho")
                keys.append("j" + vname[-1])
            elif vname[0] == 'u':  # for four-velocity
                keys.append("rho")
                keys.append("p" + vname[-1])
            else:  # for pressure tensor
                keys.append("rho")
                keys.append("j" + vname[-2])
                keys.append("p" + vname[-1])
                vtmp = "t" + vname[2:]
                if vtmp in group:
                    keys.append(vtmp)
                else:
                    keys.append("t" + vname[-1] + vname[-2])
            for key in keys:
                dset = group[key]
                if self.normal == 'x':
                    hydro[key] = dset[self.plane_index, :, :]
                elif self.normal == 'y':
                    hydro[key] = dset[:, self.plane_index, :]
                else:
                    hydro[key] = dset[:, :, self.plane_index]
        return hydro

    def read_hydro(self, vname, tindex):
        """Read hydro data from file

        Args:
            vname (string): variable name
            tindex (int): time index
        """
        if vname in self.ehydro_list:
            pmass = 1.0
            if turbulence_mixing:
                hydro = self.read_hydro_species(vname, tindex, "electronTop")
                hydro_bot = self.read_hydro_species(vname, tindex,
                                                    "electronBot")
                for var in hydro:
                    hydro[var] += hydro_bot[var]
            else:
                hydro = self.read_hydro_species(vname, tindex, "electron")
        else:
            pmass = vpic_info["mi/me"]
            if turbulence_mixing:
                hydro = self.read_hydro_species(vname, tindex, "ionTop")
                hydro_bot = self.read_hydro_species(vname, tindex, "ionBot")
                for var in hydro:
                    hydro[var] += hydro_bot[var]
            else:
                hydro = self.read_hydro_species(vname, tindex, "ion")

        if vname[0] == 'n':  # number density
            field_2d = np.abs(hydro["rho"])
        elif vname[0] == 'v':  # velocity
            field_2d = div0(hydro["j" + vname[-1]], hydro["rho"])
        elif vname[0] == 'u':  # four-velocity
            field_2d = div0(hydro["p" + vname[-1]], (pmass * np.abs(hydro["rho"])))
        else:  # pressure tensor
            vtmp = "t" + vname[2:]
            if vtmp in hydro:
                tvar = vtmp
            else:
                tvar = "t" + vname[-1] + vname[-2]
            jvar = "j" + vname[-2]
            pvar = "p" + vname[-1]
            field_2d = hydro[tvar] - div0(hydro[jvar], hydro["rho"]) * hydro[pvar]
        return field_2d

    def get_jdote(self, tindex):
        """get the diagnostics data of j.E
        """
        if hdf5_fields:
            jx = self.read_current_density("jx", tindex)
            jy = self.read_current_density("jy", tindex)
            jz = self.read_current_density("jz", tindex)
            ex = self.read_fields("ex", tindex)
            ey = self.read_fields("ey", tindex)
            ez = self.read_fields("ez", tindex)
            field_2d = jx * ex + jy * ey + jz * ez
            return field_2d, field_2d
        else:
            j2d, j3d = self.read_gda_file("jx", tindex)
            e2d, e3d = self.read_gda_file("ex", tindex)
            field_2d = j2d * e2d
            field_3d = j2d * e3d
            j2d, j3d = self.read_gda_file("jy", tindex)
            e2d, e3d = self.read_gda_file("ey", tindex)
            field_2d += j2d * e2d
            field_3d += j2d * e3d
            j2d, j3d = self.read_gda_file("jz", tindex)
            e2d, e3d = self.read_gda_file("ez", tindex)
            field_2d += j2d * e2d
            field_3d += j2d * e3d
            return field_2d, field_3d

    def electron_mixing_fraction(self, tindex):
        """get the electron mixing fraction

        (ne_bot - ne_top) / (ne_bot + ne_top)
        """
        ne_top = self.read_hydro_species("ne", tindex, "electronTop")
        ne_bot = self.read_hydro_species("ne", tindex, "electronBot")
        field_2d = (ne_bot["rho"] - ne_top["rho"]) / (ne_bot["rho"] +
                                                      ne_top["rho"])
        return field_2d, field_2d

    def read_data(self, vname, tindex):
        """Read data from file

        Args:
            vname (string): variable name
            tindex (int): time index
        """
        if self.diag_plot:
            if self.diag_var_name == "jdotE":
                self.field_2d, self.field_3d = self.get_jdote(tindex)
            elif self.diag_var_name == "emix":
                self.field_2d, self.field_3d = self.electron_mixing_fraction(
                    tindex)
        else:
            if hdf5_fields:
                if vname in self.fields_list:  # electric and magnetic fields
                    self.field_2d = self.read_fields(vname, tindex)
                elif vname in self.jlist:  # current density
                    self.field_2d = self.read_current_density(vname, tindex)
                else:  # density, velocity, momentum, pressure tensor
                    self.field_2d = self.read_hydro(vname, tindex)
            else:
                self.field_2d, self.field_3d = self.read_gda_file(
                    vname, tindex)

    def update_plot(self):
        if self.fix_cbar_range:
            dmin = float(self.cbar_min_LineEdit.text())
            dmax = float(self.cbar_max_LineEdit.text())
        else:
            if np.any(self.field_2d < 0) and np.any(self.field_2d > 0):
                dmax = min(-self.field_2d.min(), self.field_2d.max())
                dmin = -dmax
            else:
                dmin, dmax = self.field_2d.min(), self.field_2d.max()

        if self.fix_cmap:
            cmap = self.cmap_LineEdit.text()
        else:
            if dmin < 0 and dmax > 0:
                cmap = "coolwarm"
            else:
                cmap = "viridis"

        h = self.hv[0]
        v = self.hv[1]
        hmin = self.range_dist[h + "min_box"].value()
        hmax = self.range_dist[h + "max_box"].value()
        vmin = self.range_dist[v + "min_box"].value()
        vmax = self.range_dist[v + "max_box"].value()
        lhp = hmax - hmin
        lvp = vmax - vmin

        hmin_s = hmin - self.vpic_domain[h + "min"]
        hmax_s = hmax - self.vpic_domain[h + "min"]
        vmin_s = vmin - self.vpic_domain[v + "min"]
        vmax_s = vmax - self.vpic_domain[v + "min"]
        ihs = math.floor(hmin_s / self.vpic_domain["d" + h])
        ihe = math.ceil(hmax_s / self.vpic_domain["d" + h])
        ivs = math.floor(vmin_s / self.vpic_domain["d" + v])
        ive = math.ceil(vmax_s / self.vpic_domain["d" + v])

        # rescale the plot when additional panels are included
        orientation = "vertical"
        if ("Contour+" + self.hv[0].upper()) in self.plot_type:
            lhp *= 1.25
        elif self.plot_type == "Contour+" + self.hv[1].upper() + "-Slice":
            lvp *= 1.25
        denp = max(lhp / self.width_max, lvp / self.height_max)
        canvas_h = int(lhp / denp)
        canvas_v = int(lvp / denp)
        if canvas_h < self.width_max // 3:
            canvas_h = self.width_max // 3
        if canvas_v < self.height_max // 3:
            canvas_v = self.height_max // 3
        if self.plot_type == "Contour+" + self.hv[1].upper() + "-Slice":
            if canvas_v < self.height_max // 2:
                canvas_v = self.height_max // 2
        canvas_l = self.middle - canvas_h // 2
        self.canvas.ax_main.clear()
        self.plot_vLayoutWidget.setGeometry(
            QtCore.QRect(canvas_l, self.top, canvas_h, canvas_v))

        im = self.canvas.ax_main.imshow(self.field_2d[ihs:ihe, ivs:ive].T,
                                        extent=[hmin, hmax, vmin, vmax],
                                        vmin=dmin,
                                        vmax=dmax,
                                        cmap=cmap,
                                        aspect='auto',
                                        origin='lower',
                                        interpolation='none')
        self.canvas.ax_cbar.clear()
        self.canvas.fig.colorbar(im,
                                 cax=self.canvas.ax_cbar,
                                 orientation=orientation)

        # 1D plot
        if self.plot_type == "Contour+" + self.hv[0].upper() + "-Average":
            self.field_1d = np.mean(self.field_2d[ihs:ihe, :], axis=0)
        elif self.plot_type == "Contour+" + self.hv[0].upper() + "-Slice":
            hslice = self.slice_dist[h + "slice_box"].value()
            ih_slice = int((hslice - self.vpic_domain[h + "min"]) /
                           self.vpic_domain["d" + h])
            if ih_slice == self.vpic_domain["n" + h]:
                ih_slice = self.vpic_domain["n" + h] - 1
            self.field_1d = self.field_2d[ih_slice, :]
            ylim = self.canvas.ax_main.get_ylim()
            self.canvas.ax_main.plot([hslice, hslice], ylim, color='w')
            self.canvas.ax_main.set_ylim(ylim)
        elif self.plot_type == "Contour+" + self.hv[1].upper() + "-Slice":
            vslice = self.slice_dist[v + "slice_box"].value()
            iv_slice = int((vslice - self.vpic_domain[v + "min"]) /
                           self.vpic_domain["d" + v])
            if iv_slice == self.vpic_domain["n" + v]:
                iv_slice = self.vpic_domain["n" + v] - 1
            self.field_1d = self.field_2d[:, iv_slice]
            xlim = self.canvas.ax_main.get_xlim()
            self.canvas.ax_main.plot(xlim, [vslice, vslice], color='w')
            self.canvas.ax_main.set_xlim(xlim)
        if ("Contour+" + self.hv[0].upper()) in self.plot_type:
            self.canvas.ax1d.clear()
            self.canvas.ax1d.plot(self.field_1d[ivs:ive],
                                  self.vpic_domain[v + "grid"][ivs:ive])
            self.canvas.ax1d.set_ylim(self.canvas.ax_main.get_ylim())
        elif self.plot_type == "Contour+" + self.hv[1].upper() + "-Slice":
            self.canvas.ax1d.clear()
            self.canvas.ax1d.plot(self.vpic_domain[h + "grid"][ihs:ihe],
                                  self.field_1d[ihs:ihe])
            self.canvas.ax1d.set_xlim(self.canvas.ax_main.get_xlim())

        if self.plot_type == "Contour+" + self.hv[1].upper() + "-Slice":
            self.canvas.ax1d.set_xlabel(r"$" + h + "/d_e$", fontsize=12)
        else:
            self.canvas.ax_main.set_xlabel(r"$" + h + "/d_e$", fontsize=12)
        self.canvas.ax_main.set_ylabel(r"$" + v + "/d_e$", fontsize=12)

        # Trigger the canvas to update and redraw.
        self.canvas.draw()
        # plt.tight_layout()

        # save the figure
        if self.save_jpegs and self.is_animation:
            img_dir = ("./img/" + self.var_name + "/" + self.normal + "_" +
                       str(self.plane_index) + "/")
            mkdir_p(img_dir)
            fname = img_dir + self.var_name + "_" + str(self.tframe) + ".jpg"
            self.canvas.fig.savefig(fname)

    def start_animation(self):
        self.timer = QtCore.QTimer()
        self.timer.setInterval(animation_tinterval)
        self.timer.timeout.connect(self.tick_timer)
        self.timer.start()
        self.tframe_hSlider.setValue(tmin)
        self.stop_animateButton.setDisabled(False)
        self.continue_animateButton.setDisabled(False)
        self.is_animation = True
        self.auto_update_old = self.auto_update
        self.autoplot_checkBox.setChecked(True)

    def tick_timer(self):
        tframe = ((self.tframe - tmin + 1) % nt) + tmin
        if self.tframe == tmax:
            self.savejpg_checkBox.setChecked(False)
        self.tframe_hSlider.setValue(tframe)

    def stop_animation(self):
        self.timer.stop()
        self.is_animation = False
        self.autoplot_checkBox.setChecked(self.auto_update_old)

    def continue_animation(self):
        self.timer.start()
        self.is_animation = True
        self.auto_update_old = self.auto_update
        self.autoplot_checkBox.setChecked(True)


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
