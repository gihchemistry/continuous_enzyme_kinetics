# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 17:12:32 2018

@author: molp
"""

#from functools import lru_cache

from os import listdir
from os.path import dirname, join

import re
import math
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

from bokeh.io import curdoc
from bokeh.layouts import row, column, widgetbox, layout
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Div, BasicTickFormatter
from bokeh.models.widgets import DataTable, Select, TableColumn, Button, TextInput, RadioButtonGroup, RangeSlider
from bokeh.plotting import figure

from io import StringIO
import base64

fit_choice = RadioButtonGroup(labels=["Maximize Slope Magnitude", "Linear Fit", "Logarithmic Fit"], active=0, width=375)
fit_choice_dict = {}

# equations & methods for curve fitting

def linear(x, m , b):
    '''
    linear equation for raw data samples
    '''
    return m*x + b

def logarithmic(x, yo, b, to):
    '''
    logarithmic equation from Lu & Fei et. al, 2003
    '''
    return yo + b*np.log(1 + x*to)

def mmfit(x, km, vmax):
    '''
    Michaelis Menten equation
    '''
    return vmax * x / (km + x)

def icfit(x, bottom, top, IcFifty):
    '''
    IC50 equation
    '''
    return bottom + (top-bottom)/(1+(x/IcFifty))

def subtract_blank(data, sample):
    '''
    subtract blank sample slope from model y-values
    '''
    value = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(sample))[0])
    new_y = data.y - data[data.x == value].y.values[0]
    data.y = new_y
    return data

# methods for data import and analysis

def get_data(x, y):
    data = df[[x, y]]
    data.columns = ['x', 'y']
    return data

def fit_linear(x, y):
    xvalues = np.array(x)
    yvalues = np.array(y)
    popt, pcov = curve_fit(linear, xvalues, yvalues)
    xix = [min(x), max(x)] # plot line across entire range of x
    xfit = np.linspace(xix[0], xix[1], len(xvalues))
    yfit = linear(xfit, *popt)
    fit = pd.DataFrame({'xfit' : xfit, 'yfit' : yfit})
    residuals = np.array(fit.yfit) - yvalues

    return popt[0], fit, residuals

def fit_raw(x, y):
    x = np.array(x)
    y = np.array(y)
    # smooth data
    yfit = UnivariateSpline(x, y)(x)
    # find steepest part of curve
    derivative = np.abs(np.diff(yfit)/np.diff(x))
    threshold = 0.7*(np.max(derivative)-np.min(derivative))+np.min(derivative)
    ix_arr = np.where(derivative > threshold)[0]
    ix_arr = list(np.sort(np.array(ix_arr)))
    xvalues, yvalues = x[ix_arr], y[ix_arr]
    tmp_df = pd.DataFrame(data={'x' : xvalues, 'y' : yvalues})
    tmp_df = tmp_df.sort_values('x')
    xvalues, yvalues = tmp_df.x, tmp_df.y
    popt, pcov = curve_fit(linear, xvalues, yvalues)
    xix = [min(x), max(x)] # plot line across entire range of x
    xfit = np.linspace(xix[0], xix[1], len(yfit))
    yfit = linear(xfit, *popt)
    fit = pd.DataFrame({'xfit' : xfit, 'yfit' : yfit})
    residuals = np.array(yfit) - y

    return popt[0], fit, residuals

def fit_raw_logarithmic(x, y):
    x = np.array(x)
    y = np.array(y)
    n = 0
    try:
        n = eval(offset_time.value)
    except:
        pass
    x = np.array([i+n for i in x])
    popt, pcov = curve_fit(logarithmic, x, y, maxfev=100000)
    #slope = popt[1]/popt[2]*0.001
    xix = [np.min(x), np.max(x)] # plot line across entire range of x
    xfit = np.linspace(0, xix[1], len(x))
    yfit = logarithmic(xfit, *popt)
    slope = np.mean((np.diff(yfit)/np.diff(xfit))[:n+1])
    xfit = [xi-n for xi in xfit]
    fit = pd.DataFrame({'xfit' : xfit, 'yfit' : yfit})
    residuals = np.array(fit.yfit) - y

    return slope, fit, residuals

def fit_slopes(database):

    # gather concentrtions and slopes for x and y axes respectively
    x, y = [], []
    for sample in list(database):
        if any(char.isdigit() for char in str(sample)) == True:
            # linear fit and get slope
            sample_data = database[sample]
            xtmp, ytmp = sample_data[list(sample_data)[0]], sample_data[list(sample_data)[1]]
            if fit_choice_dict[sample] == 0:
                slope, linear_fit, residuals = fit_raw(xtmp, ytmp)
            elif fit_choice_dict[sample] == 2:
                slope, linear_fit, residuals = fit_raw_logarithmic(xtmp, ytmp)
            elif fit_choice_dict[sample] == 1:
                slope, linear_fit, residuals = fit_linear(xtmp, ytmp)
            y.append(np.abs(slope))
            x.append(float(re.findall(r"[-+]?\d*\.\d+|\d+", str(sample))[0]))

    # fit slopes to chosen kinetics model
    xfit = np.linspace(min(x), max(x), 100)
    if model_select.value == 'Michaelis-Menten':
        popt, pcov = curve_fit(mmfit, x, y)
        yfit = mmfit(xfit, *popt)
        model.title.text='Vmax = %.2e, Km = %.1e' % tuple(popt[::-1])
        model.title.text_font_size = '10pt'
    elif model_select.value == 'EC50/IC50':
        popt, pcov = curve_fit(icfit, x, y)
        yfit = icfit(xfit, *popt)
        model.title.text='EC50/IC50 = %.2e' % tuple([popt[2]])
        model.title.text_font_size = '10pt'
    data = pd.DataFrame({'x' : x, 'y' : y}).sort_values('x') # points for scatter
    fit = pd.DataFrame({'xfit' : xfit, 'yfit' : yfit}).sort_values('xfit')

    return data, fit

def refit_slopes(slope_data):
    x = slope_data.x
    y = slope_data.y

    # fit slopes to chosen kinetics model
    xfit = np.linspace(min(x), max(x), 100)
    if model_select.value == 'Michaelis-Menten':
        popt, pcov = curve_fit(mmfit, x, y)
        yfit = mmfit(xfit, *popt)
        model.title.text=r'Vmax = %.2e, Km = %.2e' % tuple(popt[::-1])
        model.title.text_font_size = '10pt'
    elif model_select.value == 'EC50/IC50':
        popt, pcov = curve_fit(icfit, x, y)
        yfit = icfit(xfit, *popt)
        model.title.text='EC50/IC50 = %.2e' % tuple([popt[2]])
        model.title.text_font_size = '10pt'
    data = pd.DataFrame({'x' : x, 'y' : y}).sort_values('x') # points for scatter
    fit = pd.DataFrame({'xfit' : xfit, 'yfit' : yfit}).sort_values('xfit')

    return data, fit

# bokeh methods

def update_tickers(attrname, old, new):
    update()

def update_range_slider(attrname, old, new):
    range_slider.start=df[x_sample_choice.value].values[0]
    range_slider.end=df[x_sample_choice.value].values[-1]
    range_slider.value=(df[x_sample_choice.value].values[0], df[x_sample_guess].values[-1])
    range_slider.step=step=df[x_sample_choice.value].values[1]-df[x_sample_guess].values[0]
    update()

def update_time(attrname, old, new):
    if range_slider.value[0] != float(start_time.value) or range_slider.value[1] != float(end_time.value):
        x = x_sample_choice.value
        y = sample_select.value
        data = get_data(x, y)

        #data = data.sort_values('x')
        selected_start = min(enumerate(list(data['x'].values)), key=lambda x: abs(x[1]-float(start_time.value)))[1]
        selected_end = min(enumerate(list(data['x'].values)), key=lambda x: abs(x[1]-float(end_time.value)))[1]
        range_slider.value=(selected_start, selected_end)

def update():
    # update raw plot
    x = x_sample_choice.value
    y = sample_select.value
    data = get_data(x, y)
    data = data.sort_values('x')

    selected_start = list(data['x'].values)[0]
    selected_end = list(data['x'].values)[-1]
    range_slider.value=(selected_start, selected_end)
    start_time.value=str(list(data['x'].values)[0])
    end_time.value=str(list(data['x'].values)[-1])

    source_raw.data = data[['x', 'y']].to_dict('list')
    if fit_choice_dict[sample_select.value] == 0:
        slope, raw_fit, residuals = fit_raw(data.x, data.y)
    elif fit_choice_dict[sample_select.value] == 2:
        slope, raw_fit, residuals = fit_raw_logarithmic(data.x, data.y)
    elif fit_choice_dict[sample_select.value] == 1:
        slope, raw_fit, residuals = fit_linear(data.x, data.y)
    raw_fit = raw_fit.sort_values('xfit')

    source_raw_line.data = raw_fit[['xfit', 'yfit']].to_dict('list')

    xr = np.array(raw_fit[['xfit']].values.T[0])
    yr = np.array(residuals)
    resi_df = pd.DataFrame(data={'xr' : xr, 'yr' : yr})
    source_resi.data = resi_df[['xr', 'yr']].to_dict('list')

    # update model plot & data table
    slope_data, slope_data_fit = fit_slopes(exp_database)

    if subtract_sample_choice.value != ' ':
        tmp_data = subtract_blank(slope_data, subtract_sample_choice.value)
        slope_data, slope_data_fit = refit_slopes(tmp_data)
    if transform.value != " ":
        transform_string = transform.value.replace('x', 'slope_data.y')
        tmp_data = slope_data
        tmp_data.y = eval(transform_string)
        slope_data, slope_data_fit = refit_slopes(tmp_data)

    slope_data = slope_data.sort_values('x')

    slope_data_fit = slope_data_fit.sort_values('xfit')

    source_data_table.data = slope_data[['x', 'y']].to_dict('list')

    source_model.data = slope_data[['x', 'y']].to_dict('list')

    source_model_line.data = slope_data_fit[['xfit', 'yfit']].to_dict('list')


def selection_change(attrname, old, new):
    fit_choice_dict[sample_select.value] = fit_choice.active

    # select range of current raw plot data
    x = x_sample_choice.value
    y = sample_select.value
    data = get_data(x, y)

    #data = data.sort_values('x')
    selected_start = min(range(len(list(data['x'].values))), key=lambda i: abs(list(data['x'].values)[i]-range_slider.value[0]))
    selected_end = min(range(len(list(data['x'].values))), key=lambda i: abs(list(data['x'].values)[i]-range_slider.value[1]))
    start_time.value=str(range_slider.value[0])
    end_time.value=str(range_slider.value[1])
    selected = list(range(selected_start, selected_end))
    if len(selected) > 0:
        # subset data according to selected range
        tmp_data = data.iloc[selected, :].reset_index(drop=True)
        tmp_data = tmp_data.sort_values('x')
    else:
        tmp_data = data
        tmp_data = tmp_data.sort_values('x')

        # re-fit straight line based on selection
    if fit_choice_dict[sample_select.value] == 0:
        slope, new_fit, residuals = fit_raw(tmp_data.x, tmp_data.y)
    elif fit_choice_dict[sample_select.value] == 2:
        slope, new_fit, residuals = fit_raw_logarithmic(tmp_data.x, tmp_data.y)
    elif fit_choice_dict[sample_select.value] == 1:
        slope, new_fit, residuals = fit_linear(tmp_data.x, tmp_data.y)
    new_fit = new_fit.sort_values('xfit')
    source_raw_line.data = new_fit[['xfit', 'yfit']].to_dict('list')

    xr = np.array(new_fit[['xfit']].values.T[0])
    yr = np.array(residuals)
    resi_df = pd.DataFrame(data={'xr' : xr, 'yr' : yr})
    source_resi.data = resi_df[['xr', 'yr']].to_dict('list')

    # re-fit model based on selection
    exp_database[y] = tmp_data

    slope_data, slope_data_fit = fit_slopes(exp_database)

    if subtract_sample_choice.value != ' ':
        tmp_data = subtract_blank(slope_data, subtract_sample_choice.value)
        slope_data, slope_data_fit = refit_slopes(tmp_data)
    slope_data = slope_data.sort_values('x')
    if transform.value != " ":
        transform_string = transform.value.replace('x', 'slope_data.y')
        tmp_data = slope_data
        tmp_data.y = eval(transform_string)
        slope_data, slope_data_fit = refit_slopes(tmp_data)
    source_data_table.data = slope_data[['x', 'y']].to_dict('list')
    source_model.data = slope_data[['x', 'y']].to_dict('list')
    source_model_line.data = slope_data_fit[['xfit', 'yfit']].to_dict('list')

def file_callback(attrname, old, new):
    global fit_choice
    fit_choice = RadioButtonGroup(labels=["Maximize Slope Magnitude", "Linear", "Logarithmic"], active=0, width=375)
    global fit_choice_dict
    fit_choice_dict = {}
    global file_source
    global filename
    filename=file_source.data['file_name']
    # read data file
    global output_filename
    output_filename = filename[0]+'-out.csv'
    raw_contents = file_source.data['file_contents'][0]
    prefix, b64_contents = raw_contents.split(",", 1)
    file_contents = base64.b64decode(b64_contents)
    file_contents = file_contents.decode("utf-8-sig")
    file_io = StringIO(file_contents)
    global df
    df = pd.read_csv(file_io)
    df = df.dropna(axis=1)
    sample_names = list(df)
    string_names = [str(i).strip(' ') for i in sample_names]
    df.columns = string_names
    global exp_database
    exp_database = {}
    for sample in list(df)[1:]:
        exp_database[sample] = df[[list(df)[0], sample]]
        fit_choice_dict[sample] = 0
    global model_select
    model_select = Select(title='Choose Model', value='Michaelis-Menten', options=['Michaelis-Menten', 'EC50/IC50'], width=350)
    global source_model
    source_model = ColumnDataSource(data=dict(x=[], y=[]))
    global source_model_line
    source_model_line = ColumnDataSource(data=dict(xfit=[], yfit=[]))

    tools_model = 'wheel_zoom,pan,reset,save, hover'

    hover = HoverTool(tooltips=[
        ("(x,y)", "($x, $y)")
    ])

    global model
    model = figure(title="Model Fit", x_axis_label="Concentration", y_axis_label="Rate",
                   plot_width=350, plot_height=300, tools=tools_model)
    model.yaxis.formatter = BasicTickFormatter(precision=2, use_scientific=True)

    global slope_data, slope_data_fit
    slope_data, slope_data_fit = fit_slopes(exp_database)
    slope_data = slope_data.sort_values('x')

    # set up widgets
    x_sample_guess = ' '
    blank_sample_guess = ' '
    for s in list(df):
        try:
            if 'ime' in s:
                x_sample_guess = s
        except:
            pass
        try:
            if 'lank' in s:
                blank_sample_guess = s
        except:
            pass

    file_source = ColumnDataSource({'file_contents':[], 'file_name':[]})
    file_source.on_change('data', file_callback)
    global upload_button
    upload_button = Button(label="Upload Local File", button_type="success", width=350)
    upload_button.callback = CustomJS(args=dict(file_source=file_source),
                               code=open(join(dirname(__file__), "upload.js")).read())
    global x_sample_choice
    x_sample_choice = Select(title='X Axis Column', value=x_sample_guess, options=string_names+[' '], width=350)
    global subtract_sample_choice
    subtract_sample_choice = Select(title='Select Blank Sample for Subtraction', value=blank_sample_guess, options=string_names+[' '], width=350)
    global sample_select
    sample_select = Select(title='Y Axis Sample', value=string_names[1], options=string_names+[' '], width=350)
    global transform
    transform = TextInput(value=" ", title="Enter Transform Equation", width=350)
    global offset_time
    offset_time = TextInput(value=" ", title="Enter Time Between Mixing and First Read", width=350)

    global source_data_table
    source_data_table = ColumnDataSource(slope_data)
    columns = [
            TableColumn(field="x", title="Concentration"),
            TableColumn(field="y", title="Slope (Initial Rate)"),
        ]
    global data_table
    data_table = DataTable(source=source_data_table, columns=columns, width=350, height=450, selectable=True, editable=True)
    global download_button
    download_button = Button(label="Download Table to CSV", button_type="primary", width=350)

    try:
        output_filename = file_source.data['file_name']+'-out.csv'
    except:
        output_filename = 'output.csv'
    download_button.callback = CustomJS(args=dict(source=source_data_table, file_name=output_filename),
                               code=open(join(dirname(__file__), "download.js")).read())
    global copy_button
    copy_button = Button(label="Copy Table to Clipboard", button_type="primary", width=350)
    copy_button.callback = CustomJS(args=dict(source=source_data_table),
                               code=open(join(dirname(__file__), "copy.js")).read())

    # configure plots
    global source_resi
    source_resi = ColumnDataSource(data=dict(xr=[], yr=[]))
    global resi
    resi = figure(title="Progress Curve Fit Residuals", x_axis_label="Time", y_axis_label="Residual",
             plot_width=700, plot_height=200, tools='wheel_zoom,pan,reset')
    resi.yaxis.formatter = BasicTickFormatter(precision=2, use_scientific=True)
    resi.circle('xr', 'yr', size=5, source=source_resi, color='grey', alpha=0.6)

    global source_raw
    source_raw = ColumnDataSource(data=dict(x=[], y=[]))
    global source_raw_line
    source_raw_line = ColumnDataSource(data=dict(xfit=[], yfit=[]))
    tools_raw = 'wheel_zoom,pan,reset,save'
    global raw
    raw = figure(title="Progress Curve Fit", x_axis_label="Time", y_axis_label="Signal",
                 plot_width=350, plot_height=300, tools=tools_raw)

    raw.circle('x', 'y', size=2, source=source_raw, color='gray',
                selection_color="black", alpha=0.6, nonselection_alpha=0.2, selection_alpha=0.6)
    raw.line('xfit', 'yfit', source=source_raw_line, color='red')

    source_model.data = slope_data[['x', 'y']].to_dict('list')
    source_model_line.data = slope_data_fit[['xfit', 'yfit']].to_dict('list')

    model.circle('x', 'y', size=8, source=source_model, color='grey', alpha=0.6)
    model.line('xfit', 'yfit', source=source_model_line, line_width=3, color='black', alpha=0.4)

    # update plots according to raw data selection
    global range_slider
    range_slider = RangeSlider(start=df[x_sample_guess].values[0], end=df[x_sample_guess].values[-1],
                value=(df[x_sample_guess].values[0], df[x_sample_guess].values[-1]),
                step=df[x_sample_guess].values[1]-df[x_sample_guess].values[0],
                title='X-Axis Range', width=650)
    range_slider.on_change('value', selection_change)
    start_time = TextInput(value=str(df[x_sample_guess].values[0]), title="Start Time")
    end_time = TextInput(value=str(df[x_sample_guess].values[-1]), title='End Time')

    # update plots based on ticker selections
    fit_choice.on_change('active', selection_change)
    x_sample_choice.on_change('value', update_tickers)
    subtract_sample_choice.on_change('value', update_tickers)
    sample_select.on_change('value', update_tickers)
    model_select.on_change('value', update_tickers)
    transform.on_change('value', update_tickers)
    offset_time.on_change('value', update_tickers)
    start_time.on_change('value', update_time)
    end_time.on_change('value', update_time)

    # document formatting
    desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1400)
    widgets = widgetbox(model_select, sample_select, x_sample_choice,
                        subtract_sample_choice, transform, offset_time)
    table = widgetbox(data_table)

    main_row = row(column(upload_button, fit_choice, widgets),
                    column(row(raw, model), resi, range_slider, row(start_time, end_time)),
                    column(download_button, copy_button, table))

    sizing_mode = 'scale_width'
    l = layout([
        [desc],
        [main_row]
    ], sizing_mode=sizing_mode)
    curdoc().clear()
    curdoc().add_root(l)
    curdoc().title = "Kinetics"

    update()

# read data file
data_file = join(dirname(__file__), 'test.csv')
df = pd.read_csv(data_file)
sample_names = list(df)
string_names = [str(i) for i in sample_names]
df.columns = string_names
exp_database = {}
for sample in list(df)[1:]:
    exp_database[sample] = df[[list(df)[0], sample]]
    fit_choice_dict[sample] = 0

model_select = Select(title='Choose Model', value='Michaelis-Menten', options=['Michaelis-Menten', 'EC50/IC50'], width=350)

source_model = ColumnDataSource(data=dict(x=[], y=[]))

source_model_line = ColumnDataSource(data=dict(xfit=[], yfit=[]))

tools_model = 'wheel_zoom,pan,reset,save, hover'

hover = HoverTool(tooltips=[
    ("(x,y)", "($x, $y)")
])

model = figure(title="Spline Model Fit", x_axis_label="Concentration", y_axis_label="Rate",
               plot_width=350, plot_height=300, tools=tools_model)
model.yaxis.formatter = BasicTickFormatter(precision=2, use_scientific=True)

slope_data, slope_data_fit = fit_slopes(exp_database)

slope_data = slope_data.sort_values('x')

# set up widgets
x_sample_guess = ' '
blank_sample_guess = ' '
for s in list(df):
    try:
        if 'ime' in s:
            x_sample_guess = s
    except:
        pass
    try:
        if 'lank' in s:
            blank_sample_guess = s
    except:
        pass

file_source = ColumnDataSource({'file_contents':[], 'file_name':[]})
file_source.on_change('data', file_callback)
upload_button = Button(label="Upload Local File", button_type="success", width=350)
upload_button.callback = CustomJS(args=dict(file_source=file_source),
                           code=open(join(dirname(__file__), "upload.js")).read())

x_sample_choice = Select(title='X Axis Column', value=x_sample_guess, options=string_names+[' '], width=350)
subtract_sample_choice = Select(title='Select Blank Sample for Subtraction', value=blank_sample_guess, options=string_names+[' '], width=350)
sample_select = Select(title='Y Axis Sample', value=string_names[1], options=string_names+[' '], width=350)
transform = TextInput(value=" ", title="Enter Transform Equation", width=350)
offset_time = TextInput(value=" ", title="Enter Time Between Mixing and First Read", width=350)

source_data_table = ColumnDataSource(slope_data)
columns = [
        TableColumn(field="x", title="Concentration"),
        TableColumn(field="y", title="Slope (Initial Rate)"),
    ]

data_table = DataTable(source=source_data_table, columns=columns, width=350, height=450, selectable=True, editable=True)

download_button = Button(label="Download Table to CSV", button_type="primary", width=350)
try:
    output_filename = file_source.data['file_name']+'-out.csv'
except:
    output_filename = 'output.csv'
download_button.callback = CustomJS(args=dict(source=source_data_table, file_name=output_filename),
                           code=open(join(dirname(__file__), "download.js")).read())

copy_button = Button(label="Copy Table to Clipboard", button_type="primary", width=350)
copy_button.callback = CustomJS(args=dict(source=source_data_table),
                           code=open(join(dirname(__file__), "copy.js")).read())

# configure plots
source_raw = ColumnDataSource(data=dict(x=[], y=[]))
source_resi = ColumnDataSource(data=dict(xr=[], yr=[]))

source_raw_line = ColumnDataSource(data=dict(xfit=[], yfit=[]))

tools_raw = 'wheel_zoom,pan,reset,save'

raw = figure(title="Progress Curve Fit", x_axis_label="Time", y_axis_label="Signal",
             plot_width=350, plot_height=300, tools=tools_raw)
resi = figure(title="Progress Curve Fit Residuals", x_axis_label="Time", y_axis_label="Residual",
             plot_width=700, plot_height=200, tools='wheel_zoom,pan,reset')
resi.yaxis.formatter = BasicTickFormatter(precision=2, use_scientific=True)

raw.circle('x', 'y', size=2, source=source_raw, color='gray',
            selection_color="black", alpha=0.6, nonselection_alpha=0.2, selection_alpha=0.6)
resi.circle('xr', 'yr', size=5, source=source_resi, color='grey', alpha=0.6)

raw.line('xfit', 'yfit', source=source_raw_line, color='red')

source_model.data = slope_data[['x', 'y']].to_dict('list')

source_model_line.data = slope_data_fit[['xfit', 'yfit']].to_dict('list')

model.circle('x', 'y', size=8, source=source_model, color='grey', alpha=0.6)

model.line('xfit', 'yfit', source=source_model_line, line_width=3, color='black', alpha=0.4)

# update plots according to raw data selection
range_slider = RangeSlider(start=df[x_sample_guess].values[0], end=df[x_sample_guess].values[-1],
                value=(df[x_sample_guess].values[0], df[x_sample_guess].values[-1]),
                step=df[x_sample_guess].values[1]-df[x_sample_guess].values[0],
                title='X-Axis Range', width=650)

range_slider.on_change('value', selection_change)

start_time = TextInput(value=str(df[x_sample_guess].values[0]), title="Start Time")
end_time = TextInput(value=str(df[x_sample_guess].values[-1]), title='End Time')

# update plots based on ticker selections
fit_choice.on_change('active', selection_change)
x_sample_choice.on_change('value', update_range_slider)
subtract_sample_choice.on_change('value', update_tickers)
sample_select.on_change('value', update_tickers)
model_select.on_change('value', update_tickers)
transform.on_change('value', update_tickers)
offset_time.on_change('value', update_tickers)
start_time.on_change('value', update_time)
end_time.on_change('value', update_time)

# document formatting
desc = Div(text=open(join(dirname(__file__), "description.html")).read(), width=1400)
widgets = widgetbox(model_select, sample_select, x_sample_choice,
                    subtract_sample_choice, transform, offset_time)
table = widgetbox(data_table)

main_row = row(column(upload_button, fit_choice, widgets),
                column(row(raw, model), resi, range_slider, row(start_time, end_time)),
                column(download_button, copy_button, table))

sizing_mode = 'scale_width'
l = layout([
    [desc],
    [main_row]
], sizing_mode=sizing_mode)

update()

curdoc().add_root(l)
curdoc().title = "Kinetics"
