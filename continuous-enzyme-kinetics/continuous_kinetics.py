import re
from decimal import Decimal
import math
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from uncertainties import ufloat

def linear(x, m, b):
    '''
    straight line
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

def icfit(x, bottom, top, slope, IcFifty):
    '''
    IC50 equation
    '''
    return bottom + (top-bottom)/(slope+(x/IcFifty))

def spline_fit(x, y):
    
    x, y = x.values, y.values
    spline = UnivariateSpline(x, y)(x)
    derivative = np.abs(np.diff(spline)/np.diff(x))
    threshold = 0.7*(np.max(derivative) - np.min(derivative)) + np.min(derivative)
    indices = np.where(derivative > threshold)[0]
    while len(indices) < 4:
        threshold = threshold*0.9
        indices = np.where(derivative > threshold)[0]
    xi, yi = x[indices], y[indices]
    df = pd.DataFrame(data={'x' : xi, 'y' : yi}).sort_values('x').dropna()
    xi, yi = df.x, df.y
    popt, pcov = curve_fit(linear, xi, yi)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(spline))
    yfit = linear(xfit, *popt)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(popt[0]), 'error' : perr[0],
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

def linear_fit(x, y):

    popt, pcov = curve_fit(linear, x, y)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(x))
    yfit = linear(xfit, *popt)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(popt[0]), 'error' : perr[0],
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

def logarithmic_fit(x, y):
    
    popt, pcov = curve_fit(logarithmic, x, y, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    xfit = np.linspace(np.min(x), np.max(x), len(x))
    yfit = logarithmic(xfit, *popt)
    yerr = logarithmic(xfit, *perr)
    fit_dict = { 'x' : x, 'y' : y,
                    'rate' : np.abs(np.array(np.diff(yfit)/np.diff(xfit))[0]),
                    'error' : np.abs(np.array(np.diff(yerr)/np.diff(xfit))[0]),
                    'xfit' : xfit,
                    'yfit' : yfit, 'resi' : np.array(yfit) - y }

    return fit_dict

class progress_curve(object):

    def __init__(self, dataframe):

        self.data = dataframe

    def spline(self, start, end):

        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))]
        x, y = df[df.columns[0]], df[df.columns[1]]
        spline = spline_fit(x, y)
        self.spline = spline
        
        return self.spline

    def linear(self, start, end):

        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))]
        x, y = df[df.columns[0]], df[df.columns[1]]
        linear = linear_fit(x, y)
        self.linear = linear
        
        return self.linear
    
    def logarithmic(self, start, end, offset):
        
        df = self.data[(self.data[self.data.columns[0]] >= float(start)) &
                       (self.data[self.data.columns[0]] <= float(end))]
        x, y = df[df.columns[0]], df[df.columns[1]]
        try:
            x = x + offset
        except:
            pass
        logarithmic = logarithmic_fit(x, y)
        self.logarithmic = logarithmic

        return self.logarithmic

class kinetic_model(object):

    def __init__(self, dictionary):
        self.dict = dictionary

    def model(self, subtract, transform, threshold):
        
        result = {}
        df = pd.DataFrame()
        for s in self.dict:
            if type(self.dict[s]) == progress_curve:
                x = float(re.findall(r"[-+]?\d*\.\d+|\d+", str(s))[0])
                if self.dict[s+'_fit'] == 0:
                        sdf = self.dict[s].spline
                elif self.dict[s+'_fit'] == 1:
                        sdf = self.dict[s].linear
                else:
                        sdf = self.dict[s].logarithmic
                df.at[s, 'rate'] = sdf['rate']
                df.at[s, 'error'] = sdf['error']
                df.at[s, 'x'] = x
        df = df.sort_values(['x'])
        uRates = [ufloat(ur, ue) for ur, ue in zip(df['rate'], df['error'])]
        df['uRates'] = uRates
        try:
            df['uRates'] = df['uRates'] - df.loc[subtract]['uRates']
            df['rate'] = [ur.nominal_value for ur in df['uRates']]
            df['error'] = [ur.std_dev for ur in df['uRates']]
            df = df[df.index != subtract]:wq
        except:
            pass
        try:
            x = np.array([ufloat(r, e) for r, e in zip(df['rate'], df['error'])])
            x = eval(transform)
            df['rate'] = [xi.nominal_value for xi in x]
            df['error'] = [xi.std_dev for xi in x]
        except:
            pass
        n, x, y, e = df.index.values, df['x'].values, df['rate'].values, df['error'].values
        result['n'] = n
        result['y'] = y
        result['yt'] = ['%.2E' % Decimal(str(yi)) for yi in y]
        result['e'] = ['%.2E' % Decimal(str(ei)) for ei in e]
        result['l'] = y - e
        result['u'] = y + e
        xfit = np.linspace(np.min(x), np.max(x), 100)
        if self.dict['model'] == 'Michaelis-Menten':
            result['x'] = x
            result['xfit'] = xfit
            try:
                popt_mm, pcov_mm = curve_fit(mmfit, x, y, sigma=e, absolute_sigma=True)
            except:
                popt_mm, pcov_mm = curve_fit(mmfit, x, y)
            perr_mm = np.sqrt(np.diag(pcov_mm))
            ymm = mmfit(xfit, *popt_mm)
            result['yfit'] = ymm
            result['Km'] = tuple(['%.2E' % Decimal(str(popt_mm[0])), 
                                    '%.2E' % Decimal(str(perr_mm[0]))])
            result['Vmax'] = tuple(['%.2E' % Decimal(str(popt_mm[1])), 
                                    '%.2E' % Decimal(str(perr_mm[1]))])
            result['c'] = ['grey']*len(result['x'])
            result['ct'] = ['white']*len(result['x'])
        elif self.dict['model'] == 'EC50/IC50':
            result['x'] = x
            result['xfit'] = xfit
            try:
                popt_ic, pcov_ic = curve_fit(icfit, x, y, sigma=e, absolute_sigma=True)
            except:
                popt_ic, pcov_ic = curve_fit(icfit, x, y)
            perr_ic = np.sqrt(np.abs(np.diag(pcov_ic)))
            yic = icfit(xfit, *popt_ic)
            result['yfit'] = yic
            if popt_ic[0] < popt_ic[1]:
                result['Bottom'] = np.array(['%.2E' % Decimal(str(popt_ic[0])),
                                             '%.2E' % Decimal(str(perr_ic[0]))])
                result['Top'] = np.array(['%.2E' % Decimal(str(popt_ic[1])),
                                          '%.2E' % Decimal(str(perr_ic[1]))])
            else:
                result['Top'] = np.array(['%.2E' % Decimal(str(popt_ic[0])),
                                          '%.2E' % Decimal(str(perr_ic[0]))])
                result['Bottom'] = np.array(['%.2E' % Decimal(str(popt_ic[1])),
                                             '%.2E' % Decimal(str(perr_ic[1]))])
            result['Slope'] = np.array(['%.2E' % Decimal(str(popt_ic[2])), 
                                        '%.2E' % Decimal(str(perr_ic[2]))])
            result['IcFifty'] = np.array(['%.2E' % Decimal(str(popt_ic[3])),
                                          '%.2E' % Decimal(str(perr_ic[3]))])
            result['c'] = ['grey']*len(result['x'])
            result['ct'] = ['white']*len(result['x'])
        else:
            result['x'] = np.linspace(1, len(n), len(n))
            result['xfit'] = np.linspace(1, len(n), len(n))
            result['yfit'] = np.repeat(np.mean(result['y']), len(n))
            std = np.std(result['y'])
            avg = np.mean(result['y'])
            color, colort = [], []
            for r in result['y']:
                if r >= avg + std*threshold:
                    color.append('red')
                    colort.append('#EC7063')
                elif r <= avg - std*threshold:
                    color.append('blue')
                    colort.append('#5DADE2')
                else:
                    color.append('grey')
                    colort.append('white')
            result['c'] = color
            result['ct'] = colort
        self.result = result
        
        return self.result
