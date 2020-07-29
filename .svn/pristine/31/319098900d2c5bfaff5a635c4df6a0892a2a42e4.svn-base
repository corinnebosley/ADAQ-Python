#!/usr/bin/env python
'''
Top level script for plotting wind roses.

.. note::
   windrose_plot.py imports windrose.py from adaqsandpit

In the sample data directory two sets of data which can be
used with this script can be found

  * StAthan_HourlyWind_070714to130714.csv - a csv file containing met data
  * Met_Data_traj_C1.txt - output from a NAME run containing met data

Example of use::

   windrose_plot.py StAthan_HourlyWind_070714to130714.csv

There are two other options which may be passed to windrose_plot.py

  :-n: Used to plot percentages instead of counts
  :-t: Used to add a title (e.g. -t "St Athan")

'''
from six.moves.builtins import range
from six.moves.builtins import object
from windrose import WindroseAxes
from matplotlib import pyplot as plt
import numpy.random as random
import csv, sys, re, argparse

class CCUObservationalData(object):
    """ Interprets the header of the CCU Observational CSV file and returns
    a tuple of wind speeds and direction. Units of the wind speed and
    direction can also be queried.  """

    def __init__(self, csvfile):
        csvfile = file(csvfile, 'rb')
        self._reader = csv.reader(csvfile)
        self._interpreted = False
        self._data = []

    def _getDateColumn(self, row):
        c = 0
        for i in range(len(row)):
            m = self._dateMatch(row[i])
            if m:
                return c
            c += 1
        raise AssertionError('Incorrect Format')

    def _getTimeColumn(self, row):
        c = 0
        for i in range(len(row)):
            m = self._timeMatch(row[i])
            if m:
                return c
            c += 1
        raise AssertionError('Incorrect Format')


    def _getWindSpeedColumn(self, row):
        c = 0
        for i in range(len(row)):
            m = self._speedMatch(row[i])
            if m:
                return c
            c += 1
        raise AssertionError('Incorrect Format')

    def _getWindSpeedText(self, row):
        for i in range(len(row)):
            m = self._speedMatch(row[i])
            if m:
                return m.group(0)
        raise AssertionError('Incorrect Format')

    def _getWindDirectionColumn(self, row):
        c = 0
        for i in range(len(row)):
            m = self._directionMatch(row[i])
            if m:
                return c
            c += 1
        raise AssertionError('Incorrect Format')



    def _interpretHeader(self):
        row = next(self._reader)
        self._dateColumn      = self._getDateColumn(row)
        self._timeColumn      = self._getTimeColumn(row)
        self._speedColumn     = self._getWindSpeedColumn(row)
        self._directionColumn = self._getWindDirectionColumn(row)
        self._speedText       = self._getWindSpeedText(row)
        self._interpreted     = True

    def _dateMatch(self, item):
        return re.match('(date$)', item, re.I)

    def _timeMatch(self, item):
        return re.match('(time$)', item, re.I)

    def _directionMatch(self, item):
        return re.match('(.*direction.*$)', item, re.I)

    def _speedMatch(self, item):
        return re.match('(.*speed.*$)', item, re.I)

    def _readRows(self):
        for i in self._reader:
            self._data.append(i)

    def getWindSpeed(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [float(x[self._speedColumn]) for x in self._data]

    def getWindDirection(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [float(x[self._directionColumn]) for x in self._data]

    def getDate(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [x[self._dateColumn] for x in self._data]

    def getTime(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [x[self._timeColumn] for x in self._data]


    def getWindSpeedUnits(self):
        if not self._interpreted:
            self._interpretHeader()
            self._interpreted = True

        try:
            units = re.match('(.*speed.*\((\S+){1}\).*$)', self._speedText, re.I).group(2)

        except:
            raise AssertionError('Incorrect Format (wind speed)')

        return units

    def getWindDirectionUnits(self):
        return 'deg'


class NAMENWPData(CCUObservationalData):
    """ Interprets the header of the CSV file created from NWP data (NAME)
    and returns a tuple of wind speeds and direction. Units of the wind speed
    and direction can also be queried.  """

    def _interpretHeader(self):
        for i in self._reader:
            if i[0] == ' Fields:':
                break

        row = next(self._reader)
        row = next(self._reader)
        row = next(self._reader)
        self._dateColumn      = 0
        self._timeColumn      = 0
        self._speedColumn     = self._getWindSpeedColumn(row)
        self._directionColumn = self._getWindDirectionColumn(row)
        self._interpreted     = True

        row = next(self._reader)
        row = next(self._reader)

        self._speedUnits = row[self._speedColumn].lstrip()
        self._directionUnits = row[self._directionColumn].lstrip()

        endHeader = False
        for i in self._reader:
            if i[0].lstrip() == 'T':
                endHeader = True
                break

        if not endHeader:
            raise AssertionError('Incorrect Format')

    def getWindSpeedUnits(self):
        if not self._interpreted:
            self._interpretHeader()
            self._interpreted = True

        return self._speedUnits

    def getWindDirectionUnits(self): return self._directionUnits

    def getDate(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [x[self._dateColumn].split()[0] for x in self._data]

    def getTime(self):
        if not self._interpreted:
            self._interpretHeader()

        if not self._data:
            self._readRows()

        return [x[self._timeColumn].split()[1] for x in self._data]


def new_axes(title):
    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = WindroseAxes(fig, rect, axisbg='w')
    fig.add_axes(ax)
    fig.suptitle(title)
    return ax

def set_legend(ax):
    # A borderaxespad of -3 places the legend outside of the rose and away
    # from the text for the compass direction in the bottom left corner of
    # the image.
    l = ax.legend(ws_units, borderaxespad=-3)
    plt.setp(l.get_texts(), fontsize=8)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce windrose plot from csv data obtained from either \
             (1) Observations from CCU. (2) NWP data extracted using NAME.')

    parser.add_argument('csvfile', action="store", help="Comma separated file with header")
    parser.add_argument('-t', '--title', action="store", dest="title", help="Centered title for plot")
    parser.add_argument('-n', "--normed", help="Write percentages  rather than the count of winds in the plot",
                        action="store_true")
    args = parser.parse_args()

    if args.title:
        title = args.title

    else:
        title = 'Title not specified'

    # Need to determine which type of CSV file we have.
    data = CCUObservationalData(args.csvfile)
    try:
        data.getWindSpeedUnits()

    except AssertionError:
        data = NAMENWPData(args.csvfile)

    # Read wind speed and direction and get units.
    (ws, wd) = (data.getWindSpeed(), data.getWindDirection())
    ws_units = data.getWindSpeedUnits()

    date = data.getDate()
    time = data.getTime()

    title = title + '  (' + time[0] + 'Z ' + date[0] + ' - ' + time[-1] + 'Z ' + date[-1] + ')'

    # Convert knots to metres per sec
    if re.match('kn|knots', ws_units, re.I):
        ws = [x*0.514444444444444 for x in ws]
        ws_units = 'm/s'

    bins = [0., 1., 3., 5.]


    # windrose like a stacked histogram with normed (displayed in percent) results
    ax = new_axes(title)
    ax.bar(wd, ws, bins=bins, colors=['#ffffff', '#7fff00', '#ffff00', '#ff4500'],
           normed=args.normed, opening=0.5, edgecolor='black')
    set_legend(ax)

    plt.show()





