# -*- coding: utf-8 -*-
"""
@author: cutter

Plot echogram of Sv for each channel of a raw EK60 file. 

Modified from towler's examples. 

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from echolab2.instruments import EK60
# for processed data
from echolab2.processing import processed_data, line
# for echogram fig
from matplotlib.pyplot import figure, show 
from echolab2.plotting.matplotlib import echogram


#  User: specify infile path and name
#infp = 'C:\\dev\\pyEcholab\\data\\'
infp = '..\\..\\data\\'
infn = 'ek-701061512.raw'
raw_filename = infp + infn


#  specify the ping number to plot power data profiles
ping_number = 11


#  create an instance of the EK60 object
ek60 = EK60.EK60()


#  read the .raw data file
print('Reading raw file %s' % (raw_filename))
ek60.read_raw(raw_filename)

#  print some info about the state of the EK60 object
print(ek60)

#  colormap for plotting
color = iter(cm.rainbow(np.linspace(0, 1, len(ek60.channel_ids))))

#  create a figure
fig = plt.figure(figsize=(7,7))

#  info ~ num chans
print("N channels = %d" % len(ek60.channel_ids))

#  plot power from the specified ping from all channels
for channel_id in ek60.channel_ids:

    #  get a reference to the raw_data for this channel
    raw_data = ek60.get_raw_data(channel_id=channel_id)

    print("raw_data.power shape = %s" % str(raw_data.power[:,:].shape))

    #  get a color for this channel's power plot
    c = next(color)

    #  create a sample index so we can pass it as the Y axis in the figure
    #  so we plot the data up and down
    n_samples = len(raw_data.power[ping_number,:])
    yaxis = np.arange(n_samples)

    #  plot power
    plt.plot(raw_data.power[ping_number,:], yaxis, color=c, label=channel_id)
    
    
    # plot raw power in matrix / image form 
    #hfigimpow = plt.imshow( raw_data.power )


    # convert from raw power to Sv 
    Sv = raw_data.get_Sv()
    
    
    # plot 
    #hfigeg = plt.imshow( Sv ) #NO -- ERROR "image data cannot be converted to float."
    # ...Sv is not just a matrix
    
    
    ## Echogram plot
    fig_1 = figure()
    eg = echogram.echogram(fig_1, Sv )
    titstr = 'Sv Echogram, %s' % channel_id
    eg.axes.set_title( titstr )
    
    
    show()
    #  end iteration through channel ids

# end