# import required modules
import os
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
import string
import re

def data_processing(file, platemap_plot=False,endremove=2,
	start_row = 2, start_col = 'B'):
	# read data
	df_in = pd.read_table(file, skiprows=[0,1], 
		index_col = 1, header=0)
	# detect shape of the dataframe
	plate_cols = df_in.columns.tolist()
	# if the input data has been broken down to individual images, exit the script
	if 'Image' in plate_cols[0]:
		print('Input dataframe not supported')
		exit()
	p = re.compile(r'[A-Z]+\d+')
	columns_valid = [s for s in plate_cols if p.match(s)]
	shape_x = int(re.findall(r'\d+',columns_valid[-1])[0]) - int(re.findall(r'\d+',columns_valid[0])[0]) + 1
	shape_y = ord(re.findall(r'[A-Z]',columns_valid[-1])[0]) - ord(re.findall(r'[A-Z]',columns_valid[0])[0]) + 1
	print('Output: DIP-rate dataframe of {} columns and {} rows.'.format(shape_x, shape_y))
	# Normalize data to log2 cell doublings
	df_in = df_in[columns_valid]
	df_in = np.log2(df_in.divide(df_in.loc[0]))

	# Find DIP for each column, output results to dataframe
	columns = ['DIP','intercept','r_value','p_value']
	df_in_DIP = pd.DataFrame(index=df_in.columns, columns=columns)
	if platemap_plot is True:
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		import matplotlib as mpl
		figure = plt.figure(figsize=(shape_x*0.6, shape_y*0.6))
		gs = gridspec.GridSpec(shape_y, shape_x)
		
	# For every column of data, find DIP rate
	def find_DIP(x,y):
		# find and remove outliers
		slope, intercept, r_value, p_value, std_err = linregress(x,y)
		remove_list = []
		for i in range(0,len(x)):
			y_predict = slope*x[i]+intercept
			if (y_predict-y[i]) > 0.33*(max(y)-min(y)):
				remove_list.append(i)
		for index in sorted(remove_list, reverse=True):
			del x[index]
			del y[index]
		# final dip fit (later time point biased)
		list_p = []
		for i in range(0,len(x)-endremove):
			temp_x = x[i:]
			slope, intercept, r_value, p_value, std_err = linregress(temp_x,y[i:])
			list_p.append(abs(r_value))
		best = list_p.index(max(list_p))
		slope, intercept, r_value, p_value, std_err = linregress(x[best:],y[best:])
		return slope,intercept,r_value,p_value
		
	for column in df_in:
		a = ord(column[0])-66
		b = int(column[1:])-2
		DIP_list = find_DIP(df_in.index.tolist(), df_in[column].tolist())
		if platemap_plot is True:
			axes = plt.subplot(gs[a,b])
			axes.scatter(df_in.index.tolist(), df_in[column].tolist(),
				s=1.5,marker='x')
			axes.plot(
				(0,df_in.index.tolist()[-1]),
				(DIP_list[1],df_in.index.tolist()[-1]*DIP_list[0]+DIP_list[1]),
				color='r',alpha=0.8
				)
			axes.set_ylim(-2,8)
			axes.set_xticks([], [])
			axes.set_yticks([], [])
		DIP_dict = pd.Series(dict(zip(columns, DIP_list)))
		# print(DIP_dict)
		df_in_DIP.loc[column] = DIP_dict
	# Out DIP for each well into a heatmap
	plate_df_out = pd.DataFrame(df_in_DIP.DIP.values.reshape((shape_x,shape_y)))
	plate_df_out = plate_df_out.T
	plate_df_out.index = list(string.ascii_uppercase[ord(start_col)-65:shape_y+ord(start_col)-65])
	plate_df_out.columns = range(start_row,shape_x+start_row)
	# plate_df_out.to_csv('temp.csv', header=True, index=True)
	# plate_df_out = pd.read_csv('temp.csv',header=0,index_col=0)
	# os.remove('temp.csv')
	plate_df_out = plate_df_out*24
	return plate_df_out

