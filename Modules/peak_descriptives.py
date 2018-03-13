import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


#######################################################
def calculate_distance(peaks):
	'''
	Calculates the distance between each peak if they are on the same chromosome 
	Args
		1) peaks: Pandas dataframe with at least columns Chr, Start, End
	Returns:
		1) peaks: Same dataframe, but with an additional column 'Distance to Next Peak'. 
				  Distance is from the start of the peak+1 entry - current peaks end
	'''
	
	peaks.sort_values(['Chr','Start','End'],inplace=True)
	peaks['Distance to Next Peak'] = None
	for ind,val in	peaks.iterrows():
	    if ind ==	peaks.index[-1]:
	        break
	    if	peaks.loc[ind,'Chr'] ==	peaks.loc[ind+1,'Chr']:
	       	peaks.loc[ind,'Distance to Next Peak'] = peaks.loc[ind+1,'Start'] - val['End']
	       	
    
	return peaks
#######################################################


#######################################################
def calculate_divergent_transcripts(peaks,divergent_width,is_plot=False,f_name = 'divergent_transcripts'):
	''' See if peaks are near each other and on the opposite strand
	'''
	peaks.sort_values(['Chr','Start','End'],inplace=True)
	peaks.reset_index(inplace=True,drop=True)
	peaks['has Divergent'] = False
	peaks['Associated Divergent'] = None
	peaks['Divergent Distance'] = None

	for ind,val in peaks.iterrows():
	    if ind == peaks.index[-1]:
	        break
	    if peaks.loc[ind,'Chr'] == peaks.loc[ind+1,'Chr']:
	        if peaks.loc[ind,'Strand'] == '-' and peaks.loc[ind+1,'Strand'] == '+':
	            val2 = peaks.loc[ind+1]
	            #See if center of peaks are less than divergent_width distance
	            div_dist = (val2['Start'] + (val2['End']-val2['Start'])/2) - (val['Start'] + (val['End']-val['Start'])/2)
	            if  div_dist < divergent_width:
	                peaks.loc[ind,'has Divergent'] = True
	                peaks.loc[ind,'Associated Divergent'] = ind+1
	                peaks.loc[ind,'Divergent Distance'] = div_dist
	                peaks.loc[ind+1,'has Divergent'] = True
	                peaks.loc[ind+1,'Associated Divergent'] = ind
	                peaks.loc[ind+1,'Divergent Distance'] = div_dist

	if is_plot:
		peaks[peaks['has Divergent']]['Divergent Distance'].hist()
		num_div = np.sum(peaks['has Divergent'])
		total_peaks = len(peaks)
		title = '%s divergent peaks out of %s total peaks (divergent counts both forward and reverse)' % (num_div,total_peaks)
		plt.title(title)
		plt.savefig('./Figures' + f_name + '.png')
	return peaks
#######################################################


#######################################################
def group_by_adjacent_peaks(peaks,multiple_peaks_width):
	peaks.sort_values(['Chr','Start','End'],inplace=True)
	peaks['Group TSS'] = None
	curr_group = 0
	curr_chr = peaks.iloc[0]['Chr']
	curr_end = peaks.iloc[0]['Start']
	for ind,val in peaks.iterrows():
	    if not (val['Chr'] == curr_chr and val['Start'] - curr_end < multiple_peaks_width):
	        curr_group += 1
	        curr_chr = val['Chr']
	    curr_end = val['End']
	    peaks.loc[ind,'Group TSS'] = curr_group
	return peaks
#######################################################


#######################################################
def stats_all_in_one(peaks,multiple_peaks_width,divergent_width,f_name = '', is_save = False,is_plot = False):
	peaks.sort_values(['Chr','Start','End'],inplace=True)
	peaks.reset_index(inplace=True)

	# Initialize additional columns
	peaks['Distance to Next Peak'] = None
	peaks['has Divergent'] = False
	peaks['Associated Divergent'] = None
	peaks['Divergent Distance'] = None
	peaks['Group TSS'] = None

	curr_group = 0
	curr_chr = peaks.iloc[0]['Chr']
	curr_end = peaks.iloc[0]['Start']
	for ind,val in peaks.iterrows():
	    
		## Group if they are nearby
	    if not (val['Chr'] == curr_chr and val['Start'] - curr_end < multiple_peaks_width):
	        curr_group += 1
	        curr_chr = val['Chr']
	    curr_end = val['End']
	    peaks.loc[ind,'Group TSS'] = curr_group
	    
	    if not ind ==	peaks.index[-1]:	        
		    ## Calculate distance between peaks
		    if	peaks.loc[ind,'Chr'] ==	peaks.loc[ind+1,'Chr']:
		       	peaks.loc[ind,'Distance to Next Peak'] = peaks.loc[ind+1,'Start'] - val['End']

		    ## Determine if adjacent peaks are divergent
		    if peaks.loc[ind,'Chr'] == peaks.loc[ind+1,'Chr']:
		        if peaks.loc[ind,'Strand'] == '-' and peaks.loc[ind+1,'Strand'] == '+':
		            val2 = peaks.loc[ind+1]
		            #See if center of peaks are less than divergent_width distance
		            div_dist = (val2['Start'] + (val2['End']-val2['Start'])/2) - (val['Start'] + (val['End']-val['Start'])/2)
		            if  div_dist < divergent_width:
		                peaks.loc[ind,'has Divergent'] = True
		                peaks.loc[ind,'Associated Divergent'] = ind+1
		                peaks.loc[ind,'Divergent Distance'] = div_dist
		                peaks.loc[ind+1,'has Divergent'] = True
		                peaks.loc[ind+1,'Associated Divergent'] = ind
		                peaks.loc[ind+1,'Divergent Distance'] = div_dist		
	if is_plot:
	 	#Fig 1. Divergent distance. 
	 	f = plt.figure()
	 	peaks[peaks['has Divergent']]['Divergent Distance'].hist()
		num_div = np.sum(peaks['has Divergent'])
		total_peaks = len(peaks)
		title = '%s divergent peaks out of %s total peaks (divergent counts both forward and reverse)' % (num_div,total_peaks)
		plt.title(title)
		if not f_name  == '':
			plt.savefig(f_name + '_divergent.png')
			plt.close()
		# Fig 2. Distance between peaks
		f = plt.figure()
	 	peaks['Distance to Next Peak'].hist()
		num_div = np.sum(peaks['Distance to Next Peak'])
		total_peaks = len(peaks)
		title = 'Distance to the next peak (same or opposite strand)'
		plt.title(title)
		if not f_name  == '':
			plt.savefig(f_name + '_adjacent_distance.png')
			plt.close()
	if is_save:
		if f_name == '':
			print 'No name to save in. Not saving to csv' 
		else:
			peaks.to_csv(f_name + '.csv')                
	return peaks


#######################################################
def wrap_stats_all_in_one(in_file,out_file,multiple_peaks_width,divergent_width,type_merge=None):
	''' Loads the peaks, filters based on which peaks to keep, and calls stats_all_in_one '''

	peaks = pd.read_csv(in_file,sep='\t')

	if type_merge == 'all' and 'Parent files' in peaks.columns.values:
		#try: 
		num_breaks = max(peaks['Parent files'].apply(lambda x: x.count('|')))
		filt = '^([^\|]*\|){%d}[^\|]*$' % (num_breaks)
		peaks = peaks[peaks['Parent files'].str.contains(filt)]
		print 'peaks filtered here'
	
	if len(peaks) == 0:
		print 'Peak file empty: ', in_file
		return

	stats_all_in_one(peaks,multiple_peaks_width,divergent_width,f_name = out_file, is_save = True, is_plot = True)

