### Helper functions
def list_or_single(item):
    ''' Returns list if item is single element, returns the item itself if a list '''
    try:
        len(item)
    except TypeError:
        item = [item]
    return item

########################################
def createPeakFileFromGFF(annotation_file,output_file,anno_of_interest='mRNA'):
    ''' 
    Function to create a peak file based on GFF and specific annotation of interest
	
	Input
	annotation_file: gff3 file
	output_file: File to save to
	anno_of_interest: a text term that has to be in the 'Annotation' column of the annotation file
    
    Creates a peak file of the Homer format and saves it to output_file

    '''
    #Read in annotation file
    genome_ann = pd.read_csv(annotation_file,skiprows=1830,sep='\t',header=None)
    col_names = list(genome_ann.columns.values)
    col_names[0]= 'Chr'
    col_names[1] = 'How'
    col_names[2] = 'Annotation'
    col_names[3]= 'Start'
    col_names[4] = 'End'
    col_names[6] = 'Strand'
    col_names[8] = 'ID'

    genome_ann.columns = col_names
    genome_ann.sort_values(['Chr','Start','End'],inplace=True)
    
    curr = genome_ann[genome_ann['Annotation'] == anno_of_interest]
    curr = curr[['ID','Chr','Start','End','Strand']]
    curr.to_csv(output_file,index=None,sep='\t')
    return
########################################

########################################
def peakFileToPeakFile(desired_peak_file,tag_peak_file,distance=1000,is_save=1):
    ''' 
    Function to filter the peaks in a file to only the ones that are within a certain distance
    of at least one peak in another file. 
	
	Input:
	desired_peak_file: File which contains the peaks to be filtered
	tag_peak_file: File that contains peaks to see if nearby the desired peaks
	distance: Maximum distance to consider between peaks
	is_save: To save the file or not. If on, will save peaks to file desired_peak_file+'filt'

	Output:
	updated_desired_peaks: The filtered peaks.
    '''
    desired_peaks = pd.read_csv(desired_peak_file,sep='\t')
    tag_peaks = pd.read_csv(tag_peak_file,sep='\t')
    inds_to_keep = []
    for ind,val in desired_peaks.iterrows():
        curr = tag_peaks[tag_peaks['Chr'] == val['Chr']]
        diff =  val['Start'] - curr['Start']
        if np.sum(np.abs(diff)<distance) > 0:
            inds_to_keep.append(ind)
    
    updated_desired_peaks = desired_peaks.loc[inds_to_keep]
    if is_save:
        updated_desired_peaks.to_csv(desired_peak_file + 'filt',index=None,sep='\t')
    return updated_desired_peaks
########################################

