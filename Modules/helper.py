### Helper functions
def list_or_single(item):
    ''' Returns list if item is single element, returns the item itself if a list '''
    try:
        len(item)
    except TypeError:
        item = [item]
    return item

########################################


########################################
###########
### These two functions allow for parallelization of pandas DataFrames
def multi_run_wrapper(args):    
    curr_func = args[0]
    return curr_func(args[1],*args[2])


def parallel_df(df,func,func_args = (),num_processes= 4):
    '''
    Function to parallelize a function that loops over dataframes
    
    ''' 
    import pandas as pd
    import multiprocessing
    import itertools
    
    # calculate the chunk size as an integer
    chunk_size = int(df.shape[0]/num_processes)

    # this solution was reworked from the above link.
    # will work even if the length of the dataframe is not evenly divisible by num_processes
    chunks = [df.loc[df.index[i:i + chunk_size]] for i in range(0, df.shape[0], chunk_size)]
    
    # create our pool with `num_processes` processes
    pool = multiprocessing.Pool(processes=num_processes)

    
    # apply our function to each chunk in the list
    result = pool.map(multi_run_wrapper,itertools.izip(itertools.repeat(func),chunks, itertools.repeat(func_args)))
    
    #print result[0]

    new_df = pd.DataFrame()
    for i in result:
        # since i is just a dataframe
        # we can reassign the original dataframe based on the index of each chunk
        new_df = new_df.append(i)
        #print df.loc[result[i].index,:]
        #df.loc[result[i].index,:] = result[i]
    pool.terminate()
    return new_df

########################################