## Incucyte_ZOOM_data_processing ##updated: 11/20/2019Author: Dr. Yunkai Zhang, Vanderbilt University Medical Center### This Python script processes the IncuCyte ZOOM data ###Input: the data exported from IncuCyte ZOOM, usually a txt file. *Current, this scirpt do not support 'break data down into individual images'.* Output: A dataframe of DIP-rate (DIP: Drug-induced proliferation). Usage: import the scirpt, then call *'data_processing'* function.Parameters:- file: *str*, filename.- platemap_plot: *bool, default False*. Plot for the whole plate or not.- endremove: *int, default 2*. The end point to remove.- start_row: *int, default 2*. The smallest row number.- start_col: *str, default 'B'*. The smallest col number.