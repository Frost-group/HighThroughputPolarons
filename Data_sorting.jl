include("Data_workflow.jl")
#%%
# Commands for saving data

#df_multimode = multi_mode() # Too long to calculate! Use multiprocessing unit to do so.
df_Feynman = Feynman_Data()
df_Standard = Standard_Data()
df_General = General_Data()
df_General_final = General_Comparison_Data(df_General)

#%%
df_multi = multi_mode_sorting(0)
df_multi_0_K = multi_mode_sorting(1)

#%%
df_multi_single = multi_single_comparison()