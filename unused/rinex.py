import georinex as gr
# gr.crx2rnx
import glob
# for name in glob.glob('./archive/*.crx'):
print('stopped')
nav = gr.load("/Users/sabrinaberger/DRAO00CAN_R_20222410000_15M_01S_MO.crx", fast=True)
print(nav)
nav.sel(sv='G13')
print(nav.sel(sv='G13').values())
