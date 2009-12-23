import os, sys

# Program options
errorRate = float(sys.argv[3])   # 0.04 meaning a 4% single-base error rate
insertions = False # change to True to allow insertions
deletions = False  # change to True to allow deletions
Ns = False 	   # change to True to allow Ns (unknown base)

geneFile = sys.argv[1]

digFile = open(geneFile+'.dig', 'w')
digFile.write('countsamplerandom 350 600 %d 0\n' %(int(sys.argv[2])))
digFile.close()

if len(sys.argv) > 4:
	mstFile = sys.argv[4]
	estFile = sys.argv[5]
else:
	mstFile = '.'
	estFile = '.'

insNum = 0
delNum = 0
nNum   = 0
if insertions:
	insNum = 1
if deletions:
	delNum = 1
if Ns:
	nNum = 2

# Run the pipeline
if errorRate <= 0:
	os.system('python pipeline.py '+geneFile+' '+geneFile+'.dig 1 -p -c 1')
else:
	os.system('python pipeline.py '+geneFile+' '+geneFile+'.dig 1 sbe %.6f 20 %.6f 1 1 10 %d %d %d -p -c 1' %(errorRate, errorRate, insNum, delNum, nNum))

# Clean up output
if errorRate <= 0:
	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_]/peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_].mst '+mstFile)
	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_]/estsim_'+geneFile+'_'+geneFile+'.dig[1_]_fmt.fa '+estFile)
#	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_]/peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_].out clusters.out')
	os.system('rm -rf peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_]')
else:
	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_]' %(errorRate, errorRate, insNum, delNum, nNum) + '/peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_].mst ' %(errorRate, errorRate, insNum, delNum, nNum) + mstFile)
	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_]' %(errorRate, errorRate, insNum, delNum, nNum) + '/estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_]_fmt.fa '%(errorRate, errorRate, insNum, delNum, nNum) + estFile)
#	os.system('mv peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_]' %(errorRate, errorRate, insNum, delNum, nNum) + '/peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_].out clusters.out' %(errorRate, errorRate, insNum, delNum, nNum))
	os.system('rm -rf peace_estsim_'+geneFile+'_'+geneFile+'.dig[1_sbe_%.6f_20_%.6f_1_1_10_%d_%d_%d_]/' %(errorRate, errorRate, insNum, delNum, nNum))

#os.system('rm peace_scr*.job')
os.system('rm '+geneFile+'.dig')

# MST file will be in file mstFile
# ESTs will be in file estFile
# Clusters will be in file clusters.out
