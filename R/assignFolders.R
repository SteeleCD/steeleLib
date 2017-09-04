# QOL function for setting various commmonly used directories across different PCs/clusters
assignFolders = function()
	{
	user = Sys.info()['user']
	if(user=='c.steele')
		{
		# MAC
		patholDir='/Volumes/grp13$/CI_POGB_Pathology/'
		ucbtcdsDir='/Volumes/ucbtcdsh$/'
		commonDir='/Volumes/common$/'
		grpDir='/Volumes/lgrp10$/CI_Pathology_Steele/'
		dropDir='/Users/c.steele/Dropbox/'
		} else if(user=='chris') {
		# WORK
		patholDir='/media/pathology/CI_POGB_Pathology/'
		ucbtcdsDir='/media/ucbtcds/'
		commonDir='/media/common/'
		grpDir='/media/grp11/CI_Pathology_Steele/'
		grpDir2='/run/user/1000/gvfs/smb-share:server=file03.ucl.ac.uk,share=grp11$/CI_Pathology_Steele/'
		dropDir='/home/chris/Dropbox/'
		} else {
		# LEGION
		scratchDir='/scratch/scratch/ucbtcds'
		}
	if(user%in%c('c.steele','chris'))
		{
		out = list(patholDir=patholDir,
			ucbtcdsDir=ucbtcdsDir,
			commonDir=commonDir,
			grpDir=grpDir,
			dropDir=dropDir)
		if(user=="chris") out = append(out,list(grpDir2=grpDir2))
		return(out)
		} else {
		return(list(scratchDir=scratchDir))
		}
	}

