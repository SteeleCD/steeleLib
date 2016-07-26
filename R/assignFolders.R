assignFolders = function()
	{
	user = Sys.info()['user']
	if(user=='c.steele')
		{
		# MAC
		patholDir='/Volumes/pathology$/'
		ucbtcdsDir='/Volumes/ucbtcdsh$/'
		commonDir='/Volumes/common$/'
		grpDir='/Volumes/grp11$/CI_Pathology_Steele/'
		dropDir='/Users/c.steele/Dropbox/'
		} else if(user=='chris') {
		# WORK
		patholDir='/media/pathology/'
		ucbtcdsDir='/media/ucbtcds/'
		commonDir='/media/common/'
		grpDir='/media/grp11/CI_Pathology_Steele/'
		dropDir='/home/chris/Dropbox/'
		} else {
		# LEGION
		scratchDir='/scratch/scratch/ucbtcds'
		}
	if(user%in%c('c.steele','chris'))
		{
		return(list(patholDir=patholDir,
			ucbtcdsDir=ucbtcdsDir,
			commonDir=commonDir,
			grpDir=grpDir,
			dropDir=dropDir))
		} else {
		return(list(scratchDir=scratchDir))
		}
	}

