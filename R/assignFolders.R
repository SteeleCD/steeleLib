assignFolders = function()
	{
		if(Sys.info()['sysname']=='Darwin')
		{
		patholDir='/Volumes/pathology$/'
		ucbtcdsDir='/Volumes/ucbtcdsh$/'
		commonDir='/Volumes/common$/'
		grpDir='/Volumes/grp11$/CI_Pathology_Steele/'
		dropDir='/Users/c.steele/Dropbox/'
		} else {
		patholDir='/media/pathology/'
		ucbtcdsDir='/media/ucbtcds/'
		commonDir='/media/common/'
		grpDir='/media/grp11/CI_Pathology_Steele/'
		dropDir='/home/chris/Dropbox/'
		}
	return(list(patholDir=patholDir,
			ucbtcdsDir=ucbtcdsDir,
			commonDir=commonDir,
			grpDir=grpDir,
			dropDir=dropDir))
	}

