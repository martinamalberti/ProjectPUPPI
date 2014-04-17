import os
import sys
import numpy
import glob

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir",dest="inputdir",type="string",help="Path where ntuples are located. Example: /store/cmst3/user/malberti/HIGGS/VERTEX/2012/DATA/")
parser.add_option("-w","--workdir",dest="workdir",type="string",default="EffJobs",help="Name of the directory for jobs")
parser.add_option("-o","--outputname",dest="outputname",type="string",default="testEfficiency",help="Name of the output file")
parser.add_option("-n","--njobs",dest="njobs",type="int",help="Number of jobs")
parser.add_option("-m","--maxEvents",dest="maxEvents",type="int",default=-1,help="Max number of events per job. If maxEvents=-1, all events are analyzed.")
parser.add_option("-r","--radius",dest="radius",type="float",default=0.5,help="Jet radius")
parser.add_option("-q","--queue",dest="queue",type="string",default="1nh",help="Name of the queue on lxbatch")
parser.add_option("","--checkJobs",dest="checkJobs",action="store_true",default=False,help="Checks job status")
parser.add_option("","--resubmit",dest="resubmit",action="store_true",default=False,help="Resubmit job ")
parser.add_option("-j","--job",dest="job",type="int",help="Job number (for resubmission). You can resubmit one job at time for now.")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true",help="Do not submit jobs")
parser.add_option("-e","--executable",dest="executable",type="string",default="MyVertexAnalysis.exe",help="Name of the executable. Default is: MyVertexAnalysis.exe")

(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'


def makeFilesList(indir,wdir):
    list = []
    command = ('%s find -f %s | grep root > %s/list.txt' % (eos,indir,wdir))
    #print command
    os.system(command)
    file = open('%s/list.txt'%wdir, 'r')
    list =[line.replace('/eos/cms/','root://eoscms//eos/cms/').replace('\n','') for line in file ]
    print 'Found %d files' %len(list)
    #print list             
    return list



def writeJobs(indir,wdir,njobs,maxevents,analysis,r,output):
    #---------------------------------------------
    # --- prepare the list of files to be analyzed
    #---------------------------------------------
    listoffiles = []
    listoffiles = makeFilesList(indir,wdir)

    #---------------------------------------------
    # --- now split the jobs
    #---------------------------------------------
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        os.mkdir(jobdir)
        
        #--- prepare the list of files for each job
        f = open('%s/input_%d.txt'%(jobdir,job), 'w')
        sublist = [file for i,file in enumerate(listoffiles) if (i%njobs==job)]
        for fname in sublist:
            f.write('%s \n'%fname)

        #--- prepare the jobs scripts
        jobscript = open('%s/sub_%d.sh'%(jobdir,job),'w')
        jobscript.write('cd %s \n'%jobdir)
#        jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc462 \n')
        jobscript.write('export SCRAM_ARCH=slc6_amd64_gcc472 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:%s/../../Bacon/BaconAna/DataFormatsOffline/ \n'%wdir)
        jobscript.write('cd - \n')
        jobscript.write('cp %s/../../%s ./ \n'%(wdir,analysis))
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/sub_%d.run \n'%(jobdir,job))
        jobscript.write('\t ./%s %d %s/input_%d.txt %s_%d %f'%(analysis,maxevents,jobdir,job,output,job,r))
        jobscript.write(') then \n')
#        jobscript.write('\t mv ./outtree_%s_%d.root %s \n'%(output,job,jobdir))
        jobscript.write('\t cp ./outtree_%s_%d.root %s \n'%(output,job,jobdir))
        jobscript.write('\t touch %s/sub_%d.done \n'%(jobdir,job))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/sub_%d.fail \n'%(jobdir,job))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/sub_%d.sh'%(jobdir,job))
   

def submitJobs(wdir, njobs, queue):
    for job in range(njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        jobname = '%s/sub_%d.sh'%(jobdir,job)
        print 'bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname )
        os.system('bsub -q %s -o %s/sub_%d.log %s'%(queue,jobdir,job,jobname ))
        


def checkJobs(wdir,output, queue):
    jobs = glob.glob( '%s/JOB_*/sub*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)
    
    listdone = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.done' % (wdir,j,j))]
    print 'Total number of DONE jobs: %s ' % len(listdone)
    print '  %s' %listdone
    for j in listdone:
        f = '%s/JOB_%d/sub_%d.run'%(wdir,j,j)
        #print 'rm %s/JOB_%d/sub_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
            
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
    print '   %s' %listrun
            
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/sub_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed
                        
    if (len(listdone) == len(jobs)):
        print "All jobs successful! Merging output files..."
        os.chdir(wdir)
        #        os.system('ls JOB*/*.root')
        os.system('hadd -f outtree_%s.root JOB_*/*.root'%output)
    else:
        for j in listfailed:
            print 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(queue,wdir,j,j,wdir,j,j )
                    


#-----------------------------------------
#--- MAIN
#-----------------------------------------

path = os.getcwd()
workingdir = path+'/'+options.workdir

if not options.checkJobs and not options.resubmit:

    # -- write jobs scripts
    # N.B can add an option for the radius 
    os.mkdir(workingdir)
    writeJobs(options.inputdir, workingdir, options.njobs, options.maxEvents,options.executable, options.radius,options.outputname )

    # -- submit jobs 
    if not options.dryRun:
        submitJobs(workingdir, options.njobs, options.queue)
        
elif options.resubmit and options.job >-1 :
    print 'Resubmitting job %d ' %options.job
    resubcmd = 'bsub -q %s -o %s/JOB_%d/sub_%d.log %s/JOB_%d/sub_%d.sh'%(options.queue,workingdir,options.job,options.job,wdir,options.job,options.job )
    #print resubcmd
    os.system(resubcmd)

elif options.checkJobs:
    checkJobs(workingdir,options.outputname, options.queue)
    
