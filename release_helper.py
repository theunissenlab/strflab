import string
import os
import sys
import re
import shutil


def list_files(fileList, rootDir, origRootDir, excludeExpr):

   for fname in os.listdir(rootDir):

      fullName = os.path.join(rootDir, fname)

      if re.search(excludeExpr, fname) is None:
      
         if os.path.isdir(fullName):
            list_files(fileList, fullName, origRootDir, excludeExpr)
         else:
            sindx = len(origRootDir)
            fileList.append(fullName[sindx+1:])


def print_usage():

   print 'This script helps with releases of STRFLAB. Make sure you run it from the STRFLAB trunk.'
   print '\nUsages:'
   print '\tpython release_helper.py detect'
   print '\tpython release_helper.py create <release_dir>'
   print '\tpython release_helper.py tag <release_dir>'


def copy_files_relative(sourceDir, relativeFilePaths, destDir):

      for rfile in relativeFilePaths:

         sourceName = os.path.join(sourceDir, rfile)

         if not os.path.isdir(sourceName):

            #make sure directories exist
            [relPath, fileName] = os.path.split(rfile)
            dirName = os.path.join(destDir, relPath)
            try:
               os.makedirs(dirName)
            except OSError:
               pass

            outputName = os.path.join(destDir, rfile)
            shutil.copy2(sourceName, outputName)



if __name__ == '__main__':

   if len(sys.argv) < 2:
      print_usage()
      exit()

   cmd = sys.argv[1]

   if cmd == 'detect':

      releaseFileName = 'strflab_release_includes.txt'
      if not os.path.exists(releaseFileName):
         print 'Cannot find strflab_release_includes.txt, are you running this script from the trunk?'
         exit()

      f = open(releaseFileName)   
      releaseFiles = [l.strip() for l in f.readlines()]

      wdir = os.getcwd()

      print '\nThese files are listed in strflab_release_includes.txt, but no longer exist! Make sure to remove them:'
      for rfile in releaseFiles:
         fullPath = os.path.join(wdir, rfile)
         if not os.path.exists(fullPath):
            print rfile

      print '\nThese files are not currently found in strflab_release_includes.txt, you may want to add some of them:'
      currentFiles = []
      list_files(currentFiles, wdir, wdir, '\.svn|.m~')
      for cfile in currentFiles:         
         if cfile not in releaseFiles:      
            print cfile


   elif cmd == 'create':

      if len(sys.argv) < 3:
         print 'Creation requires the uncreated release directory as the last argument.'
         exit()

      releaseDir = sys.argv[2]

      releaseFileName = 'strflab_release_includes.txt'
      if not os.path.exists(releaseFileName):
         print 'Cannot find strflab_release_includes.txt, are you running this script from the trunk?'
         exit()

      f = open(releaseFileName)   
      releaseFiles = [l.strip() for l in f.readlines()]

      wdir = os.getcwd()

      try:
         os.mkdir(releaseDir)
      except OSError:
         print 'Could not make release directory %s, it probably already exists...' % releaseDir

      copy_files_relative(wdir, releaseFiles, releaseDir)


   else:
      print 'Invalid command, valid commands are detect or create'
      exit()

   
