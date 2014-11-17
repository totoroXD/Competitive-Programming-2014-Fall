import sys,os
import cStringIO
import path
output = cStringIO.StringIO()


codes = open('codebook.md', 'w')
for filename in path.path:
	f = open(filename, 'r')
	# print f.read()
	output.write('#'+filename+'\n``` cpp\n')
	output.write(f.read())
	output.write('\n```\n')

codes.write(output.getvalue())