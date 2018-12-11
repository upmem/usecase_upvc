class GenomeParser(object):
	def __init__(self, filename):
		self.filename = filename
		with open(self.filename, 'r') as ref:
			self.headerlen = len(ref.readline())
			self.readlen = len(ref.readline().rstrip())
		#print("Read length of {}: {}".format(self.filename, str(self.readlen)))

	def getSeek(self,pos):
		pos-=1
		lineNumber = int(pos / self.readlen)
		offset = pos % self.readlen
		return self.headerlen + (self.readlen + 1) * lineNumber + offset

	def getNuc(self,line):
		with open(self.filename, 'r') as ref:
			if line == "\n":
				pass
			else:
				try:
					pos = int(line)
					ref.seek(self.getSeek(pos))
					return ref.read(1)
				except ValueError:
					(startPos, endPos) = line.split('-')
					(startPos, endPos) = (int(startPos), int(endPos))
					output = ""
					for i in range(startPos, endPos+1):
						ref.seek(self.getSeek(i))
						output += ref.read(1)
					return output

