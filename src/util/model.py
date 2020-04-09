class Model:

    def __init__(self):
        self.fileName = None
        self.fileContent = ''

    def is_valid(self, fname):
        try:
            file = open(fname, 'r')
            file.close()
            return True
        except FileNotFoundError:
            return False

    def set_filename(self, fname):
        if self.is_valid(fname):
            self.fileName = fname
            self.fileContents = open(fname, 'r').read()

    def get_filename(self):
        return self.fileName

    def get_filecontents(self):
        return self.fileContents

    def write_doc(self, text):
        if self.is_valid(self.fileName):
            fileName = self.fileName + ".bak"
            file = open(fileName, 'w')
            file.write(text)
            file.close()
